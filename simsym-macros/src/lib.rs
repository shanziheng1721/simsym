use proc_macro::TokenStream;
use proc_macro2::{Span, TokenTree};
use quote::quote;
use syn::Ident;

fn crate_ident() -> Ident {
    Ident::new("simsym", Span::call_site())
}

#[proc_macro]
pub fn expr(input: TokenStream) -> TokenStream {
    let tokens: Vec<TokenTree> = proc_macro2::TokenStream::from(input).into_iter().collect();
    let mut pos = 0;
    match parse_expr(&tokens, &mut pos) {
        Ok(ts) if pos == tokens.len() => ts.into(),
        Ok(_) => syn::Error::new(Span::call_site(), "trailing tokens in expr!")
            .to_compile_error()
            .into(),
        Err(e) => e.to_compile_error().into(),
    }
}

fn parse_expr(tokens: &[TokenTree], pos: &mut usize) -> syn::Result<proc_macro2::TokenStream> {
    parse_additive(tokens, pos)
}

fn parse_additive(tokens: &[TokenTree], pos: &mut usize) -> syn::Result<proc_macro2::TokenStream> {
    let root = crate_ident();
    let mut left = parse_multiplicative(tokens, pos)?;
    while let Some(op) = peek_op(tokens, *pos) {
        if op != '+' && op != '-' {
            break;
        }
        *pos += 1;
        let right = parse_multiplicative(tokens, pos)?;
        left = if op == '+' {
            quote! { #root::expr::add(#left, #right) }
        } else {
            quote! { #root::expr::sub(#left, #right) }
        };
    }
    Ok(left)
}

fn parse_multiplicative(tokens: &[TokenTree], pos: &mut usize) -> syn::Result<proc_macro2::TokenStream> {
    let root = crate_ident();
    let mut left = parse_power(tokens, pos)?;
    while let Some(op) = peek_op(tokens, *pos) {
        if op != '*' && op != '/' {
            break;
        }
        *pos += 1;
        let right = parse_power(tokens, pos)?;
        left = if op == '*' {
            quote! { #root::expr::mul(#left, #right) }
        } else {
            quote! { #root::expr::div(#left, #right) }
        };
    }
    Ok(left)
}

fn parse_power(tokens: &[TokenTree], pos: &mut usize) -> syn::Result<proc_macro2::TokenStream> {
    let root = crate_ident();
    let mut left = parse_unary(tokens, pos)?;
    if peek_op(tokens, *pos) == Some('^') {
        *pos += 1;
        let right = parse_power(tokens, pos)?;
        left = quote! { #root::expr::pow(#left, #right) };
    }
    Ok(left)
}

fn parse_unary(tokens: &[TokenTree], pos: &mut usize) -> syn::Result<proc_macro2::TokenStream> {
    let root = crate_ident();
    if peek_op(tokens, *pos) == Some('-') {
        *pos += 1;
        let inner = parse_unary(tokens, pos)?;
        return Ok(quote! { #root::expr::neg(#inner) });
    }
    parse_atom(tokens, pos)
}

fn parse_atom(tokens: &[TokenTree], pos: &mut usize) -> syn::Result<proc_macro2::TokenStream> {
    let root = crate_ident();
    let Some(tok) = tokens.get(*pos) else {
        return Err(syn::Error::new(Span::call_site(), "unexpected end of expr!"));
    };
    match tok {
        TokenTree::Group(g) if g.delimiter() == proc_macro2::Delimiter::Parenthesis => {
            *pos += 1;
            let inner: Vec<TokenTree> = g.stream().into_iter().collect();
            let mut p = 0;
            let e = parse_expr(&inner, &mut p)?;
            if p != inner.len() {
                return Err(syn::Error::new(g.span(), "trailing tokens in parentheses"));
            }
            Ok(e)
        }
        TokenTree::Ident(id) => {
            let name = id.to_string();
            if matches!(
                name.as_str(),
                "sin" | "cos" | "tan" | "cot" | "sec" | "csc"
                    | "asin" | "acos" | "atan" | "acot" | "asec" | "acsc"
                    | "sinh" | "cosh" | "tanh" | "coth" | "sech" | "csch"
                    | "asinh" | "acosh" | "atanh" | "acoth" | "asech" | "acsch"
                    | "exp" | "ln"
            ) {
                *pos += 1;
                let args = parse_paren_args(tokens, pos)?;
                let fname = syn::Ident::new(&name, id.span());
                return Ok(quote! { #root::expr::#fname(#args) });
            }
            // `e^x` is exp(x), not the symbol e raised to x
            if name == "e" && peek_op(tokens, *pos) == Some('^') {
                *pos += 1; // ^
                let exp = parse_power(tokens, pos)?;
                return Ok(quote! { #root::expr::exp(#exp) });
            }
            *pos += 1;
            Ok(quote! { #root::expr::var(#root::symbol(#name)) })
        }
        TokenTree::Literal(lit) => {
            *pos += 1;
            let s = lit.to_string();
            if s.contains('.') {
                let v: f64 = s.parse().map_err(|_| lit_err(lit))?;
                let n = (v * 1_000_000.0).round() as i64;
                Ok(quote! { #root::expr::const_(#root::rational(#n, 1_000_000i64)) })
            } else {
                let v: i64 = s.parse().map_err(|_| lit_err(lit))?;
                Ok(quote! { #root::expr::const_(#root::rational(#v, 1i64)) })
            }
        }
        _ => Err(syn::Error::new(tok.span(), "expected atom in expr!")),
    }
}

fn parse_paren_args(tokens: &[TokenTree], pos: &mut usize) -> syn::Result<proc_macro2::TokenStream> {
    let Some(TokenTree::Group(g)) = tokens.get(*pos) else {
        return Err(syn::Error::new(Span::call_site(), "expected '('"));
    };
    if g.delimiter() != proc_macro2::Delimiter::Parenthesis {
        return Err(syn::Error::new(g.span(), "expected '('"));
    }
    *pos += 1;
    let inner: Vec<TokenTree> = g.stream().into_iter().collect();
    let mut p = 0;
    let e = parse_expr(&inner, &mut p)?;
    if p != inner.len() {
        return Err(syn::Error::new(g.span(), "trailing tokens in function call"));
    }
    Ok(e)
}

fn peek_op(tokens: &[TokenTree], pos: usize) -> Option<char> {
    match tokens.get(pos)? {
        TokenTree::Punct(p) if p.spacing() == proc_macro2::Spacing::Alone => {
            Some(p.as_char())
        }
        _ => None,
    }
}

fn lit_err(lit: &proc_macro2::Literal) -> syn::Error {
    syn::Error::new(lit.span(), "invalid numeric literal")
}
