use std::collections::HashMap;
use std::fmt;
use std::sync::{Mutex, OnceLock};

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug)]
pub struct Symbol(u32);

struct InternTable {
    names: Vec<&'static str>,
    index: HashMap<String, u32>,
}

fn table() -> &'static Mutex<InternTable> {
    static TABLE: OnceLock<Mutex<InternTable>> = OnceLock::new();
    TABLE.get_or_init(|| {
        Mutex::new(InternTable {
            names: Vec::new(),
            index: HashMap::new(),
        })
    })
}

/// Create or reuse an interned symbol by name.
pub fn symbol(name: &str) -> Symbol {
    let mut guard = table().lock().expect("symbol intern table poisoned");
    if let Some(&id) = guard.index.get(name) {
        return Symbol(id);
    }
    let id = guard.names.len() as u32;
    let leaked: &'static str = Box::leak(name.to_owned().into_boxed_str());
    guard.names.push(leaked);
    guard.index.insert(name.to_owned(), id);
    Symbol(id)
}

impl Symbol {
    pub fn name(self) -> &'static str {
        let guard = table().lock().expect("symbol intern table poisoned");
        guard.names[self.0 as usize]
    }
}

impl fmt::Display for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.name())
    }
}
