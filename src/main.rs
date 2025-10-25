mod args;
mod long;
mod short;
use crate::args::CommandParse;
use crate::args::Commands;
use crate::long::longread;
use crate::short::shortread;
use clap::Parser;
use figlet_rs::FIGfont;

/*
 Author Gaurav Sablok,
 Email: codeprog@icloud.com
*/

fn main() {
    let fontgenerate = FIGfont::standard().unwrap();
    let repgenerate = fontgenerate.convert("bactiPAN");
    println!("{}", repgenerate.unwrap());
    let argsparse = CommandParse::parse();
    match &argsparse.command {
        Commands::Selectoption { optionvalue } => {
            let value = optionvalue;
            if optionvalue.to_string() == "short" {
                let _ = shortread().unwrap();
            } else if optionvalue.to_string() == "long" {
                let _ = longread().unwrap();
            }
        }
    }
}
