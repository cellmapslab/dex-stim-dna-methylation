//
//  ContentView.swift
//  Shared
//
//  Created by Hryhorzhevska, Anastasiia on 25.12.20.
//

import SwiftUI

struct ContentView: View {
    @Binding var document: dex_stim_dna_methylationDocument

    var body: some View {
        TextEditor(text: $document.text)
    }
}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView(document: .constant(dex_stim_dna_methylationDocument()))
    }
}
