//
//  dex_stim_dna_methylationApp.swift
//  Shared
//
//  Created by Hryhorzhevska, Anastasiia on 25.12.20.
//

import SwiftUI

@main
struct dex_stim_dna_methylationApp: App {
    var body: some Scene {
        DocumentGroup(newDocument: dex_stim_dna_methylationDocument()) { file in
            ContentView(document: file.$document)
        }
    }
}
