% Geometrische Jacobi-Matrix für Segment Nr. 1 (0=Basis) von
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% Jg [3x6]
%   Geometrische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S6RRRRRR10V2_jacobig_1_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S6RRRRRR10V2_jacobia_transl_1_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S6RRRRRR10V2_jacobig_rot_1_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
