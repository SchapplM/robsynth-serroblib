% Geometrische Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% Jg [3x5]
%   Geometrische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S5RRRRR2_jacobig_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S5RRRRR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S5RRRRR2_jacobig_rot_5_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
