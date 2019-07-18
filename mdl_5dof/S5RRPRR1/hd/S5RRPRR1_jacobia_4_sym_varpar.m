% Analytische Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% Ja [6x5]
%   Analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = S5RRPRR1_jacobia_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)

Ja_transl = S5RRPRR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Ja_rot = S5RRPRR1_jacobia_rot_4_sym_varpar(qJ, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
