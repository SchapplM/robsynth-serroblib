% Analytische Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S4PRPP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
%
% Output:
% Ja [6x4]
%   Analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = S4PRPP1_jacobia_0_sym_varpar(qJ, r_i_i_C, ...
  pkin)

Ja_transl = S4PRPP1_jacobia_transl_0_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Ja_rot = S4PRPP1_jacobia_rot_0_sym_varpar(qJ, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
