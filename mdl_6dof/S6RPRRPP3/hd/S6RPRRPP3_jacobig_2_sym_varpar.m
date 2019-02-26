% Geometrische Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% Jg [3x6]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S6RPRRPP3_jacobig_2_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S6RPRRPP3_jacobia_transl_2_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S6RPRRPP3_jacobig_rot_2_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
