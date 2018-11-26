% Geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Jg [3x6]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:45
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg = S6RRRPRP6_jacobig_5_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S6RRRPRP6_jacobia_transl_5_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S6RRRPRP6_jacobig_rot_5_floatb_twist_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
