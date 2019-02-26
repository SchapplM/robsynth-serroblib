% Geometrischen Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S6RRRPRR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Jg [3x6]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:07
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S6RRRPRR8_jacobia_0_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S6RRRPRR8_jacobia_transl_0_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S6RRRPRR8_jacobig_rot_0_floatb_twist_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
