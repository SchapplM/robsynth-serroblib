% Geometrischen Jacobi-Matrix für Segment Nr. 1 (0=Basis) von
% S4PRPP4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
%
% Output:
% Jg [3x4]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:09
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg = S4PRPP4_jacobia_1_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S4PRPP4_jacobia_transl_1_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S4PRPP4_jacobig_rot_1_floatb_twist_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];