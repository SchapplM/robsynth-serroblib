% Geometrischen Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S4PPPR3
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
%   pkin=[a2,a3,a4,d4,theta3]';
%
% Output:
% Jg [3x4]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:56
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg = S4PPPR3_jacobia_0_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S4PPPR3_jacobia_transl_0_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S4PPPR3_jacobig_rot_0_floatb_twist_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];