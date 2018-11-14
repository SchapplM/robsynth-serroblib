% Geometrischen Jacobi-Matrix für Segment Nr. 2 (0=Basis) von
% S3PRP2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
%
% Output:
% Jg [3x3]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg = S3PRP2_jacobia_2_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S3PRP2_jacobia_transl_2_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S3PRP2_jacobig_rot_2_floatb_twist_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
