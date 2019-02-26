% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_transl [3x7]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-26 21:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S7RRRRRRR1_jacobia_transl_3_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_transl_3_floatb_twist_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobia_transl_3_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_transl_3_floatb_twist_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-26 21:20:52
% EndTime: 2018-11-26 21:20:53
% DurationCPUTime: 0.08s
% Computational Cost: add. (26->15), mult. (70->31), div. (0->0), fcn. (78->6), ass. (0->17)
t15 = pkin(2) + r_i_i_C(3);
t6 = sin(qJ(2));
t13 = t15 * t6;
t7 = sin(qJ(1));
t9 = cos(qJ(2));
t16 = t7 * t9;
t10 = cos(qJ(1));
t14 = t10 * t9;
t5 = sin(qJ(3));
t8 = cos(qJ(3));
t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5;
t11 = -t12 * t6 - t15 * t9;
t4 = t8 * t14 - t7 * t5;
t3 = -t5 * t14 - t7 * t8;
t2 = -t10 * t5 - t8 * t16;
t1 = -t10 * t8 + t5 * t16;
t17 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * t13, t11 * t10, t3 * r_i_i_C(1) - t4 * r_i_i_C(2), 0, 0, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t10 * t13, t11 * t7, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2), 0, 0, 0, 0; 0, t12 * t9 - t13 (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6, 0, 0, 0, 0;];
Ja_transl  = t17;
