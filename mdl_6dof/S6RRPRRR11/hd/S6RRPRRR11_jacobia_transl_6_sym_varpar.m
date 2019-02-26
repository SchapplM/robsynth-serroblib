% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:59:48
% EndTime: 2019-02-26 21:59:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (201->39), mult. (188->50), div. (0->0), fcn. (207->10), ass. (0->31)
t22 = sin(qJ(2));
t24 = cos(qJ(2));
t29 = pkin(2) + r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t28 = t29 * t24;
t21 = qJ(4) + qJ(5);
t16 = sin(t21);
t10 = pkin(5) * t16 + sin(qJ(4)) * pkin(4);
t30 = qJ(3) + t10;
t38 = -t30 * t22 - pkin(1) - t28;
t17 = cos(t21);
t13 = pkin(5) * t17;
t11 = t13 + cos(qJ(4)) * pkin(4);
t36 = pkin(7) + pkin(3) + t11;
t18 = qJ(6) + t21;
t14 = sin(t18);
t15 = cos(t18);
t23 = sin(qJ(1));
t25 = cos(qJ(1));
t31 = t25 * t22;
t5 = -t23 * t14 + t15 * t31;
t6 = t14 * t31 + t23 * t15;
t35 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t32 = t23 * t22;
t7 = t25 * t14 + t15 * t32;
t8 = -t14 * t32 + t25 * t15;
t34 = t7 * r_i_i_C(1) + t8 * r_i_i_C(2);
t33 = r_i_i_C(1) * t15;
t27 = r_i_i_C(1) * t14 + r_i_i_C(2) * t15 + t30;
t26 = -t29 * t22 + t27 * t24;
t12 = t24 * t14 * r_i_i_C(2);
t1 = [t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t38 * t23 + t36 * t25, t26 * t25, t31, -t23 * t10 + t11 * t31 + t35 (-t16 * t23 + t17 * t31) * pkin(5) + t35, t35; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t36 * t23 - t38 * t25, t26 * t23, t32, t25 * t10 + t11 * t32 + t34 (t16 * t25 + t17 * t32) * pkin(5) + t34, t34; 0, t27 * t22 + t28, -t24, t12 + (-t11 - t33) * t24, t12 + (-t33 - t13) * t24, -t24 * t33 + t12;];
Ja_transl  = t1;
