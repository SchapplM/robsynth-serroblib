% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR8
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:16
% EndTime: 2019-02-26 21:58:16
% DurationCPUTime: 0.14s
% Computational Cost: add. (279->40), mult. (177->54), div. (0->0), fcn. (196->12), ass. (0->33)
t26 = cos(qJ(2));
t24 = sin(qJ(2));
t37 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8) + qJ(3);
t31 = t37 * t24;
t23 = pkin(11) + qJ(4);
t21 = qJ(5) + t23;
t18 = cos(t21);
t12 = pkin(5) * t18 + pkin(4) * cos(t23);
t9 = cos(pkin(11)) * pkin(3) + pkin(2) + t12;
t42 = t26 * t9 + pkin(1) + t31;
t19 = qJ(6) + t21;
t14 = sin(t19);
t15 = cos(t19);
t27 = cos(qJ(1));
t33 = t27 * t15;
t25 = sin(qJ(1));
t36 = t25 * t26;
t5 = t14 * t36 + t33;
t34 = t27 * t14;
t6 = -t15 * t36 + t34;
t41 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t25 * t15 - t26 * t34;
t8 = t25 * t14 + t26 * t33;
t40 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t17 = sin(t21);
t39 = pkin(5) * t17;
t11 = -pkin(4) * sin(t23) - t39;
t38 = pkin(7) + sin(pkin(11)) * pkin(3) - t11;
t35 = t26 * t27;
t30 = -r_i_i_C(1) * t14 - r_i_i_C(2) * t15;
t29 = r_i_i_C(1) * t15 - r_i_i_C(2) * t14 + t9;
t28 = -t29 * t24 + t37 * t26;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t42 * t25 + t38 * t27, t28 * t27, t27 * t24, t11 * t35 + t25 * t12 + t40 (-t17 * t35 + t18 * t25) * pkin(5) + t40, t40; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t38 * t25 + t42 * t27, t28 * t25, t25 * t24, t11 * t36 - t27 * t12 + t41 (-t17 * t36 - t18 * t27) * pkin(5) + t41, t41; 0, t29 * t26 + t31, -t26 (t11 + t30) * t24 (t30 - t39) * t24, t30 * t24;];
Ja_transl  = t1;
