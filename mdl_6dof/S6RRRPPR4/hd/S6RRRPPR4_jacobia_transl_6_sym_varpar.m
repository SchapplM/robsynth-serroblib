% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPPR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:08
% EndTime: 2019-02-26 22:05:09
% DurationCPUTime: 0.19s
% Computational Cost: add. (246->48), mult. (314->72), div. (0->0), fcn. (380->10), ass. (0->36)
t22 = cos(qJ(3));
t12 = t22 * pkin(3) + pkin(2);
t23 = cos(qJ(2));
t19 = sin(qJ(2));
t34 = -r_i_i_C(3) - pkin(9) + qJ(4) + pkin(8);
t29 = t34 * t19;
t41 = t23 * t12 + pkin(1) + t29;
t15 = qJ(3) + pkin(10);
t13 = sin(t15);
t14 = cos(t15);
t21 = cos(qJ(6));
t17 = sin(qJ(6));
t39 = pkin(4) + pkin(5);
t30 = t17 * r_i_i_C(2) - t39;
t26 = t21 * r_i_i_C(1) - t30;
t31 = -t17 * r_i_i_C(1) - qJ(5);
t27 = t21 * r_i_i_C(2) - t31;
t40 = t27 * t13 + t26 * t14 + t12;
t18 = sin(qJ(3));
t38 = t18 * pkin(3);
t24 = cos(qJ(1));
t35 = t24 * t14;
t20 = sin(qJ(1));
t37 = t20 * t23;
t6 = t13 * t37 + t35;
t3 = t6 * t21;
t36 = t24 * t13;
t33 = pkin(7) + t38;
t28 = (-(-t13 * t21 + t14 * t17) * r_i_i_C(1) - (t13 * t17 + t14 * t21) * r_i_i_C(2)) * t19;
t25 = -t40 * t19 + t34 * t23;
t9 = t20 * t13 + t23 * t35;
t8 = -t20 * t14 + t23 * t36;
t7 = t14 * t37 - t36;
t2 = t8 * t17 + t9 * t21;
t1 = -t9 * t17 + t8 * t21;
t4 = [-t3 * r_i_i_C(2) - t41 * t20 + t33 * t24 - t26 * t7 + t31 * t6, t25 * t24, t27 * t9 - t26 * t8 + (-t24 * t23 * t18 + t20 * t22) * pkin(3), t24 * t19, t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t8 * qJ(5) + t33 * t20 + t41 * t24 + t39 * t9, t25 * t20, -t3 * r_i_i_C(1) + t27 * t7 + t30 * t6 + (-t18 * t37 - t24 * t22) * pkin(3), t20 * t19, t6 (-t7 * t17 + t3) * r_i_i_C(1) + (-t6 * t17 - t7 * t21) * r_i_i_C(2); 0, t40 * t23 + t29 (qJ(5) * t14 - t39 * t13 - t38) * t19 - t28, -t23, t19 * t13, t28;];
Ja_transl  = t4;
