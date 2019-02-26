% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR13_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR13_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:49
% EndTime: 2019-02-26 21:44:50
% DurationCPUTime: 0.17s
% Computational Cost: add. (232->51), mult. (456->80), div. (0->0), fcn. (578->12), ass. (0->37)
t44 = r_i_i_C(3) + pkin(10) + qJ(5);
t24 = sin(pkin(6));
t27 = sin(qJ(2));
t43 = t24 * t27;
t28 = sin(qJ(1));
t42 = t24 * t28;
t30 = cos(qJ(2));
t41 = t24 * t30;
t31 = cos(qJ(1));
t40 = t24 * t31;
t39 = cos(pkin(6));
t38 = t24 * (pkin(3) + pkin(8));
t37 = t28 * t39;
t36 = t31 * t39;
t35 = sin(pkin(11)) * pkin(5) + pkin(2) + pkin(9);
t19 = cos(pkin(11)) * pkin(5) + pkin(4);
t22 = pkin(11) + qJ(6);
t20 = sin(t22);
t21 = cos(t22);
t34 = t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t19;
t13 = t28 * t27 - t30 * t36;
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t7 = -t13 * t26 + t29 * t40;
t5 = t13 * t29 + t26 * t40;
t33 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + t35;
t32 = t34 * t26 - t44 * t29 + qJ(3);
t16 = -t27 * t37 + t31 * t30;
t15 = t31 * t27 + t30 * t37;
t14 = t27 * t36 + t28 * t30;
t12 = -t26 * t41 + t39 * t29;
t11 = t39 * t26 + t29 * t41;
t4 = t15 * t26 + t29 * t42;
t3 = -t15 * t29 + t26 * t42;
t2 = t16 * t20 + t4 * t21;
t1 = t16 * t21 - t4 * t20;
t6 = [-t28 * pkin(1) - t13 * qJ(3) - t33 * t14 + t31 * t38 + t34 * t7 + t44 * t5, -t33 * t15 + t32 * t16, t15, -t34 * t3 + t44 * t4, t3, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t31 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * qJ(3) + t35 * t16 + t4 * t19 + t28 * t38 + t44 * t3, -t33 * t13 + t32 * t14, t13, t34 * t5 - t44 * t7, -t5 (t14 * t21 + t7 * t20) * r_i_i_C(1) + (-t14 * t20 + t7 * t21) * r_i_i_C(2); 0 (t32 * t27 + t33 * t30) * t24, -t41, -t34 * t11 + t44 * t12, t11 (-t12 * t20 + t21 * t43) * r_i_i_C(1) + (-t12 * t21 - t20 * t43) * r_i_i_C(2);];
Ja_transl  = t6;
