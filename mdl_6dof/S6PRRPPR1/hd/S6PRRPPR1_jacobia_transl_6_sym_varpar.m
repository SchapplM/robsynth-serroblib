% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:58:13
% EndTime: 2019-02-26 19:58:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (234->45), mult. (334->76), div. (0->0), fcn. (420->14), ass. (0->36)
t21 = qJ(3) + pkin(11);
t17 = sin(t21);
t19 = cos(t21);
t30 = cos(qJ(3));
t20 = pkin(12) + qJ(6);
t16 = sin(t20);
t18 = cos(t20);
t34 = t18 * r_i_i_C(1) - t16 * r_i_i_C(2) + cos(pkin(12)) * pkin(5) + pkin(4);
t43 = r_i_i_C(3) + pkin(9) + qJ(5);
t44 = t30 * pkin(3) + t43 * t17 + t34 * t19 + pkin(2);
t23 = sin(pkin(10));
t24 = sin(pkin(6));
t42 = t23 * t24;
t25 = cos(pkin(10));
t41 = t24 * t25;
t29 = sin(qJ(2));
t40 = t24 * t29;
t39 = t24 * t30;
t31 = cos(qJ(2));
t38 = t24 * t31;
t37 = cos(pkin(6));
t36 = t29 * t37;
t35 = t31 * t37;
t33 = sin(pkin(12)) * pkin(5) + t16 * r_i_i_C(1) + t18 * r_i_i_C(2) + qJ(4) + pkin(8);
t28 = sin(qJ(3));
t10 = -t23 * t36 + t25 * t31;
t9 = t23 * t35 + t25 * t29;
t8 = t23 * t31 + t25 * t36;
t7 = t23 * t29 - t25 * t35;
t6 = t37 * t17 + t19 * t40;
t5 = t17 * t40 - t37 * t19;
t4 = t10 * t19 + t17 * t42;
t3 = t10 * t17 - t19 * t42;
t2 = -t17 * t41 + t8 * t19;
t1 = t8 * t17 + t19 * t41;
t11 = [0, t33 * t10 - t44 * t9, t43 * t4 + (-t10 * t28 + t23 * t39) * pkin(3) - t34 * t3, t9, t3 (-t4 * t16 + t9 * t18) * r_i_i_C(1) + (-t9 * t16 - t4 * t18) * r_i_i_C(2); 0, t33 * t8 - t44 * t7, t43 * t2 + (-t25 * t39 - t28 * t8) * pkin(3) - t34 * t1, t7, t1 (-t2 * t16 + t7 * t18) * r_i_i_C(1) + (-t7 * t16 - t2 * t18) * r_i_i_C(2); 1 (t33 * t29 + t44 * t31) * t24, t43 * t6 + (-t28 * t40 + t37 * t30) * pkin(3) - t34 * t5, -t38, t5 (-t6 * t16 - t18 * t38) * r_i_i_C(1) + (t16 * t38 - t6 * t18) * r_i_i_C(2);];
Ja_transl  = t11;
