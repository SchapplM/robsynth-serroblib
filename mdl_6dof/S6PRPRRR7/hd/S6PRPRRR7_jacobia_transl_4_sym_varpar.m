% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRR7_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_transl_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:28
% EndTime: 2019-02-26 19:57:29
% DurationCPUTime: 0.25s
% Computational Cost: add. (167->63), mult. (489->125), div. (0->0), fcn. (635->14), ass. (0->50)
t57 = pkin(10) + r_i_i_C(3);
t21 = sin(pkin(14));
t28 = cos(pkin(7));
t56 = t21 * t28;
t22 = sin(pkin(8));
t23 = sin(pkin(7));
t55 = t22 * t23;
t27 = cos(pkin(8));
t54 = t23 * t27;
t24 = sin(pkin(6));
t26 = cos(pkin(13));
t53 = t24 * t26;
t25 = cos(pkin(14));
t52 = t25 * t28;
t30 = sin(qJ(2));
t51 = t28 * t30;
t32 = cos(qJ(2));
t50 = t28 * t32;
t49 = t23 * qJ(3);
t48 = cos(pkin(6));
t47 = sin(pkin(13));
t46 = t24 * t47;
t45 = t26 * t48;
t44 = t48 * t23;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t43 = t29 * r_i_i_C(1) + t31 * r_i_i_C(2);
t16 = -t47 * t30 + t32 * t45;
t11 = -t16 * t23 - t28 * t53;
t17 = t30 * t45 + t47 * t32;
t34 = t16 * t28 - t23 * t53;
t42 = (-t17 * t21 + t34 * t25) * t27 + t11 * t22;
t39 = t48 * t47;
t18 = -t26 * t30 - t32 * t39;
t12 = -t18 * t23 + t28 * t46;
t19 = t26 * t32 - t30 * t39;
t33 = t18 * t28 + t23 * t46;
t41 = t12 * t22 + t27 * (-t19 * t21 + t33 * t25);
t15 = -t24 * t32 * t23 + t48 * t28;
t40 = t15 * t22 + t27 * (t25 * t44 + (-t21 * t30 + t25 * t50) * t24);
t5 = -t16 * t21 - t17 * t52;
t37 = t17 * t55 + t27 * t5;
t7 = -t18 * t21 - t19 * t52;
t35 = t19 * t55 + t27 * t7;
t10 = t24 * t30 * t25 + (t24 * t50 + t44) * t21;
t8 = t18 * t25 - t19 * t56;
t6 = t16 * t25 - t17 * t56;
t4 = t19 * t25 + t33 * t21;
t2 = t17 * t25 + t34 * t21;
t1 = [0 (t35 * t29 + t8 * t31) * r_i_i_C(1) + (-t8 * t29 + t35 * t31) * r_i_i_C(2) + t8 * pkin(3) + t18 * pkin(2) + t19 * t49 + t57 * (t19 * t54 - t7 * t22) t12 (-t4 * t29 + t41 * t31) * r_i_i_C(1) + (-t41 * t29 - t4 * t31) * r_i_i_C(2), 0, 0; 0 (t37 * t29 + t6 * t31) * r_i_i_C(1) + (-t6 * t29 + t37 * t31) * r_i_i_C(2) + t6 * pkin(3) + t16 * pkin(2) + t17 * t49 + t57 * (t17 * t54 - t5 * t22) t11 (-t2 * t29 + t42 * t31) * r_i_i_C(1) + (-t2 * t31 - t42 * t29) * r_i_i_C(2), 0, 0; 1 ((t31 * r_i_i_C(1) - t29 * r_i_i_C(2) + pkin(3)) * (-t21 * t51 + t25 * t32) + (-t57 * t22 + t43 * t27) * (-t21 * t32 - t25 * t51) + t32 * pkin(2) + (t43 * t22 + t57 * t27 + qJ(3)) * t30 * t23) * t24, t15 (-t10 * t29 + t40 * t31) * r_i_i_C(1) + (-t10 * t31 - t40 * t29) * r_i_i_C(2), 0, 0;];
Ja_transl  = t1;
