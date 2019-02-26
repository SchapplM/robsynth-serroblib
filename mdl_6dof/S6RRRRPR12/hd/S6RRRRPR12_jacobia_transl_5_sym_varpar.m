% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR12_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR12_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:58
% EndTime: 2019-02-26 22:36:58
% DurationCPUTime: 0.26s
% Computational Cost: add. (338->74), mult. (792->125), div. (0->0), fcn. (1017->14), ass. (0->51)
t41 = sin(qJ(4));
t69 = t41 * pkin(4) + pkin(10);
t45 = cos(qJ(4));
t33 = t45 * pkin(4) + pkin(3);
t36 = qJ(4) + pkin(13);
t34 = sin(t36);
t35 = cos(t36);
t50 = t35 * r_i_i_C(1) - t34 * r_i_i_C(2) + t33;
t68 = t34 * r_i_i_C(1) + t35 * r_i_i_C(2) + t69;
t66 = r_i_i_C(3) + qJ(5) + pkin(11);
t43 = sin(qJ(2));
t44 = sin(qJ(1));
t47 = cos(qJ(2));
t48 = cos(qJ(1));
t56 = cos(pkin(6));
t52 = t48 * t56;
t24 = t43 * t52 + t44 * t47;
t42 = sin(qJ(3));
t65 = t24 * t42;
t38 = sin(pkin(6));
t64 = t38 * t44;
t63 = t38 * t48;
t39 = cos(pkin(7));
t62 = t39 * t42;
t46 = cos(qJ(3));
t61 = t39 * t46;
t60 = t42 * t43;
t59 = t42 * t47;
t58 = t43 * t46;
t57 = t46 * t47;
t37 = sin(pkin(7));
t55 = t37 * t64;
t54 = t37 * t63;
t53 = t44 * t56;
t51 = t56 * t37;
t23 = t44 * t43 - t47 * t52;
t15 = -t23 * t37 + t39 * t63;
t25 = -t48 * t43 - t47 * t53;
t17 = -t25 * t37 + t39 * t64;
t6 = t23 * t62 - t24 * t46 + t42 * t54;
t49 = t68 * t37;
t26 = -t43 * t53 + t48 * t47;
t22 = -t38 * t47 * t37 + t56 * t39;
t14 = t42 * t51 + (t39 * t59 + t58) * t38;
t13 = -t46 * t51 + (-t39 * t57 + t60) * t38;
t8 = t26 * t46 + (t25 * t39 + t55) * t42;
t7 = -t25 * t61 + t26 * t42 - t46 * t55;
t3 = t23 * t61 + t46 * t54 + t65;
t2 = t17 * t34 + t8 * t35;
t1 = t17 * t35 - t8 * t34;
t4 = [-t24 * pkin(2) - t44 * pkin(1) + pkin(9) * t63 + t66 * (-t65 + (-t23 * t39 - t54) * t46) + t50 * t6 + t68 * t15, t25 * pkin(2) + t50 * (t25 * t46 - t26 * t62) + t26 * t49 + t66 * (t25 * t42 + t26 * t61) -t50 * t7 + t66 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (t17 * t45 - t41 * t8) * pkin(4), t7, 0; t48 * pkin(1) + t26 * pkin(2) + pkin(9) * t64 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t69 * t17 + t8 * t33 + t66 * t7, -t23 * pkin(2) + t66 * (-t23 * t42 + t24 * t61) + t50 * (-t23 * t46 - t24 * t62) + t24 * t49, -t50 * t3 - t6 * t66 (-t15 * t35 + t6 * t34) * r_i_i_C(1) + (t15 * t34 + t6 * t35) * r_i_i_C(2) + (-t15 * t45 + t41 * t6) * pkin(4), t3, 0; 0 (t66 * (t39 * t58 + t59) + t50 * (-t39 * t60 + t57) + t47 * pkin(2) + t43 * t49) * t38, -t50 * t13 + t66 * t14 (-t14 * t34 + t22 * t35) * r_i_i_C(1) + (-t14 * t35 - t22 * t34) * r_i_i_C(2) + (-t14 * t41 + t22 * t45) * pkin(4), t13, 0;];
Ja_transl  = t4;
