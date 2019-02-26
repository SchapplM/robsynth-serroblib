% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:03
% EndTime: 2019-02-26 20:14:04
% DurationCPUTime: 0.23s
% Computational Cost: add. (356->75), mult. (1003->139), div. (0->0), fcn. (1308->14), ass. (0->59)
t35 = sin(pkin(7));
t71 = pkin(9) * t35;
t36 = sin(pkin(6));
t70 = t35 * t36;
t40 = cos(pkin(6));
t69 = t35 * t40;
t41 = sin(qJ(4));
t68 = t35 * t41;
t43 = sin(qJ(2));
t67 = t35 * t43;
t44 = cos(qJ(4));
t66 = t35 * t44;
t39 = cos(pkin(7));
t65 = t36 * t39;
t42 = sin(qJ(3));
t64 = t39 * t42;
t45 = cos(qJ(3));
t63 = t39 * t45;
t62 = t40 * t43;
t46 = cos(qJ(2));
t61 = t40 * t46;
t60 = t42 * t43;
t59 = t42 * t46;
t58 = t43 * t45;
t57 = t45 * t46;
t56 = r_i_i_C(3) + qJ(5);
t55 = t36 * t67;
t33 = sin(pkin(13));
t37 = cos(pkin(13));
t54 = r_i_i_C(1) * t37 - r_i_i_C(2) * t33 + pkin(4);
t53 = t33 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(10);
t34 = sin(pkin(12));
t38 = cos(pkin(12));
t28 = -t34 * t43 + t38 * t61;
t52 = -t28 * t35 - t38 * t65;
t51 = t28 * t39 - t38 * t70;
t30 = -t34 * t61 - t38 * t43;
t50 = -t30 * t35 + t34 * t65;
t49 = t30 * t39 + t34 * t70;
t48 = t40 * t39 - t46 * t70;
t47 = t56 * t41 + t54 * t44 + pkin(3);
t31 = -t34 * t62 + t38 * t46;
t29 = t34 * t46 + t38 * t62;
t26 = (-t39 * t60 + t57) * t36;
t25 = (t39 * t58 + t59) * t36;
t24 = t42 * t69 + (t39 * t59 + t58) * t36;
t20 = t26 * t44 + t41 * t55;
t18 = t30 * t45 - t31 * t64;
t17 = t30 * t42 + t31 * t63;
t16 = t28 * t45 - t29 * t64;
t15 = t28 * t42 + t29 * t63;
t13 = t24 * t41 - t48 * t44;
t12 = t31 * t45 + t49 * t42;
t10 = t29 * t45 + t51 * t42;
t8 = t18 * t44 + t31 * t68;
t6 = t16 * t44 + t29 * t68;
t3 = t12 * t41 - t50 * t44;
t1 = t10 * t41 - t52 * t44;
t2 = [0 (t17 * t33 + t8 * t37) * r_i_i_C(1) + (t17 * t37 - t8 * t33) * r_i_i_C(2) + t8 * pkin(4) + t18 * pkin(3) + t17 * pkin(10) + t30 * pkin(2) + t31 * t71 + t56 * (t18 * t41 - t31 * t66) t53 * t12 + t47 * (-t31 * t42 + t49 * t45) t56 * (t12 * t44 + t50 * t41) - t54 * t3, t3, 0; 0 (t15 * t33 + t6 * t37) * r_i_i_C(1) + (t15 * t37 - t6 * t33) * r_i_i_C(2) + t6 * pkin(4) + t16 * pkin(3) + t15 * pkin(10) + t28 * pkin(2) + t29 * t71 + t56 * (t16 * t41 - t29 * t66) t53 * t10 + t47 * (-t29 * t42 + t51 * t45) t56 * (t10 * t44 + t52 * t41) - t54 * t1, t1, 0; 1 (t20 * t37 + t25 * t33) * r_i_i_C(1) + (-t20 * t33 + t25 * t37) * r_i_i_C(2) + t20 * pkin(4) + t26 * pkin(3) + t25 * pkin(10) + (t46 * pkin(2) + pkin(9) * t67) * t36 + t56 * (t26 * t41 - t44 * t55) t53 * t24 + t47 * (t45 * t69 + (t39 * t57 - t60) * t36) t56 * (t24 * t44 + t48 * t41) - t54 * t13, t13, 0;];
Ja_transl  = t2;
