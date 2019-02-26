% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6PRPRRR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobia_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:29
% EndTime: 2019-02-26 19:57:29
% DurationCPUTime: 0.31s
% Computational Cost: add. (466->90), mult. (1356->175), div. (0->0), fcn. (1783->16), ass. (0->73)
t79 = r_i_i_C(3) + pkin(11);
t39 = sin(pkin(14));
t42 = sin(pkin(6));
t43 = cos(pkin(14));
t52 = cos(qJ(2));
t46 = cos(pkin(7));
t49 = sin(qJ(2));
t71 = t46 * t49;
t31 = (-t39 * t52 - t43 * t71) * t42;
t40 = sin(pkin(8));
t78 = t31 * t40;
t77 = t39 * t46;
t41 = sin(pkin(7));
t76 = t40 * t41;
t45 = cos(pkin(8));
t75 = t41 * t45;
t44 = cos(pkin(13));
t74 = t42 * t44;
t73 = t42 * t49;
t72 = t43 * t46;
t70 = t46 * t52;
t69 = t41 * qJ(3);
t68 = cos(pkin(6));
t67 = sin(pkin(13));
t66 = t41 * t73;
t65 = t42 * t67;
t64 = t44 * t68;
t63 = t68 * t41;
t35 = t49 * t64 + t67 * t52;
t34 = -t67 * t49 + t52 * t64;
t55 = t34 * t46 - t41 * t74;
t18 = -t35 * t39 + t55 * t43;
t29 = -t34 * t41 - t46 * t74;
t62 = t18 * t45 + t29 * t40;
t59 = t68 * t67;
t37 = t44 * t52 - t49 * t59;
t36 = -t44 * t49 - t52 * t59;
t53 = t36 * t46 + t41 * t65;
t20 = -t37 * t39 + t53 * t43;
t30 = -t36 * t41 + t46 * t65;
t61 = t20 * t45 + t30 * t40;
t27 = t43 * t63 + (-t39 * t49 + t43 * t70) * t42;
t33 = -t42 * t52 * t41 + t68 * t46;
t60 = t27 * t45 + t33 * t40;
t47 = sin(qJ(5));
t50 = cos(qJ(5));
t58 = t50 * r_i_i_C(1) - t47 * r_i_i_C(2) + pkin(4);
t22 = -t34 * t39 - t35 * t72;
t13 = -t22 * t40 + t35 * t75;
t57 = t22 * t45 + t35 * t76;
t24 = -t36 * t39 - t37 * t72;
t14 = -t24 * t40 + t37 * t75;
t56 = t24 * t45 + t37 * t76;
t54 = t31 * t45 + t40 * t66;
t51 = cos(qJ(4));
t48 = sin(qJ(4));
t32 = (-t39 * t71 + t43 * t52) * t42;
t28 = t43 * t73 + (t42 * t70 + t63) * t39;
t26 = t45 * t66 - t78;
t25 = t36 * t43 - t37 * t77;
t23 = t34 * t43 - t35 * t77;
t21 = t37 * t43 + t53 * t39;
t19 = t35 * t43 + t55 * t39;
t17 = -t27 * t40 + t33 * t45;
t16 = t32 * t51 + t54 * t48;
t12 = -t20 * t40 + t30 * t45;
t11 = -t18 * t40 + t29 * t45;
t10 = t28 * t51 + t60 * t48;
t8 = t25 * t51 + t56 * t48;
t6 = t23 * t51 + t57 * t48;
t4 = t21 * t51 + t61 * t48;
t2 = t19 * t51 + t62 * t48;
t1 = [0 (t14 * t47 + t8 * t50) * r_i_i_C(1) + (t14 * t50 - t8 * t47) * r_i_i_C(2) + t8 * pkin(4) + t25 * pkin(3) + t36 * pkin(2) + t37 * t69 + t79 * (t25 * t48 - t56 * t51) + t14 * pkin(10), t30, t79 * t4 + t58 * (-t21 * t48 + t61 * t51) (t12 * t50 - t4 * t47) * r_i_i_C(1) + (-t12 * t47 - t4 * t50) * r_i_i_C(2), 0; 0 (t13 * t47 + t6 * t50) * r_i_i_C(1) + (t13 * t50 - t6 * t47) * r_i_i_C(2) + t6 * pkin(4) + t23 * pkin(3) + t34 * pkin(2) + t35 * t69 + t79 * (t23 * t48 - t57 * t51) + t13 * pkin(10), t29, t79 * t2 + t58 * (-t19 * t48 + t62 * t51) (t11 * t50 - t2 * t47) * r_i_i_C(1) + (-t11 * t47 - t2 * t50) * r_i_i_C(2), 0; 1 (t16 * t50 + t26 * t47) * r_i_i_C(1) + (-t16 * t47 + t26 * t50) * r_i_i_C(2) + t16 * pkin(4) + t32 * pkin(3) - pkin(10) * t78 + t79 * (t32 * t48 - t54 * t51) + (t52 * pkin(2) + (pkin(10) * t45 + qJ(3)) * t49 * t41) * t42, t33, t79 * t10 + t58 * (-t28 * t48 + t60 * t51) (-t10 * t47 + t17 * t50) * r_i_i_C(1) + (-t10 * t50 - t17 * t47) * r_i_i_C(2), 0;];
Ja_transl  = t1;
