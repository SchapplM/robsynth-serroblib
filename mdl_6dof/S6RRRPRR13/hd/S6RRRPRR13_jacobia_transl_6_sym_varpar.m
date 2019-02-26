% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR13_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR13_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:03
% EndTime: 2019-02-26 22:23:04
% DurationCPUTime: 0.50s
% Computational Cost: add. (679->108), mult. (1554->181), div. (0->0), fcn. (2029->16), ass. (0->73)
t56 = sin(qJ(2));
t58 = cos(qJ(2));
t59 = cos(qJ(1));
t77 = cos(pkin(6));
t70 = t59 * t77;
t85 = sin(qJ(1));
t37 = t85 * t56 - t58 * t70;
t38 = t56 * t70 + t85 * t58;
t55 = sin(qJ(3));
t52 = cos(pkin(7));
t86 = cos(qJ(3));
t73 = t52 * t86;
t50 = sin(pkin(7));
t51 = sin(pkin(6));
t81 = t51 * t59;
t75 = t50 * t81;
t15 = t37 * t73 + t38 * t55 + t86 * t75;
t80 = t52 * t55;
t16 = -t37 * t80 + t38 * t86 - t55 * t75;
t30 = -t37 * t50 + t52 * t81;
t48 = pkin(13) + qJ(5);
t46 = sin(t48);
t47 = cos(t48);
t4 = t16 * t47 - t30 * t46;
t54 = sin(qJ(6));
t57 = cos(qJ(6));
t94 = -t15 * t57 + t4 * t54;
t93 = -t15 * t54 - t4 * t57;
t90 = pkin(10) + pkin(4) * sin(pkin(13));
t45 = cos(pkin(13)) * pkin(4) + pkin(3);
t88 = r_i_i_C(3) + pkin(12);
t89 = t88 * t46 + t45;
t84 = t46 * t50;
t83 = t47 * t50;
t82 = t50 * t51;
t79 = t55 * t56;
t78 = t55 * t58;
t76 = t56 * t82;
t74 = t51 * t85;
t72 = t86 * t56;
t71 = t86 * t58;
t69 = t77 * t50;
t68 = t50 * t74;
t67 = t90 * t50;
t66 = t77 * t85;
t65 = t57 * r_i_i_C(1) - t54 * r_i_i_C(2) + pkin(5);
t53 = -pkin(11) - qJ(4);
t64 = t54 * r_i_i_C(1) + t57 * r_i_i_C(2) - t53;
t63 = t59 * t56 + t58 * t66;
t62 = t63 * t86;
t61 = -t65 * t47 - t89;
t60 = t63 * t50 + t52 * t74;
t39 = -t56 * t66 + t59 * t58;
t36 = t77 * t52 - t58 * t82;
t35 = (-t52 * t79 + t71) * t51;
t34 = (t52 * t72 + t78) * t51;
t29 = t55 * t69 + (t52 * t78 + t72) * t51;
t28 = -t86 * t69 + (-t52 * t71 + t79) * t51;
t26 = -t39 * t80 - t62;
t25 = t39 * t73 - t63 * t55;
t24 = -t37 * t86 - t38 * t80;
t23 = -t37 * t55 + t38 * t73;
t22 = t35 * t47 + t46 * t76;
t20 = t39 * t86 + (-t63 * t52 + t68) * t55;
t19 = t39 * t55 + t52 * t62 - t86 * t68;
t14 = t29 * t47 + t36 * t46;
t12 = t26 * t47 + t39 * t84;
t10 = t24 * t47 + t38 * t84;
t8 = t20 * t47 + t60 * t46;
t7 = t20 * t46 - t60 * t47;
t2 = t19 * t54 + t8 * t57;
t1 = t19 * t57 - t8 * t54;
t3 = [t93 * r_i_i_C(1) + t94 * r_i_i_C(2) - t4 * pkin(5) + t15 * t53 - t38 * pkin(2) - t85 * pkin(1) + pkin(9) * t81 - t89 * t16 + (-t88 * t47 + t90) * t30 (t12 * t57 + t25 * t54) * r_i_i_C(1) + (-t12 * t54 + t25 * t57) * r_i_i_C(2) + t12 * pkin(5) + t26 * t45 - t25 * t53 - t63 * pkin(2) + t39 * t67 + t88 * (t26 * t46 - t39 * t83) t61 * t19 + t64 * t20, t19, -t65 * t7 + t88 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t59 * pkin(1) + t39 * pkin(2) + t8 * pkin(5) + pkin(9) * t74 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t19 * t53 + t20 * t45 + t60 * t90 + t88 * t7 (t10 * t57 + t23 * t54) * r_i_i_C(1) + (-t10 * t54 + t23 * t57) * r_i_i_C(2) + t10 * pkin(5) + t24 * t45 - t23 * t53 - t37 * pkin(2) + t88 * (t24 * t46 - t38 * t83) + t38 * t67, t61 * t15 + t64 * t16, t15, t88 * t4 + t65 * (-t16 * t46 - t30 * t47) -t94 * r_i_i_C(1) + t93 * r_i_i_C(2); 0 (t22 * t57 + t34 * t54) * r_i_i_C(1) + (-t22 * t54 + t34 * t57) * r_i_i_C(2) + t22 * pkin(5) + t35 * t45 - t34 * t53 + t88 * (t35 * t46 - t47 * t76) + (t58 * pkin(2) + t56 * t67) * t51, t61 * t28 + t64 * t29, t28, t88 * t14 + t65 * (-t29 * t46 + t36 * t47) (-t14 * t54 + t28 * t57) * r_i_i_C(1) + (-t14 * t57 - t28 * t54) * r_i_i_C(2);];
Ja_transl  = t3;
