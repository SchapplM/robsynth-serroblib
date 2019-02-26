% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RRRRPR12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:58
% EndTime: 2019-02-26 22:36:59
% DurationCPUTime: 0.52s
% Computational Cost: add. (700->114), mult. (1612->190), div. (0->0), fcn. (2102->16), ass. (0->75)
t56 = sin(qJ(2));
t59 = cos(qJ(2));
t60 = cos(qJ(1));
t78 = cos(pkin(6));
t71 = t60 * t78;
t86 = sin(qJ(1));
t37 = t56 * t86 - t59 * t71;
t38 = t56 * t71 + t59 * t86;
t55 = sin(qJ(3));
t51 = cos(pkin(7));
t87 = cos(qJ(3));
t74 = t51 * t87;
t49 = sin(pkin(7));
t50 = sin(pkin(6));
t82 = t50 * t60;
t76 = t49 * t82;
t15 = t37 * t74 + t38 * t55 + t76 * t87;
t81 = t51 * t55;
t16 = -t37 * t81 + t38 * t87 - t55 * t76;
t30 = -t37 * t49 + t51 * t82;
t48 = qJ(4) + pkin(13);
t46 = sin(t48);
t47 = cos(t48);
t4 = t16 * t47 - t30 * t46;
t53 = sin(qJ(6));
t57 = cos(qJ(6));
t95 = -t15 * t57 + t4 * t53;
t94 = -t15 * t53 - t4 * t57;
t54 = sin(qJ(4));
t91 = pkin(4) * t54 + pkin(10);
t58 = cos(qJ(4));
t45 = pkin(4) * t58 + pkin(3);
t89 = r_i_i_C(3) + pkin(12);
t90 = t46 * t89 + t45;
t85 = t46 * t49;
t84 = t47 * t49;
t83 = t49 * t50;
t80 = t55 * t56;
t79 = t55 * t59;
t77 = t56 * t83;
t75 = t50 * t86;
t73 = t87 * t56;
t72 = t87 * t59;
t70 = t78 * t49;
t69 = t49 * t75;
t68 = t91 * t49;
t67 = t78 * t86;
t66 = r_i_i_C(1) * t57 - r_i_i_C(2) * t53 + pkin(5);
t52 = -qJ(5) - pkin(11);
t65 = r_i_i_C(1) * t53 + r_i_i_C(2) * t57 - t52;
t64 = t56 * t60 + t59 * t67;
t63 = t64 * t87;
t62 = -t47 * t66 - t90;
t61 = t49 * t64 + t51 * t75;
t39 = -t56 * t67 + t59 * t60;
t36 = t51 * t78 - t59 * t83;
t35 = (-t51 * t80 + t72) * t50;
t34 = (t51 * t73 + t79) * t50;
t29 = t55 * t70 + (t51 * t79 + t73) * t50;
t28 = -t87 * t70 + (-t51 * t72 + t80) * t50;
t26 = -t39 * t81 - t63;
t25 = t39 * t74 - t55 * t64;
t24 = -t37 * t87 - t38 * t81;
t23 = -t37 * t55 + t38 * t74;
t22 = t35 * t47 + t46 * t77;
t20 = t39 * t87 + (-t51 * t64 + t69) * t55;
t19 = t39 * t55 + t51 * t63 - t69 * t87;
t14 = t29 * t47 + t36 * t46;
t12 = t26 * t47 + t39 * t85;
t10 = t24 * t47 + t38 * t85;
t8 = t20 * t47 + t46 * t61;
t7 = t20 * t46 - t47 * t61;
t2 = t19 * t53 + t57 * t8;
t1 = t19 * t57 - t53 * t8;
t3 = [t94 * r_i_i_C(1) + t95 * r_i_i_C(2) - t4 * pkin(5) + t15 * t52 - t38 * pkin(2) - t86 * pkin(1) + pkin(9) * t82 - t90 * t16 + (-t47 * t89 + t91) * t30 (t12 * t57 + t25 * t53) * r_i_i_C(1) + (-t12 * t53 + t25 * t57) * r_i_i_C(2) + t12 * pkin(5) + t26 * t45 - t25 * t52 - t64 * pkin(2) + t39 * t68 + t89 * (t26 * t46 - t39 * t84) t19 * t62 + t20 * t65, t89 * t8 + (-t20 * t54 + t58 * t61) * pkin(4) - t66 * t7, t19, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t60 * pkin(1) + t39 * pkin(2) + t8 * pkin(5) + pkin(9) * t75 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t19 * t52 + t20 * t45 + t61 * t91 + t89 * t7 (t10 * t57 + t23 * t53) * r_i_i_C(1) + (-t10 * t53 + t23 * t57) * r_i_i_C(2) + t10 * pkin(5) + t24 * t45 - t23 * t52 - t37 * pkin(2) + t89 * (t24 * t46 - t38 * t84) + t38 * t68, t15 * t62 + t16 * t65, t89 * t4 + (-t16 * t54 - t30 * t58) * pkin(4) + t66 * (-t16 * t46 - t30 * t47) t15, -r_i_i_C(1) * t95 + r_i_i_C(2) * t94; 0 (t22 * t57 + t34 * t53) * r_i_i_C(1) + (-t22 * t53 + t34 * t57) * r_i_i_C(2) + t22 * pkin(5) + t35 * t45 - t34 * t52 + t89 * (t35 * t46 - t47 * t77) + (t59 * pkin(2) + t56 * t68) * t50, t28 * t62 + t29 * t65, t89 * t14 + (-t29 * t54 + t36 * t58) * pkin(4) + t66 * (-t29 * t46 + t36 * t47) t28 (-t14 * t53 + t28 * t57) * r_i_i_C(1) + (-t14 * t57 - t28 * t53) * r_i_i_C(2);];
Ja_transl  = t3;
