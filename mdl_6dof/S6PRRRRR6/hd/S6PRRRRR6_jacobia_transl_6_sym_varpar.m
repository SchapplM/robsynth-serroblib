% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_jacobia_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:21:51
% EndTime: 2019-02-26 20:21:52
% DurationCPUTime: 0.63s
% Computational Cost: add. (1404->161), mult. (4083->287), div. (0->0), fcn. (5397->18), ass. (0->106)
t122 = r_i_i_C(3) + pkin(13);
t64 = sin(pkin(8));
t121 = pkin(11) * t64;
t65 = sin(pkin(7));
t120 = t65 * pkin(10);
t119 = cos(qJ(4));
t73 = sin(qJ(2));
t76 = cos(qJ(3));
t107 = t73 * t76;
t72 = sin(qJ(3));
t77 = cos(qJ(2));
t108 = t72 * t77;
t66 = sin(pkin(6));
t68 = cos(pkin(7));
t56 = (-t68 * t107 - t108) * t66;
t118 = t56 * t64;
t70 = sin(qJ(5));
t117 = t64 * t70;
t75 = cos(qJ(5));
t116 = t64 * t75;
t115 = t65 * t64;
t67 = cos(pkin(8));
t114 = t65 * t67;
t113 = t65 * t73;
t71 = sin(qJ(4));
t112 = t67 * t71;
t111 = t68 * t72;
t110 = t68 * t76;
t109 = t72 * t73;
t106 = t76 * t77;
t105 = cos(pkin(6));
t104 = cos(pkin(14));
t103 = sin(pkin(14));
t102 = t66 * t113;
t101 = t67 * t119;
t100 = t66 * t104;
t99 = t66 * t103;
t98 = t105 * t65;
t97 = t119 * t115;
t96 = t65 * t100;
t95 = t105 * t104;
t94 = t105 * t103;
t69 = sin(qJ(6));
t74 = cos(qJ(6));
t93 = t74 * r_i_i_C(1) - t69 * r_i_i_C(2) + pkin(5);
t92 = t69 * r_i_i_C(1) + t74 * r_i_i_C(2) + pkin(12);
t58 = -t103 * t73 + t77 * t95;
t59 = t103 * t77 + t73 * t95;
t46 = -t59 * t110 - t58 * t72;
t37 = t59 * t114 - t46 * t64;
t60 = -t104 * t73 - t77 * t94;
t61 = t104 * t77 - t73 * t94;
t48 = -t61 * t110 - t60 * t72;
t38 = t61 * t114 - t48 * t64;
t91 = -t68 * t100 - t58 * t65;
t90 = -t60 * t65 + t68 * t99;
t89 = -t66 * t77 * t65 + t105 * t68;
t88 = t60 * t68 + t65 * t99;
t87 = t91 * t64;
t86 = t90 * t64;
t85 = t89 * t64;
t84 = -t122 * t70 - t93 * t75 - pkin(4);
t83 = t59 * t72 - (t58 * t68 - t96) * t76;
t82 = t61 * t72 - t88 * t76;
t81 = t76 * t98 + (t68 * t106 - t109) * t66;
t80 = t83 * t119;
t79 = t82 * t119;
t78 = t81 * t119;
t57 = (-t68 * t109 + t106) * t66;
t54 = t72 * t98 + (t68 * t108 + t107) * t66;
t50 = t67 * t102 - t118;
t49 = -t61 * t111 + t60 * t76;
t47 = -t59 * t111 + t58 * t76;
t45 = t61 * t76 + t88 * t72;
t44 = t58 * t111 + t59 * t76 - t72 * t96;
t43 = -t81 * t64 + t89 * t67;
t40 = t57 * t119 + (t64 * t102 + t56 * t67) * t71;
t39 = -t66 * t73 * t97 - t56 * t101 + t57 * t71;
t36 = -t54 * t112 + t78;
t35 = t54 * t101 + t81 * t71;
t34 = t82 * t64 + t90 * t67;
t33 = t83 * t64 + t91 * t67;
t32 = t54 * t119 + (t81 * t67 + t85) * t71;
t31 = -t119 * t85 + t54 * t71 - t67 * t78;
t30 = t54 * t117 + t36 * t75;
t28 = t40 * t75 + t50 * t70;
t26 = -t45 * t112 - t79;
t25 = t45 * t101 - t82 * t71;
t24 = -t44 * t112 - t80;
t23 = t44 * t101 - t83 * t71;
t22 = t49 * t119 + (t61 * t115 + t48 * t67) * t71;
t21 = -t48 * t101 + t49 * t71 - t61 * t97;
t20 = t47 * t119 + (t59 * t115 + t46 * t67) * t71;
t19 = -t46 * t101 + t47 * t71 - t59 * t97;
t18 = t45 * t119 + (-t82 * t67 + t86) * t71;
t17 = -t119 * t86 + t45 * t71 + t67 * t79;
t16 = t44 * t119 + (-t83 * t67 + t87) * t71;
t15 = -t119 * t87 + t44 * t71 + t67 * t80;
t14 = t32 * t75 + t43 * t70;
t12 = t45 * t117 + t26 * t75;
t10 = t44 * t117 + t24 * t75;
t8 = t22 * t75 + t38 * t70;
t6 = t20 * t75 + t37 * t70;
t4 = t18 * t75 + t34 * t70;
t2 = t16 * t75 + t33 * t70;
t1 = [0 (t21 * t69 + t8 * t74) * r_i_i_C(1) + (t21 * t74 - t8 * t69) * r_i_i_C(2) + t8 * pkin(5) + t22 * pkin(4) + t21 * pkin(12) + t49 * pkin(3) + t60 * pkin(2) + t61 * t120 + t122 * (t22 * t70 - t38 * t75) + t38 * pkin(11) (t12 * t74 + t25 * t69) * r_i_i_C(1) + (-t12 * t69 + t25 * t74) * r_i_i_C(2) + t12 * pkin(5) + t26 * pkin(4) + t25 * pkin(12) - t82 * pkin(3) + t45 * t121 + t122 * (-t45 * t116 + t26 * t70) t84 * t17 + t92 * t18, t122 * t4 + t93 * (-t18 * t70 + t34 * t75) (t17 * t74 - t4 * t69) * r_i_i_C(1) + (-t17 * t69 - t4 * t74) * r_i_i_C(2); 0 (t19 * t69 + t6 * t74) * r_i_i_C(1) + (t19 * t74 - t6 * t69) * r_i_i_C(2) + t6 * pkin(5) + t20 * pkin(4) + t19 * pkin(12) + t47 * pkin(3) + t58 * pkin(2) + t59 * t120 + t122 * (t20 * t70 - t37 * t75) + t37 * pkin(11) (t10 * t74 + t23 * t69) * r_i_i_C(1) + (-t10 * t69 + t23 * t74) * r_i_i_C(2) + t10 * pkin(5) + t24 * pkin(4) + t23 * pkin(12) - t83 * pkin(3) + t44 * t121 + t122 * (-t44 * t116 + t24 * t70) t84 * t15 + t92 * t16, t122 * t2 + t93 * (-t16 * t70 + t33 * t75) (t15 * t74 - t2 * t69) * r_i_i_C(1) + (-t15 * t69 - t2 * t74) * r_i_i_C(2); 1 (t28 * t74 + t39 * t69) * r_i_i_C(1) + (-t28 * t69 + t39 * t74) * r_i_i_C(2) + t28 * pkin(5) + t40 * pkin(4) + t39 * pkin(12) + t57 * pkin(3) - pkin(11) * t118 + t122 * (t40 * t70 - t50 * t75) + (t77 * pkin(2) + (pkin(11) * t67 + pkin(10)) * t113) * t66 (t30 * t74 + t35 * t69) * r_i_i_C(1) + (-t30 * t69 + t35 * t74) * r_i_i_C(2) + t30 * pkin(5) + t36 * pkin(4) + t35 * pkin(12) + t81 * pkin(3) + t54 * t121 + t122 * (-t54 * t116 + t36 * t70) t84 * t31 + t92 * t32, t122 * t14 + t93 * (-t32 * t70 + t43 * t75) (-t14 * t69 + t31 * t74) * r_i_i_C(1) + (-t14 * t74 - t31 * t69) * r_i_i_C(2);];
Ja_transl  = t1;
