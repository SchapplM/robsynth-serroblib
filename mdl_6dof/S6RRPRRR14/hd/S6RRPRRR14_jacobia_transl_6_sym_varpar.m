% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:55:44
% EndTime: 2019-02-26 22:55:45
% DurationCPUTime: 0.96s
% Computational Cost: add. (1439->146), mult. (4127->251), div. (0->0), fcn. (5475->18), ass. (0->98)
t118 = cos(qJ(4));
t116 = sin(qJ(2));
t117 = sin(qJ(1));
t71 = cos(qJ(1));
t107 = cos(pkin(6));
t70 = cos(qJ(2));
t97 = t70 * t107;
t56 = t117 * t116 - t71 * t97;
t92 = t107 * t116;
t57 = t117 * t70 + t71 * t92;
t60 = sin(pkin(14));
t105 = cos(pkin(14));
t106 = cos(pkin(7));
t90 = t106 * t105;
t62 = sin(pkin(7));
t63 = sin(pkin(6));
t98 = t63 * t105;
t93 = t62 * t98;
t40 = t56 * t90 + t57 * t60 + t71 * t93;
t61 = sin(pkin(8));
t64 = cos(pkin(8));
t99 = t63 * t106;
t87 = t56 * t62 - t71 * t99;
t120 = t40 * t64 - t87 * t61;
t109 = t63 * t71;
t41 = (t106 * t56 + t62 * t109) * t60 - t57 * t105;
t67 = sin(qJ(4));
t18 = t41 * t118 + t120 * t67;
t30 = t40 * t61 + t87 * t64;
t66 = sin(qJ(5));
t69 = cos(qJ(5));
t6 = t18 * t69 - t30 * t66;
t65 = sin(qJ(6));
t129 = t6 * t65;
t68 = cos(qJ(6));
t128 = t6 * t68;
t127 = t18 * t66 + t30 * t69;
t124 = t41 * t67;
t58 = -t117 * t92 + t71 * t70;
t83 = t71 * t116 + t117 * t97;
t80 = t83 * t105;
t74 = t106 * t80 - t117 * t93 + t58 * t60;
t78 = t117 * t99 + t83 * t62;
t122 = -t78 * t61 + t74 * t64;
t96 = t107 * t62;
t77 = t105 * t96 + (-t116 * t60 + t70 * t90) * t63;
t111 = t62 * t63;
t84 = t107 * t106 - t70 * t111;
t121 = t84 * t61 + t77 * t64;
t119 = r_i_i_C(3) + pkin(13);
t54 = (-t116 * t90 - t60 * t70) * t63;
t113 = t54 * t61;
t112 = t62 * t61;
t110 = t62 * t64;
t108 = t62 * qJ(3);
t104 = t61 * t118;
t103 = t63 * t117;
t102 = t64 * t118;
t100 = t60 * t106;
t95 = t62 * t104;
t94 = t116 * t111;
t91 = t61 * t94;
t89 = t68 * r_i_i_C(1) - t65 * r_i_i_C(2) + pkin(5);
t88 = t65 * r_i_i_C(1) + t68 * r_i_i_C(2) + pkin(12);
t43 = t56 * t60 - t57 * t90;
t34 = t57 * t110 - t43 * t61;
t45 = -t58 * t90 + t83 * t60;
t35 = t58 * t110 - t45 * t61;
t81 = -t119 * t66 - t89 * t69 - pkin(4);
t72 = t74 * t61 + t78 * t64;
t55 = (-t116 * t100 + t105 * t70) * t63;
t51 = t116 * t98 + (t70 * t99 + t96) * t60;
t47 = t64 * t94 - t113;
t46 = -t58 * t100 - t80;
t44 = -t57 * t100 - t56 * t105;
t42 = t58 * t105 + (t62 * t103 - t83 * t106) * t60;
t38 = -t77 * t61 + t84 * t64;
t33 = t55 * t118 + (t54 * t64 + t91) * t67;
t32 = -t54 * t102 - t118 * t91 + t55 * t67;
t28 = t51 * t118 + t121 * t67;
t27 = -t121 * t118 + t51 * t67;
t26 = t33 * t69 + t47 * t66;
t24 = t46 * t118 + (t58 * t112 + t45 * t64) * t67;
t23 = -t45 * t102 + t46 * t67 - t58 * t95;
t22 = t44 * t118 + (t57 * t112 + t43 * t64) * t67;
t21 = -t43 * t102 + t44 * t67 - t57 * t95;
t20 = t42 * t118 - t122 * t67;
t19 = t122 * t118 + t42 * t67;
t17 = -t40 * t102 + t104 * t87 + t124;
t15 = t120 * t118 - t124;
t14 = t28 * t69 + t38 * t66;
t12 = t24 * t69 + t35 * t66;
t10 = t22 * t69 + t34 * t66;
t8 = t20 * t69 + t72 * t66;
t7 = t20 * t66 - t72 * t69;
t2 = t19 * t65 + t8 * t68;
t1 = t19 * t68 - t8 * t65;
t3 = [(t17 * t65 + t128) * r_i_i_C(1) + (t17 * t68 - t129) * r_i_i_C(2) + t6 * pkin(5) + t18 * pkin(4) + t17 * pkin(12) + t41 * pkin(3) - t57 * pkin(2) - t117 * pkin(1) + pkin(10) * t109 + t119 * t127 - t87 * qJ(3) - t30 * pkin(11) (t12 * t68 + t23 * t65) * r_i_i_C(1) + (-t12 * t65 + t23 * t68) * r_i_i_C(2) + t12 * pkin(5) + t24 * pkin(4) + t23 * pkin(12) + t46 * pkin(3) - t83 * pkin(2) + t58 * t108 + t119 * (t24 * t66 - t35 * t69) + t35 * pkin(11), t78, t81 * t19 + t88 * t20, t119 * t8 - t89 * t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t71 * pkin(1) + t58 * pkin(2) + t42 * pkin(3) + t20 * pkin(4) + t8 * pkin(5) + pkin(10) * t103 + t72 * pkin(11) + t19 * pkin(12) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t78 * qJ(3) + t119 * t7 (t10 * t68 + t21 * t65) * r_i_i_C(1) + (-t10 * t65 + t21 * t68) * r_i_i_C(2) + t10 * pkin(5) + t22 * pkin(4) + t21 * pkin(12) + t44 * pkin(3) - t56 * pkin(2) + t57 * t108 + t119 * (t22 * t66 - t34 * t69) + t34 * pkin(11), t87, t81 * t15 - t18 * t88, -t119 * t6 + t89 * t127 (t15 * t68 + t129) * r_i_i_C(1) + (-t15 * t65 + t128) * r_i_i_C(2); 0 (t26 * t68 + t32 * t65) * r_i_i_C(1) + (-t26 * t65 + t32 * t68) * r_i_i_C(2) + t26 * pkin(5) + t33 * pkin(4) + t32 * pkin(12) + t55 * pkin(3) - pkin(11) * t113 + t119 * (t33 * t66 - t47 * t69) + (t70 * pkin(2) + (pkin(11) * t64 + qJ(3)) * t62 * t116) * t63, t84, t81 * t27 + t88 * t28, t119 * t14 + t89 * (-t28 * t66 + t38 * t69) (-t14 * t65 + t27 * t68) * r_i_i_C(1) + (-t14 * t68 - t27 * t65) * r_i_i_C(2);];
Ja_transl  = t3;
