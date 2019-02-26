% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_6_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:52:56
% EndTime: 2019-02-26 22:52:57
% DurationCPUTime: 1.31s
% Computational Cost: add. (1756->187), mult. (5055->316), div. (0->0), fcn. (6693->18), ass. (0->115)
t76 = sin(pkin(6));
t88 = cos(qJ(1));
t128 = t76 * t88;
t120 = cos(pkin(6));
t113 = t88 * t120;
t138 = sin(qJ(1));
t83 = sin(qJ(2));
t87 = cos(qJ(2));
t68 = -t87 * t113 + t138 * t83;
t75 = sin(pkin(7));
t78 = cos(pkin(7));
t64 = t78 * t128 - t68 * t75;
t74 = sin(pkin(8));
t137 = t64 * t74;
t139 = cos(qJ(4));
t118 = t75 * t128;
t69 = t83 * t113 + t138 * t87;
t82 = sin(qJ(3));
t86 = cos(qJ(3));
t143 = (t68 * t78 + t118) * t86 + t69 * t82;
t126 = t78 * t82;
t53 = t82 * t118 + t68 * t126 - t69 * t86;
t77 = cos(pkin(8));
t81 = sin(qJ(4));
t22 = (t143 * t77 + t137) * t81 + t53 * t139;
t40 = t143 * t74 - t64 * t77;
t80 = sin(qJ(5));
t85 = cos(qJ(5));
t6 = t22 * t85 - t40 * t80;
t79 = sin(qJ(6));
t151 = t6 * t79;
t84 = cos(qJ(6));
t150 = t6 * t84;
t149 = t22 * t80 + t40 * t85;
t146 = t53 * t81;
t142 = r_i_i_C(3) + pkin(14);
t110 = t120 * t138;
t70 = -t83 * t110 + t87 * t88;
t104 = t87 * t110 + t88 * t83;
t116 = t76 * t138;
t96 = -t104 * t78 + t75 * t116;
t92 = t70 * t82 - t96 * t86;
t97 = t104 * t75 + t78 * t116;
t89 = t92 * t74 + t97 * t77;
t141 = pkin(11) * t75;
t140 = pkin(12) * t74;
t122 = t83 * t86;
t123 = t82 * t87;
t66 = (-t78 * t122 - t123) * t76;
t135 = t66 * t74;
t133 = t74 * t75;
t132 = t74 * t80;
t131 = t74 * t85;
t130 = t75 * t77;
t129 = t75 * t83;
t127 = t77 * t81;
t125 = t78 * t86;
t124 = t82 * t83;
t121 = t86 * t87;
t119 = t76 * t129;
t117 = t74 * t139;
t115 = t77 * t139;
t114 = t75 * t120;
t112 = t75 * t117;
t109 = r_i_i_C(1) * t84 - r_i_i_C(2) * t79 + pkin(5);
t108 = r_i_i_C(1) * t79 + r_i_i_C(2) * t84 + pkin(13);
t55 = -t69 * t125 + t68 * t82;
t46 = t69 * t130 - t55 * t74;
t57 = t104 * t82 - t70 * t125;
t47 = t70 * t130 - t57 * t74;
t105 = -t76 * t87 * t75 + t120 * t78;
t103 = t105 * t74;
t100 = -t109 * t85 - t142 * t80 - pkin(4);
t99 = t86 * t114 + (t78 * t121 - t124) * t76;
t98 = t143 * t139;
t95 = t99 * t139;
t94 = t97 * t74;
t90 = t92 * t139;
t67 = (-t78 * t124 + t121) * t76;
t63 = t82 * t114 + (t78 * t123 + t122) * t76;
t59 = t77 * t119 - t135;
t58 = -t104 * t86 - t70 * t126;
t56 = -t69 * t126 - t68 * t86;
t54 = t70 * t86 + t96 * t82;
t50 = t105 * t77 - t99 * t74;
t45 = t67 * t139 + (t74 * t119 + t66 * t77) * t81;
t44 = -t76 * t83 * t112 - t66 * t115 + t67 * t81;
t43 = -t63 * t127 + t95;
t42 = t63 * t115 + t99 * t81;
t38 = t63 * t139 + (t99 * t77 + t103) * t81;
t37 = -t139 * t103 + t63 * t81 - t77 * t95;
t36 = t63 * t132 + t43 * t85;
t34 = t45 * t85 + t59 * t80;
t32 = -t54 * t127 - t90;
t31 = t54 * t115 - t92 * t81;
t30 = t127 * t53 - t98;
t29 = -t115 * t53 - t143 * t81;
t28 = t58 * t139 + (t70 * t133 + t57 * t77) * t81;
t27 = -t70 * t112 - t57 * t115 + t58 * t81;
t26 = t56 * t139 + (t69 * t133 + t55 * t77) * t81;
t25 = -t69 * t112 - t55 * t115 + t56 * t81;
t24 = t54 * t139 + (-t92 * t77 + t94) * t81;
t23 = -t139 * t94 + t54 * t81 + t77 * t90;
t21 = -t115 * t143 - t64 * t117 + t146;
t19 = t137 * t139 + t77 * t98 - t146;
t18 = t38 * t85 + t50 * t80;
t16 = t54 * t132 + t32 * t85;
t14 = -t132 * t53 + t30 * t85;
t12 = t28 * t85 + t47 * t80;
t10 = t26 * t85 + t46 * t80;
t8 = t24 * t85 + t89 * t80;
t7 = t24 * t80 - t89 * t85;
t2 = t23 * t79 + t8 * t84;
t1 = t23 * t84 - t79 * t8;
t3 = [(t21 * t79 + t150) * r_i_i_C(1) + (t21 * t84 - t151) * r_i_i_C(2) + t6 * pkin(5) + t22 * pkin(4) + t21 * pkin(13) + t53 * pkin(3) - t69 * pkin(2) - t138 * pkin(1) + pkin(10) * t128 + t142 * t149 - t40 * pkin(12) + t64 * pkin(11) (t12 * t84 + t27 * t79) * r_i_i_C(1) + (-t12 * t79 + t27 * t84) * r_i_i_C(2) + t12 * pkin(5) + t28 * pkin(4) + t27 * pkin(13) + t58 * pkin(3) - t104 * pkin(2) + t70 * t141 + t142 * (t28 * t80 - t47 * t85) + t47 * pkin(12) (t16 * t84 + t31 * t79) * r_i_i_C(1) + (-t16 * t79 + t31 * t84) * r_i_i_C(2) + t16 * pkin(5) + t32 * pkin(4) + t31 * pkin(13) - t92 * pkin(3) + t54 * t140 + t142 * (-t54 * t131 + t32 * t80) t100 * t23 + t108 * t24, -t109 * t7 + t142 * t8, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t88 * pkin(1) + t70 * pkin(2) + t54 * pkin(3) + t24 * pkin(4) + t8 * pkin(5) + pkin(10) * t116 + t97 * pkin(11) + t89 * pkin(12) + t23 * pkin(13) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t142 * t7 (t10 * t84 + t25 * t79) * r_i_i_C(1) + (-t10 * t79 + t25 * t84) * r_i_i_C(2) + t10 * pkin(5) + t26 * pkin(4) + t25 * pkin(13) + t56 * pkin(3) - t68 * pkin(2) + t69 * t141 + t142 * (t26 * t80 - t46 * t85) + t46 * pkin(12) (t14 * t84 + t29 * t79) * r_i_i_C(1) + (-t14 * t79 + t29 * t84) * r_i_i_C(2) + t14 * pkin(5) + t30 * pkin(4) + t29 * pkin(13) - t143 * pkin(3) - t53 * t140 + t142 * (t131 * t53 + t30 * t80) t100 * t19 - t108 * t22, t109 * t149 - t142 * t6 (t19 * t84 + t151) * r_i_i_C(1) + (-t19 * t79 + t150) * r_i_i_C(2); 0 (t34 * t84 + t44 * t79) * r_i_i_C(1) + (-t34 * t79 + t44 * t84) * r_i_i_C(2) + t34 * pkin(5) + t45 * pkin(4) + t44 * pkin(13) + t67 * pkin(3) - pkin(12) * t135 + t142 * (t45 * t80 - t59 * t85) + (t87 * pkin(2) + (pkin(12) * t77 + pkin(11)) * t129) * t76 (t36 * t84 + t42 * t79) * r_i_i_C(1) + (-t36 * t79 + t42 * t84) * r_i_i_C(2) + t36 * pkin(5) + t43 * pkin(4) + t42 * pkin(13) + t99 * pkin(3) + t63 * t140 + t142 * (-t63 * t131 + t43 * t80) t100 * t37 + t108 * t38, t142 * t18 + t109 * (-t38 * t80 + t50 * t85) (-t18 * t79 + t37 * t84) * r_i_i_C(1) + (-t18 * t84 - t37 * t79) * r_i_i_C(2);];
Ja_transl  = t3;
