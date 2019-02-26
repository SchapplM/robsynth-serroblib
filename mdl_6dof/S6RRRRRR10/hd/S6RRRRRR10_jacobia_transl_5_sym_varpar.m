% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR10_jacobia_transl_5_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_transl_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10_jacobia_transl_5_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_transl_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:17
% EndTime: 2018-11-23 11:27:18
% DurationCPUTime: 0.81s
% Computational Cost: add. (3186->160), mult. (3192->244), div. (0->0), fcn. (3166->28), ass. (0->99)
t100 = cos(qJ(1));
t87 = sin(pkin(6));
t121 = t100 * t87;
t119 = pkin(6) + qJ(2);
t108 = cos(t119) / 0.2e1;
t120 = pkin(6) - qJ(2);
t114 = cos(t120);
t101 = t114 / 0.2e1 + t108;
t94 = sin(qJ(2));
t95 = sin(qJ(1));
t56 = -t100 * t101 + t94 * t95;
t105 = sin(t119) / 0.2e1;
t111 = sin(t120);
t68 = t105 - t111 / 0.2e1;
t99 = cos(qJ(2));
t57 = t100 * t68 + t95 * t99;
t117 = pkin(7) + qJ(3);
t104 = sin(t117) / 0.2e1;
t118 = pkin(7) - qJ(3);
t110 = sin(t118);
t65 = t104 + t110 / 0.2e1;
t107 = cos(t117) / 0.2e1;
t113 = cos(t118);
t72 = t113 / 0.2e1 + t107;
t93 = sin(qJ(3));
t36 = t121 * t65 + t56 * t72 + t57 * t93;
t86 = sin(pkin(7));
t89 = cos(pkin(7));
t52 = t121 * t89 - t56 * t86;
t85 = sin(pkin(8));
t88 = cos(pkin(8));
t26 = t36 * t85 - t52 * t88;
t66 = t104 - t110 / 0.2e1;
t71 = t107 - t113 / 0.2e1;
t98 = cos(qJ(3));
t34 = t121 * t71 - t56 * t66 + t57 * t98;
t115 = pkin(8) + qJ(4);
t103 = sin(t115) / 0.2e1;
t116 = pkin(8) - qJ(4);
t109 = sin(t116);
t64 = t103 - t109 / 0.2e1;
t106 = cos(t115) / 0.2e1;
t112 = cos(t116);
t69 = t106 - t112 / 0.2e1;
t97 = cos(qJ(4));
t4 = t34 * t97 - t36 * t64 + t52 * t69;
t91 = sin(qJ(5));
t96 = cos(qJ(5));
t134 = t26 * t96 - t4 * t91;
t133 = -t26 * t91 - t4 * t96;
t63 = t103 + t109 / 0.2e1;
t70 = t112 / 0.2e1 + t106;
t92 = sin(qJ(4));
t132 = -t34 * t92 - t36 * t70 - t52 * t63;
t131 = r_i_i_C(3) + pkin(13);
t130 = t85 * pkin(12);
t129 = t86 * pkin(11);
t128 = t63 * t86;
t127 = t69 * t86;
t73 = t108 - t114 / 0.2e1;
t126 = t73 * t86;
t125 = t85 * t91;
t124 = t85 * t96;
t123 = t86 * t88;
t122 = t87 * t95;
t59 = -t100 * t94 - t101 * t95;
t60 = t100 * t99 - t68 * t95;
t38 = t122 * t65 + t59 * t72 - t60 * t93;
t54 = t122 * t89 - t59 * t86;
t28 = -t38 * t85 + t54 * t88;
t102 = r_i_i_C(1) * t96 - r_i_i_C(2) * t91 + pkin(4);
t42 = t56 * t93 - t57 * t72;
t29 = t123 * t57 - t42 * t85;
t44 = -t59 * t93 - t60 * t72;
t30 = t123 * t60 - t44 * t85;
t67 = t105 + t111 / 0.2e1;
t50 = -t67 * t93 + t72 * t73;
t41 = -t123 * t73 - t50 * t85;
t40 = t122 * t71 - t59 * t66 - t60 * t98;
t9 = t38 * t64 - t40 * t97 - t54 * t69;
t90 = cos(pkin(6));
t46 = t65 * t90 + t67 * t72 + t73 * t93;
t47 = t66 * t67 - t71 * t90 - t73 * t98;
t55 = -t67 * t86 + t89 * t90;
t20 = t46 * t64 + t47 * t97 - t55 * t69;
t51 = t66 * t73 + t67 * t98;
t45 = t59 * t98 - t60 * t66;
t43 = -t56 * t98 - t57 * t66;
t32 = -t46 * t85 + t55 * t88;
t25 = t126 * t69 + t50 * t64 + t51 * t97;
t23 = t46 * t97 - t47 * t64;
t18 = t38 * t97 + t40 * t64;
t16 = -t34 * t64 - t36 * t97;
t14 = -t127 * t60 + t44 * t64 + t45 * t97;
t12 = -t127 * t57 + t42 * t64 + t43 * t97;
t8 = -t38 * t70 - t40 * t92 - t54 * t63;
t2 = t28 * t91 + t9 * t96;
t1 = t28 * t96 - t9 * t91;
t3 = [-t95 * pkin(1) - t57 * pkin(2) - t34 * pkin(3) - pkin(4) * t4 + pkin(10) * t121 + t52 * pkin(11) - pkin(12) * t26 + r_i_i_C(1) * t133 - r_i_i_C(2) * t134 + t131 * t132 (t14 * t96 + t30 * t91) * r_i_i_C(1) + (-t14 * t91 + t30 * t96) * r_i_i_C(2) + t14 * pkin(4) + t45 * pkin(3) + t59 * pkin(2) + t60 * t129 + t131 * (-t128 * t60 - t44 * t70 + t45 * t92) + t30 * pkin(12) (-t125 * t40 + t18 * t96) * r_i_i_C(1) + (-t124 * t40 - t18 * t91) * r_i_i_C(2) + t18 * pkin(4) + t38 * pkin(3) - t40 * t130 + t131 * (t38 * t92 - t40 * t70) -t102 * t8 + t131 * t9, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t100 * pkin(1) + t60 * pkin(2) - pkin(3) * t40 + t9 * pkin(4) + pkin(10) * t122 + pkin(11) * t54 + pkin(12) * t28 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t131 * t8 (t12 * t96 + t29 * t91) * r_i_i_C(1) + (-t12 * t91 + t29 * t96) * r_i_i_C(2) + t12 * pkin(4) + t43 * pkin(3) - t56 * pkin(2) + t57 * t129 + t131 * (-t128 * t57 - t42 * t70 + t43 * t92) + t29 * pkin(12) (t125 * t34 + t16 * t96) * r_i_i_C(1) + (t124 * t34 - t16 * t91) * r_i_i_C(2) + t16 * pkin(4) - t36 * pkin(3) + t34 * t130 + t131 * (t34 * t70 - t36 * t92) t102 * t132 + t131 * t4, r_i_i_C(1) * t134 + r_i_i_C(2) * t133, 0; 0 (t25 * t96 + t41 * t91) * r_i_i_C(1) + (-t25 * t91 + t41 * t96) * r_i_i_C(2) + t25 * pkin(4) + t51 * pkin(3) + t67 * pkin(2) - pkin(11) * t126 + t131 * (t126 * t63 - t50 * t70 + t51 * t92) + t41 * pkin(12) (t125 * t47 + t23 * t96) * r_i_i_C(1) + (t124 * t47 - t23 * t91) * r_i_i_C(2) + t23 * pkin(4) + t46 * pkin(3) + t47 * t130 + t131 * (t46 * t92 + t47 * t70) t131 * t20 + t102 * (t46 * t70 - t47 * t92 + t55 * t63) (-t20 * t91 + t32 * t96) * r_i_i_C(1) + (-t20 * t96 - t32 * t91) * r_i_i_C(2), 0;];
Ja_transl  = t3;
