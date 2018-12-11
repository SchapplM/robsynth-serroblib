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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR14_jacobia_transl_6_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobia_transl_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobia_transl_6_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobia_transl_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:20
% EndTime: 2018-12-10 18:38:22
% DurationCPUTime: 1.08s
% Computational Cost: add. (5345->164), mult. (5429->244), div. (0->0), fcn. (5475->30), ass. (0->110)
t122 = pkin(8) + qJ(4);
t109 = sin(t122) / 0.2e1;
t123 = pkin(8) - qJ(4);
t115 = sin(t123);
t100 = t109 + t115 / 0.2e1;
t111 = cos(t122) / 0.2e1;
t117 = cos(t123);
t101 = t117 / 0.2e1 + t111;
t127 = cos(pkin(7));
t77 = sin(pkin(6));
t119 = t77 * t127;
t124 = pkin(6) + qJ(2);
t112 = cos(t124) / 0.2e1;
t125 = pkin(6) - qJ(2);
t118 = cos(t125);
t102 = t118 / 0.2e1 + t112;
t136 = sin(qJ(2));
t83 = sin(qJ(1));
t88 = cos(qJ(1));
t59 = -t88 * t102 + t83 * t136;
t76 = sin(pkin(7));
t104 = -t88 * t119 + t59 * t76;
t126 = sin(pkin(14));
t110 = sin(t124) / 0.2e1;
t116 = sin(t125);
t68 = t110 - t116 / 0.2e1;
t87 = cos(qJ(2));
t60 = t68 * t88 + t83 * t87;
t120 = pkin(7) + pkin(14);
t107 = sin(t120) / 0.2e1;
t121 = pkin(7) - pkin(14);
t113 = sin(t121);
t98 = t107 + t113 / 0.2e1;
t95 = t77 * t98;
t108 = cos(t121) / 0.2e1;
t114 = cos(t120);
t99 = t108 + t114 / 0.2e1;
t43 = t60 * t126 + t59 * t99 + t88 * t95;
t130 = t77 * t88;
t64 = t107 - t113 / 0.2e1;
t65 = t108 - t114 / 0.2e1;
t78 = cos(pkin(14));
t44 = t65 * t130 + t59 * t64 - t60 * t78;
t82 = sin(qJ(4));
t18 = t104 * t100 - t43 * t101 + t44 * t82;
t66 = t109 - t115 / 0.2e1;
t69 = t111 - t117 / 0.2e1;
t86 = cos(qJ(4));
t19 = t104 * t69 + t43 * t66 + t44 * t86;
t75 = sin(pkin(8));
t79 = cos(pkin(8));
t35 = t104 * t79 + t43 * t75;
t81 = sin(qJ(5));
t85 = cos(qJ(5));
t6 = t19 * t85 - t35 * t81;
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t142 = -t18 * t84 + t6 * t80;
t141 = t18 * t80 + t6 * t84;
t138 = t19 * t81 + t35 * t85;
t137 = r_i_i_C(3) + pkin(13);
t133 = t76 * t69;
t132 = t76 * t79;
t131 = t77 * t83;
t129 = t76 * qJ(3);
t128 = cos(pkin(6));
t62 = t83 * t68 - t87 * t88;
t106 = t84 * r_i_i_C(1) - t80 * r_i_i_C(2) + pkin(5);
t105 = t80 * r_i_i_C(1) + t84 * r_i_i_C(2) + pkin(12);
t48 = t59 * t126 - t60 * t99;
t37 = t132 * t60 - t48 * t75;
t93 = t83 * t102 + t88 * t136;
t50 = t93 * t126 + t62 * t99;
t38 = -t62 * t132 - t50 * t75;
t67 = t110 + t116 / 0.2e1;
t70 = t112 - t118 / 0.2e1;
t55 = -t67 * t126 + t70 * t99;
t47 = -t70 * t132 - t55 * t75;
t103 = -t128 * t127 + t67 * t76;
t96 = t76 * t100;
t94 = -t106 * t85 - t137 * t81 - pkin(4);
t92 = -t83 * t119 - t93 * t76;
t91 = t70 * t126 + t128 * t98 + t67 * t99;
t52 = t128 * t65 + t67 * t64 - t70 * t78;
t30 = t103 * t69 + t52 * t86 + t91 * t66;
t90 = -t62 * t126 - t83 * t95 + t93 * t99;
t89 = t90 * t75 - t92 * t79;
t45 = t65 * t131 - t62 * t78 - t93 * t64;
t21 = t45 * t86 - t90 * t66 + t92 * t69;
t56 = t64 * t70 + t67 * t78;
t51 = t62 * t64 - t93 * t78;
t49 = -t59 * t78 - t60 * t64;
t41 = -t103 * t79 - t91 * t75;
t34 = t70 * t133 + t55 * t66 + t56 * t86;
t33 = -t55 * t101 + t56 * t82 + t70 * t96;
t29 = t103 * t100 - t91 * t101 + t52 * t82;
t28 = t62 * t133 + t50 * t66 + t51 * t86;
t27 = -t50 * t101 + t51 * t82 + t62 * t96;
t26 = -t133 * t60 + t48 * t66 + t49 * t86;
t25 = -t48 * t101 + t49 * t82 - t60 * t96;
t24 = t34 * t85 + t47 * t81;
t20 = t92 * t100 + t90 * t101 + t45 * t82;
t14 = t30 * t85 + t41 * t81;
t12 = t28 * t85 + t38 * t81;
t10 = t26 * t85 + t37 * t81;
t8 = t21 * t85 + t89 * t81;
t7 = t21 * t81 - t89 * t85;
t2 = t20 * t80 + t8 * t84;
t1 = t20 * t84 - t8 * t80;
t3 = [-t83 * pkin(1) - t60 * pkin(2) + t44 * pkin(3) + t19 * pkin(4) + t6 * pkin(5) + pkin(10) * t130 - t35 * pkin(11) + t18 * pkin(12) + t141 * r_i_i_C(1) - t142 * r_i_i_C(2) - t104 * qJ(3) + t137 * t138 (t12 * t84 + t27 * t80) * r_i_i_C(1) + (-t12 * t80 + t27 * t84) * r_i_i_C(2) + t12 * pkin(5) + t28 * pkin(4) + t27 * pkin(12) + t51 * pkin(3) - t93 * pkin(2) - t62 * t129 + t137 * (t28 * t81 - t38 * t85) + t38 * pkin(11), -t92, t105 * t21 + t94 * t20, -t106 * t7 + t137 * t8, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t88 * pkin(1) - t62 * pkin(2) + t45 * pkin(3) + t21 * pkin(4) + t8 * pkin(5) + pkin(10) * t131 + t89 * pkin(11) + t20 * pkin(12) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t92 * qJ(3) + t137 * t7 (t10 * t84 + t25 * t80) * r_i_i_C(1) + (-t10 * t80 + t25 * t84) * r_i_i_C(2) + t10 * pkin(5) + t26 * pkin(4) + t25 * pkin(12) + t49 * pkin(3) - t59 * pkin(2) + t60 * t129 + t137 * (t26 * t81 - t37 * t85) + t37 * pkin(11), t104, -t105 * t19 - t18 * t94, t106 * t138 - t137 * t6, t142 * r_i_i_C(1) + t141 * r_i_i_C(2); 0 (t24 * t84 + t33 * t80) * r_i_i_C(1) + (-t24 * t80 + t33 * t84) * r_i_i_C(2) + t24 * pkin(5) + t34 * pkin(4) + t33 * pkin(12) + t56 * pkin(3) + t67 * pkin(2) - t70 * t129 + t137 * (t34 * t81 - t47 * t85) + t47 * pkin(11), -t103, t105 * t30 + t94 * t29, t137 * t14 + t106 * (-t30 * t81 + t41 * t85) (-t14 * t80 + t29 * t84) * r_i_i_C(1) + (-t14 * t84 - t29 * t80) * r_i_i_C(2);];
Ja_transl  = t3;
