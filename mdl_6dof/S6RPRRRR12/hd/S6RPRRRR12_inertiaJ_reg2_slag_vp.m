% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_inertiaJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 06:43:27
% EndTime: 2019-05-06 06:43:44
% DurationCPUTime: 4.80s
% Computational Cost: add. (12691->339), mult. (35394->691), div. (0->0), fcn. (41596->16), ass. (0->166)
t108 = sin(qJ(4));
t112 = cos(qJ(4));
t103 = cos(pkin(8));
t101 = sin(pkin(6));
t105 = cos(pkin(6));
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t102 = cos(pkin(14));
t104 = cos(pkin(7));
t142 = t102 * t104;
t100 = sin(pkin(7));
t145 = t100 * t109;
t98 = sin(pkin(14));
t51 = t105 * t145 + (t109 * t142 + t113 * t98) * t101;
t178 = pkin(11) * t51;
t144 = t100 * t113;
t148 = t98 * t101;
t177 = pkin(1) * t105;
t86 = t102 * t177;
t52 = t105 * pkin(2) + t86 + (-pkin(10) * t104 - qJ(2)) * t148;
t167 = t104 * t52;
t116 = t100 * t105 + t101 * t142;
t143 = t101 * t102;
t72 = qJ(2) * t143 + t98 * t177;
t49 = t116 * pkin(10) + t72;
t61 = (-pkin(10) * t100 * t98 - pkin(2) * t102 - pkin(1)) * t101;
t27 = -t109 * t49 + t113 * t167 + t61 * t144;
t69 = -t100 * t143 + t105 * t104;
t22 = t69 * pkin(3) - t103 * t178 + t27;
t35 = -t100 * t52 + t104 * t61;
t50 = -t109 * t148 + t116 * t113;
t99 = sin(pkin(8));
t24 = -t50 * pkin(3) - t99 * t178 + t35;
t123 = t103 * t22 + t24 * t99;
t122 = t103 * t50 + t69 * t99;
t28 = t113 * t49 + (t100 * t61 + t167) * t109;
t18 = t122 * pkin(11) + t28;
t11 = -t108 * t18 + t123 * t112;
t107 = sin(qJ(5));
t111 = cos(qJ(5));
t32 = t122 * t108 + t51 * t112;
t37 = t69 * t103 - t50 * t99;
t19 = t32 * t107 - t37 * t111;
t193 = t19 ^ 2;
t140 = t103 * t112;
t146 = t99 * t112;
t30 = t51 * t108 - t50 * t140 - t69 * t146;
t192 = t30 ^ 2;
t139 = t103 * t113;
t147 = t99 * t108;
t55 = t104 * t147 + (t108 * t139 + t109 * t112) * t100;
t70 = t103 * t104 - t99 * t144;
t38 = t107 * t55 - t111 * t70;
t191 = t38 ^ 2;
t53 = -t112 * t100 * t139 - t104 * t146 + t108 * t145;
t190 = t53 ^ 2;
t73 = -t111 * t103 + t107 * t147;
t189 = t73 ^ 2;
t188 = -0.2e1 * t19;
t187 = -0.2e1 * t30;
t106 = sin(qJ(6));
t110 = cos(qJ(6));
t75 = t107 * t103 + t111 * t147;
t57 = t106 * t75 + t110 * t146;
t186 = -0.2e1 * t57;
t185 = 0.2e1 * t69;
t184 = -0.2e1 * t75;
t183 = 0.2e1 * t99;
t182 = 0.2e1 * t101;
t181 = -0.2e1 * t107;
t180 = 0.2e1 * t111;
t179 = pkin(3) * t99;
t176 = pkin(3) * t108;
t175 = pkin(5) * t110;
t174 = pkin(12) * t106;
t173 = t107 * pkin(12);
t172 = t19 * t73;
t13 = t103 * t24 - t99 * t22;
t10 = t30 * pkin(4) - t32 * pkin(12) + t13;
t12 = t123 * t108 + t112 * t18;
t9 = t37 * pkin(12) + t12;
t5 = t111 * t10 - t107 * t9;
t3 = -t30 * pkin(5) - t5;
t171 = t3 * t106;
t170 = t3 * t110;
t94 = t106 ^ 2;
t96 = t110 ^ 2;
t169 = t94 + t96;
t93 = t101 ^ 2;
t168 = t102 * t93;
t166 = t106 * t19;
t165 = t106 * t73;
t164 = t107 * t30;
t163 = t110 * t19;
t162 = t110 * t73;
t161 = t111 * t30;
t91 = t99 ^ 2;
t160 = t112 * t91;
t21 = t37 * t107 + t32 * t111;
t14 = t21 * t106 - t30 * t110;
t159 = t14 * t110;
t16 = t30 * t106 + t21 * t110;
t158 = t16 * t106;
t157 = t19 * t111;
t156 = t21 * t107;
t155 = t38 * t107;
t132 = pkin(11) * t146;
t67 = t132 + (pkin(12) + t176) * t103;
t68 = (-pkin(4) * t112 - pkin(12) * t108 - pkin(3)) * t99;
t44 = -t107 * t67 + t111 * t68;
t42 = pkin(5) * t146 - t44;
t154 = t42 * t106;
t153 = t42 * t110;
t152 = t57 * t110;
t59 = -t106 * t146 + t110 * t75;
t151 = t59 * t106;
t150 = t73 * t111;
t149 = t75 * t107;
t141 = t103 * t108;
t138 = t106 * t107;
t137 = t106 * t110;
t136 = t106 * t111;
t135 = t110 * t107;
t134 = t110 * t111;
t133 = 0.2e1 * t146;
t131 = t105 * t182;
t130 = t107 * t180;
t129 = t30 * t146;
t128 = t107 * t146;
t127 = t111 * t146;
t126 = t106 * t135;
t6 = t107 * t10 + t111 * t9;
t4 = t30 * pkin(13) + t6;
t8 = -t37 * pkin(4) - t11;
t7 = t19 * pkin(5) - t21 * pkin(13) + t8;
t1 = -t106 * t4 + t110 * t7;
t2 = t106 * t7 + t110 * t4;
t125 = -t1 * t106 + t2 * t110;
t124 = -t5 * t107 + t6 * t111;
t87 = pkin(11) * t147;
t66 = t87 + (-pkin(3) * t112 - pkin(4)) * t103;
t41 = t73 * pkin(5) - t75 * pkin(13) + t66;
t45 = t107 * t68 + t111 * t67;
t43 = -pkin(13) * t146 + t45;
t25 = -t106 * t43 + t110 * t41;
t26 = t106 * t41 + t110 * t43;
t121 = -t25 * t106 + t26 * t110;
t40 = t107 * t70 + t111 * t55;
t33 = -t106 * t40 + t110 * t53;
t34 = t106 * t53 + t110 * t40;
t120 = -t33 * t106 + t34 * t110;
t79 = -t111 * pkin(5) - t107 * pkin(13) - pkin(4);
t63 = -pkin(12) * t136 + t110 * t79;
t64 = pkin(12) * t134 + t106 * t79;
t119 = -t63 * t106 + t64 * t110;
t118 = t40 * t111 + t155;
t117 = -t44 * t107 + t45 * t111;
t115 = pkin(12) ^ 2;
t97 = t111 ^ 2;
t95 = t107 ^ 2;
t90 = t95 * t115;
t88 = t91 * t112 ^ 2;
t77 = pkin(3) * t141 + t132;
t76 = pkin(3) * t140 - t87;
t71 = -qJ(2) * t148 + t86;
t15 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t93 * t98 ^ 2, 0.2e1 * t98 * t168, t98 * t131, t93 * t102 ^ 2, t102 * t131, t105 ^ 2, 0.2e1 * pkin(1) * t168 + 0.2e1 * t71 * t105, -0.2e1 * t93 * pkin(1) * t98 - 0.2e1 * t72 * t105 (t102 * t72 - t71 * t98) * t182, t93 * pkin(1) ^ 2 + t71 ^ 2 + t72 ^ 2, t51 ^ 2, 0.2e1 * t51 * t50, t51 * t185, t50 ^ 2, t50 * t185, t69 ^ 2, 0.2e1 * t27 * t69 - 0.2e1 * t35 * t50, -0.2e1 * t28 * t69 + 0.2e1 * t35 * t51, -0.2e1 * t27 * t51 + 0.2e1 * t28 * t50, t27 ^ 2 + t28 ^ 2 + t35 ^ 2, t32 ^ 2, t32 * t187, 0.2e1 * t37 * t32, t192, t37 * t187, t37 ^ 2, 0.2e1 * t11 * t37 + 0.2e1 * t13 * t30, -0.2e1 * t12 * t37 + 0.2e1 * t13 * t32, -0.2e1 * t11 * t32 - 0.2e1 * t12 * t30, t11 ^ 2 + t12 ^ 2 + t13 ^ 2, t21 ^ 2, t21 * t188, 0.2e1 * t21 * t30, t193, t19 * t187, t192, 0.2e1 * t8 * t19 + 0.2e1 * t5 * t30, 0.2e1 * t8 * t21 - 0.2e1 * t6 * t30, -0.2e1 * t6 * t19 - 0.2e1 * t5 * t21, t5 ^ 2 + t6 ^ 2 + t8 ^ 2, t16 ^ 2, -0.2e1 * t16 * t14, 0.2e1 * t16 * t19, t14 ^ 2, t14 * t188, t193, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t14, 0.2e1 * t3 * t16 - 0.2e1 * t2 * t19, -0.2e1 * t1 * t16 - 0.2e1 * t2 * t14, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, t148, 0, -t101 * pkin(1), 0, 0, 0, 0, 0, 0, -t104 * t50 + t69 * t144, t104 * t51 - t69 * t145 (t109 * t50 - t113 * t51) * t100, t35 * t104 + (t109 * t28 + t113 * t27) * t100, 0, 0, 0, 0, 0, 0, t70 * t30 - t53 * t37, t70 * t32 - t55 * t37, -t55 * t30 + t53 * t32, -t11 * t53 + t12 * t55 + t13 * t70, 0, 0, 0, 0, 0, 0, t53 * t19 - t38 * t30, t53 * t21 - t40 * t30, -t40 * t19 + t38 * t21, -t5 * t38 + t6 * t40 + t8 * t53, 0, 0, 0, 0, 0, 0, t38 * t14 + t33 * t19, t38 * t16 - t34 * t19, -t34 * t14 - t33 * t16, t1 * t33 + t2 * t34 + t3 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 ^ 2 + (t109 ^ 2 + t113 ^ 2) * t100 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 ^ 2 + t70 ^ 2 + t190, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 ^ 2 + t190 + t191, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 ^ 2 + t34 ^ 2 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, t50, t69, t27, -t28, 0, 0, t32 * t147 (-t108 * t30 + t112 * t32) * t99, t103 * t32 + t37 * t147, -t129, -t103 * t30 + t37 * t146, t37 * t103, t11 * t103 + t76 * t37 + (-pkin(3) * t30 - t112 * t13) * t99, -t12 * t103 - t77 * t37 + (-pkin(3) * t32 + t108 * t13) * t99, -t77 * t30 - t76 * t32 + (-t108 * t11 + t112 * t12) * t99, t11 * t76 + t12 * t77 - t13 * t179, t21 * t75, -t75 * t19 - t21 * t73, -t146 * t21 + t75 * t30, t172, t146 * t19 - t73 * t30, -t129, -t146 * t5 + t66 * t19 + t44 * t30 + t8 * t73, t146 * t6 + t66 * t21 - t45 * t30 + t8 * t75, -t45 * t19 - t44 * t21 - t5 * t75 - t6 * t73, t5 * t44 + t6 * t45 + t8 * t66, t16 * t59, -t59 * t14 - t16 * t57, t16 * t73 + t59 * t19, t14 * t57, -t14 * t73 - t57 * t19, t172, t1 * t73 + t42 * t14 + t25 * t19 + t3 * t57, t42 * t16 - t26 * t19 - t2 * t73 + t3 * t59, -t1 * t59 - t26 * t14 - t25 * t16 - t2 * t57, t1 * t25 + t2 * t26 + t3 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, -t145, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t103 - t70 * t146, -t55 * t103 + t70 * t147 (t108 * t53 + t112 * t55) * t99, -t70 * t179 - t53 * t76 + t55 * t77, 0, 0, 0, 0, 0, 0, t146 * t38 + t53 * t73, t146 * t40 + t53 * t75, t38 * t75 - t40 * t73, -t38 * t44 + t40 * t45 + t53 * t66, 0, 0, 0, 0, 0, 0, t33 * t73 + t38 * t57, -t34 * t73 + t38 * t59, -t33 * t59 - t34 * t57, t33 * t25 + t34 * t26 + t38 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t91 * t108 ^ 2, 0.2e1 * t108 * t160, t141 * t183, t88, t103 * t133, t103 ^ 2, 0.2e1 * pkin(3) * t160 + 0.2e1 * t76 * t103, -0.2e1 * t77 * t103 - 0.2e1 * t91 * t176 (-t108 * t76 + t112 * t77) * t183, t91 * pkin(3) ^ 2 + t76 ^ 2 + t77 ^ 2, t75 ^ 2, t73 * t184, t146 * t184, t189, t73 * t133, t88, -0.2e1 * t146 * t44 + 0.2e1 * t66 * t73, 0.2e1 * t146 * t45 + 0.2e1 * t66 * t75, -0.2e1 * t44 * t75 - 0.2e1 * t45 * t73, t44 ^ 2 + t45 ^ 2 + t66 ^ 2, t59 ^ 2, t59 * t186, 0.2e1 * t59 * t73, t57 ^ 2, t73 * t186, t189, 0.2e1 * t25 * t73 + 0.2e1 * t42 * t57, -0.2e1 * t26 * t73 + 0.2e1 * t42 * t59, -0.2e1 * t25 * t59 - 0.2e1 * t26 * t57, t25 ^ 2 + t26 ^ 2 + t42 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t30, t37, t11, -t12, 0, 0, t156, -t107 * t19 + t21 * t111, t164, -t157, t161, 0, -pkin(4) * t19 - pkin(12) * t164 - t8 * t111, -pkin(4) * t21 - pkin(12) * t161 + t8 * t107 (t156 - t157) * pkin(12) + t124, -t8 * pkin(4) + pkin(12) * t124, t16 * t135 (-t158 - t159) * t107, -t16 * t111 + t135 * t19, t14 * t138, t14 * t111 - t138 * t19, -t157, -t1 * t111 + t63 * t19 + (pkin(12) * t14 + t171) * t107, t2 * t111 - t64 * t19 + (pkin(12) * t16 + t170) * t107, -t64 * t14 - t63 * t16 + (-t1 * t110 - t106 * t2) * t107, t1 * t63 + t173 * t3 + t2 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t55, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * t111, t53 * t107, t118, -t53 * pkin(4) + pkin(12) * t118, 0, 0, 0, 0, 0, 0, -t33 * t111 + t138 * t38, t34 * t111 + t135 * t38 (-t106 * t34 - t110 * t33) * t107, pkin(12) * t155 + t33 * t63 + t34 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, t146, t103, t76, -t77, 0, 0, t149, -t107 * t73 + t75 * t111, -t128, -t150, -t127, 0, -pkin(4) * t73 + pkin(12) * t128 - t66 * t111, -pkin(4) * t75 + pkin(12) * t127 + t66 * t107 (t149 - t150) * pkin(12) + t117, -t66 * pkin(4) + pkin(12) * t117, t59 * t135 (-t151 - t152) * t107, -t59 * t111 + t135 * t73, t57 * t138, t57 * t111 - t138 * t73, -t150, -t25 * t111 + t63 * t73 + (pkin(12) * t57 + t154) * t107, t26 * t111 - t64 * t73 + (pkin(12) * t59 + t153) * t107, -t64 * t57 - t63 * t59 + (-t106 * t26 - t110 * t25) * t107, t173 * t42 + t25 * t63 + t26 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t95, t130, 0, t97, 0, 0, pkin(4) * t180, pkin(4) * t181, 0.2e1 * (t95 + t97) * pkin(12), pkin(4) ^ 2 + t97 * t115 + t90, t96 * t95, -0.2e1 * t95 * t137, t134 * t181, t94 * t95, t106 * t130, t97, -0.2e1 * t63 * t111 + 0.2e1 * t174 * t95, 0.2e1 * t95 * pkin(12) * t110 + 0.2e1 * t64 * t111, 0.2e1 * (-t106 * t64 - t110 * t63) * t107, t63 ^ 2 + t64 ^ 2 + t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, t30, t5, -t6, 0, 0, t158, -t106 * t14 + t16 * t110, t166, -t159, t163, 0, -pkin(5) * t14 - pkin(13) * t166 - t170, -pkin(5) * t16 - pkin(13) * t163 + t171 (t158 - t159) * pkin(13) + t125, -t3 * pkin(5) + pkin(13) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, -t40, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t110, t38 * t106, t120, -t38 * pkin(5) + pkin(13) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, -t73, -t146, t44, -t45, 0, 0, t151, -t106 * t57 + t59 * t110, t165, -t152, t162, 0, -pkin(5) * t57 - pkin(13) * t165 - t153, -pkin(5) * t59 - pkin(13) * t162 + t154 (t151 - t152) * pkin(13) + t121, -t42 * pkin(5) + pkin(13) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, 0, t111, 0, -t173, -t111 * pkin(12), 0, 0, t126 (-t94 + t96) * t107, -t136, -t126, -t134, 0, -pkin(12) * t135 + (-pkin(5) * t107 + pkin(13) * t111) * t106, pkin(13) * t134 + (t174 - t175) * t107, t119, -pkin(5) * t173 + pkin(13) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t94, 0.2e1 * t137, 0, t96, 0, 0, 0.2e1 * t175, -0.2e1 * pkin(5) * t106, 0.2e1 * t169 * pkin(13), pkin(13) ^ 2 * t169 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, t19, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t57, t73, t25, -t26, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, 0, -t138, -t111, t63, -t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, 0, t110, 0, -t106 * pkin(13), -t110 * pkin(13), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t15;
