% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:41
% EndTime: 2019-12-05 16:17:44
% DurationCPUTime: 1.04s
% Computational Cost: add. (1709->224), mult. (2393->309), div. (0->0), fcn. (1432->10), ass. (0->145)
t87 = sin(qJ(3));
t126 = qJDD(2) * t87;
t89 = cos(qJ(3));
t134 = qJD(3) * t89;
t76 = qJDD(2) + qJDD(3);
t137 = t76 * qJ(4);
t80 = qJD(2) + qJD(3);
t25 = t137 + t80 * qJD(4) + (qJD(2) * t134 + t126) * pkin(2);
t84 = sin(pkin(9));
t85 = cos(pkin(9));
t17 = -t85 * qJDD(1) + t84 * t25;
t165 = t17 * t84;
t18 = t84 * qJDD(1) + t85 * t25;
t100 = t18 * t85 + t165;
t136 = qJD(2) * t89;
t101 = -pkin(2) * t136 + qJD(4);
t79 = pkin(8) + qJ(2);
t72 = sin(t79);
t169 = pkin(2) * t72;
t74 = qJ(3) + t79;
t65 = sin(t74);
t59 = g(1) * t65;
t66 = cos(t74);
t58 = g(2) * t66;
t168 = t76 * pkin(3);
t167 = t85 * pkin(4);
t166 = t89 * pkin(2);
t140 = pkin(2) * qJD(2);
t122 = t87 * t140;
t43 = t80 * qJ(4) + t122;
t35 = -t85 * qJD(1) + t84 * t43;
t164 = t35 * t84;
t155 = t85 * t76;
t51 = -qJDD(5) + t155;
t163 = t51 * t85;
t55 = pkin(2) * t134 + qJD(4);
t162 = t55 * t80;
t161 = t65 * t85;
t160 = t66 * t84;
t86 = sin(qJ(5));
t159 = t76 * t86;
t88 = cos(qJ(5));
t158 = t76 * t88;
t75 = t80 ^ 2;
t77 = t84 ^ 2;
t157 = t77 * t75;
t156 = t84 * t76;
t154 = t85 * t80;
t153 = t85 * t86;
t152 = t85 * t88;
t151 = t85 * t89;
t150 = t86 * t88;
t139 = qJ(4) * t85;
t45 = -t84 * pkin(7) - pkin(3) - t167;
t33 = -t86 * t139 + t88 * t45;
t131 = qJD(5) * t33;
t133 = qJD(4) * t85;
t149 = t88 * t133 + t131 - (t88 * t151 + t86 * t87) * t140;
t34 = t88 * t139 + t86 * t45;
t130 = qJD(5) * t34;
t148 = -t86 * t133 - t130 - (-t86 * t151 + t87 * t88) * t140;
t147 = g(2) * t160 - t84 * t59;
t146 = t66 * pkin(3) + t65 * qJ(4);
t145 = g(1) * t66 + g(2) * t65;
t144 = t59 - t58;
t143 = -qJD(3) * t122 + qJDD(2) * t166;
t78 = t85 ^ 2;
t142 = t77 + t78;
t81 = t86 ^ 2;
t82 = t88 ^ 2;
t141 = t81 - t82;
t23 = t45 * t80 + t101;
t36 = t84 * qJD(1) + t85 * t43;
t8 = t86 * t23 + t88 * t36;
t138 = qJD(5) * t8;
t135 = qJD(3) * t87;
t40 = t45 - t166;
t67 = t87 * pkin(2) + qJ(4);
t21 = -t67 * t153 + t88 * t40;
t132 = qJD(5) * t21;
t129 = qJD(5) * t86;
t128 = qJD(5) * t88;
t127 = qJ(4) * qJD(5);
t118 = t84 * t129;
t7 = t88 * t23 - t86 * t36;
t125 = t7 * t118 - t147;
t116 = qJDD(4) - t143;
t37 = t116 - t168;
t124 = t37 * t84 + t147;
t123 = pkin(2) * t135;
t120 = t80 * t135;
t119 = t80 * t128;
t52 = -qJD(5) + t154;
t117 = t52 * t129;
t54 = t66 * qJ(4);
t115 = -t65 * pkin(3) + t54;
t114 = -t37 - t58;
t113 = -t17 * t85 - g(3);
t112 = t142 * t76;
t111 = t51 - t155;
t110 = t51 + t155;
t109 = t80 * (-qJD(5) - t52);
t108 = t150 * t157;
t107 = pkin(7) * t160 + t66 * t167 + t146;
t106 = -t143 - t144;
t105 = t80 * t122;
t104 = t86 * t119;
t73 = cos(t79);
t103 = g(1) * t72 - g(2) * t73;
t102 = t7 * t86 - t8 * t88;
t99 = t36 * t85 + t164;
t98 = -t145 + t100;
t22 = t67 * t152 + t86 * t40;
t97 = g(3) * t84 - qJD(5) * t23 - t18;
t96 = -qJD(5) * t36 - t80 * t164;
t95 = -t105 - t168;
t68 = -pkin(3) - t166;
t94 = pkin(2) * t120 + t68 * t76;
t93 = -t52 ^ 2 - t157;
t16 = t45 * t76 + t116;
t13 = t88 * t16;
t3 = -t86 * t18 + t13 - t138;
t30 = -t65 * t152 + t66 * t86;
t32 = t66 * t152 + t65 * t86;
t92 = -g(1) * t30 - g(2) * t32 + t128 * t164 + t86 * t165 - t3 * t85;
t2 = t7 * qJD(5) + t86 * t16 + t88 * t18;
t29 = t65 * t153 + t66 * t88;
t31 = -t66 * t153 + t65 * t88;
t91 = -g(1) * t29 - g(2) * t31 - t35 * t118 + t88 * t165 + t2 * t85;
t83 = qJDD(1) - g(3);
t64 = pkin(2) * t73;
t62 = t78 * t76;
t61 = t77 * t76;
t50 = g(1) * t161;
t44 = 0.2e1 * t84 * t155;
t42 = -t80 * pkin(3) + t101;
t41 = t118 * t154;
t27 = (t76 * t82 - 0.2e1 * t104) * t77;
t26 = (t76 * t81 + 0.2e1 * t104) * t77;
t14 = 0.2e1 * (t141 * t80 * qJD(5) - t76 * t150) * t77;
t10 = -t22 * qJD(5) + t88 * t123 - t55 * t153;
t9 = t86 * t123 + t55 * t152 + t132;
t6 = (t110 * t86 + (t52 + t154) * t128) * t84;
t5 = t41 + (-t110 * t88 + t117) * t84;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t84 + t113, 0, 0, 0, 0, 0, 0, (t111 * t86 + (t52 - t154) * t128) * t84, t41 + (t111 * t88 - t117) * t84, 0, (t2 * t88 - t3 * t86 + (-t7 * t88 - t8 * t86) * qJD(5)) * t84 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t103, g(1) * t73 + g(2) * t72, 0, 0, 0, 0, 0, 0, 0, t76, (t76 * t89 - t120) * pkin(2) - t106, ((-qJDD(2) - t76) * t87 + (-qJD(2) - t80) * t134) * pkin(2) + t145, 0, (t103 + (t87 ^ 2 + t89 ^ 2) * qJDD(2) * pkin(2)) * pkin(2), t61, t44, 0, t62, 0, 0, t50 + (t114 - t94) * t85, t94 * t84 + t124, t67 * t112 + t142 * t162 + t98, t37 * t68 + t42 * t123 - g(1) * (t115 - t169) - g(2) * (t64 + t146) + t100 * t67 + t99 * t55, t27, t14, t5, t26, t6, t163, -t10 * t52 - t21 * t51 + (t86 * t162 + (t119 + t159) * t67) * t77 + t92, t22 * t51 + t9 * t52 + (t88 * t162 + (-t129 * t80 + t158) * t67) * t77 + t91, ((-t22 * t76 - t2 + (-t9 + t132) * t80) * t86 + (-t10 * t80 - t21 * t76 - t3 + (-t22 * t80 - t8) * qJD(5)) * t88) * t84 + t125, t2 * t22 + t8 * t9 + t3 * t21 + t7 * t10 - g(1) * (t54 - t169) - g(2) * (t64 + t107) + (t17 * t67 + t35 * t55) * t84 - t45 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t105 - t106, (-t126 + (-qJD(3) + t80) * t136) * pkin(2) + t145, 0, 0, t61, t44, 0, t62, 0, 0, t50 + (t114 - t95) * t85, t95 * t84 + t124, t101 * t80 * t142 + qJ(4) * t112 + t98, -t37 * pkin(3) - g(1) * t115 - g(2) * t146 + t99 * qJD(4) + t100 * qJ(4) + (-t42 * t87 - t99 * t89) * t140, t27, t14, t5, t26, t6, t163, -t33 * t51 - t148 * t52 + (t86 * t137 + (t101 * t86 + t88 * t127) * t80) * t77 + t92, t34 * t51 + t149 * t52 + (t88 * t137 + (t101 * t88 - t127 * t86) * t80) * t77 + t91, ((-t34 * t76 - t2) * t86 + (-t33 * t76 - t138 - t3) * t88 + ((-t130 - t148) * t88 + (t131 - t149) * t86) * t80) * t84 + t125, t2 * t34 + t3 * t33 - g(1) * (-pkin(4) * t161 + t115) - g(2) * t107 + t149 * t8 + t148 * t7 + (pkin(7) * t59 + t17 * qJ(4) + t101 * t35) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t156, -t142 * t75, -t99 * t80 - t144 + t37, 0, 0, 0, 0, 0, 0, -t88 * t51 + t86 * t93, t86 * t51 + t88 * t93, (-t81 - t82) * t156, t2 * t86 + t3 * t88 - t102 * qJD(5) + (t102 * t85 - t164) * t80 - t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t141 * t157, (t86 * t109 + t158) * t84, -t108, (t88 * t109 - t159) * t84, -t51, -g(1) * t31 + g(2) * t29 - t8 * t52 + t86 * t97 + t88 * t96 + t13, g(1) * t32 - g(2) * t30 - t7 * t52 + t97 * t88 + (-t16 - t96) * t86, 0, 0;];
tau_reg = t1;
