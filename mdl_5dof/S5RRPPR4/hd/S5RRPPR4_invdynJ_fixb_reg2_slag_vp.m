% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:53
% EndTime: 2019-12-31 19:27:55
% DurationCPUTime: 0.97s
% Computational Cost: add. (1826->215), mult. (2230->253), div. (0->0), fcn. (1199->10), ass. (0->135)
t90 = qJ(1) + qJ(2);
t81 = sin(t90);
t82 = cos(t90);
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t40 = -t81 * t91 - t82 * t92;
t41 = -t81 * t92 + t82 * t91;
t120 = -g(1) * t41 + g(2) * t40;
t150 = pkin(1) * qJD(1);
t97 = cos(qJ(2));
t133 = t97 * t150;
t117 = qJD(3) - t133;
t86 = qJD(1) + qJD(2);
t99 = -pkin(2) - pkin(3);
t37 = t99 * t86 + t117;
t94 = sin(qJ(2));
t134 = t94 * t150;
t46 = t86 * qJ(3) + t134;
t16 = t91 * t37 + t92 * t46;
t14 = -t86 * pkin(7) + t16;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t11 = t96 * qJD(4) - t93 * t14;
t144 = t11 * qJD(5);
t130 = qJD(2) * t150;
t149 = pkin(1) * qJDD(1);
t153 = -t94 * t130 + t97 * t149;
t129 = -qJDD(3) + t153;
t85 = qJDD(1) + qJDD(2);
t21 = t99 * t85 - t129;
t156 = t97 * t130 + t94 * t149;
t76 = t85 * qJ(3);
t77 = t86 * qJD(3);
t22 = t76 + t77 + t156;
t8 = t91 * t21 + t92 * t22;
t6 = -t85 * pkin(7) + t8;
t2 = t93 * qJDD(4) + t96 * t6 + t144;
t159 = t96 * t14;
t12 = t93 * qJD(4) + t159;
t79 = t96 * qJDD(4);
t3 = -t12 * qJD(5) - t93 * t6 + t79;
t173 = t2 * t96 - t3 * t93;
t172 = -g(1) * t82 - g(2) * t81;
t39 = (t91 * t94 + t92 * t97) * t150;
t123 = t92 * qJD(3) - t39;
t171 = t123 * t86;
t38 = t91 * t133 - t92 * t134;
t170 = (t91 * qJD(3) - t38) * t86;
t84 = t86 ^ 2;
t166 = t85 * pkin(2);
t95 = sin(qJ(1));
t165 = t95 * pkin(1);
t15 = t92 * t37 - t91 * t46;
t13 = t86 * pkin(4) - t15;
t164 = t13 * t86;
t148 = qJD(2) * t94;
t135 = pkin(1) * t148;
t147 = qJD(2) * t97;
t59 = pkin(1) * t147 + qJD(3);
t30 = -t92 * t135 + t91 * t59;
t163 = t30 * t86;
t31 = t91 * t135 + t92 * t59;
t162 = t31 * t86;
t161 = t92 * t84;
t160 = t93 * t85;
t158 = t96 * t85;
t157 = t120 * t93;
t73 = -t97 * pkin(1) - pkin(2);
t62 = -pkin(3) + t73;
t65 = t94 * pkin(1) + qJ(3);
t29 = t91 * t62 + t92 * t65;
t155 = t82 * pkin(2) + t81 * qJ(3);
t154 = g(1) * t81 - g(2) * t82;
t48 = t92 * qJ(3) + t91 * t99;
t88 = t93 ^ 2;
t89 = t96 ^ 2;
t152 = t88 - t89;
t151 = t88 + t89;
t146 = qJD(5) * t86;
t145 = qJD(5) * t93;
t142 = qJDD(4) + g(3);
t141 = qJDD(5) * t93;
t140 = qJDD(5) * t96;
t100 = qJD(5) ^ 2;
t139 = t100 + t84;
t138 = t93 * t84 * t96;
t67 = t82 * pkin(3);
t137 = t67 + t155;
t98 = cos(qJ(1));
t83 = t98 * pkin(1);
t136 = t83 + t155;
t132 = t86 * t148;
t131 = t96 * t146;
t64 = t82 * qJ(3);
t128 = -t81 * pkin(2) + t64;
t127 = t151 * t85;
t7 = t92 * t21 - t91 * t22;
t126 = t153 + t154;
t125 = t156 + t172;
t122 = t93 * t131;
t121 = t99 * t81 + t64;
t119 = g(1) * t40 + g(2) * t41;
t118 = g(1) * t95 - g(2) * t98;
t116 = t11 * t93 - t12 * t96;
t115 = t15 * t91 - t16 * t92;
t28 = t92 * t62 - t91 * t65;
t114 = -qJDD(3) + t126;
t47 = -t91 * qJ(3) + t92 * t99;
t44 = pkin(4) - t47;
t45 = -pkin(7) + t48;
t113 = t100 * t45 - t44 * t85;
t32 = -t129 - t166;
t112 = -t40 * pkin(4) + t41 * pkin(7) + t137;
t111 = t119 + t8;
t23 = pkin(4) - t28;
t24 = -pkin(7) + t29;
t110 = t100 * t24 - t23 * t85 - t163;
t109 = t41 * pkin(4) + t40 * pkin(7) + t121;
t108 = t86 * t133 - t125;
t107 = t12 * t145 + t96 * t144 - t119 - t173;
t106 = -qJDD(5) * t24 + (-t23 * t86 - t13 - t31) * qJD(5);
t105 = t120 + t170;
t104 = -qJD(4) * qJD(5) - t119 + t164 - t6;
t103 = (-t11 * t96 - t12 * t93) * qJD(5) + t173;
t102 = -qJDD(5) * t45 + (-t44 * t86 - t123 - t13) * qJD(5);
t53 = t86 * t134;
t50 = -t100 * t93 + t140;
t49 = -t100 * t96 - t141;
t43 = -t86 * pkin(2) + t117;
t34 = t88 * t85 + 0.2e1 * t122;
t33 = t89 * t85 - 0.2e1 * t122;
t18 = -0.2e1 * t152 * t146 + 0.2e1 * t93 * t158;
t5 = t85 * pkin(4) - t7;
t4 = t5 * t96;
t1 = [0, 0, 0, 0, 0, qJDD(1), t118, g(1) * t98 + g(2) * t95, 0, 0, 0, 0, 0, 0, 0, t85, (t85 * t97 - t132) * pkin(1) + t126, (-t86 * t147 - t85 * t94) * pkin(1) - t125, 0, (t118 + (t94 ^ 2 + t97 ^ 2) * t149) * pkin(1), 0, 0, 0, t85, 0, 0, -pkin(1) * t132 + (pkin(2) - t73) * t85 + t114, 0, t59 * t86 + t65 * t85 + t172 + t22, t22 * t65 + t46 * t59 + t32 * t73 + t43 * t135 - g(1) * (t128 - t165) - g(2) * t136, 0, 0, 0, 0, 0, t85, -t28 * t85 + t120 + t163 - t7, t29 * t85 + t111 + t162, 0, t8 * t29 + t16 * t31 + t7 * t28 - t15 * t30 - g(1) * (t121 - t165) - g(2) * (t67 + t136), t34, t18, t49, t33, -t50, 0, t4 + t106 * t93 + (-t110 + t120) * t96, t106 * t96 + (t110 - t5) * t93 - t157, -t24 * t127 - t151 * t162 + t107, t5 * t23 + t13 * t30 - g(1) * (t109 - t165) - g(2) * (t112 + t83) - t116 * t31 + t103 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t53 + t126, t108, 0, 0, 0, 0, 0, t85, 0, 0, t114 + t53 + 0.2e1 * t166, 0, -t108 + 0.2e1 * t76 + 0.2e1 * t77, t22 * qJ(3) + t46 * qJD(3) - t32 * pkin(2) - g(1) * t128 - g(2) * t155 + (-t43 * t94 - t46 * t97) * t150, 0, 0, 0, 0, 0, t85, -t47 * t85 + t105 - t7, t48 * t85 + t111 + t171, 0, -g(1) * t121 - g(2) * t137 - t115 * qJD(3) + t15 * t38 - t16 * t39 + t7 * t47 + t8 * t48, t34, t18, t49, t33, -t50, 0, t4 + t102 * t93 + (t105 - t113) * t96, t102 * t96 + (t113 - t5 - t170) * t93 - t157, -t45 * t127 - t151 * t171 + t107, t5 * t44 - t13 * t38 - g(1) * t109 - g(2) * t112 + t116 * t39 + (-t116 * t92 + t13 * t91) * qJD(3) + t103 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, 0, -t84, -t46 * t86 - t154 + t32, 0, 0, 0, 0, 0, 0, -t91 * t84 - t92 * t85, t91 * t85 - t161, 0, t115 * t86 + t7 * t92 + t8 * t91 - t154, 0, 0, 0, 0, 0, 0, (0.2e1 * t86 * t145 - t158) * t92 + (-t139 * t96 - t141) * t91, (0.2e1 * t131 + t160) * t92 + (t139 * t93 - t140) * t91, -t91 * t127 + t151 * t161, (t116 * t86 - t5) * t92 + (t103 - t164) * t91 - t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, 0, 0, 0, 0, 0, 0, t50, t49, 0, -t116 * qJD(5) + t2 * t93 + t3 * t96 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, t152 * t84, -t160, t138, -t158, qJDD(5), g(3) * t96 + t79 + (t12 - t159) * qJD(5) + t104 * t93, t144 + (qJD(5) * t14 - t142) * t93 + t104 * t96, 0, 0;];
tau_reg = t1;
