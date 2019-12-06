% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:40:01
% EndTime: 2019-12-05 17:40:05
% DurationCPUTime: 1.02s
% Computational Cost: add. (1208->187), mult. (2435->243), div. (0->0), fcn. (1848->12), ass. (0->117)
t100 = sin(qJ(5));
t103 = cos(qJ(5));
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t97 = sin(pkin(8));
t98 = cos(pkin(8));
t57 = t101 * t98 + t104 * t97;
t47 = t57 * qJD(1);
t143 = qJD(1) * t97;
t133 = t101 * t143;
t145 = t104 * t98;
t135 = qJD(1) * t145;
t49 = -t133 + t135;
t18 = t100 * t49 + t103 * t47;
t93 = qJD(4) + qJD(5);
t156 = t18 * t93;
t139 = qJD(5) * t103;
t140 = qJD(5) * t100;
t137 = t98 * qJDD(1);
t138 = t97 * qJDD(1);
t122 = -t101 * t138 + t104 * t137;
t50 = t57 * qJD(4);
t25 = -qJD(1) * t50 + t122;
t113 = -qJD(4) * t133 + t57 * qJDD(1);
t141 = qJD(4) * t104;
t134 = t98 * t141;
t26 = qJD(1) * t134 + t113;
t4 = -t100 * t26 + t103 * t25 - t47 * t139 - t140 * t49;
t174 = t4 + t156;
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t168 = g(1) * t102 - g(2) * t105;
t118 = -t100 * t47 + t103 * t49;
t173 = t118 * t18;
t157 = t118 * t93;
t5 = qJD(5) * t118 + t100 * t25 + t103 * t26;
t172 = -t5 + t157;
t58 = -t101 * t97 + t145;
t171 = t118 ^ 2 - t18 ^ 2;
t99 = -pkin(1) - qJ(3);
t67 = t99 * qJD(1) + qJD(2);
t131 = -pkin(6) * qJD(1) + t67;
t41 = t131 * t97;
t42 = t131 * t98;
t116 = -t101 * t42 - t104 * t41;
t13 = -pkin(7) * t47 - t116;
t161 = -qJD(1) * qJD(3) + qJDD(1) * t99;
t59 = qJDD(2) + t161;
t127 = -pkin(6) * qJDD(1) + t59;
t31 = t127 * t97;
t32 = t127 * t98;
t129 = -t101 * t31 + t104 * t32;
t2 = qJDD(4) * pkin(4) - pkin(7) * t25 + qJD(4) * t116 + t129;
t79 = qJD(1) * qJ(2) + qJD(3);
t62 = pkin(3) * t143 + t79;
t29 = pkin(4) * t47 + t62;
t92 = pkin(8) + qJ(4);
t82 = qJ(5) + t92;
t75 = sin(t82);
t76 = cos(t82);
t170 = t29 * t18 + g(3) * t76 + t13 * t140 + (-t13 * t93 - t2) * t100 + t168 * t75;
t150 = t97 ^ 2 + t98 ^ 2;
t94 = qJDD(1) * qJ(2);
t95 = qJD(1) * qJD(2);
t167 = t94 + t95;
t166 = t150 * t67;
t165 = -t101 * t41 + t104 * t42;
t28 = -t100 * t57 + t103 * t58;
t124 = g(1) * t105 + g(2) * t102;
t65 = qJDD(3) + t167;
t164 = t65 - t124;
t27 = t100 * t58 + t103 * t57;
t142 = qJD(4) * t101;
t51 = -t142 * t97 + t134;
t6 = -qJD(5) * t27 - t100 * t51 - t103 * t50;
t89 = qJDD(4) + qJDD(5);
t163 = t28 * t89 + t6 * t93;
t12 = -pkin(7) * t49 + t165;
t11 = qJD(4) * pkin(4) + t12;
t146 = t103 * t13;
t121 = -t100 * t11 - t146;
t117 = t101 * t32 + t104 * t31;
t3 = -pkin(7) * t26 + t165 * qJD(4) + t117;
t162 = g(3) * t75 + qJD(5) * t121 - t100 * t3 + t103 * t2 - t29 * t118 - t168 * t76;
t160 = 0.2e1 * t95;
t158 = -pkin(6) + t99;
t154 = -t50 * qJD(4) + t58 * qJDD(4);
t60 = t158 * t97;
t61 = t158 * t98;
t153 = t101 * t61 + t104 * t60;
t152 = t105 * pkin(1) + t102 * qJ(2);
t74 = t97 * pkin(3) + qJ(2);
t144 = pkin(1) * qJDD(1);
t132 = t150 * t59;
t128 = -t101 * t60 + t104 * t61;
t53 = pkin(3) * t138 + t65;
t126 = qJDD(2) - t144;
t7 = t28 * qJD(5) - t100 * t50 + t103 * t51;
t125 = -t27 * t89 - t7 * t93;
t15 = -pkin(7) * t58 + t128;
t16 = -pkin(7) * t57 + t153;
t120 = -t100 * t16 + t103 * t15;
t119 = t100 * t15 + t103 * t16;
t115 = -qJD(4) * t51 - qJDD(4) * t57;
t110 = t164 + t167;
t109 = -qJD(3) * t57 + t61 * t141 - t142 * t60;
t108 = -t58 * qJD(3) - qJD(4) * t153;
t106 = qJD(1) ^ 2;
t85 = t105 * qJ(2);
t81 = cos(t92);
t80 = sin(t92);
t36 = pkin(4) * t51 + qJD(2);
t34 = pkin(4) * t57 + t74;
t14 = pkin(4) * t26 + t53;
t10 = pkin(7) * t50 + t108;
t9 = -pkin(7) * t51 + t109;
t1 = [qJDD(1), t168, t124, qJDD(2) - 0.2e1 * t144 - t168, -t124 + 0.2e1 * t94 + t160, -t126 * pkin(1) - g(1) * (-pkin(1) * t102 + t85) - g(2) * t152 + (t94 + t160) * qJ(2), t110 * t97, t110 * t98, t168 + t150 * (-t161 - t59), t65 * qJ(2) + t79 * qJD(2) - g(1) * (t99 * t102 + t85) - g(2) * (qJ(3) * t105 + t152) + t99 * t132 - qJD(3) * t166, t25 * t58 - t49 * t50, -t25 * t57 - t26 * t58 + t47 * t50 - t49 * t51, t154, t115, 0, qJD(2) * t47 + qJD(4) * t108 + qJDD(4) * t128 - t124 * t80 + t74 * t26 + t62 * t51 + t53 * t57, qJD(2) * t49 - qJD(4) * t109 - qJDD(4) * t153 - t124 * t81 + t74 * t25 - t62 * t50 + t53 * t58, t118 * t6 + t28 * t4, -t118 * t7 - t18 * t6 - t27 * t4 - t28 * t5, t163, t125, 0, t36 * t18 + t34 * t5 + t14 * t27 + t29 * t7 + (-qJD(5) * t119 + t10 * t103 - t100 * t9) * t93 + t120 * t89 - t124 * t75, t36 * t118 + t34 * t4 + t14 * t28 + t29 * t6 - (qJD(5) * t120 + t10 * t100 + t103 * t9) * t93 - t119 * t89 - t124 * t76; 0, 0, 0, qJDD(1), -t106, -qJ(2) * t106 + t126 - t168, -t106 * t97, -t106 * t98, -t150 * qJDD(1), -qJD(1) * t79 + t132 - t168, 0, 0, 0, 0, 0, -qJD(1) * t47 + t154, -qJD(1) * t49 + t115, 0, 0, 0, 0, 0, -qJD(1) * t18 + t163, -qJD(1) * t118 + t125; 0, 0, 0, 0, 0, 0, t138, t137, -t150 * t106, qJD(1) * t166 + t164, 0, 0, 0, 0, 0, (t49 + t135) * qJD(4) + t113, -0.2e1 * qJD(4) * t47 + t122, 0, 0, 0, 0, 0, t5 + t157, t4 - t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t47, -t47 ^ 2 + t49 ^ 2, t122, (t49 - t135) * qJD(4) - t113, qJDD(4), g(3) * t80 - t168 * t81 - t62 * t49 + t129, g(3) * t81 + t168 * t80 + t62 * t47 - t117, t173, t171, t174, t172, t89, -(-t100 * t12 - t146) * t93 + (t103 * t89 - t140 * t93 - t18 * t49) * pkin(4) + t162, (-qJD(5) * t11 + t12 * t93 - t3) * t103 + (-t100 * t89 - t118 * t49 - t139 * t93) * pkin(4) + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t171, t174, t172, t89, -t121 * t93 + t162, (-t3 + (-qJD(5) + t93) * t11) * t103 + t170;];
tau_reg = t1;
