% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:39:22
% EndTime: 2020-01-03 11:39:28
% DurationCPUTime: 1.21s
% Computational Cost: add. (1460->199), mult. (3181->277), div. (0->0), fcn. (2292->14), ass. (0->129)
t105 = qJD(3) + qJD(5);
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t109 = sin(pkin(9));
t111 = cos(pkin(9));
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t80 = -t109 * t115 + t111 * t118;
t73 = t80 * qJD(1);
t65 = t117 * t73;
t81 = t109 * t118 + t111 * t115;
t75 = t81 * qJD(1);
t35 = -t114 * t75 + t65;
t156 = t35 * t105;
t152 = qJD(5) * t114;
t74 = t81 * qJD(3);
t42 = -qJD(1) * t74 + qJDD(1) * t80;
t150 = qJD(1) * qJD(3);
t146 = t118 * t150;
t147 = t115 * t150;
t43 = qJDD(1) * t81 - t109 * t147 + t111 * t146;
t7 = qJD(5) * t65 + t114 * t42 + t117 * t43 - t152 * t75;
t174 = t7 - t156;
t135 = t114 * t73 + t117 * t75;
t173 = t135 * t35;
t110 = sin(pkin(8));
t95 = pkin(1) * t110 + pkin(6);
t155 = qJ(4) + t95;
t157 = t135 * t105;
t8 = qJD(5) * t135 + t114 * t43 - t117 * t42;
t172 = -t8 + t157;
t106 = qJ(1) + pkin(8);
t100 = cos(t106);
t99 = sin(t106);
t142 = g(2) * t99 - g(3) * t100;
t171 = t135 ^ 2 - t35 ^ 2;
t165 = t73 * pkin(7);
t143 = t155 * qJD(1);
t59 = t115 * qJD(2) + t118 * t143;
t158 = t111 * t59;
t159 = qJD(3) * pkin(3);
t58 = qJD(2) * t118 - t115 * t143;
t53 = t58 + t159;
t24 = t109 * t53 + t158;
t11 = t24 + t165;
t112 = cos(pkin(8));
t160 = t112 * pkin(1);
t98 = pkin(3) * t118 + pkin(2);
t133 = -t98 - t160;
t71 = qJD(1) * t133 + qJD(4);
t44 = -t73 * pkin(4) + t71;
t102 = qJ(3) + pkin(9) + qJ(5);
t93 = sin(t102);
t94 = cos(t102);
t170 = g(1) * t93 + t11 * t152 + t142 * t94 - t44 * t35;
t168 = qJD(5) - t105;
t101 = t118 * qJDD(2);
t86 = t95 * qJDD(1);
t126 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t86;
t134 = t143 * qJD(3);
t21 = qJDD(3) * pkin(3) - t115 * t126 - t118 * t134 + t101;
t25 = (qJDD(2) - t134) * t115 + t126 * t118;
t4 = -t109 * t25 + t111 * t21;
t2 = qJDD(3) * pkin(4) - pkin(7) * t43 + t4;
t5 = t109 * t21 + t111 * t25;
t3 = pkin(7) * t42 + t5;
t167 = -g(1) * t94 - t114 * t3 + t117 * t2 - t44 * t135 + t142 * t93;
t125 = pkin(3) * t147 + qJDD(1) * t133 + qJDD(4);
t164 = t75 * pkin(7);
t163 = pkin(3) * t109;
t162 = g(1) * t118;
t49 = t109 * t59;
t27 = t111 * t58 - t49;
t144 = qJD(3) * t155;
t62 = t118 * qJD(4) - t115 * t144;
t63 = -t115 * qJD(4) - t118 * t144;
t29 = t109 * t63 + t111 * t62;
t78 = t155 * t115;
t79 = t155 * t118;
t41 = -t109 * t78 + t111 * t79;
t97 = -pkin(2) - t160;
t89 = qJD(1) * t97;
t154 = qJDD(2) - g(1);
t107 = t115 ^ 2;
t153 = -t118 ^ 2 + t107;
t149 = t118 * qJDD(1);
t148 = t115 * t159;
t23 = t111 * t53 - t49;
t26 = -t109 * t58 - t158;
t28 = -t109 * t62 + t111 * t63;
t40 = -t109 * t79 - t111 * t78;
t141 = g(2) * t100 + g(3) * t99;
t116 = sin(qJ(1));
t119 = cos(qJ(1));
t140 = -g(2) * t119 - g(3) * t116;
t10 = qJD(3) * pkin(4) - t164 + t23;
t139 = -t114 * t10 - t117 * t11;
t104 = qJDD(3) + qJDD(5);
t45 = t114 * t81 - t117 * t80;
t77 = t80 * qJD(3);
t12 = -qJD(5) * t45 - t114 * t74 + t117 * t77;
t46 = t114 * t80 + t117 * t81;
t138 = t104 * t46 + t105 * t12;
t30 = -pkin(7) * t81 + t40;
t31 = pkin(7) * t80 + t41;
t137 = -t114 * t31 + t117 * t30;
t136 = t114 * t30 + t117 * t31;
t96 = pkin(3) * t111 + pkin(4);
t131 = t114 * t96 + t117 * t163;
t130 = -t114 * t163 + t117 * t96;
t129 = -qJD(1) * t89 + t142 - t86;
t128 = 0.2e1 * qJD(3) * t89 - qJDD(3) * t95;
t120 = qJD(3) ^ 2;
t123 = 0.2e1 * qJDD(1) * t97 + t120 * t95 + t141;
t121 = qJD(1) ^ 2;
t113 = -qJ(4) - pkin(6);
t85 = qJDD(3) * t118 - t115 * t120;
t84 = qJDD(3) * t115 + t118 * t120;
t61 = pkin(4) * t74 + t148;
t60 = pkin(3) * qJD(1) * t115 + pkin(4) * t75;
t57 = -pkin(4) * t80 + t133;
t22 = -t42 * pkin(4) + t125;
t17 = -pkin(7) * t74 + t29;
t16 = -pkin(7) * t77 + t28;
t15 = t27 - t164;
t14 = t26 - t165;
t13 = qJD(5) * t46 + t114 * t77 + t117 * t74;
t6 = -t104 * t45 - t105 * t13;
t1 = [qJDD(1), t140, g(2) * t116 - g(3) * t119, (t140 + (t110 ^ 2 + t112 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t107 + 0.2e1 * t115 * t146, 0.2e1 * t115 * t149 - 0.2e1 * t150 * t153, t84, t85, 0, t115 * t128 - t118 * t123, t115 * t123 + t118 * t128, -t23 * t77 - t24 * t74 - t28 * t75 + t29 * t73 - t4 * t81 - t40 * t43 + t41 * t42 + t5 * t80 - t142, t5 * t41 + t24 * t29 + t4 * t40 + t23 * t28 + t71 * t148 - g(2) * (pkin(1) * t119 + t100 * t98 - t113 * t99) - g(3) * (pkin(1) * t116 + t100 * t113 + t98 * t99) + t125 * t133, t12 * t135 + t46 * t7, t12 * t35 - t13 * t135 - t45 * t7 - t46 * t8, t138, t6, 0, -t61 * t35 + t57 * t8 + t22 * t45 + t44 * t13 + (-qJD(5) * t136 - t114 * t17 + t117 * t16) * t105 + t137 * t104 - t141 * t94, t61 * t135 + t57 * t7 + t22 * t46 + t44 * t12 - (qJD(5) * t137 + t114 * t16 + t117 * t17) * t105 - t136 * t104 + t141 * t93; 0, 0, 0, t154, 0, 0, 0, 0, 0, t85, -t84, t42 * t81 - t43 * t80 + t73 * t77 + t74 * t75, -t23 * t74 + t24 * t77 + t4 * t80 + t5 * t81 - g(1), 0, 0, 0, 0, 0, t6, -t138; 0, 0, 0, 0, -t115 * t121 * t118, t153 * t121, t115 * qJDD(1), t149, qJDD(3), t115 * t129 + t101 - t162, -t115 * t154 + t118 * t129, (t24 + t26) * t75 + (t23 - t27) * t73 + (t109 * t42 - t111 * t43) * pkin(3), -t23 * t26 - t24 * t27 + (-t162 + t109 * t5 + t111 * t4 + (-qJD(1) * t71 + t142) * t115) * pkin(3), -t173, t171, t174, t172, t104, t130 * t104 + t60 * t35 - (-t114 * t15 + t117 * t14) * t105 + (-t105 * t131 + t139) * qJD(5) + t167, -t131 * t104 - t117 * t3 - t114 * t2 - t60 * t135 + (t114 * t14 + t117 * t15) * t105 + (-t117 * t10 - t105 * t130) * qJD(5) + t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 ^ 2 - t75 ^ 2, t23 * t75 - t24 * t73 + t125 + t141, 0, 0, 0, 0, 0, t8 + t157, t7 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t171, t174, t172, t104, t139 * t168 + t167, (-t105 * t11 - t2) * t114 + (-t10 * t168 - t3) * t117 + t170;];
tau_reg = t1;
