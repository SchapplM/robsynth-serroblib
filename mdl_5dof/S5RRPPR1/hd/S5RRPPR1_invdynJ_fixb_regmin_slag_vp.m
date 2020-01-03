% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:56:01
% EndTime: 2020-01-03 11:56:03
% DurationCPUTime: 0.74s
% Computational Cost: add. (1086->174), mult. (1749->229), div. (0->0), fcn. (1162->16), ass. (0->126)
t107 = qJ(1) + qJ(2);
t95 = pkin(8) + t107;
t81 = sin(t95);
t82 = cos(t95);
t134 = g(2) * t82 + g(3) * t81;
t102 = qJDD(1) + qJDD(2);
t109 = sin(pkin(8));
t111 = cos(pkin(8));
t116 = cos(qJ(2));
t162 = pkin(1) * qJD(2);
t147 = qJD(1) * t162;
t113 = sin(qJ(2));
t153 = qJDD(1) * t113;
t173 = pkin(1) * t153 + t116 * t147;
t165 = t116 * pkin(1);
t88 = qJDD(1) * t165;
t41 = pkin(2) * t102 - t113 * t147 + t88;
t21 = -t109 * t173 + t111 * t41;
t120 = qJDD(4) - t21;
t16 = -t102 * pkin(3) + t120;
t123 = -t134 - t16;
t110 = cos(pkin(9));
t115 = cos(qJ(5));
t156 = t115 * t110;
t108 = sin(pkin(9));
t112 = sin(qJ(5));
t159 = t112 * t108;
t52 = -t156 + t159;
t53 = t115 * t108 + t112 * t110;
t106 = qJD(1) + qJD(2);
t39 = t53 * t106;
t133 = g(2) * t81 - g(3) * t82;
t154 = t108 ^ 2 + t110 ^ 2;
t172 = t106 * t154;
t163 = pkin(1) * qJD(1);
t148 = t116 * t163;
t149 = t113 * t163;
t70 = t109 * t149;
t46 = t111 * t148 - t70;
t155 = qJD(4) - t46;
t168 = pkin(1) * t113;
t167 = t110 * pkin(4);
t166 = t111 * pkin(2);
t22 = t109 * t41 + t111 * t173;
t14 = qJ(4) * t102 + qJD(4) * t106 + t22;
t9 = t108 * qJDD(3) + t110 * t14;
t56 = pkin(2) * t106 + t148;
t71 = t111 * t149;
t34 = t109 * t56 + t71;
t160 = t111 * t113;
t87 = pkin(2) + t165;
t164 = pkin(1) * t160 + t109 * t87;
t161 = t110 * t102;
t152 = t123 * t108;
t144 = t106 * t156;
t151 = qJD(5) * t144 + t102 * t53;
t97 = cos(t107);
t86 = pkin(2) * t97;
t150 = t82 * pkin(3) + t81 * qJ(4) + t86;
t77 = t109 * t168;
t146 = t106 * t159;
t143 = -pkin(3) - t167;
t96 = sin(t107);
t142 = g(2) * t96 - g(3) * t97;
t10 = t102 * t143 + t120;
t33 = t111 * t56 - t70;
t130 = qJD(4) - t33;
t25 = t106 * t143 + t130;
t47 = t52 * qJD(5);
t105 = pkin(9) + qJ(5);
t93 = sin(t105);
t141 = t10 * t53 + t134 * t93 - t25 * t47;
t140 = t111 * t87 - t77;
t139 = t102 * t154;
t138 = qJD(1) * (-qJD(2) + t106);
t137 = qJD(2) * (-qJD(1) - t106);
t43 = -pkin(3) - t140;
t85 = pkin(2) * t96;
t135 = t81 * pkin(3) - qJ(4) * t82 + t85;
t132 = -g(2) * t97 - g(3) * t96;
t131 = t52 * t102;
t45 = (t109 * t116 + t160) * t162;
t129 = t102 * t43 + t106 * t45;
t44 = t109 * t148 + t71;
t83 = -pkin(3) - t166;
t128 = t102 * t83 - t106 * t44;
t42 = qJ(4) + t164;
t31 = (-pkin(7) - t42) * t108;
t98 = t110 * pkin(7);
t32 = t110 * t42 + t98;
t127 = t112 * t32 - t115 * t31;
t126 = t112 * t31 + t115 * t32;
t79 = pkin(2) * t109 + qJ(4);
t49 = (-pkin(7) - t79) * t108;
t50 = t110 * t79 + t98;
t125 = t112 * t50 - t115 * t49;
t124 = t112 * t49 + t115 * t50;
t122 = t132 + t88;
t121 = t111 * t116 * t162 - qJD(2) * t77;
t90 = t110 * qJDD(3);
t8 = -t108 * t14 + t90;
t119 = -t8 * t108 + t110 * t9 - t133;
t48 = t53 * qJD(5);
t94 = cos(t105);
t118 = t10 * t52 - t134 * t94 + t25 * t48;
t117 = cos(qJ(1));
t114 = sin(qJ(1));
t100 = t117 * pkin(1);
t99 = t114 * pkin(1);
t58 = t143 - t166;
t40 = qJD(4) + t121;
t37 = -t144 + t146;
t35 = t43 - t167;
t30 = qJ(4) * t106 + t34;
t29 = -pkin(3) * t106 + t130;
t27 = -qJD(5) * t48 - qJDD(5) * t52;
t26 = -qJD(5) * t47 + qJDD(5) * t53;
t24 = qJD(3) * t108 + t110 * t30;
t23 = qJD(3) * t110 - t108 * t30;
t20 = t106 * t48 + t131;
t19 = -qJD(5) * t146 + t151;
t4 = pkin(7) * t161 + t9;
t3 = t90 + (-pkin(7) * t102 - t14) * t108;
t2 = t19 * t53 - t39 * t47;
t1 = -t19 * t52 - t20 * t53 + t37 * t47 - t39 * t48;
t5 = [qJDD(1), -g(2) * t117 - g(3) * t114, g(2) * t114 - g(3) * t117, t102, (t102 * t116 + t113 * t137) * pkin(1) + t122, ((-qJDD(1) - t102) * t113 + t116 * t137) * pkin(1) + t142, t22 * t164 + t34 * t121 + t21 * t140 - t33 * t45 - g(2) * (t86 + t100) - g(3) * (t85 + t99), (t123 - t129) * t110, t108 * t129 - t152, t139 * t42 + t172 * t40 + t119, t16 * t43 + t29 * t45 - g(2) * (t100 + t150) - g(3) * (t135 + t99) + (t24 * t40 + t42 * t9) * t110 + (-t23 * t40 - t42 * t8) * t108, t2, t1, t26, t27, 0, t45 * t37 + t35 * t20 - t127 * qJDD(5) + (-qJD(5) * t126 - t40 * t53) * qJD(5) + t118, t45 * t39 + t35 * t19 - t126 * qJDD(5) + (qJD(5) * t127 + t40 * t52) * qJD(5) + t141; 0, 0, 0, t102, t138 * t168 + t122, (t116 * t138 - t153) * pkin(1) + t142, t33 * t44 - t34 * t46 + (t109 * t22 + t111 * t21 + t132) * pkin(2), (t123 - t128) * t110, t108 * t128 - t152, t79 * t139 + t155 * t172 + t119, t16 * t83 - t29 * t44 - g(2) * t150 - g(3) * t135 + (t155 * t24 + t9 * t79) * t110 + (-t155 * t23 - t8 * t79) * t108, t2, t1, t26, t27, 0, t58 * t20 - t125 * qJDD(5) - t44 * t37 + (-qJD(5) * t124 - t155 * t53) * qJD(5) + t118, t58 * t19 - t124 * qJDD(5) - t44 * t39 + (qJD(5) * t125 + t155 * t52) * qJD(5) + t141; 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, t108 * t9 + t110 * t8 - g(1), 0, 0, 0, 0, 0, t27, -t26; 0, 0, 0, 0, 0, 0, 0, -t161, t108 * t102, -t154 * t106 ^ 2, (t108 * t23 - t110 * t24) * t106 - t123, 0, 0, 0, 0, 0, 0.2e1 * t39 * qJD(5) + t131, (-t37 - t146) * qJD(5) + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, (t37 - t146) * qJD(5) + t151, -t131, qJDD(5), -g(1) * t94 - t112 * t4 + t115 * t3 + t133 * t93 - t25 * t39, g(1) * t93 - t112 * t3 - t115 * t4 + t133 * t94 + t25 * t37;];
tau_reg = t5;
