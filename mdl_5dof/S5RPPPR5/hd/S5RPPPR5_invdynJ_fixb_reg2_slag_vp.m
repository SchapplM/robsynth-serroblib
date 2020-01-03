% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:35
% EndTime: 2019-12-31 17:46:36
% DurationCPUTime: 0.82s
% Computational Cost: add. (1389->206), mult. (2469->260), div. (0->0), fcn. (1648->10), ass. (0->122)
t104 = cos(qJ(5));
t103 = sin(qJ(5));
t97 = sin(pkin(8));
t148 = t103 * t97;
t99 = cos(pkin(8));
t46 = -t104 * t99 + t148;
t98 = sin(pkin(7));
t37 = t46 * t98;
t48 = t103 * t99 + t104 * t97;
t161 = t48 * qJD(1);
t160 = qJDD(1) * pkin(3) + qJDD(4);
t159 = t161 ^ 2;
t105 = -pkin(1) - pkin(2);
t100 = cos(pkin(7));
t54 = t100 * qJ(2) + t98 * t105;
t50 = -qJ(4) + t54;
t158 = pkin(6) - t50;
t157 = cos(qJ(1));
t156 = sin(qJ(1));
t144 = qJD(1) * t99;
t130 = t104 * t144;
t131 = qJD(1) * t148;
t42 = -t130 + t131;
t155 = t42 * t161;
t135 = t100 * qJDD(1);
t57 = t105 * qJDD(1) + qJDD(2);
t152 = qJ(2) * t135 + t98 * t57;
t138 = qJD(1) * qJD(2);
t68 = t100 * t138;
t26 = t68 + t152;
t20 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4) + t26;
t13 = t97 * qJDD(3) + t99 * t20;
t141 = qJD(1) * t100;
t58 = t105 * qJD(1) + qJD(2);
t39 = qJ(2) * t141 + t98 * t58;
t29 = -qJD(1) * qJ(4) + t39;
t22 = t97 * qJD(3) + t99 * t29;
t154 = qJD(5) * t37 + t100 * t161;
t36 = t48 * t98;
t153 = qJD(5) * t36 - t46 * t141;
t151 = t157 * pkin(1) + t156 * qJ(2);
t150 = g(1) * t156 - g(2) * t157;
t93 = t97 ^ 2;
t94 = t99 ^ 2;
t149 = t93 + t94;
t146 = pkin(1) * qJDD(1);
t145 = qJD(1) * t98;
t143 = qJD(2) * t98;
t106 = qJD(1) ^ 2;
t142 = t100 * t106;
t140 = t99 * qJDD(1);
t139 = qJ(2) * qJDD(1);
t137 = qJDD(1) * t103;
t136 = qJDD(1) * t104;
t134 = -qJD(5) * t130 - t97 * t136 - t99 * t137;
t133 = t157 * pkin(2) + t151;
t132 = 0.2e1 * t138;
t66 = t98 * t138;
t79 = t99 * qJDD(3);
t10 = t79 + (pkin(6) * qJDD(1) - t20) * t97;
t11 = -pkin(6) * t140 + t13;
t129 = t104 * t10 - t103 * t11;
t128 = t100 * t57 - t98 * t139;
t38 = -qJ(2) * t145 + t100 * t58;
t53 = -t98 * qJ(2) + t100 * t105;
t127 = qJDD(1) * t149;
t126 = qJDD(2) - t146;
t51 = pkin(3) - t53;
t25 = t128 - t66;
t125 = -t156 * pkin(1) + t157 * qJ(2);
t47 = -t157 * t100 - t156 * t98;
t49 = -t156 * t100 + t157 * t98;
t124 = g(1) * t49 - g(2) * t47;
t123 = g(1) * t47 + g(2) * t49;
t122 = -t99 * t136 + t97 * t137;
t12 = -t97 * t20 + t79;
t121 = -t12 * t97 + t13 * t99;
t18 = qJD(5) * t131 + t134;
t45 = t48 * qJD(5);
t120 = -t161 * t45 + t18 * t46;
t19 = qJD(1) * t45 + t122;
t44 = t46 * qJD(5);
t119 = -t48 * t19 + t44 * t42;
t81 = t99 * qJD(3);
t118 = (-t97 * t29 + t81) * t97 - t22 * t99;
t117 = t39 * t100 - t38 * t98;
t116 = t103 * t10 + t104 * t11;
t14 = t81 + (pkin(6) * qJD(1) - t29) * t97;
t15 = -pkin(6) * t144 + t22;
t3 = -t103 * t15 + t104 * t14;
t4 = t103 * t14 + t104 * t15;
t34 = t158 * t97;
t35 = t158 * t99;
t8 = t103 * t35 + t104 * t34;
t9 = t103 * t34 - t104 * t35;
t115 = t44 * qJD(5) - t48 * qJDD(5);
t114 = t45 * qJD(5) + t46 * qJDD(5);
t113 = t98 * t106 + t135;
t28 = qJD(1) * pkin(3) + qJD(4) - t38;
t23 = -t25 + t160;
t111 = g(1) * t157 + g(2) * t156;
t110 = -t124 - t128;
t109 = -t156 * pkin(2) + t125;
t108 = t110 + t66 + t160;
t107 = qJDD(1) * t51 - t124 + t23 + t66;
t102 = -pkin(6) - qJ(4);
t95 = pkin(8) + qJ(5);
t87 = t99 * pkin(4);
t83 = cos(t95);
t82 = sin(t95);
t76 = pkin(4) * t140;
t75 = t87 + pkin(3);
t64 = t100 * qJD(2) - qJD(4);
t41 = t42 ^ 2;
t40 = t51 + t87;
t24 = pkin(4) * t144 + t28;
t16 = t23 + t76;
t6 = -t9 * qJD(5) - t48 * t64;
t5 = t8 * qJD(5) - t46 * t64;
t2 = -t4 * qJD(5) + t129;
t1 = t3 * qJD(5) + t116;
t7 = [0, 0, 0, 0, 0, qJDD(1), t150, t111, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t146 + t150, 0, -t111 + t132 + 0.2e1 * t139, -t126 * pkin(1) - g(1) * t125 - g(2) * t151 + (t132 + t139) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -t53 * qJDD(1) + t110 + 0.2e1 * t66, t54 * qJDD(1) + t123 + t152 + 0.2e1 * t68, 0, -g(1) * t109 - g(2) * t133 + t117 * qJD(2) + t25 * t53 + t26 * t54, t93 * qJDD(1), 0.2e1 * t97 * t140, 0, t94 * qJDD(1), 0, 0, t107 * t99, -t107 * t97, -t149 * t64 * qJD(1) - t50 * t127 - t121 - t123, t23 * t51 + t28 * t143 - g(1) * (t49 * pkin(3) + t47 * qJ(4) + t109) - g(2) * (-t47 * pkin(3) + t49 * qJ(4) + t133) - t118 * t64 + t121 * t50, -t161 * t44 - t18 * t48, t119 + t120, t115, t19 * t46 + t42 * t45, t114, 0, t6 * qJD(5) + t8 * qJDD(5) - t124 * t83 - t42 * t143 - t16 * t46 - t40 * t19 - t24 * t45, -t5 * qJD(5) - t9 * qJDD(5) + t124 * t82 - t143 * t161 - t16 * t48 + t40 * t18 + t24 * t44, t1 * t46 + t161 * t6 - t8 * t18 + t9 * t19 + t2 * t48 - t3 * t44 + t4 * t45 + t5 * t42 - t123, t1 * t9 + t4 * t5 + t2 * t8 + t3 * t6 + t16 * t40 + t24 * t143 - g(1) * (-t47 * t102 + t49 * t75 + t109) - g(2) * (-t49 * t102 - t47 * t75 + t133); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t106, -t106 * qJ(2) + t126 - t150, 0, 0, 0, 0, 0, 0, -t113, t98 * qJDD(1) - t142, 0, -t117 * qJD(1) + t25 * t100 + t26 * t98 - t150, 0, 0, 0, 0, 0, 0, -t113 * t99, t113 * t97, -t98 * t127 + t149 * t142, -t23 * t100 + t121 * t98 + (t118 * t100 - t28 * t98) * qJD(1) - t150, 0, 0, 0, 0, 0, 0, t154 * qJD(5) - t36 * qJDD(5) + t100 * t19 + t42 * t145, t153 * qJD(5) + t37 * qJDD(5) - t100 * t18 + t145 * t161, -t153 * t42 + t154 * t161 + t36 * t18 - t37 * t19, -t1 * t37 - t16 * t100 - t24 * t145 - t153 * t4 + t154 * t3 - t2 * t36 - t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t99 + t13 * t97 + g(3), 0, 0, 0, 0, 0, 0, -t114, t115, -t119 + t120, t1 * t48 - t2 * t46 - t3 * t45 - t4 * t44 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t140, -t97 * qJDD(1), -t149 * t106, -t118 * qJD(1) + t108, 0, 0, 0, 0, 0, 0, -0.2e1 * t161 * qJD(5) - t122, (t42 + t131) * qJD(5) + t134, -t41 - t159, -t161 * t3 - t4 * t42 + t108 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t41 + t159, (-t42 + t131) * qJD(5) + t134, -t155, t122, qJDD(5), g(3) * t83 - t123 * t82 + t161 * t24 + t129, -g(3) * t82 - t123 * t83 - t24 * t42 - t116, 0, 0;];
tau_reg = t7;
