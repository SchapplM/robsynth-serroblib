% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRPRP5
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:44
% EndTime: 2019-12-05 15:38:49
% DurationCPUTime: 1.17s
% Computational Cost: add. (1138->208), mult. (2501->261), div. (0->0), fcn. (1847->10), ass. (0->124)
t130 = qJD(1) * qJD(2);
t86 = cos(qJ(2));
t131 = t86 * qJDD(1);
t85 = sin(qJ(2));
t110 = t85 * t130 + qJDD(3) - t131;
t141 = qJDD(2) * pkin(2);
t45 = t110 - t141;
t80 = sin(pkin(7));
t82 = cos(pkin(7));
t112 = g(1) * t82 + g(2) * t80;
t75 = g(3) * t86;
t98 = t112 * t85 - t75;
t167 = -t45 + t98;
t155 = cos(qJ(4));
t81 = cos(pkin(8));
t125 = t155 * t81;
t79 = sin(pkin(8));
t84 = sin(qJ(4));
t148 = t84 * t79;
t103 = t125 - t148;
t101 = t103 * t86;
t147 = pkin(6) + qJ(3);
t57 = t147 * t79;
t58 = t147 * t81;
t104 = -t155 * t57 - t84 * t58;
t146 = -qJD(1) * t101 + t103 * qJD(3) + t104 * qJD(4);
t27 = t155 * t58 - t84 * t57;
t78 = pkin(8) + qJ(4);
t73 = sin(t78);
t166 = t146 * qJD(4) + t27 * qJDD(4) + t98 * t73;
t137 = t86 * qJD(1);
t114 = qJD(3) - t137;
t164 = t112 * t86;
t54 = t155 * t79 + t84 * t81;
t162 = t54 * qJD(2);
t51 = t54 * qJD(4);
t161 = (t112 + t130) * t85 + t141 - t45 - t75;
t159 = t162 ^ 2;
t156 = g(3) * t85;
t113 = qJD(2) * t125;
t126 = qJD(2) * t148;
t46 = -t113 + t126;
t154 = t162 * t46;
t74 = cos(t78);
t153 = t80 * t74;
t152 = t80 * t86;
t151 = t82 * t74;
t150 = t82 * t86;
t138 = t85 * qJD(1);
t62 = qJD(2) * qJ(3) + t138;
t120 = pkin(6) * qJD(2) + t62;
t40 = t120 * t81;
t149 = t84 * t40;
t145 = t27 * qJD(4) + t114 * t54;
t144 = t79 ^ 2 + t81 ^ 2;
t143 = qJD(2) * t85;
t142 = qJD(4) * t84;
t140 = qJDD(4) * pkin(4);
t39 = t120 * t79;
t13 = t155 * t40 - t84 * t39;
t139 = t13 * qJD(4);
t12 = -t155 * t39 - t149;
t136 = qJD(5) - t12;
t135 = qJDD(1) - g(3);
t133 = t79 * qJDD(2);
t132 = t81 * qJDD(2);
t129 = qJDD(4) * qJ(5);
t121 = qJD(4) * t155;
t41 = qJDD(2) * qJ(3) + t85 * qJDD(1) + (qJD(3) + t137) * qJD(2);
t116 = pkin(6) * qJDD(2) + t41;
t24 = t116 * t79;
t25 = t116 * t81;
t128 = -t39 * t121 + t155 * t25 - t84 * t24;
t118 = qJDD(2) * t155;
t127 = qJD(4) * t113 + t79 * t118 + t84 * t132;
t71 = t81 * pkin(3) + pkin(2);
t124 = t144 * t41;
t123 = t144 * t86;
t122 = t12 + t149;
t50 = -t81 * t121 + t79 * t142;
t6 = t51 * pkin(4) + t50 * qJ(5) - t54 * qJD(5);
t119 = -t6 + t138;
t117 = t40 * t121 - t39 * t142 + t155 * t24 + t84 * t25;
t115 = t144 * qJDD(2);
t111 = -t81 * t118 + t84 * t133;
t87 = qJD(2) ^ 2;
t109 = t86 * qJDD(2) - t87 * t85;
t108 = pkin(4) * t74 + qJ(5) * t73 + t71;
t43 = t103 * t85;
t99 = -t156 - t164;
t15 = qJD(4) * t43 + t162 * t86;
t19 = qJD(2) * t51 + t111;
t42 = t54 * t85;
t97 = -t15 * qJD(4) - t42 * qJDD(4) + t46 * t143 - t86 * t19;
t35 = t73 * t152 + t151;
t37 = t73 * t150 - t153;
t96 = g(1) * t37 + g(2) * t35 + t73 * t156 - t117;
t95 = t114 * t144;
t52 = -t71 * qJD(2) + t114;
t36 = t74 * t152 - t82 * t73;
t38 = t74 * t150 + t80 * t73;
t94 = g(1) * t38 + g(2) * t36 + t74 * t156 - t128;
t31 = -t71 * qJDD(2) + t110;
t14 = qJD(2) * t101 - t85 * t51;
t18 = qJD(4) * t126 - t127;
t93 = t14 * qJD(4) + t43 * qJDD(4) - t143 * t162 - t86 * t18;
t92 = -t145 * qJD(4) + t104 * qJDD(4) - t74 * t75 + (g(1) * t151 + g(2) * t153) * t85;
t7 = t46 * pkin(4) - qJ(5) * t162 + t52;
t91 = t162 * t7 + qJDD(5) - t96;
t90 = t124 + t99;
t89 = t19 * pkin(4) + t18 * qJ(5) + t31;
t88 = 0.2e1 * t162 * qJD(4) + t111;
t59 = -qJD(2) * pkin(2) + t114;
t44 = t46 ^ 2;
t17 = -pkin(4) * t103 - t54 * qJ(5) - t71;
t16 = pkin(4) * t162 + t46 * qJ(5);
t9 = qJD(4) * qJ(5) + t13;
t8 = -qJD(4) * pkin(4) + t136;
t5 = (t46 - t126) * qJD(4) + t127;
t4 = (t46 + t126) * qJD(4) - t127;
t3 = -qJD(5) * t162 + t89;
t2 = qJDD(5) + t117 - t140;
t1 = t129 + (qJD(5) - t149) * qJD(4) + t128;
t10 = [t135, 0, t109, -qJDD(2) * t85 - t87 * t86, t109 * t81, -t109 * t79, t85 * t115 + t87 * t123, -t45 * t86 - g(3) + t85 * t124 + (t62 * t123 + t59 * t85) * qJD(2), 0, 0, 0, 0, 0, t97, -t93, t97, -t14 * t46 + t15 * t162 - t42 * t18 - t43 * t19, t93, t1 * t43 + t9 * t14 + t7 * t143 + t8 * t15 + t2 * t42 - t3 * t86 - g(3); 0, qJDD(2), t131 + t98, -t135 * t85 + t164, t161 * t81, -t161 * t79, qJ(3) * t115 + t95 * qJD(2) + t90, t167 * pkin(2) + t90 * qJ(3) - t59 * t138 + t95 * t62, -t162 * t50 - t18 * t54, -t103 * t18 - t162 * t51 - t54 * t19 + t50 * t46, -t50 * qJD(4) + t54 * qJDD(4), -t51 * qJD(4) + qJDD(4) * t103, 0, -t103 * t31 - t46 * t138 - t71 * t19 + t52 * t51 + t92, -t138 * t162 + t71 * t18 + t31 * t54 - t52 * t50 - t166, -t103 * t3 - t119 * t46 + t17 * t19 + t7 * t51 + t92, t1 * t103 + t104 * t18 + t145 * t162 - t146 * t46 - t27 * t19 + t2 * t54 - t8 * t50 - t9 * t51 + t99, t119 * t162 + t17 * t18 - t3 * t54 + t7 * t50 + t166, t1 * t27 + t3 * t17 - t2 * t104 + t7 * t6 + t146 * t9 + t145 * t8 + (-g(3) * t108 - t112 * t147) * t86 + (-g(3) * t147 - t7 * qJD(1) + t112 * t108) * t85; 0, 0, 0, 0, -t132, t133, -t144 * t87, -t144 * t62 * qJD(2) - t167, 0, 0, 0, 0, 0, t88, -t4, t88, -t44 - t159, t4, t9 * t46 + (-qJD(5) - t8) * t162 + t89 - t98; 0, 0, 0, 0, 0, 0, 0, 0, t154, -t44 + t159, t5, -t111, qJDD(4), -t162 * t52 + t139 + t96, t122 * qJD(4) + t52 * t46 + t94, -t16 * t46 + t139 + 0.2e1 * t140 - t91, pkin(4) * t18 - t19 * qJ(5) + (-t13 + t9) * t162 + (t8 - t136) * t46, 0.2e1 * t129 + t16 * t162 - t7 * t46 + (0.2e1 * qJD(5) - t122) * qJD(4) - t94, t1 * qJ(5) - t2 * pkin(4) - t7 * t16 - t8 * t13 - g(1) * (-t37 * pkin(4) + t38 * qJ(5)) - g(2) * (-t35 * pkin(4) + t36 * qJ(5)) + t136 * t9 - (-pkin(4) * t73 + qJ(5) * t74) * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t154, t5, -qJD(4) ^ 2 - t159, -t9 * qJD(4) - t140 + t91;];
tau_reg = t10;
