% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:45
% EndTime: 2019-12-31 17:53:48
% DurationCPUTime: 1.07s
% Computational Cost: add. (953->223), mult. (2173->265), div. (0->0), fcn. (1547->6), ass. (0->117)
t136 = (qJ(2) * qJDD(1));
t135 = (qJD(1) * qJD(2));
t128 = 2 * t135;
t160 = (t128 + t136) * qJ(2);
t159 = t135 + t136;
t94 = sin(pkin(7));
t81 = t94 * qJ(3);
t95 = cos(pkin(7));
t50 = -t95 * pkin(2) - pkin(1) - t81;
t44 = t95 * pkin(3) - t50;
t100 = cos(qJ(1));
t98 = sin(qJ(1));
t147 = g(1) * t100 + g(2) * t98;
t97 = sin(qJ(4));
t99 = cos(qJ(4));
t47 = t94 * t97 + t95 * t99;
t40 = t47 * qJD(1);
t137 = t95 * qJDD(1);
t134 = qJD(1) * qJD(3);
t69 = t94 * t134;
t80 = t94 * qJDD(1);
t158 = -pkin(2) * t137 - qJ(3) * t80 - t69;
t144 = qJD(1) * t95;
t130 = t97 * t144;
t112 = qJD(4) * t130 - t47 * qJDD(1);
t157 = t40 ^ 2;
t145 = qJD(1) * t94;
t129 = t99 * t145;
t42 = t129 - t130;
t156 = t42 ^ 2;
t155 = g(1) * t98;
t89 = g(2) * t100;
t154 = t42 * t40;
t153 = t94 * t99;
t152 = t95 * t98;
t150 = -pkin(6) + qJ(2);
t53 = t150 * t95;
t49 = qJD(1) * t53;
t151 = t97 * t49;
t92 = t95 ^ 2;
t149 = t160 * t92;
t148 = t100 * pkin(1) + t98 * qJ(2);
t146 = t100 * t95;
t143 = qJD(4) * t97;
t142 = qJD(4) * t99;
t93 = qJDD(1) * pkin(1);
t141 = qJDD(4) * pkin(4);
t57 = qJ(2) * t145 + qJD(3);
t46 = -pkin(6) * t145 + t57;
t17 = t97 * t46 + t99 * t49;
t140 = t17 * qJD(4);
t139 = t94 * qJD(3);
t16 = t99 * t46 - t151;
t138 = qJD(5) - t16;
t79 = qJDD(2) - t93;
t132 = qJDD(4) * qJ(5);
t45 = t159 * t94 + qJDD(3);
t28 = -pkin(6) * t80 + t45;
t29 = (t150 * qJDD(1) + t135) * t95;
t131 = t46 * t142 + t97 * t28 + t99 * t29;
t85 = t100 * qJ(2);
t127 = -t98 * pkin(1) + t85;
t126 = -t89 + t155;
t125 = t16 + t151;
t102 = qJD(1) ^ 2;
t91 = t94 ^ 2;
t51 = (-t91 - t92) * t102;
t124 = t49 * t142 + t46 * t143 - t99 * t28 + t97 * t29;
t123 = -t147 + (t128 + 2 * t136) * t92;
t122 = pkin(2) * t146 + t100 * t81 + t148;
t37 = -qJD(1) * pkin(1) - pkin(2) * t144 - qJ(3) * t145 + qJD(2);
t119 = t97 * t137 - t99 * t80;
t25 = pkin(3) * t144 - t37;
t48 = -t95 * t97 + t153;
t118 = -t47 * pkin(4) + t48 * qJ(5);
t52 = t150 * t94;
t117 = t99 * t52 - t97 * t53;
t21 = t97 * t52 + t99 * t53;
t24 = t79 + t158;
t116 = -t126 + t79;
t19 = pkin(3) * t137 - t24;
t115 = -t79 + t93 - t89;
t114 = -qJDD(1) * t50 - t24 - t89;
t38 = t47 * qJD(4);
t113 = t159 * t91;
t101 = qJD(4) ^ 2;
t111 = qJDD(4) * t97 + t101 * t99 + t42 * t145;
t33 = t47 * t98;
t35 = t47 * t100;
t110 = g(1) * t35 + g(2) * t33 + g(3) * t48 - t131;
t9 = -t48 * qJD(2) + t21 * qJD(4);
t109 = g(1) * t33 - g(2) * t35 - t9 * qJD(4) + qJDD(4) * t117;
t32 = t97 * t152 - t98 * t153;
t34 = -t100 * t153 + t97 * t146;
t8 = t47 * qJD(2) + t117 * qJD(4);
t108 = g(1) * t32 - g(2) * t34 + t8 * qJD(4) + t21 * qJDD(4);
t107 = g(1) * t34 + g(2) * t32 + g(3) * t47 - t124;
t14 = qJD(1) * t38 + t119;
t15 = qJD(4) * t129 - t112;
t106 = -t15 * pkin(4) - t14 * qJ(5) - t19;
t4 = t40 * pkin(4) - t42 * qJ(5) + t25;
t105 = t4 * t42 + qJDD(5) - t107;
t104 = (-t42 - t129) * qJD(4) + t112;
t86 = g(3) * t95;
t66 = g(1) * t152;
t39 = t94 * t142 - t95 * t143;
t18 = qJDD(4) * t99 - t101 * t97 - t40 * t145;
t13 = t42 * pkin(4) + t40 * qJ(5);
t12 = qJD(4) * qJ(5) + t17;
t11 = -qJD(4) * pkin(4) + t138;
t10 = -t118 + t44;
t6 = 0.2e1 * t40 * qJD(4) + t119;
t5 = t39 * pkin(4) + t38 * qJ(5) - t48 * qJD(5) + t139;
t3 = qJDD(5) + t124 - t141;
t2 = t132 + (qJD(5) - t151) * qJD(4) + t131;
t1 = -t42 * qJD(5) - t106;
t7 = [qJDD(1), t126, t147, t115 * t95 + t66, (-t115 - t155) * t94, 0.2e1 * t113 + t123, -t79 * pkin(1) - g(1) * t127 - g(2) * t148 + t160 * t91 + t149, t66 + (t114 + t69) * t95, t45 * t94 + t113 + t123, t91 * t134 + (t114 + t155) * t94, t24 * t50 - g(1) * (-pkin(2) * t152 + t127) - g(2) * t122 + (t45 * qJ(2) + qJ(3) * t155 + t57 * qJD(2) - t37 * qJD(3)) * t94 + t149, -t14 * t48 - t42 * t38, t14 * t47 - t48 * t15 + t38 * t40 - t42 * t39, -t38 * qJD(4) + t48 * qJDD(4), -t39 * qJD(4) - t47 * qJDD(4), 0, t40 * t139 + t44 * t15 + t19 * t47 + t25 * t39 + t109, t42 * t139 - t44 * t14 + t19 * t48 - t25 * t38 - t108, t1 * t47 + t10 * t15 + t4 * t39 + t5 * t40 + t109, -t11 * t38 + t117 * t14 - t12 * t39 - t21 * t15 - t2 * t47 + t3 * t48 - t8 * t40 + t9 * t42 + t147, -t1 * t48 + t10 * t14 + t4 * t38 - t5 * t42 + t108, t2 * t21 + t12 * t8 + t1 * t10 + t4 * t5 - t3 * t117 + t11 * t9 - g(1) * (-t33 * pkin(4) - t100 * pkin(6) - t32 * qJ(5) + t85) - g(2) * (pkin(3) * t146 + t35 * pkin(4) + t34 * qJ(5) + t122) + (g(2) * pkin(6) + g(1) * t44) * t98; 0, 0, 0, -t137, t80, t51, qJ(2) * t51 + t116, -t137, t51, -t80, -t92 * t102 * qJ(2) - t57 * t145 + t116 + t158, 0, 0, 0, 0, 0, t104, t6, t104, t156 + t157, -t6, -t12 * t40 + (qJD(5) + t11) * t42 + t106 - t126; 0, 0, 0, 0, 0, 0, 0, -t94 * t102 * t95, t80, -t91 * t102, t86 + (qJD(1) * t37 - t147) * t94 + t45, 0, 0, 0, 0, 0, t18, -t111, t18, t99 * t14 - t97 * t15 + (-t40 * t99 + t42 * t97) * qJD(4), t111, t2 * t97 - t3 * t99 + t86 + (t11 * t97 + t12 * t99) * qJD(4) + (-qJD(1) * t4 - t147) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, t156 - t157, -t119, (t42 - t129) * qJD(4) + t112, qJDD(4), -t25 * t42 + t107 + t140, t125 * qJD(4) + t25 * t40 + t110, -t13 * t40 - t105 + t140 + 0.2e1 * t141, pkin(4) * t14 - t15 * qJ(5) + (t12 - t17) * t42 + (t11 - t138) * t40, 0.2e1 * t132 + t13 * t42 - t4 * t40 + (0.2e1 * qJD(5) - t125) * qJD(4) - t110, t2 * qJ(5) - t3 * pkin(4) - t4 * t13 - t11 * t17 - g(1) * (-t34 * pkin(4) + t35 * qJ(5)) - g(2) * (-t32 * pkin(4) + t33 * qJ(5)) - g(3) * t118 + t138 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) + t154, -t119, -t101 - t156, -t12 * qJD(4) + t105 - t141;];
tau_reg = t7;
