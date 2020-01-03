% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPPR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:47:53
% EndTime: 2019-12-31 17:47:57
% DurationCPUTime: 1.36s
% Computational Cost: add. (845->218), mult. (1920->315), div. (0->0), fcn. (1371->8), ass. (0->134)
t140 = (qJ(2) * qJDD(1));
t172 = pkin(3) + qJ(2);
t139 = (qJD(1) * qJD(2));
t171 = (2 * t139 + t140) * qJ(2);
t170 = t139 + t140;
t96 = sin(qJ(1));
t168 = g(1) * t96;
t98 = cos(qJ(1));
t86 = g(2) * t98;
t130 = -t86 + t168;
t94 = cos(pkin(7));
t90 = t94 ^ 2;
t127 = t90 * t139;
t141 = t94 * qJDD(1);
t30 = pkin(3) * t141 + t170 * t94 + qJDD(4);
t54 = t172 * t94;
t169 = (qJDD(1) * t54 + t30) * t94 + t127;
t148 = qJD(1) * t94;
t93 = cos(pkin(8));
t52 = t93 * t148 + qJD(5);
t97 = cos(qJ(5));
t167 = t52 * t97;
t91 = sin(pkin(8));
t92 = sin(pkin(7));
t166 = t91 * t92;
t99 = qJD(1) ^ 2;
t165 = t91 * t99;
t164 = t92 * t93;
t95 = sin(qJ(5));
t163 = t92 * t95;
t162 = t92 * t97;
t161 = t92 * t99;
t160 = t93 * t94;
t159 = t93 * t99;
t158 = t94 * t95;
t157 = t94 * t96;
t156 = t94 * t97;
t155 = t94 * t98;
t150 = t92 * qJ(3);
t125 = -pkin(1) - t150;
t42 = (-pkin(2) - qJ(4)) * t94 + t125;
t50 = -t92 * qJD(3) - t94 * qJD(4);
t16 = t50 * qJD(1) + t42 * qJDD(1) + qJDD(2);
t43 = t170 * t92 + qJDD(3);
t79 = t92 * qJDD(1);
t29 = pkin(3) * t79 + t43;
t4 = t93 * t16 + t91 * t29;
t26 = t42 * qJD(1) + qJD(2);
t149 = qJD(1) * t92;
t143 = qJ(2) * qJD(1);
t58 = t92 * t143 + qJD(3);
t44 = pkin(3) * t149 + t58;
t8 = t93 * t26 + t91 * t44;
t53 = t172 * t92;
t18 = t93 * t42 + t91 * t53;
t154 = t171 * t90;
t153 = t98 * pkin(1) + t96 * qJ(2);
t152 = -t91 ^ 2 - t93 ^ 2;
t88 = t92 ^ 2;
t151 = t88 + t90;
t147 = qJD(2) * t92;
t146 = qJD(2) * t94;
t145 = qJD(5) * t52;
t144 = qJDD(1) * pkin(1);
t48 = -t94 * pkin(2) + t125;
t142 = qJDD(1) * t48;
t138 = qJD(1) * qJD(3);
t136 = t91 * t156;
t135 = t91 * t148;
t134 = t95 * t149;
t45 = pkin(3) * t148 + t94 * t143 + qJD(4);
t133 = t91 * t79;
t132 = t93 * t79;
t83 = t98 * qJ(2);
t131 = -t96 * pkin(1) + t83;
t129 = t152 * t99;
t49 = t151 * t99;
t128 = qJD(5) * t149;
t126 = t92 * t138;
t124 = pkin(2) * t155 + t98 * t150 + t153;
t78 = qJDD(2) - t144;
t121 = qJD(5) * t135;
t120 = g(1) * t98 + g(2) * t96;
t111 = (pkin(4) * t93 + pkin(6) * t91) * t94;
t13 = qJDD(1) * t111 + t30;
t2 = pkin(6) * t79 + t4;
t119 = t95 * t13 + t97 * t2;
t19 = qJD(1) * t111 + t45;
t6 = pkin(6) * t149 + t8;
t118 = t97 * t19 - t95 * t6;
t117 = -t95 * t19 - t97 * t6;
t15 = t92 * pkin(6) + t18;
t23 = t111 + t54;
t116 = t97 * t15 + t95 * t23;
t115 = -t95 * t15 + t97 * t23;
t3 = -t91 * t16 + t93 * t29;
t7 = -t91 * t26 + t93 * t44;
t17 = -t91 * t42 + t93 * t53;
t41 = -t96 * t166 + t98 * t93;
t114 = t96 * t156 + t41 * t95;
t113 = t95 * t157 - t41 * t97;
t39 = t136 - t163;
t38 = t91 * t158 + t162;
t112 = t144 - t78 - t86;
t51 = t93 * t141 + qJDD(5);
t110 = -t97 * t145 - t95 * t51;
t109 = t95 * t145 - t97 * t51;
t103 = qJDD(2) + t142;
t25 = t103 - t126;
t108 = -t142 - t25 - t86;
t27 = -t93 * t147 + t91 * t50;
t107 = -qJD(1) * t27 + qJDD(1) * t17 + t3;
t28 = t91 * t147 + t93 * t50;
t106 = -t28 * qJD(1) - t18 * qJDD(1) - t4;
t105 = 0.2e1 * t90 * t140 - t120 + 0.2e1 * t127;
t11 = -qJDD(1) * t136 + t97 * t128 + (t121 + t79) * t95;
t104 = t170 * t88;
t102 = t117 * qJD(5) + t97 * t13 - t95 * t2;
t101 = (-t7 * t91 + t8 * t93) * qJD(1) - t120;
t12 = t38 * qJDD(1) + t97 * t121 - t95 * t128;
t84 = g(3) * t94;
t67 = g(1) * t157;
t40 = t98 * t166 + t96 * t93;
t37 = t48 * qJD(1) + qJD(2);
t35 = t97 * t135 - t134;
t34 = t39 * qJD(5);
t33 = t38 * qJD(5);
t32 = t38 * qJD(1);
t21 = t95 * t155 + t40 * t97;
t20 = t97 * t155 - t40 * t95;
t14 = -t92 * pkin(4) - t17;
t5 = -pkin(4) * t149 - t7;
t1 = -pkin(4) * t79 - t3;
t9 = [qJDD(1), t130, t120, t112 * t94 + t67, (-t112 - t168) * t92, 0.2e1 * t104 + t105, -t78 * pkin(1) - g(1) * t131 - g(2) * t153 + t171 * t88 + t154, t43 * t92 + t104 + t105, -t67 + (-t108 - t126) * t94, t88 * t138 + (t108 + t168) * t92, t25 * t48 - g(1) * (-pkin(2) * t157 + t131) - g(2) * t124 + (t43 * qJ(2) + qJ(3) * t168 + t58 * qJD(2) - t37 * qJD(3)) * t92 + t154, -g(1) * t41 - g(2) * t40 + t107 * t92 + t169 * t93, (t130 * t93 + t106) * t92 + (t120 - t169) * t91, t67 + (t106 * t93 + t107 * t91 - t86) * t94, t4 * t18 + t8 * t28 + t3 * t17 - t7 * t27 + t30 * t54 + t45 * t146 - g(1) * (t98 * pkin(3) + t83) - g(2) * (qJ(4) * t155 + t124) + (-g(2) * pkin(3) - g(1) * t42) * t96, -t11 * t39 - t35 * t33, t11 * t38 - t39 * t12 + t33 * t32 - t35 * t34, t11 * t160 + t33 * t52 - t39 * t51, t12 * t160 + t34 * t52 + t38 * t51, t51 * t160, (-qJD(5) * t116 + t146 * t97 - t95 * t28) * t52 + t115 * t51 + t102 * t160 - t27 * t32 - t14 * t12 - t1 * t38 - t5 * t34 + g(1) * t113 - g(2) * t21, -(qJD(5) * t115 + t146 * t95 + t97 * t28) * t52 - t116 * t51 - (qJD(5) * t118 + t119) * t160 - t27 * t35 + t14 * t11 - t1 * t39 + t5 * t33 + g(1) * t114 - g(2) * t20; 0, 0, 0, -t141, t79, -t49, -qJ(2) * t49 - t130 + t78, -t49, t141, -t79, -t90 * t99 * qJ(2) + (-qJD(3) - t58) * t149 + t103 - t130, -t93 * t49 - t133, t151 * t165 - t132, t152 * t141, -t3 * t91 + t4 * t93 + (-t45 * t94 + (-t7 * t93 - t8 * t91) * t92) * qJD(1) - t130, 0, 0, 0, 0, 0, -t91 * t12 + t110 * t93 + (-(-t91 * t163 + t156) * t52 - t32 * t164) * qJD(1), t91 * t11 + t109 * t93 + ((t91 * t162 + t158) * t52 - t35 * t164) * qJD(1); 0, 0, 0, 0, 0, 0, 0, t79, t94 * t161, -t88 * t99, t84 + (qJD(1) * t37 - t120) * t92 + t43, -t88 * t165 + t132, -t88 * t159 - t133, t94 * t92 * t129, t101 * t92 + t3 * t93 + t4 * t91 + t84, 0, 0, 0, 0, 0, (-t134 * t52 + t12) * t93 + (-t32 * t149 + t110) * t91, (-t149 * t167 - t11) * t93 + (-t35 * t149 + t109) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (qJDD(1) * t93 - t91 * t161) * t94, (-qJDD(1) * t91 - t92 * t159) * t94, t90 * t129, -g(3) * t92 + t101 * t94 + t30, 0, 0, 0, 0, 0, (-t52 * t93 * t95 - t32 * t91) * t148 - t109, (-t93 * t167 - t35 * t91) * t148 + t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t32, -t32 ^ 2 + t35 ^ 2, -t32 * t52 + t11, -t35 * t52 + t12, t51, -g(1) * t20 - g(2) * t114 - g(3) * t38 - t117 * t52 + t5 * t35 + t102, g(1) * t21 + g(2) * t113 - g(3) * t39 - t5 * t32 - t119 + (-qJD(5) + t52) * t118;];
tau_reg = t9;
