% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPR2
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:38
% EndTime: 2020-01-03 11:57:41
% DurationCPUTime: 0.98s
% Computational Cost: add. (1213->197), mult. (1950->280), div. (0->0), fcn. (1204->14), ass. (0->137)
t139 = pkin(1) * qJD(2);
t125 = qJD(1) * t139;
t92 = sin(qJ(2));
t134 = qJDD(1) * t92;
t95 = cos(qJ(2));
t169 = pkin(1) * t134 + t95 * t125;
t160 = t95 * pkin(1);
t71 = qJDD(1) * t160;
t80 = qJDD(1) + qJDD(2);
t34 = t80 * pkin(2) - t92 * t125 + t71;
t88 = sin(pkin(8));
t90 = cos(pkin(8));
t17 = -t169 * t88 + t90 * t34;
t104 = qJDD(4) - t17;
t14 = -t80 * pkin(3) + t104;
t86 = qJ(1) + qJ(2);
t74 = pkin(8) + t86;
t64 = sin(t74);
t65 = cos(t74);
t168 = g(2) * t65 + g(3) * t64 + t14;
t18 = t169 * t90 + t88 * t34;
t83 = qJD(1) + qJD(2);
t12 = t80 * qJ(4) + t83 * qJD(4) + t18;
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t8 = -t89 * qJDD(3) + t87 * t12;
t162 = t8 * t87;
t9 = t87 * qJDD(3) + t89 * t12;
t167 = t9 * t89 + t162;
t140 = pkin(1) * qJD(1);
t129 = t95 * t140;
t130 = t92 * t140;
t54 = t88 * t130;
t39 = t90 * t129 - t54;
t136 = qJD(4) - t39;
t150 = t89 * t83;
t51 = -qJD(5) + t150;
t135 = qJD(5) + t51;
t43 = t83 * pkin(2) + t129;
t55 = t90 * t130;
t27 = t88 * t43 + t55;
t23 = t83 * qJ(4) + t27;
t19 = -t89 * qJD(3) + t87 * t23;
t159 = t19 * t87;
t20 = t87 * qJD(3) + t89 * t23;
t166 = -t135 * t20 - t83 * t159;
t165 = pkin(1) * t92;
t161 = t90 * pkin(2);
t61 = t88 * t165;
t108 = t90 * t95 * t139 - qJD(2) * t61;
t33 = qJD(4) + t108;
t158 = t33 * t83;
t91 = sin(qJ(5));
t157 = t33 * t91;
t151 = t89 * t80;
t48 = -qJDD(5) + t151;
t156 = t48 * t89;
t155 = t80 * t91;
t94 = cos(qJ(5));
t154 = t80 * t94;
t79 = t83 ^ 2;
t81 = t87 ^ 2;
t153 = t81 * t79;
t152 = t87 * t91;
t149 = t89 * t91;
t148 = t89 * t94;
t147 = t90 * t92;
t146 = t91 * t48;
t145 = t91 * t94;
t144 = t94 * t48;
t70 = pkin(2) + t160;
t143 = pkin(1) * t147 + t88 * t70;
t142 = t89 ^ 2 + t81;
t85 = t94 ^ 2;
t141 = t91 ^ 2 - t85;
t138 = qJD(5) * t91;
t137 = qJD(5) * t94;
t132 = t168 * t87;
t76 = cos(t86);
t69 = pkin(2) * t76;
t131 = t65 * pkin(3) + t64 * qJ(4) + t69;
t128 = t83 * t137;
t127 = t51 * t138;
t75 = sin(t86);
t124 = g(2) * t75 - g(3) * t76;
t123 = t142 * t80;
t26 = t90 * t43 - t54;
t122 = t48 - t151;
t121 = t48 + t151;
t120 = t90 * t70 - t61;
t119 = t83 * t135;
t118 = t136 * t94;
t117 = qJD(1) * (-qJD(2) + t83);
t116 = qJD(2) * (-qJD(1) - t83);
t68 = pkin(2) * t75;
t114 = t64 * pkin(3) - t65 * qJ(4) + t68;
t113 = -g(2) * t76 - g(3) * t75;
t112 = qJD(4) - t26;
t111 = t20 * t89 + t159;
t36 = -pkin(3) - t120;
t38 = (t88 * t95 + t147) * t139;
t110 = t36 * t80 + t38 * t83;
t37 = t88 * t129 + t55;
t66 = -pkin(3) - t161;
t109 = -t37 * t83 + t66 * t80;
t107 = -t89 * pkin(4) - t87 * pkin(7) - pkin(3);
t105 = t113 + t71;
t40 = t107 - t161;
t63 = t88 * pkin(2) + qJ(4);
t103 = -t63 * t149 + t94 * t40;
t15 = t107 * t83 + t112;
t28 = -t64 * t149 - t65 * t94;
t30 = t65 * t149 - t64 * t94;
t7 = t107 * t80 + t104;
t102 = g(2) * t30 - g(3) * t28 + (t91 * t7 + t94 * t9 + (t94 * t15 - t91 * t20) * qJD(5)) * t89 + t94 * t162;
t29 = t64 * t148 - t65 * t91;
t31 = t65 * t148 + t64 * t91;
t101 = -g(2) * t31 - g(3) * t29 + t137 * t159 + t8 * t152;
t100 = -g(2) * t64 + g(3) * t65 + t167;
t99 = g(1) * t87 - t135 * t15 - t9;
t98 = t136 * t91 + t63 * t137;
t97 = -t51 ^ 2 - t153;
t96 = cos(qJ(1));
t93 = sin(qJ(1));
t78 = t96 * pkin(1);
t77 = t93 * pkin(1);
t42 = t87 * t138 * t150;
t35 = qJ(4) + t143;
t25 = (-0.2e1 * t91 * t128 + t80 * t85) * t81;
t24 = t107 - t120;
t22 = -t83 * pkin(3) + t112;
t21 = 0.2e1 * (t141 * t83 * qJD(5) - t80 * t145) * t81;
t11 = (t121 * t91 + (t51 + t150) * t137) * t87;
t10 = t42 + (-t121 * t94 + t127) * t87;
t5 = t94 * t7;
t2 = -t91 * t9 + t5 + (-t91 * t15 - t94 * t20) * qJD(5);
t1 = [qJDD(1), -g(2) * t96 - g(3) * t93, g(2) * t93 - g(3) * t96, t80, (t92 * t116 + t80 * t95) * pkin(1) + t105, ((-qJDD(1) - t80) * t92 + t95 * t116) * pkin(1) + t124, t18 * t143 + t27 * t108 + t17 * t120 - t26 * t38 - g(2) * (t69 + t78) - g(3) * (t68 + t77), (-t168 - t110) * t89, t110 * t87 + t132, t35 * t123 + t142 * t158 + t100, t14 * t36 + t22 * t38 - g(2) * (t78 + t131) - g(3) * (t114 + t77) + t167 * t35 + t111 * t33, t25, t21, t10, t11, t156, -(-t24 * t138 + t94 * t38) * t51 - t24 * t144 + (-(-t35 * t137 - t157) * t51 + t35 * t146 - t2) * t89 + (t83 * t157 + (t128 + t155) * t35) * t81 + t101, (t33 * t148 + t91 * t38) * t51 + (t35 * t148 + t91 * t24) * t48 + (t35 * t80 + t158) * t94 * t81 + (t94 * t24 * t51 + (-t159 + (-t51 * t89 - t81 * t83) * t35) * t91) * qJD(5) + t102; 0, 0, 0, t80, t117 * t165 + t105, (t95 * t117 - t134) * pkin(1) + t124, t26 * t37 - t27 * t39 + (t17 * t90 + t18 * t88 + t113) * pkin(2), (-t168 - t109) * t89, t109 * t87 + t132, t136 * t83 * t142 + t63 * t123 + t100, t14 * t66 - t22 * t37 - g(2) * t131 - g(3) * t114 + (t136 * t20 + t9 * t63) * t89 + (t136 * t19 + t8 * t63) * t87, t25, t21, t10, t11, t156, -t103 * t48 - t2 * t89 + (t40 * t138 + t94 * t37 + t98 * t89) * t51 + (t63 * t155 + t98 * t83) * t81 + t101, (t63 * t148 + t91 * t40) * t48 + (t89 * t118 - t91 * t37) * t51 + (t103 * t51 - t19 * t152) * qJD(5) + (t63 * t154 + (-t63 * t138 + t118) * t83) * t81 + t102; 0, 0, 0, 0, 0, 0, qJDD(3) - g(1), 0, 0, 0, -t8 * t89 + t9 * t87 - g(1), 0, 0, 0, 0, 0, (t122 * t91 + (t51 - t150) * t137) * t87, t42 + (t122 * t94 - t127) * t87; 0, 0, 0, 0, 0, 0, 0, -t151, t87 * t80, -t142 * t79, -t111 * t83 + t168, 0, 0, 0, 0, 0, t97 * t91 - t144, t97 * t94 + t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145 * t153, -t141 * t153, (-t91 * t119 + t154) * t87, (-t94 * t119 - t155) * t87, -t48, -g(2) * t28 - g(3) * t30 + t166 * t94 + t99 * t91 + t5, g(2) * t29 - g(3) * t31 + t99 * t94 + (-t166 - t7) * t91;];
tau_reg = t1;
