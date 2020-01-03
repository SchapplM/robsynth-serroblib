% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% tau_reg [5x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:44
% EndTime: 2019-12-31 18:16:47
% DurationCPUTime: 0.94s
% Computational Cost: add. (772->209), mult. (1306->236), div. (0->0), fcn. (601->4), ass. (0->132)
t124 = qJ(5) * qJDD(1);
t162 = qJD(1) * qJD(5) + t124;
t128 = qJ(5) * qJD(1);
t82 = -pkin(1) - pkin(6);
t38 = t82 * qJD(1) + qJD(2);
t79 = cos(qJ(3));
t18 = (t38 + t128) * t79;
t130 = qJD(4) - t18;
t81 = -pkin(3) - pkin(4);
t10 = qJD(3) * t81 + t130;
t80 = cos(qJ(1));
t69 = g(2) * t80;
t78 = sin(qJ(1));
t70 = g(1) * t78;
t144 = t70 - t69;
t109 = qJD(3) * pkin(3) - qJD(4);
t151 = t79 * t38;
t20 = -t109 - t151;
t77 = sin(qJ(3));
t31 = t77 * t38;
t74 = qJD(3) * qJ(4);
t22 = t31 + t74;
t27 = qJD(3) * t151;
t37 = t82 * qJDD(1) + qJDD(2);
t29 = t77 * t37;
t72 = qJDD(3) * qJ(4);
t73 = qJD(3) * qJD(4);
t6 = t27 + t29 + t72 + t73;
t135 = qJD(3) * t77;
t26 = t38 * t135;
t30 = t79 * t37;
t117 = t30 - qJDD(4) - t26;
t134 = qJDD(3) * pkin(3);
t7 = -t117 - t134;
t85 = t6 * t77 - t7 * t79 + (t20 * t77 + t22 * t79) * qJD(3);
t161 = t85 - t144;
t120 = t81 * t77;
t102 = -qJ(2) + t120;
t61 = t79 * qJ(4);
t25 = t61 + t102;
t160 = qJD(1) * t25;
t155 = t77 * pkin(3);
t116 = qJ(2) + t155;
t34 = t116 - t61;
t21 = qJD(1) * t34;
t143 = g(1) * t80 + g(2) * t78;
t137 = qJ(5) + t82;
t158 = t81 * qJDD(3);
t75 = t77 ^ 2;
t76 = t79 ^ 2;
t141 = t75 + t76;
t35 = t141 * qJDD(1);
t123 = qJD(1) * qJD(3);
t114 = qJ(5) * t123;
t157 = t79 * t114 + t162 * t77;
t150 = t79 * t80;
t152 = t78 * t79;
t156 = g(1) * t152 - g(2) * t150 - g(3) * t77;
t139 = qJ(4) * t77;
t95 = t79 * t81 - t139;
t154 = t77 * t78;
t153 = t77 * t80;
t83 = qJD(3) ^ 2;
t149 = t82 * t83;
t148 = pkin(3) * t152 + t78 * t139;
t147 = g(1) * t150 + g(2) * t152;
t131 = t79 * qJD(4);
t58 = t79 * qJDD(1);
t146 = -qJ(4) * t58 - qJD(1) * t131;
t145 = t80 * pkin(1) + t78 * qJ(2);
t142 = t75 - t76;
t84 = qJD(1) ^ 2;
t140 = t83 + t84;
t138 = t84 * qJ(2);
t136 = pkin(1) * qJDD(1);
t57 = t77 * t128;
t13 = t57 + t22;
t133 = t13 * qJD(3);
t132 = t21 * qJD(1);
t12 = qJD(5) + t160;
t129 = qJD(5) + t12;
t127 = qJDD(3) * t77;
t126 = qJDD(3) * t82;
t125 = t77 * qJDD(1);
t121 = qJDD(1) * qJ(2);
t119 = 0.2e1 * qJD(1) * qJD(2);
t118 = -0.2e1 * t123;
t115 = qJ(2) * t123;
t113 = pkin(3) * t154 + t80 * pkin(6) + t145;
t112 = -t30 + t156;
t111 = -t12 - t160;
t110 = qJD(3) * t137;
t107 = t79 * t118;
t106 = qJDD(2) - t136;
t105 = -g(2) * t153 + g(3) * t79 - t29;
t101 = pkin(3) * t79 + t139;
t99 = -qJDD(4) - t112;
t62 = t80 * qJ(2);
t98 = pkin(3) * t153 - t80 * t61 + t62;
t97 = 0.2e1 * t21 * qJD(3);
t96 = t119 + 0.2e1 * t121;
t87 = t95 * qJD(3) - qJD(2);
t11 = t87 + t131;
t86 = t102 * qJDD(1) + qJDD(5) - t146;
t2 = t87 * qJD(1) + t86;
t94 = qJD(1) * t11 + qJDD(1) * t25 + t2;
t40 = t77 * t114;
t92 = -t117 + t40 + t158;
t3 = -t162 * t79 + t92;
t4 = t6 + t157;
t93 = t10 * t135 + t4 * t77 - t144 + (t133 - t3) * t79;
t91 = t101 * qJD(3) + qJD(2);
t90 = t96 - t149;
t16 = t91 - t131;
t5 = t91 * qJD(1) + t116 * qJDD(1) + t146;
t89 = -qJD(1) * t16 - qJDD(1) * t34 + t149 - t5;
t88 = -g(1) * t154 - t105 + 0.2e1 * t72 + 0.2e1 * t73;
t59 = qJDD(3) * t79;
t43 = t79 * t84 * t77;
t42 = t79 * t126;
t39 = -t76 * t84 - t83;
t36 = -qJDD(3) + t43;
t33 = t137 * t79;
t32 = t137 * t77;
t28 = t101 * qJD(1);
t24 = -t140 * t77 + t59;
t23 = t140 * t79 + t127;
t19 = t95 * qJD(1);
t17 = t31 + t57;
t15 = t77 * qJD(5) + t79 * t110;
t14 = -t79 * qJD(5) + t77 * t110;
t1 = [qJDD(1), t144, t143, qJDD(2) - 0.2e1 * t136 - t144, t96 - t143, -t106 * pkin(1) - g(1) * (-t78 * pkin(1) + t62) - g(2) * t145 + (t119 + t121) * qJ(2), t76 * qJDD(1) + t77 * t107, 0.2e1 * t142 * t123 - 0.2e1 * t77 * t58, -t83 * t77 + t59, -t83 * t79 - t127, 0, 0.2e1 * t79 * t115 + t42 + (-t143 + t90) * t77, (-0.2e1 * t115 - t126) * t77 + t90 * t79 - t147, t42 + t79 * t97 + (-t143 - t89) * t77, -t35 * t82 - t161, (t97 + t126) * t77 + t89 * t79 + t147, t5 * t34 + t21 * t16 - g(1) * (t82 * t78 + t98) - g(2) * (-t78 * t61 + t113) + t85 * t82, t33 * qJDD(3) + (t111 * t79 - t14) * qJD(3) + (-t143 - t94) * t77, t32 * qJDD(3) + t94 * t79 + (t111 * t77 + t15) * qJD(3) + t147, (t32 * t77 + t33 * t79) * qJDD(1) + (-t14 * t79 + t15 * t77 + (t32 * t79 - t33 * t77) * qJD(3)) * qJD(1) + t93, t4 * t32 + t13 * t15 - t3 * t33 + t10 * t14 + t2 * t25 + t12 * t11 - g(1) * (pkin(4) * t153 + t98) - g(2) * (-t80 * qJ(5) + t113) + (-g(1) * t137 - g(2) * (t77 * pkin(4) - t61)) * t78; 0, 0, 0, qJDD(1), -t84, t106 - t138 - t144, 0, 0, 0, 0, 0, t24, -t23, t24, -t35, t23, -t132 + t161, t24, t23, t35, t12 * qJD(1) + t93; 0, 0, 0, 0, 0, 0, t43, -t142 * t84, t58, -t125, qJDD(3), -t79 * t138 - t112, (t138 + t70) * t77 + t105, 0.2e1 * t134 + (-t21 * t79 - t28 * t77) * qJD(1) + t99, -t101 * qJDD(1) + ((t22 - t74) * t79 + (t109 + t20) * t77) * qJD(1), (-t21 * t77 + t28 * t79) * qJD(1) + t88, t6 * qJ(4) - t7 * pkin(3) - t21 * t28 - t20 * t31 - g(1) * t148 - g(3) * (t61 - t155) + t101 * t69 + (qJD(4) - t151) * t22, t79 * t124 + t17 * qJD(3) - t26 - t40 - 0.2e1 * t158 + (t129 * t79 + t19 * t77) * qJD(1) + t99, -t18 * qJD(3) + t27 + (t12 * t77 - t19 * t79) * qJD(1) + t88 + t157, -t95 * qJDD(1) + (-t13 + t17 + t74) * t79 * qJD(1), t4 * qJ(4) + t3 * t81 - t10 * t17 - t12 * t19 - g(1) * (pkin(4) * t152 + t148) - g(3) * (t61 + t120) + t130 * t13 - t95 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t58, t39, -t22 * qJD(3) + t79 * t132 + t156 + t7, t36, t39, -t58, -t133 + (-t129 * qJD(1) - t124) * t79 + t92 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107 - t125, t77 * t118 + t58, -t141 * t84, (t10 * t79 - t13 * t77 + t87) * qJD(1) + t86 + t143;];
tau_reg = t1;
