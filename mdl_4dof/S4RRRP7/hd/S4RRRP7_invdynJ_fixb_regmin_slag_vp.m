% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tau_reg [4x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S4RRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:06
% EndTime: 2019-12-31 17:21:10
% DurationCPUTime: 1.67s
% Computational Cost: add. (1333->274), mult. (3020->363), div. (0->0), fcn. (1869->6), ass. (0->128)
t78 = cos(qJ(2));
t131 = qJD(1) * t78;
t167 = qJD(3) - t131;
t75 = sin(qJ(2));
t124 = qJD(3) * t75;
t166 = qJD(1) * t124 - qJDD(2);
t118 = t75 * qJDD(1);
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t12 = ((qJD(3) + t131) * qJD(2) + t118) * t74 + t166 * t77;
t76 = sin(qJ(1));
t79 = cos(qJ(1));
t104 = g(1) * t79 + g(2) * t76;
t165 = t104 * t75;
t117 = qJD(1) * qJD(2);
t70 = t78 * qJDD(1);
t164 = -t75 * t117 + t70;
t120 = t77 * qJD(2);
t132 = qJD(1) * t75;
t43 = t74 * t132 - t120;
t128 = qJD(2) * t74;
t45 = t77 * t132 + t128;
t53 = -qJD(2) * pkin(2) + pkin(5) * t132;
t10 = pkin(3) * t43 - qJ(4) * t45 + t53;
t41 = qJDD(3) - t164;
t159 = pkin(6) * t41;
t163 = -t10 * t167 + t159;
t161 = t45 ^ 2;
t160 = pkin(3) * t41;
t156 = g(3) * t75;
t155 = g(3) * t78;
t98 = pkin(2) * t78 + pkin(6) * t75 + pkin(1);
t38 = t98 * qJD(1);
t67 = pkin(5) * t131;
t54 = qJD(2) * pkin(6) + t67;
t16 = -t38 * t74 + t54 * t77;
t9 = qJ(4) * t167 + t16;
t154 = t167 * t9;
t111 = t78 * t117;
t11 = -qJD(3) * t120 + (-t111 - t118) * t77 + t166 * t74;
t153 = t11 * t74;
t152 = t16 * t167;
t151 = t43 * t167;
t150 = t45 * t43;
t149 = t45 * t167;
t148 = t45 * t77;
t147 = t74 * t75;
t146 = t74 * t78;
t145 = t75 * t77;
t144 = t76 * t77;
t143 = t77 * t41;
t142 = t77 * t78;
t141 = t77 * t79;
t140 = t79 * t74;
t100 = pkin(3) * t74 - qJ(4) * t77;
t139 = -qJD(4) * t74 + t167 * t100 - t67;
t123 = qJD(3) * t77;
t107 = pkin(2) * t75 - pkin(6) * t78;
t48 = t107 * qJD(2);
t138 = -t123 * t98 + t74 * t48;
t137 = (g(1) * t141 + g(2) * t144) * t75;
t136 = pkin(5) * t142 - t74 * t98;
t72 = t75 ^ 2;
t135 = -t78 ^ 2 + t72;
t134 = pkin(6) * qJD(3);
t133 = qJ(4) * t41;
t130 = qJD(2) * t43;
t129 = qJD(2) * t45;
t127 = qJD(2) * t75;
t126 = qJD(2) * t78;
t125 = qJD(3) * t74;
t122 = qJD(3) * t78;
t121 = t53 * qJD(3);
t15 = -t38 * t77 - t54 * t74;
t119 = qJD(4) - t15;
t116 = t167 * t128;
t115 = t167 * t120;
t114 = t167 * t125;
t113 = pkin(5) * t74 + pkin(3);
t20 = qJD(1) * t48 - t98 * qJDD(1);
t30 = t164 * pkin(5) + qJDD(2) * pkin(6);
t109 = t54 * t123 - t38 * t125 - t77 * t20 + t74 * t30;
t32 = t76 * t146 + t141;
t34 = t78 * t140 - t144;
t106 = -g(1) * t32 + g(2) * t34;
t33 = t76 * t142 - t140;
t35 = t78 * t141 + t74 * t76;
t105 = g(1) * t33 - g(2) * t35;
t103 = g(1) * t76 - g(2) * t79;
t8 = -pkin(3) * t167 + t119;
t102 = -t74 * t9 + t77 * t8;
t65 = pkin(5) * t118;
t31 = -qJDD(2) * pkin(2) + pkin(5) * t111 + t65;
t101 = pkin(3) * t77 + qJ(4) * t74;
t99 = t121 - t159;
t97 = pkin(2) + t101;
t96 = pkin(5) + t100;
t95 = t134 * t167 + t155;
t93 = -t123 * t167 - t41 * t74;
t92 = -t114 + t143;
t91 = t125 * t98 + t48 * t77;
t90 = -t38 * t123 - t54 * t125 + t74 * t20 + t77 * t30;
t89 = -0.2e1 * pkin(1) * t117 - pkin(5) * qJDD(2);
t2 = pkin(3) * t12 + qJ(4) * t11 - qJD(4) * t45 + t31;
t88 = -t2 - t95;
t81 = qJD(1) ^ 2;
t87 = pkin(1) * t81 + t104;
t86 = -t104 * t78 - t156;
t85 = g(1) * t34 + g(2) * t32 + g(3) * t147 - t109;
t80 = qJD(2) ^ 2;
t84 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t80 + t103;
t83 = t10 * t45 + qJDD(4) - t85;
t82 = -g(1) * t35 - g(2) * t33 - g(3) * t145 + t90;
t47 = t107 * qJD(1);
t36 = t74 * t47;
t28 = t96 * t75;
t22 = t113 * t78 + t77 * t98;
t21 = -qJ(4) * t78 + t136;
t19 = pkin(3) * t45 + qJ(4) * t43;
t18 = -t113 * t132 - t47 * t77;
t17 = t36 + (-pkin(5) * t77 + qJ(4)) * t132;
t7 = (t101 * qJD(3) - qJD(4) * t77) * t75 + t96 * t126;
t6 = -pkin(3) * t127 + (t77 * t122 - t74 * t127) * pkin(5) - t91;
t5 = -t11 + t151;
t4 = qJ(4) * t127 - qJD(4) * t78 + (-t75 * t120 - t74 * t122) * pkin(5) + t138;
t3 = qJDD(4) + t109 - t160;
t1 = qJD(4) * t167 + t133 + t90;
t13 = [qJDD(1), t103, t104, qJDD(1) * t72 + 0.2e1 * t75 * t111, -0.2e1 * t135 * t117 + 0.2e1 * t75 * t70, qJDD(2) * t75 + t78 * t80, qJDD(2) * t78 - t75 * t80, 0, t89 * t75 + t84 * t78, -t84 * t75 + t89 * t78, -t11 * t145 + (t78 * t120 - t74 * t124) * t45, (-t43 * t77 - t45 * t74) * t126 + (t153 - t12 * t77 + (t43 * t74 - t148) * qJD(3)) * t75, (t11 + t115) * t78 + (t92 + t129) * t75, (t12 - t116) * t78 + (t93 - t130) * t75, t127 * t167 - t41 * t78, t91 * t167 - t98 * t143 + (t53 * t128 + (t93 + t130) * pkin(5) + t109) * t78 + (t77 * t121 + t15 * qJD(2) + t31 * t74 + (t12 + t116) * pkin(5)) * t75 + t105, -t138 * t167 - t136 * t41 + (t53 * t120 + (t114 + t129) * pkin(5) + t90) * t78 + (-t74 * t121 - qJD(2) * t16 + t31 * t77 + (-t11 + t115) * pkin(5)) * t75 + t106, t12 * t28 - t22 * t41 + t43 * t7 - t167 * t6 + (t10 * t128 + t3) * t78 + (-qJD(2) * t8 + t10 * t123 + t2 * t74) * t75 + t105, -t11 * t22 - t12 * t21 - t4 * t43 + t45 * t6 + t102 * t126 + (-t1 * t74 + t3 * t77 + (-t74 * t8 - t77 * t9) * qJD(3) + t103) * t75, t11 * t28 + t21 * t41 + t4 * t167 - t45 * t7 + (-t10 * t120 - t1) * t78 + (qJD(2) * t9 + t10 * t125 - t2 * t77) * t75 - t106, t1 * t21 + t9 * t4 + t2 * t28 + t10 * t7 + t3 * t22 + t8 * t6 - g(1) * (-pkin(3) * t33 - qJ(4) * t32) - g(2) * (pkin(3) * t35 + qJ(4) * t34) + (-g(1) * pkin(5) - g(2) * t98) * t79 + (-g(2) * pkin(5) + g(1) * t98) * t76; 0, 0, 0, -t75 * t81 * t78, t135 * t81, t118, t70, qJDD(2), t87 * t75 - t155 - t65, t156 + (-pkin(5) * qJDD(1) + t87) * t78, t148 * t167 - t153, (-t11 - t151) * t77 + (-t12 - t149) * t74, (-t142 * t167 - t45 * t75) * qJD(1) - t93, (t146 * t167 + t43 * t75) * qJD(1) + t92, -t167 * t132, -pkin(2) * t12 + t99 * t74 + (-t155 - t31 - (t47 + t134) * t167) * t77 + (-t53 * t146 - t15 * t75 + (-t147 * t167 - t43 * t78) * pkin(5)) * qJD(1) + t137, pkin(2) * t11 + t36 * t167 + t99 * t77 + (-t53 * t142 + t16 * t75 + (-t145 * t167 - t45 * t78) * pkin(5)) * qJD(1) + (t31 + t95 - t165) * t74, -t12 * t97 + t8 * t132 + t139 * t43 - t163 * t74 + t167 * t18 + t88 * t77 + t137, t17 * t43 - t18 * t45 + (t1 + t167 * t8 + (qJD(3) * t45 - t12) * pkin(6)) * t77 + (t3 - t154 + (qJD(3) * t43 - t11) * pkin(6)) * t74 + t86, -t9 * t132 - t11 * t97 - t17 * t167 - t139 * t45 + t163 * t77 + (t88 + t165) * t74, -t9 * t17 - t8 * t18 + t139 * t10 + (qJD(3) * t102 + t1 * t77 + t3 * t74 + t86) * pkin(6) + (-t2 - t155 + t165) * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, -t43 ^ 2 + t161, t5, -t12 + t149, t41, -t45 * t53 + t152 + t85, t15 * t167 + t43 * t53 - t82, -t19 * t43 + t152 + 0.2e1 * t160 - t83, pkin(3) * t11 - qJ(4) * t12 + (-t16 + t9) * t45 + (t8 - t119) * t43, 0.2e1 * t133 - t10 * t43 + t19 * t45 - (-0.2e1 * qJD(4) + t15) * t167 + t82, t1 * qJ(4) - t3 * pkin(3) - t10 * t19 - t8 * t16 - g(1) * (-pkin(3) * t34 + qJ(4) * t35) - g(2) * (-pkin(3) * t32 + qJ(4) * t33) + t119 * t9 + t100 * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 + t150, t5, -t167 ^ 2 - t161, -t154 + t83 - t160;];
tau_reg = t13;
