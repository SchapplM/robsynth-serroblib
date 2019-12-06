% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:22
% EndTime: 2019-12-05 17:38:24
% DurationCPUTime: 1.12s
% Computational Cost: add. (1303->217), mult. (2319->253), div. (0->0), fcn. (1335->8), ass. (0->126)
t86 = sin(qJ(1));
t89 = cos(qJ(1));
t144 = g(1) * t89 + g(2) * t86;
t83 = pkin(1) + qJ(3);
t160 = qJD(1) * t83;
t47 = -qJD(2) + t160;
t164 = -t47 * qJD(1) - t144;
t128 = qJD(1) * qJD(4);
t85 = sin(qJ(4));
t120 = t85 * t128;
t88 = cos(qJ(4));
t130 = t88 * qJDD(1);
t163 = t120 - t130;
t162 = g(1) * t86 - g(2) * t89;
t78 = t85 ^ 2;
t79 = t88 ^ 2;
t143 = t78 + t79;
t57 = qJ(2) * qJD(1) + qJD(3);
t46 = -pkin(6) * qJD(1) + t57;
t161 = t143 * t46;
t84 = sin(qJ(5));
t87 = cos(qJ(5));
t35 = t84 * t88 + t87 * t85;
t31 = t35 * qJD(1);
t73 = qJD(4) + qJD(5);
t119 = t88 * t128;
t131 = t85 * qJDD(1);
t100 = t119 + t131;
t127 = qJD(3) * qJD(1);
t81 = qJDD(1) * pkin(1);
t129 = t81 - qJDD(2);
t74 = qJDD(1) * qJ(3);
t38 = t74 + t127 + t129;
t22 = pkin(4) * t100 + t38;
t139 = qJD(4) * t85;
t75 = qJDD(1) * qJ(2);
t76 = qJD(1) * qJD(2);
t122 = qJDD(3) + t75 + t76;
t37 = -pkin(6) * qJDD(1) + t122;
t30 = t88 * t37;
t13 = qJDD(4) * pkin(4) + pkin(7) * t163 - t46 * t139 + t30;
t137 = qJD(5) * t84;
t138 = qJD(4) * t88;
t16 = -pkin(7) * t100 + t138 * t46 + t85 * t37;
t140 = qJD(1) * t88;
t28 = -pkin(7) * t140 + t88 * t46;
t24 = qJD(4) * pkin(4) + t28;
t141 = qJD(1) * t85;
t27 = -pkin(7) * t141 + t85 * t46;
t1 = (qJD(5) * t24 + t16) * t87 + t84 * t13 - t27 * t137;
t82 = -pkin(6) + qJ(2);
t159 = qJD(4) * (qJD(2) + t47 + t160) + qJDD(4) * t82;
t158 = t31 ^ 2;
t124 = t84 * t141;
t33 = t140 * t87 - t124;
t157 = t33 ^ 2;
t66 = 0.2e1 * t76;
t156 = g(3) * t85;
t155 = t85 * pkin(4);
t154 = pkin(7) - t82;
t152 = t33 * t31;
t151 = t33 * t73;
t92 = qJD(1) ^ 2;
t150 = t78 * t92;
t149 = t84 * t27;
t148 = t87 * t27;
t18 = t73 * t35;
t36 = -t84 * t85 + t87 * t88;
t72 = qJDD(4) + qJDD(5);
t147 = -t18 * t73 + t36 * t72;
t146 = t89 * pkin(1) + t86 * qJ(2);
t91 = qJD(4) ^ 2;
t142 = -t91 - t92;
t53 = t83 + t155;
t34 = qJD(1) * t53 - qJD(2);
t136 = t34 * qJD(1);
t133 = qJDD(4) * t85;
t132 = t83 * qJDD(1);
t126 = t88 * t92 * t85;
t125 = t89 * qJ(3) + t146;
t41 = t154 * t88;
t123 = qJDD(2) - t162;
t121 = t143 * t37;
t117 = 0.2e1 * t75 + t66 - t144;
t115 = t73 * t88;
t114 = -t81 + t123;
t113 = t85 * t119;
t109 = -t87 * t130 + t131 * t84;
t7 = qJD(1) * t18 + t109;
t108 = -t33 * t18 - t36 * t7;
t19 = t115 * t87 - t137 * t85 - t139 * t84;
t101 = -qJD(5) * t124 - t163 * t84;
t8 = (qJD(1) * t115 + t131) * t87 + t101;
t107 = t19 * t31 + t35 * t8;
t64 = t89 * qJ(2);
t106 = -t83 * t86 + t64;
t105 = -t19 * t73 - t35 * t72;
t12 = t84 * t24 + t148;
t40 = t154 * t85;
t21 = -t87 * t40 - t84 * t41;
t20 = t84 * t40 - t87 * t41;
t103 = -t74 + t114;
t102 = t47 * qJD(3) + t38 * t83;
t2 = -qJD(5) * t12 + t87 * t13 - t84 * t16;
t11 = t87 * t24 - t149;
t97 = t1 * t35 - t11 * t18 + t12 * t19 + t2 * t36 - t144;
t96 = -t82 * t91 + t127 + t132 + t162 + t38;
t80 = qJ(4) + qJ(5);
t60 = sin(t80);
t61 = cos(t80);
t95 = g(3) * t61 + t144 * t60 + t34 * t31 - t1;
t94 = (-t140 * t73 - t131) * t87 - t101;
t93 = g(3) * t60 - t144 * t61 - t34 * t33 + t2;
t90 = -pkin(7) - pkin(6);
t59 = t79 * t92;
t58 = qJDD(4) * t88;
t48 = pkin(4) * t138 + qJD(3);
t26 = t85 * qJD(2) - qJD(4) * t41;
t25 = t88 * qJD(2) + t154 * t139;
t15 = t87 * t28 - t149;
t14 = -t84 * t28 - t148;
t10 = t157 - t158;
t6 = t94 + t151;
t4 = -qJD(5) * t21 + t87 * t25 - t84 * t26;
t3 = qJD(5) * t20 + t84 * t25 + t87 * t26;
t5 = [0, 0, 0, 0, 0, qJDD(1), t162, t144, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t81 + t123, t117, t129 * pkin(1) - g(1) * (-t86 * pkin(1) + t64) - g(2) * t146 + (t75 + t66) * qJ(2), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t117, -t103 + 0.2e1 * t127 + t132, -g(1) * t106 - g(2) * t125 + qJ(2) * t122 + t57 * qJD(2) + t102, t79 * qJDD(1) - 0.2e1 * t113, -0.2e1 * t85 * t130 + 0.2e1 * (t78 - t79) * t128, -t91 * t85 + t58, t78 * qJDD(1) + 0.2e1 * t113, -t91 * t88 - t133, 0, t159 * t88 + t96 * t85, -t159 * t85 + t96 * t88, t144 + t143 * (-qJDD(1) * t82 - t37 - t76), -g(1) * (-t89 * pkin(6) + t106) - g(2) * (-t86 * pkin(6) + t125) + t82 * t121 + qJD(2) * t161 + t102, t108, t18 * t31 - t33 * t19 + t7 * t35 - t36 * t8, t147, t107, t105, 0, t162 * t60 + t34 * t19 + t20 * t72 + t22 * t35 + t48 * t31 + t4 * t73 + t53 * t8, t162 * t61 - t34 * t18 - t21 * t72 + t22 * t36 - t3 * t73 + t48 * t33 - t53 * t7, t20 * t7 - t21 * t8 - t3 * t31 - t4 * t33 - t97, t1 * t21 + t12 * t3 + t2 * t20 + t11 * t4 + t22 * t53 + t34 * t48 - g(1) * (t89 * t90 + t64) - g(2) * (t89 * t155 + t125) + (g(1) * t53 - g(2) * t90) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t92, -t92 * qJ(2) + t114, 0, 0, 0, 0, 0, 0, 0, -t92, -qJDD(1), (-qJD(3) - t57) * qJD(1) + t103, 0, 0, 0, 0, 0, 0, -0.2e1 * t119 - t131, 0.2e1 * t120 - t130, t59 + t150, (-qJD(3) - t161) * qJD(1) + t103, 0, 0, 0, 0, 0, 0, t94 - t151, t31 * t73 + t7, t157 + t158, -t11 * t33 - t12 * t31 - t162 - t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t92, t164 + t122, 0, 0, 0, 0, 0, 0, t142 * t85 + t58, t142 * t88 - t133, -t143 * qJDD(1), t121 + t164, 0, 0, 0, 0, 0, 0, -qJD(1) * t31 + t147, -qJD(1) * t33 + t105, -t107 - t108, t97 - t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t59 - t150, t130, -t126, -t131, qJDD(4), t164 * t88 + t156 + t30, g(3) * t88 + (-t37 - t164) * t85, 0, 0, t152, t10, -t109, -t152, t6, t72, -t14 * t73 + (-t137 * t73 - t140 * t31 + t72 * t87) * pkin(4) + t93, t15 * t73 + (-qJD(5) * t73 * t87 - t140 * t33 - t72 * t84) * pkin(4) + t95, (t12 + t14) * t33 + (-t11 + t15) * t31 + (t7 * t87 - t8 * t84 + (-t31 * t87 + t33 * t84) * qJD(5)) * pkin(4), -t11 * t14 - t12 * t15 + (t156 + t1 * t84 + t2 * t87 + (-t11 * t84 + t12 * t87) * qJD(5) + (-t144 - t136) * t88) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t10, -t109, -t152, t6, t72, t12 * t73 + t93, t11 * t73 + t95, 0, 0;];
tau_reg = t5;
