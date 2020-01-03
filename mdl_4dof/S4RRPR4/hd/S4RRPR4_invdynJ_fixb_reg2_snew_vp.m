% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPR4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:33
% EndTime: 2019-12-31 17:02:36
% DurationCPUTime: 0.71s
% Computational Cost: add. (2490->130), mult. (3603->196), div. (0->0), fcn. (2274->8), ass. (0->102)
t100 = (qJD(1) + qJD(2));
t104 = sin(qJ(2));
t107 = cos(qJ(2));
t105 = sin(qJ(1));
t108 = cos(qJ(1));
t118 = g(1) * t105 - g(2) * t108;
t114 = qJDD(1) * pkin(1) + t118;
t119 = g(1) * t108 + g(2) * t105;
t84 = -qJD(1) ^ 2 * pkin(1) - t119;
t63 = t104 * t114 + t107 * t84;
t97 = t100 ^ 2;
t98 = qJDD(1) + qJDD(2);
t121 = -pkin(2) * t97 + qJ(3) * t98 + (2 * qJD(3) * t100) + t63;
t103 = sin(qJ(4));
t101 = sin(pkin(7));
t102 = cos(pkin(7));
t106 = cos(qJ(4));
t149 = -t101 * t103 + t102 * t106;
t72 = t149 * t100;
t133 = t102 * t103;
t116 = t101 * t106 + t133;
t74 = t116 * t100;
t61 = t74 * t72;
t146 = qJDD(4) + t61;
t151 = t103 * t146;
t150 = t106 * t146;
t120 = -t101 * g(3) + t102 * t121;
t144 = t102 * g(3);
t27 = t102 * t120 + t101 * (t101 * t121 + t144);
t148 = t102 * t97;
t110 = t101 ^ 2;
t99 = t102 ^ 2;
t147 = t110 + t99;
t83 = t147 * t97;
t70 = t72 ^ 2;
t71 = t74 ^ 2;
t113 = (pkin(3) * t148 - pkin(6) * t98 - t121) * t101;
t140 = t102 * t98;
t143 = t99 * t97;
t34 = -pkin(3) * t143 + pkin(6) * t140 + t120;
t17 = t103 * t34 - t106 * (t113 - t144);
t18 = -g(3) * t133 + t103 * t113 + t106 * t34;
t6 = t103 * t18 - t106 * t17;
t145 = t101 * t6;
t62 = -t104 * t84 + t107 * t114;
t51 = -t98 * pkin(2) - t97 * qJ(3) + qJDD(3) - t62;
t142 = -pkin(2) * t51 + qJ(3) * t27;
t141 = t101 * t98;
t38 = -pkin(3) * t140 + t51 + (-t110 * t97 - t143) * pkin(6);
t139 = t103 * t38;
t54 = qJDD(4) - t61;
t138 = t103 * t54;
t137 = t106 * t38;
t136 = t106 * t54;
t135 = t107 * t98;
t131 = t72 * qJD(4);
t130 = t74 * qJD(4);
t129 = qJD(4) * t103;
t128 = qJD(4) * t106;
t78 = t147 * t148;
t127 = pkin(2) * t140 - qJ(3) * t78 - t102 * t51;
t109 = qJD(4) ^ 2;
t52 = -t109 - t70;
t29 = t103 * t52 + t150;
t30 = t106 * t52 - t151;
t11 = -t101 * t29 + t102 * t30;
t39 = t149 * t98;
t56 = -t39 + 0.2e1 * t130;
t125 = t101 * (-pkin(6) * t29 + t139) + t102 * (-pkin(3) * t56 + pkin(6) * t30 - t137) - pkin(2) * t56 + qJ(3) * t11;
t69 = t116 * t98;
t31 = t103 * t39 - t106 * t69;
t32 = t103 * t69 + t106 * t39;
t14 = -t101 * t31 + t102 * t32;
t43 = -t70 - t71;
t7 = t103 * t17 + t106 * t18;
t124 = t101 * (-pkin(6) * t31 - t6) + t102 * (-pkin(3) * t43 + pkin(6) * t32 + t7) - pkin(2) * t43 + qJ(3) * t14;
t66 = -t71 - t109;
t35 = t106 * t66 - t138;
t36 = -t103 * t66 - t136;
t22 = -t101 * t35 + t102 * t36;
t58 = t69 + 0.2e1 * t131;
t123 = t101 * (-pkin(6) * t35 + t137) + t102 * (-pkin(3) * t58 + pkin(6) * t36 + t139) - pkin(2) * t58 + qJ(3) * t22;
t93 = t110 * t98;
t94 = t99 * t98;
t81 = t94 + t93;
t122 = pkin(2) * t83 + qJ(3) * t81 + t27;
t77 = t101 * t83;
t117 = -pkin(2) * t141 + qJ(3) * t77 + t101 * t51;
t2 = t102 * t7 - t145;
t115 = -pkin(6) * t145 + qJ(3) * t2 - pkin(2) * t38 + t102 * (-pkin(3) * t38 + pkin(6) * t7);
t85 = 0.2e1 * t101 * t140;
t65 = -t71 + t109;
t64 = t70 - t109;
t59 = t69 + t131;
t57 = t39 - t130;
t28 = (t101 * (t103 * t74 + t106 * t72) + t102 * (t103 * t72 - t106 * t74)) * qJD(4);
t24 = t101 * (t106 * t59 - t129 * t74) + t102 * (t103 * t59 + t128 * t74);
t23 = t101 * (-t103 * t57 - t128 * t72) + t102 * (t106 * t57 - t129 * t72);
t21 = t101 * (-t103 * t65 + t150) + t102 * (t106 * t65 + t151);
t20 = t101 * (t106 * t64 - t138) + t102 * (t103 * t64 + t136);
t13 = t101 * (-t103 * t58 - t106 * t56) + t102 * (-t103 * t56 + t106 * t58);
t1 = [0, 0, 0, 0, 0, qJDD(1), t118, t119, 0, 0, 0, 0, 0, 0, 0, t98, pkin(1) * (-t104 * t97 + t135) + t62, pkin(1) * (-t104 * t98 - t107 * t97) - t63, 0, pkin(1) * (t104 * t63 + t107 * t62), t93, t85, 0, t94, 0, 0, pkin(1) * (t102 * t135 - t104 * t78) + t127, pkin(1) * (-t101 * t135 + t104 * t77) + t117, pkin(1) * (t104 * t81 + t107 * t83) + t122, pkin(1) * (t104 * t27 - t107 * t51) + t142, t24, t13, t21, t23, t20, t28, pkin(1) * (t104 * t11 - t107 * t56) + t125, pkin(1) * (t104 * t22 - t107 * t58) + t123, pkin(1) * (t104 * t14 - t107 * t43) + t124, pkin(1) * (t104 * t2 - t107 * t38) + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t62, -t63, 0, 0, t93, t85, 0, t94, 0, 0, t127, t117, t122, t142, t24, t13, t21, t23, t20, t28, t125, t123, t124, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t141, -t83, t51, 0, 0, 0, 0, 0, 0, t56, t58, t43, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t71 - t70, t69, t61, t39, qJDD(4), -t17, -t18, 0, 0;];
tauJ_reg = t1;
