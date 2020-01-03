% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRR2
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:21
% DurationCPUTime: 0.81s
% Computational Cost: add. (3486->165), mult. (4598->222), div. (0->0), fcn. (2877->8), ass. (0->120)
t121 = sin(qJ(4));
t115 = qJDD(3) + qJDD(4);
t118 = qJD(1) + qJD(2);
t125 = cos(qJ(4));
t126 = cos(qJ(3));
t122 = sin(qJ(3));
t150 = t118 * t122;
t84 = -t125 * t126 * t118 + t121 * t150;
t86 = (t126 * t121 + t122 * t125) * t118;
t65 = t86 * t84;
t164 = -t65 + t115;
t167 = t121 * t164;
t166 = t125 * t164;
t146 = qJD(3) * t118;
t141 = t126 * t146;
t116 = qJDD(1) + qJDD(2);
t148 = t122 * t116;
t93 = t141 + t148;
t108 = t126 * t116;
t142 = t122 * t146;
t94 = t108 - t142;
t49 = -t84 * qJD(4) + t121 * t94 + t125 * t93;
t117 = qJD(3) + qJD(4);
t81 = t117 * t84;
t165 = t49 - t81;
t114 = t118 ^ 2;
t124 = sin(qJ(1));
t128 = cos(qJ(1));
t136 = t128 * g(1) + t124 * g(2);
t100 = -qJD(1) ^ 2 * pkin(1) - t136;
t123 = sin(qJ(2));
t127 = cos(qJ(2));
t135 = t124 * g(1) - t128 * g(2);
t133 = qJDD(1) * pkin(1) + t135;
t68 = t127 * t100 + t123 * t133;
t63 = -t114 * pkin(2) + t116 * pkin(6) + t68;
t157 = t122 * t63;
t53 = t126 * g(3) + t157;
t54 = -t122 * g(3) + t126 * t63;
t30 = t122 * t53 + t126 * t54;
t153 = t114 * t122;
t162 = t93 * pkin(7);
t163 = (pkin(3) * t153 + pkin(7) * t146 - g(3)) * t126 + qJDD(3) * pkin(3) - t157 - t162;
t82 = t84 ^ 2;
t83 = t86 ^ 2;
t113 = t117 ^ 2;
t103 = qJD(3) * pkin(3) - pkin(7) * t150;
t120 = t126 ^ 2;
t110 = t120 * t114;
t34 = -pkin(3) * t110 + t94 * pkin(7) - qJD(3) * t103 + t54;
t19 = t121 * t34 - t125 * t163;
t156 = t125 * t34;
t20 = t163 * t121 + t156;
t6 = t121 * t20 - t125 * t19;
t161 = t122 * t6;
t67 = -t123 * t100 + t127 * t133;
t62 = -t116 * pkin(2) - t114 * pkin(6) - t67;
t160 = -pkin(2) * t62 + pkin(6) * t30;
t38 = -t94 * pkin(3) - pkin(7) * t110 + t103 * t150 + t62;
t159 = t121 * t38;
t59 = t65 + t115;
t158 = t121 * t59;
t155 = t125 * t38;
t154 = t125 * t59;
t152 = t117 * t121;
t151 = t117 * t125;
t105 = t126 * t153;
t101 = qJDD(3) + t105;
t149 = t122 * t101;
t147 = t126 * (qJDD(3) - t105);
t119 = t122 ^ 2;
t109 = t119 * t114;
t129 = qJD(3) ^ 2;
t76 = -t147 - t122 * (-t109 - t129);
t92 = 0.2e1 * t141 + t148;
t145 = -pkin(2) * t92 + pkin(6) * t76 + t122 * t62;
t75 = t126 * (-t110 - t129) - t149;
t95 = t108 - 0.2e1 * t142;
t144 = pkin(2) * t95 + pkin(6) * t75 - t126 * t62;
t140 = t121 * t93 - t125 * t94;
t131 = (-qJD(4) + t117) * t86 - t140;
t44 = t49 + t81;
t21 = t121 * t131 - t125 * t44;
t22 = t121 * t44 + t125 * t131;
t10 = -t122 * t21 + t126 * t22;
t50 = -t82 - t83;
t7 = t121 * t19 + t125 * t20;
t143 = t122 * (-pkin(7) * t21 - t6) + t126 * (-pkin(3) * t50 + pkin(7) * t22 + t7) - pkin(2) * t50 + pkin(6) * t10;
t57 = -t113 - t82;
t31 = t121 * t57 + t166;
t32 = t125 * t57 - t167;
t14 = -t122 * t31 + t126 * t32;
t39 = (qJD(4) + t117) * t86 + t140;
t139 = t122 * (-pkin(7) * t31 + t159) + t126 * (-pkin(3) * t39 + pkin(7) * t32 - t155) - pkin(2) * t39 + pkin(6) * t14;
t77 = -t83 - t113;
t45 = t125 * t77 - t158;
t46 = -t121 * t77 - t154;
t25 = -t122 * t45 + t126 * t46;
t138 = t122 * (-pkin(7) * t45 + t155) + t126 * (-pkin(3) * t165 + pkin(7) * t46 + t159) - pkin(2) * t165 + pkin(6) * t25;
t98 = (t119 + t120) * t116;
t99 = t109 + t110;
t137 = pkin(2) * t99 + pkin(6) * t98 + t30;
t2 = t126 * t7 - t161;
t134 = pkin(6) * t2 - pkin(7) * t161 - pkin(2) * t38 + t126 * (-pkin(3) * t38 + pkin(7) * t7);
t79 = -t83 + t113;
t78 = t82 - t113;
t74 = t149 + t126 * (-t109 + t129);
t73 = t122 * (t110 - t129) + t147;
t70 = (t93 + t141) * t122;
t69 = (t94 - t142) * t126;
t66 = t122 * t95 + t126 * t92;
t64 = t83 - t82;
t48 = -t86 * qJD(4) - t140;
t28 = (t122 * (t121 * t86 - t125 * t84) + t126 * (-t121 * t84 - t125 * t86)) * t117;
t27 = t122 * (t125 * t78 - t158) + t126 * (t121 * t78 + t154);
t26 = t122 * (-t121 * t79 + t166) + t126 * (t125 * t79 + t167);
t16 = t122 * (t125 * t49 - t86 * t152) + t126 * (t121 * t49 + t86 * t151);
t15 = t122 * (-t121 * t48 + t84 * t151) + t126 * (t125 * t48 + t84 * t152);
t9 = t122 * (-t121 * t165 - t125 * t39) + t126 * (-t121 * t39 + t125 * t165);
t1 = [0, 0, 0, 0, 0, qJDD(1), t135, t136, 0, 0, 0, 0, 0, 0, 0, t116, pkin(1) * (-t123 * t114 + t127 * t116) + t67, pkin(1) * (-t127 * t114 - t123 * t116) - t68, 0, pkin(1) * (t123 * t68 + t127 * t67), t70, t66, t74, t69, t73, 0, pkin(1) * (t123 * t75 + t127 * t95) + t144, pkin(1) * (t123 * t76 - t127 * t92) + t145, pkin(1) * (t123 * t98 + t127 * t99) + t137, pkin(1) * (t123 * t30 - t127 * t62) + t160, t16, t9, t26, t15, t27, t28, pkin(1) * (t123 * t14 - t127 * t39) + t139, pkin(1) * (t123 * t25 - t127 * t165) + t138, pkin(1) * (t123 * t10 - t127 * t50) + t143, pkin(1) * (t123 * t2 - t127 * t38) + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t67, -t68, 0, 0, t70, t66, t74, t69, t73, 0, t144, t145, t137, t160, t16, t9, t26, t15, t27, t28, t139, t138, t143, t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, t109 - t110, t148, t105, t108, qJDD(3), -t53, -t54, 0, 0, t65, t64, t44, -t65, t131, t115, pkin(3) * t31 - t19, -t156 - t121 * (pkin(7) * t141 - t162 - t53) + (-t121 * t101 + t45) * pkin(3), pkin(3) * t21, pkin(3) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t64, t44, -t65, t131, t115, -t19, -t20, 0, 0;];
tauJ_reg = t1;
