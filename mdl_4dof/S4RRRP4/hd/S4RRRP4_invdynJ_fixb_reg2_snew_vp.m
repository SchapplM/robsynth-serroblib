% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRRP4
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
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:41
% EndTime: 2019-12-31 17:15:45
% DurationCPUTime: 1.13s
% Computational Cost: add. (2500->199), mult. (5488->246), div. (0->0), fcn. (3545->6), ass. (0->127)
t107 = qJDD(2) + qJDD(3);
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t115 = cos(qJ(2));
t112 = sin(qJ(2));
t140 = qJD(1) * t112;
t85 = -t114 * t115 * qJD(1) + t111 * t140;
t87 = (t115 * t111 + t112 * t114) * qJD(1);
t67 = t87 * t85;
t170 = t67 - t107;
t171 = t170 * pkin(3);
t169 = t111 * t170;
t168 = t114 * t170;
t110 = t115 ^ 2;
t117 = qJD(1) ^ 2;
t113 = sin(qJ(1));
t159 = cos(qJ(1));
t131 = t113 * g(1) - t159 * g(2);
t124 = qJDD(1) * pkin(1) + t131;
t125 = qJD(2) * pkin(2) - pkin(6) * t140;
t102 = t115 * qJDD(1);
t138 = qJD(1) * qJD(2);
t133 = t112 * t138;
t94 = t102 - t133;
t58 = t94 * pkin(2) + (pkin(6) * t110 + pkin(5)) * t117 - t125 * t140 + t124;
t108 = qJD(2) + qJD(3);
t74 = t108 * pkin(3) - t87 * qJ(4);
t118 = -t87 * t74 - qJDD(4) + t58;
t101 = t112 * qJDD(1);
t132 = t115 * t138;
t93 = t101 + t132;
t130 = t111 * t93 - t114 * t94;
t55 = -t87 * qJD(3) - t130;
t167 = -t55 * pkin(3) - t118;
t141 = t112 * t117;
t126 = t159 * g(1) + t113 * g(2);
t144 = qJDD(1) * pkin(5);
t89 = -t117 * pkin(1) - t126 + t144;
t153 = t112 * t89;
t51 = qJDD(2) * pkin(2) - t93 * pkin(6) - t153 + (pkin(2) * t141 + pkin(6) * t138 - g(3)) * t115;
t104 = t110 * t117;
t73 = -t112 * g(3) + t115 * t89;
t52 = -pkin(2) * t104 + t94 * pkin(6) - qJD(2) * t125 + t73;
t27 = t111 * t52 - t114 * t51;
t78 = t108 * t85;
t165 = qJ(4) * t78 + 0.2e1 * qJD(4) * t87 + t171 + t27;
t28 = t111 * t51 + t114 * t52;
t83 = t85 ^ 2;
t13 = -t83 * pkin(3) + t55 * qJ(4) - 0.2e1 * qJD(4) * t85 - t108 * t74 + t28;
t84 = t87 ^ 2;
t106 = t108 ^ 2;
t121 = (-qJD(3) + t108) * t87 - t130;
t127 = t111 * t94 + t114 * t93;
t56 = -t85 * qJD(3) + t127;
t45 = t56 + t78;
t22 = t111 * t121 - t114 * t45;
t164 = pkin(6) * t22;
t61 = -t106 - t83;
t33 = t111 * t61 - t168;
t163 = pkin(6) * t33;
t63 = t67 + t107;
t155 = t111 * t63;
t71 = -t84 - t106;
t46 = t114 * t71 - t155;
t162 = pkin(6) * t46;
t23 = t111 * t45 + t114 * t121;
t57 = -t83 - t84;
t160 = pkin(5) * (-t112 * t22 + t115 * t23) - pkin(1) * t57;
t34 = t114 * t61 + t169;
t139 = qJD(3) + t108;
t40 = t139 * t87 + t130;
t158 = pkin(5) * (-t112 * t33 + t115 * t34) - pkin(1) * t40;
t43 = -t139 * t85 + t127;
t150 = t114 * t63;
t47 = -t111 * t71 - t150;
t157 = pkin(5) * (-t112 * t46 + t115 * t47) - pkin(1) * t43;
t156 = t111 * t58;
t11 = t111 * t28 - t114 * t27;
t154 = t112 * t11;
t99 = t115 * t141;
t152 = t112 * (qJDD(2) + t99);
t151 = t114 * t58;
t149 = t115 * (qJDD(2) - t99);
t148 = qJ(4) * t111;
t147 = qJ(4) * t114;
t143 = t108 * t111;
t142 = t108 * t114;
t137 = -pkin(2) * t57 + pkin(6) * t23;
t136 = -pkin(2) * t40 + pkin(6) * t34;
t135 = -pkin(2) * t43 + pkin(6) * t47;
t12 = t111 * t27 + t114 * t28;
t72 = t115 * g(3) + t153;
t129 = t112 * t72 + t115 * t73;
t123 = pkin(3) * t71 - t13;
t10 = -t56 * qJ(4) - t165;
t119 = t10 - t171;
t116 = qJD(2) ^ 2;
t109 = t112 ^ 2;
t103 = t109 * t117;
t95 = t102 - 0.2e1 * t133;
t92 = t101 + 0.2e1 * t132;
t88 = t117 * pkin(5) + t124;
t76 = -t84 + t106;
t75 = t83 - t106;
t65 = t84 - t83;
t44 = t56 - t78;
t39 = pkin(2) * t46;
t35 = pkin(3) * t45;
t32 = pkin(2) * t33;
t30 = (t112 * (t111 * t87 - t114 * t85) + t115 * (-t111 * t85 - t114 * t87)) * t108;
t29 = -pkin(3) * t43 - qJ(4) * t63;
t25 = t112 * (t114 * t75 - t155) + t115 * (t111 * t75 + t150);
t24 = t112 * (-t111 * t76 - t168) + t115 * (t114 * t76 - t169);
t20 = pkin(2) * t22;
t18 = t83 * qJ(4) - t167;
t17 = t112 * (t114 * t56 - t87 * t143) + t115 * (t111 * t56 + t87 * t142);
t16 = t112 * (-t111 * t55 + t85 * t142) + t115 * (t114 * t55 + t85 * t143);
t14 = (-t71 - t83) * qJ(4) + t167;
t9 = pkin(3) * t10;
t8 = (t61 + t83) * qJ(4) + (-t40 + t55) * pkin(3) + t118;
t7 = t112 * (-t111 * t44 - t114 * t40) + t115 * (-t111 * t40 + t114 * t44);
t5 = (t45 + t56) * qJ(4) + t165;
t4 = -pkin(3) * t57 + qJ(4) * t121 + t13;
t3 = pkin(3) * t18 + qJ(4) * t13;
t2 = -t111 * t10 + t114 * t13;
t1 = t114 * t10 + t111 * t13;
t6 = [0, 0, 0, 0, 0, qJDD(1), t131, t126, 0, 0, (t93 + t132) * t112, t112 * t95 + t115 * t92, t152 + t115 * (-t103 + t116), (t94 - t133) * t115, t112 * (t104 - t116) + t149, 0, t115 * t88 + pkin(1) * t95 + pkin(5) * (t115 * (-t104 - t116) - t152), -t112 * t88 - pkin(1) * t92 + pkin(5) * (-t149 - t112 * (-t103 - t116)), pkin(1) * (t103 + t104) + (t109 + t110) * t144 + t129, pkin(1) * t88 + pkin(5) * t129, t17, t7, t24, t16, t25, t30, t112 * (-t156 - t163) + t115 * (t136 + t151) + t158, t112 * (-t151 - t162) + t115 * (t135 - t156) + t157, t112 * (-t11 - t164) + t115 * (t12 + t137) + t160, -pkin(6) * t154 + t115 * (pkin(2) * t58 + pkin(6) * t12) + pkin(1) * t58 + pkin(5) * (t115 * t12 - t154), t17, t7, t24, t16, t25, t30, t112 * (-t111 * t8 + t147 * t170 - t163) + t115 * (t114 * t8 + t148 * t170 + t136) + t158, t112 * (-t111 * t29 + t114 * t14 - t162) + t115 * (t111 * t14 + t114 * t29 + t135) + t157, t112 * (-t111 * t4 + t114 * t5 - t164) + t115 * (t111 * t5 + t114 * t4 + t137) + t160, t112 * (-pkin(6) * t1 - t10 * t147 - t111 * t3) + t115 * (pkin(2) * t18 + pkin(6) * t2 - t10 * t148 + t114 * t3) + pkin(1) * t18 + pkin(5) * (-t112 * t1 + t115 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t103 - t104, t101, t99, t102, qJDD(2), -t72, -t73, 0, 0, t67, t65, t45, -t67, t121, t107, -t27 + t32, t39 - t28, t20, pkin(2) * t11, t67, t65, t45, -t67, t121, t107, t119 + t32, t39 + t123, -t35 + t20, pkin(2) * t1 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t65, t45, -t67, t121, t107, -t27, -t28, 0, 0, t67, t65, t45, -t67, t121, t107, t119, t123, -t35, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t43, t57, -t18;];
tauJ_reg = t6;
