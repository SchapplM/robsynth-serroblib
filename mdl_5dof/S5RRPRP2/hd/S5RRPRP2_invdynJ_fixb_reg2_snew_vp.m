% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:46
% DurationCPUTime: 0.92s
% Computational Cost: add. (2595->146), mult. (3643->185), div. (0->0), fcn. (1959->8), ass. (0->100)
t107 = sin(pkin(8));
t108 = cos(pkin(8));
t110 = sin(qJ(2));
t113 = cos(qJ(2));
t109 = sin(qJ(4));
t112 = cos(qJ(4));
t101 = qJD(1) + qJD(2);
t99 = t101 ^ 2;
t89 = t112 * t99 * t109;
t84 = qJDD(4) - t89;
t148 = t112 * t84;
t115 = qJD(4) ^ 2;
t104 = t109 ^ 2;
t156 = t104 * t99;
t85 = t115 + t156;
t59 = -t109 * t85 + t148;
t141 = qJD(4) * t101;
t100 = qJDD(1) + qJDD(2);
t145 = t109 * t100;
t71 = 0.2e1 * t112 * t141 + t145;
t37 = t107 * t59 + t108 * t71;
t167 = pkin(1) * (t110 * (-t107 * t71 + t108 * t59) + t113 * t37);
t166 = pkin(2) * t37 + pkin(3) * t71 + pkin(7) * t59;
t161 = qJ(5) * t71;
t111 = sin(qJ(1));
t114 = cos(qJ(1));
t133 = t111 * g(1) - t114 * g(2);
t81 = qJDD(1) * pkin(1) + t133;
t123 = t114 * g(1) + t111 * g(2);
t82 = -qJD(1) ^ 2 * pkin(1) - t123;
t47 = -t110 * t82 + t113 * t81;
t118 = t100 * pkin(2) + t47;
t48 = t110 * t81 + t113 * t82;
t46 = -t99 * pkin(2) + t48;
t32 = t107 * t118 + t108 * t46;
t29 = -t99 * pkin(3) + t100 * pkin(7) + t32;
t143 = -g(3) + qJDD(3);
t94 = t112 * t143;
t22 = t109 * t29 - t94;
t23 = t109 * t143 + t112 * t29;
t8 = t109 * t22 + t112 * t23;
t105 = t112 ^ 2;
t155 = t105 * t99;
t54 = t148 + t109 * (-t115 + t155);
t146 = t109 * qJ(5);
t158 = t112 * pkin(4);
t122 = -t146 - t158;
t157 = t99 * t122;
t18 = -qJDD(4) * pkin(4) - t115 * qJ(5) + t109 * (t29 + t157) + qJDD(5) - t94;
t160 = 2 * qJD(5);
t83 = qJDD(4) + t89;
t152 = t109 * t83;
t149 = t112 * t71;
t147 = t101 * t109;
t144 = t112 * t100;
t142 = t104 + t105;
t135 = t109 * t141;
t131 = t107 * t46 - t108 * t118;
t28 = -t100 * pkin(3) - t99 * pkin(7) + t131;
t117 = t28 - (-t135 + t144) * pkin(4) - t161;
t15 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t147 + t117;
t124 = qJDD(4) * qJ(5) + qJD(4) * t160 + t112 * t157 + t23;
t17 = -t115 * pkin(4) + t124;
t6 = t109 * t18 + t112 * t17;
t2 = t107 * t6 - t108 * t15;
t140 = pkin(2) * t2 - pkin(3) * t15 + pkin(7) * t6;
t4 = t107 * t8 - t108 * t28;
t139 = pkin(2) * t4 - pkin(3) * t28 + pkin(7) * t8;
t87 = -t115 - t155;
t58 = t112 * t87 - t152;
t72 = -0.2e1 * t135 + t144;
t36 = t107 * t58 + t108 * t72;
t138 = pkin(2) * t36 + pkin(3) * t72 + pkin(7) * t58;
t74 = -t107 * t100 - t108 * t99;
t137 = pkin(2) * t74 - t32;
t76 = t142 * t100;
t79 = t142 * t99;
t45 = t107 * t76 + t108 * t79;
t136 = pkin(2) * t45 + pkin(3) * t79 + pkin(7) * t76;
t130 = t109 * t28 - t166;
t129 = -t112 * t28 + t138;
t75 = t108 * t100 - t107 * t99;
t127 = pkin(2) * t75 - t131;
t126 = t109 * (qJ(5) * t79 + t18) + t112 * ((-t115 + t79) * pkin(4) + t124) + t136;
t125 = t136 + t8;
t42 = t109 * t72 + t149;
t121 = t109 * t84 + t112 * t85;
t116 = t147 * t160 - t117;
t120 = pkin(4) * t149 + t109 * (-pkin(4) * t135 + t116 + t161) + t166;
t119 = t146 * t72 + t112 * ((t72 - t135) * pkin(4) + t116) + t138;
t80 = (t104 - t105) * t99;
t57 = t152 + t112 * (t115 - t156);
t55 = t109 * t87 + t112 * t83;
t50 = t71 * t109;
t49 = t72 * t112;
t30 = pkin(1) * (t110 * (-t107 * t79 + t108 * t76) + t113 * t45);
t19 = pkin(1) * (t110 * (-t107 * t72 + t108 * t58) + t113 * t36);
t12 = t107 * t32 - t108 * t131;
t11 = pkin(2) * t12;
t1 = [0, 0, 0, 0, 0, qJDD(1), t133, t123, 0, 0, 0, 0, 0, 0, 0, t100, pkin(1) * (t113 * t100 - t110 * t99) + t47, pkin(1) * (-t110 * t100 - t113 * t99) - t48, 0, pkin(1) * (t110 * t48 + t113 * t47), 0, 0, 0, 0, 0, t100, pkin(1) * (t110 * t74 + t113 * t75) + t127, pkin(1) * (-t110 * t75 + t113 * t74) + t137, 0, pkin(1) * (t110 * (t107 * t131 + t108 * t32) + t113 * t12) + t11, t50, t42, t57, t49, t54, 0, t19 + t129, t130 - t167, t30 + t125, pkin(1) * (t110 * (t107 * t28 + t108 * t8) + t113 * t4) + t139, t50, t57, -t42, 0, -t54, t49, t119 + t19, t30 + t126, t120 + t167, pkin(1) * (t110 * (t107 * t15 + t108 * t6) + t113 * t2) - t15 * t146 - t15 * t158 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t47, -t48, 0, 0, 0, 0, 0, 0, 0, t100, t127, t137, 0, t11, t50, t42, t57, t49, t54, 0, t129, t130, t125, t139, t50, t57, -t42, 0, -t54, t49, t119, t126, t120, t122 * t15 + t140; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, 0, 0, 0, 0, 0, 0, t55, -t121, 0, t109 * t23 - t112 * t22, 0, 0, 0, 0, 0, 0, t55, 0, t121, t109 * t17 - t112 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t89, t80, t145, t89, t144, qJDD(4), -t22, -t23, 0, 0, -t89, t145, -t80, qJDD(4), -t144, t89, pkin(4) * t83 + qJ(5) * t87 - t18, (-pkin(4) * t109 + qJ(5) * t112) * t100, qJ(5) * t84 + (-t115 + t85) * pkin(4) + t124, -pkin(4) * t18 + qJ(5) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t145, -t85, t18;];
tauJ_reg = t1;
