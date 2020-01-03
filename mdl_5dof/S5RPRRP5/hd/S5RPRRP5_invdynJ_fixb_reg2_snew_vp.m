% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:58
% EndTime: 2019-12-31 18:41:01
% DurationCPUTime: 0.92s
% Computational Cost: add. (2185->147), mult. (3267->185), div. (0->0), fcn. (1752->8), ass. (0->97)
t100 = cos(pkin(8));
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t92 = qJD(1) + qJD(3);
t90 = t92 ^ 2;
t80 = t104 * t90 * t101;
t75 = qJDD(4) - t80;
t132 = t104 * t75;
t107 = qJD(4) ^ 2;
t95 = t101 ^ 2;
t145 = t95 * t90;
t76 = t107 + t145;
t50 = -t101 * t76 + t132;
t129 = qJD(4) * t92;
t91 = qJDD(1) + qJDD(3);
t136 = t101 * t91;
t60 = 0.2e1 * t104 * t129 + t136;
t32 = t102 * t50 + t105 * t60;
t99 = sin(pkin(8));
t160 = pkin(1) * (t99 * (-t102 * t60 + t105 * t50) + t100 * t32) + pkin(2) * t32;
t156 = pkin(3) * t60 + pkin(7) * t50;
t154 = qJ(5) * t60;
t128 = t101 * qJ(5);
t147 = t104 * pkin(4);
t116 = -t128 - t147;
t146 = t116 * t90;
t151 = 2 * qJD(5);
t103 = sin(qJ(1));
t106 = cos(qJ(1));
t123 = t103 * g(1) - t106 * g(2);
t72 = qJDD(1) * pkin(1) + t123;
t108 = qJD(1) ^ 2;
t118 = t106 * g(1) + t103 * g(2);
t73 = -t108 * pkin(1) - t118;
t117 = t100 * t72 - t99 * t73;
t111 = qJDD(1) * pkin(2) + t117;
t141 = t100 * t73 + t99 * t72;
t39 = -t108 * pkin(2) + t141;
t29 = t102 * t111 + t105 * t39;
t26 = -t90 * pkin(3) + t91 * pkin(7) + t29;
t97 = -g(3) + qJDD(2);
t20 = t101 * t97 + t104 * t26;
t119 = qJDD(4) * qJ(5) + (qJD(4) * t151) + t104 * t146 + t20;
t85 = t104 * t97;
t15 = -qJDD(4) * pkin(4) - t107 * qJ(5) + (t26 + t146) * t101 + qJDD(5) - t85;
t96 = t104 ^ 2;
t140 = t95 + t96;
t70 = t140 * t90;
t153 = t104 * ((-t107 + t70) * pkin(4) + t119) + t101 * (qJ(5) * t70 + t15);
t74 = qJDD(4) + t80;
t139 = t101 * t74;
t144 = t96 * t90;
t78 = -t107 - t144;
t49 = t104 * t78 - t139;
t125 = t101 * t129;
t131 = t104 * t91;
t61 = -0.2e1 * t125 + t131;
t31 = t102 * t49 + t105 * t61;
t152 = pkin(1) * (t99 * (-t102 * t61 + t105 * t49) + t100 * t31) + pkin(2) * t31;
t19 = t101 * t26 - t85;
t6 = t101 * t19 + t104 * t20;
t45 = t132 + t101 * (-t107 + t144);
t28 = -t102 * t39 + t105 * t111;
t25 = -t91 * pkin(3) - t90 * pkin(7) - t28;
t149 = -pkin(3) * t25 + pkin(7) * t6;
t110 = t25 - (-t125 + t131) * pkin(4) - t154;
t135 = t101 * t92;
t12 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t135 + t110;
t14 = -t107 * pkin(4) + t119;
t4 = t101 * t15 + t104 * t14;
t148 = -pkin(3) * t12 + pkin(7) * t4;
t143 = pkin(3) * t61 + pkin(7) * t49;
t65 = t140 * t91;
t142 = pkin(3) * t70 + pkin(7) * t65;
t133 = t104 * t60;
t127 = t101 * t25 - t156;
t126 = -t104 * t25 + t143;
t38 = t102 * t65 + t105 * t70;
t121 = pkin(1) * (t99 * (-t102 * t70 + t105 * t65) + t100 * t38) + pkin(2) * t38 + t142;
t34 = t101 * t61 + t133;
t115 = t101 * t75 + t104 * t76;
t68 = -t102 * t91 - t105 * t90;
t114 = t102 * t90 - t105 * t91;
t109 = t135 * t151 - t110;
t113 = pkin(4) * t133 + t101 * (-pkin(4) * t125 + t109 + t154) + t156;
t112 = t61 * t128 + t143 + t104 * ((t61 - t125) * pkin(4) + t109);
t71 = (t95 - t96) * t90;
t48 = t139 + t104 * (t107 - t145);
t46 = t101 * t78 + t104 * t74;
t41 = t60 * t101;
t40 = t61 * t104;
t9 = t102 * t29 + t105 * t28;
t2 = t102 * t6 - t105 * t25;
t1 = t102 * t4 - t105 * t12;
t3 = [0, 0, 0, 0, 0, qJDD(1), t123, t118, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t100 * qJDD(1) - t99 * t108) + t117, pkin(1) * (-t99 * qJDD(1) - t100 * t108) - t141, 0, pkin(1) * (t100 * t117 + t99 * t141), 0, 0, 0, 0, 0, t91, pkin(1) * (-t100 * t114 + t99 * t68) - pkin(2) * t114 + t28, pkin(1) * (t100 * t68 + t99 * t114) + pkin(2) * t68 - t29, 0, pkin(1) * (t99 * (-t102 * t28 + t105 * t29) + t100 * t9) + pkin(2) * t9, t41, t34, t48, t40, t45, 0, t126 + t152, t127 - t160, t121 + t6, pkin(1) * (t99 * (t102 * t25 + t105 * t6) + t100 * t2) + pkin(2) * t2 + t149, t41, t48, -t34, 0, -t45, t40, t112 + t152, t121 + t153, t113 + t160, pkin(1) * (t99 * (t102 * t12 + t105 * t4) + t100 * t1) + pkin(2) * t1 - t12 * t128 - t12 * t147 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, t46, -t115, 0, t101 * t20 - t104 * t19, 0, 0, 0, 0, 0, 0, t46, 0, t115, t101 * t14 - t104 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t28, -t29, 0, 0, t41, t34, t48, t40, t45, 0, t126, t127, t142 + t6, t149, t41, t48, -t34, 0, -t45, t40, t112, t142 + t153, t113, t116 * t12 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t71, t136, t80, t131, qJDD(4), -t19, -t20, 0, 0, -t80, t136, -t71, qJDD(4), -t131, t80, pkin(4) * t74 + qJ(5) * t78 - t15, (-pkin(4) * t101 + qJ(5) * t104) * t91, qJ(5) * t75 + (-t107 + t76) * pkin(4) + t119, -pkin(4) * t15 + qJ(5) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t136, -t76, t15;];
tauJ_reg = t3;
