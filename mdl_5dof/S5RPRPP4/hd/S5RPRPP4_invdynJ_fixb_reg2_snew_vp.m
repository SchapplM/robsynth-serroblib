% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:56
% DurationCPUTime: 2.17s
% Computational Cost: add. (2917->216), mult. (6483->254), div. (0->0), fcn. (3971->6), ass. (0->141)
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t119 = cos(pkin(7));
t118 = sin(pkin(7));
t158 = t118 * t122;
t95 = (-t119 * t120 - t158) * qJD(1);
t152 = qJD(1) * t122;
t97 = -qJD(1) * t118 * t120 + t119 * t152;
t177 = t97 * t95;
t192 = qJDD(3) - t177;
t172 = t118 * t192;
t124 = qJD(3) ^ 2;
t94 = t97 ^ 2;
t188 = -t94 - t124;
t29 = -t119 * t188 + t172;
t167 = t119 * t192;
t31 = t118 * t188 + t167;
t17 = t120 * t31 + t122 * t29;
t180 = pkin(6) + pkin(1);
t225 = t180 * t17;
t224 = pkin(3) * t29;
t223 = qJ(4) * t29;
t222 = qJ(4) * t31;
t182 = t95 ^ 2;
t79 = t182 - t124;
t221 = t120 * (t118 * t79 + t167) + t122 * (-t119 * t79 + t172);
t187 = -t94 - t182;
t149 = qJD(1) * qJD(3);
t143 = t122 * t149;
t147 = t120 * qJDD(1);
t101 = -t143 - t147;
t108 = t122 * qJDD(1);
t144 = t120 * t149;
t102 = t108 - t144;
t139 = -t101 * t119 + t102 * t118;
t90 = qJD(3) * t97;
t48 = t139 - t90;
t164 = qJD(3) * t95;
t69 = t101 * t118 + t102 * t119;
t51 = t69 - t164;
t199 = t118 * t51 - t119 * t48;
t200 = -t118 * t48 - t119 * t51;
t209 = t120 * t199 + t122 * t200;
t218 = qJ(2) * t187 - t180 * t209;
t47 = t139 + t90;
t50 = t69 + t164;
t217 = t120 * (-t118 * t47 + t119 * t50) + t122 * (t118 * t50 + t119 * t47);
t216 = pkin(3) * t200;
t215 = t50 * qJ(5);
t214 = qJ(4) * t200;
t191 = qJDD(3) + t177;
t171 = t118 * t191;
t185 = -t182 - t124;
t189 = t119 * t185 - t171;
t53 = t119 * t191;
t190 = t118 * t185 + t53;
t198 = t120 * t189 + t122 * t190;
t211 = qJ(2) * t47 - t180 * t198;
t210 = -pkin(3) * t187 + qJ(4) * t199;
t80 = -t94 + t124;
t208 = -t120 * (t119 * t80 + t171) + t122 * (-t118 * t80 + t53);
t203 = qJ(4) * t189;
t202 = qJ(4) * t190;
t162 = qJD(4) * t97;
t201 = pkin(3) * t190 - 0.2e1 * t162;
t197 = 0.2e1 * qJD(4);
t186 = t94 - t182;
t125 = qJD(1) ^ 2;
t155 = t122 * t125;
t121 = sin(qJ(1));
t123 = cos(qJ(1));
t141 = t121 * g(1) - g(2) * t123;
t136 = qJDD(2) - t141;
t154 = t125 * qJ(2);
t131 = t136 - t154;
t78 = -qJDD(1) * t180 + t131;
t166 = t122 * t78;
t127 = qJDD(3) * pkin(3) - t102 * qJ(4) + t166 + (-pkin(3) * t155 - qJ(4) * t149 + g(3)) * t120;
t133 = qJD(3) * pkin(3) - qJ(4) * t152;
t115 = t120 ^ 2;
t160 = t115 * t125;
t71 = t122 * g(3) - t120 * t78;
t40 = -pkin(3) * t160 + t101 * qJ(4) - qJD(3) * t133 - t71;
t20 = t118 * t127 + t119 * t40 + t197 * t95;
t184 = pkin(4) * t139 - t215;
t113 = qJDD(1) * qJ(2);
t137 = t123 * g(1) + t121 * g(2);
t134 = -t113 + t137;
t183 = -t101 * pkin(3) - (qJ(4) * t115 + t180) * t125 + t133 * t152 + qJDD(4) - t134;
t181 = 2 * qJD(5);
t179 = pkin(4) * t119;
t140 = t118 * t40 - t119 * t127;
t19 = t140 + 0.2e1 * t162;
t5 = t118 * t20 - t119 * t19;
t178 = t122 * t5;
t148 = qJD(2) * qJD(1);
t146 = -0.2e1 * t148;
t41 = t146 - t183;
t175 = t118 * t41;
t169 = t119 * t41;
t165 = qJ(5) * t119;
t161 = qJDD(1) * pkin(1);
t116 = t122 ^ 2;
t159 = t116 * t125;
t145 = t120 * t155;
t157 = t120 * (qJDD(3) + t145);
t156 = t122 * (qJDD(3) - t145);
t153 = t115 + t116;
t151 = qJD(3) * t118;
t150 = qJD(3) * t119;
t142 = -qJ(5) * t118 - pkin(3);
t6 = t118 * t19 + t119 * t20;
t54 = -pkin(4) * t95 - qJ(5) * t97;
t135 = qJDD(3) * qJ(5) + qJD(3) * t181 + t54 * t95 + t20;
t13 = -pkin(4) * t124 + t135;
t132 = -qJDD(3) * pkin(4) - qJ(5) * t124 + qJDD(5) + t140;
t14 = (t197 + t54) * t97 + t132;
t3 = t118 * t13 - t119 * t14;
t1 = t120 * (t118 * t14 + t119 * t13) + t122 * t3;
t70 = t120 * g(3) + t166;
t34 = -t120 * t71 + t122 * t70;
t129 = t122 * (t118 * t139 - t150 * t95) - t120 * (-t119 * t139 - t151 * t95);
t77 = t97 * t151;
t128 = t122 * t77 + (t122 * t119 * t95 - t120 * (t118 * t95 - t119 * t97)) * qJD(3);
t126 = t181 * t97 - t184 + t41;
t110 = 0.2e1 * t148;
t104 = t153 * qJDD(1);
t103 = t108 - 0.2e1 * t144;
t100 = 0.2e1 * t143 + t147;
t91 = -t131 + t161;
t76 = t125 * t180 + t134 + t146;
t73 = -t157 + t122 * (-t124 - t159);
t72 = t120 * (-t124 - t160) + t156;
t23 = t122 * (t119 * t69 - t77) - t120 * (t118 * t69 + t150 * t97);
t15 = t110 + (pkin(4) * qJD(3) - (2 * qJD(5))) * t97 + t183 + t184;
t10 = t126 + (-t47 - t90) * pkin(4);
t9 = -pkin(4) * t90 + t126 + t215;
t8 = -qJ(5) * t187 + t14;
t7 = (-t124 - t187) * pkin(4) + t135;
t2 = t120 * t6 + t178;
t4 = [0, 0, 0, 0, 0, qJDD(1), t141, t137, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t136 - 0.2e1 * t161, t110 + 0.2e1 * t113 - t137, pkin(1) * t91 + qJ(2) * (-t125 * pkin(1) + t110 - t134), (t102 - t144) * t122, -t100 * t122 - t103 * t120, t156 - t120 * (t124 - t159), (-t101 + t143) * t120, t122 * (-t124 + t160) - t157, 0, qJ(2) * t100 - t120 * t76 - t180 * t72, qJ(2) * t103 - t122 * t76 - t180 * t73, t104 * t180 - t153 * t154 - t34, -qJ(2) * t76 - t180 * t34, t23, -t217, t208, t129, -t221, t128, t122 * (-t175 - t202) - t120 * (-pkin(3) * t47 + t169 + t203) + t211, t122 * (-t169 + t223) - t120 * (-pkin(3) * t50 - t175 - t222) + qJ(2) * t50 + t225, t122 * (-t5 - t214) - t120 * (t210 + t6) + t218, -qJ(4) * t178 - t120 * (pkin(3) * t41 + qJ(4) * t6) - qJ(2) * t41 - t180 * t2, t23, t208, t217, t128, t221, t129, t122 * (-t10 * t118 - t165 * t47 - t202) - t120 * (t119 * t10 + t142 * t47 + t203) + t211, t122 * (-t118 * t7 + t119 * t8 - t214) - t120 * (t118 * t8 + t119 * t7 + t210) + t218, t122 * (t119 * t9 - t223) - t120 * (t118 * t9 + t222) - (pkin(4) * t158 - t120 * (-pkin(3) - t179) + qJ(2)) * t50 - t225, (t122 * (pkin(4) * t118 - t165) - t120 * (t142 - t179) + qJ(2)) * t15 + (-t180 - qJ(4)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t125, -t91, 0, 0, 0, 0, 0, 0, t72, t73, -t104, t34, 0, 0, 0, 0, 0, 0, t198, -t17, t209, t2, 0, 0, 0, 0, 0, 0, t198, t209, t17, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, (-t115 + t116) * t125, t108, -t145, -t147, qJDD(3), t70, t71, 0, 0, -t177, t186, t51, t177, -t48, qJDD(3), -t140 + t201, -t20 - t224, t216, pkin(3) * t5, -t177, t51, -t186, qJDD(3), t48, t177, pkin(4) * t191 + qJ(5) * t185 - t54 * t97 - t132 + t201, -pkin(4) * t51 - qJ(5) * t48 + t216, t224 + qJ(5) * t192 + (-t124 - t188) * pkin(4) + t135, pkin(3) * t3 - pkin(4) * t14 + qJ(5) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t50, t187, -t41, 0, 0, 0, 0, 0, 0, t47, t187, -t50, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, t51, t188, t14;];
tauJ_reg = t4;
