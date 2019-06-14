% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 16:43
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPPR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:42:45
% EndTime: 2019-05-05 16:42:53
% DurationCPUTime: 2.08s
% Computational Cost: add. (4353->271), mult. (8944->319), div. (0->0), fcn. (4606->8), ass. (0->158)
t134 = sin(pkin(9));
t135 = cos(pkin(9));
t145 = qJD(3) ^ 2;
t140 = sin(qJ(3));
t128 = t140 ^ 2;
t146 = qJD(1) ^ 2;
t200 = t128 * t146;
t108 = t145 + t200;
t143 = cos(qJ(3));
t223 = t143 * t146;
t113 = t140 * t223;
t105 = -t113 + qJDD(3);
t194 = t143 * t105;
t64 = -t140 * t108 + t194;
t188 = qJD(1) * qJD(3);
t174 = t143 * t188;
t186 = t140 * qJDD(1);
t94 = 0.2e1 * t174 + t186;
t183 = pkin(1) * (t134 * t64 + t135 * t94) + pkin(2) * t94 + pkin(7) * t64;
t95 = t174 + t186;
t228 = t95 + t174;
t139 = sin(qJ(6));
t142 = cos(qJ(6));
t191 = qJD(1) * t143;
t88 = -t142 * qJD(3) + t139 * t191;
t89 = t139 * qJD(3) + t142 * t191;
t55 = t88 * t89;
t82 = qJDD(6) + t95;
t225 = -t55 + t82;
t227 = t139 * t225;
t226 = t142 * t225;
t104 = t113 + qJDD(3);
t224 = pkin(3) * t108 + qJ(4) * t105;
t176 = pkin(1) * t134 + pkin(7);
t222 = t176 - qJ(5);
t129 = t143 ^ 2;
t199 = t129 * t146;
t62 = t194 + t140 * (-t145 + t199);
t190 = t140 * qJD(1);
t103 = -qJD(3) * pkin(4) - qJ(5) * t190;
t221 = -t103 * t190 - qJDD(5);
t122 = t140 * t188;
t185 = t143 * qJDD(1);
t96 = -t122 + t185;
t217 = t96 * pkin(4);
t220 = -t217 + t221;
t130 = -g(3) + qJDD(2);
t121 = t143 * t130;
t171 = qJDD(3) * pkin(3) + t145 * qJ(4) - qJDD(4) + t121;
t141 = sin(qJ(1));
t144 = cos(qJ(1));
t168 = t141 * g(1) - t144 * g(2);
t157 = qJDD(1) * pkin(1) + t168;
t169 = t144 * g(1) + t141 * g(2);
t92 = -t146 * pkin(1) - t169;
t215 = t134 * t157 + t135 * t92;
t47 = -t146 * pkin(2) + qJDD(1) * pkin(7) + t215;
t198 = t140 * qJ(4);
t165 = -t143 * pkin(3) - t198;
t91 = t165 * qJD(1);
t172 = qJD(1) * t91 + t47;
t28 = t172 * t140 - t171;
t80 = t88 ^ 2;
t81 = t89 ^ 2;
t114 = qJD(6) + t190;
t111 = t114 ^ 2;
t219 = pkin(3) + pkin(4);
t218 = pkin(4) + pkin(8);
t216 = pkin(5) + qJ(4);
t214 = t114 * t88;
t170 = pkin(5) * t140 + pkin(8) * t143;
t187 = qJD(4) * qJD(3);
t124 = 0.2e1 * t187;
t38 = t140 * t130 + t143 * t47;
t167 = t145 * pkin(3) - qJDD(3) * qJ(4) - t91 * t191 - t38;
t179 = 0.2e1 * qJD(5) * qJD(1);
t149 = pkin(4) * t199 - qJD(3) * t103 + t143 * t179 + t167;
t19 = -t96 * qJ(5) + t124 - t149;
t14 = qJDD(3) * pkin(5) - t145 * pkin(8) - t170 * t223 + t19;
t212 = t139 * t14;
t45 = t55 + t82;
t211 = t139 * t45;
t161 = t142 * qJDD(3) - t139 * t96;
t154 = (qJD(6) - t114) * t89 - t161;
t49 = t88 * qJD(6) - t139 * qJDD(3) - t142 * t96;
t35 = t49 - t214;
t18 = t139 * t35 + t142 * t154;
t210 = t140 * t18;
t209 = t140 * t47;
t208 = t142 * t14;
t207 = t142 * t45;
t206 = t143 * t94;
t192 = t128 + t129;
t101 = t192 * t146;
t205 = qJ(4) * t101;
t109 = t145 + t199;
t204 = qJ(4) * t109;
t203 = qJ(4) * t143;
t202 = t114 * t139;
t201 = t114 * t142;
t197 = t140 * t104;
t193 = -0.2e1 * qJD(5) + t91;
t189 = pkin(3) + t218;
t184 = -t81 - t111;
t163 = t143 * t109 + t197;
t97 = -0.2e1 * t122 + t185;
t182 = pkin(1) * (-t134 * t163 + t135 * t97) - pkin(7) * t163 + pkin(2) * t97;
t100 = t192 * qJDD(1);
t181 = pkin(1) * (t134 * t100 + t135 * t101) + pkin(7) * t100 + pkin(2) * t101;
t180 = t140 * t55;
t178 = qJ(5) * t199;
t177 = pkin(1) * t135 + pkin(2);
t175 = qJ(5) * qJD(3) * t143;
t119 = 0.2e1 * qJD(4) * t190;
t173 = -t134 * t92 + t135 * t157;
t46 = -qJDD(1) * pkin(2) - t146 * pkin(7) - t173;
t156 = -t96 * pkin(3) - t228 * qJ(4) + t46;
t152 = t156 + t221;
t10 = -t178 + t95 * pkin(5) + t119 + t218 * t96 + (pkin(5) * t143 + (-pkin(3) - pkin(8)) * t140) * t188 - t152;
t158 = -pkin(4) * t104 - t95 * qJ(5) - t171;
t153 = t158 + t209;
t15 = -t145 * pkin(5) - qJDD(3) * pkin(8) + (t175 + (-t170 * qJD(1) + t193) * t140) * qJD(1) + t153;
t5 = -t142 * t10 + t139 * t15;
t37 = -t121 + t209;
t21 = t140 * t37 + t143 * t38;
t6 = t139 * t10 + t142 * t15;
t2 = t139 * t6 - t142 * t5;
t3 = t139 * t5 + t142 * t6;
t52 = t140 * t97 + t206;
t63 = t143 * t104 - t140 * t109;
t60 = t140 * t105 + t143 * t108;
t162 = qJD(1) * (pkin(3) * qJD(3) - 0.2e1 * qJD(4));
t25 = t124 - t167;
t160 = -t143 * t219 - t177;
t159 = t49 + t214;
t155 = t119 - t156;
t150 = t140 * t216 + t143 * t189 + t177;
t148 = t140 * t162 + t156;
t147 = -pkin(3) * t122 + qJ(4) * t94 + t155;
t20 = (t193 * t140 + t175) * qJD(1) + t153;
t102 = (t128 - t129) * t146;
t69 = -t81 + t111;
t68 = t80 - t111;
t61 = t228 * t140;
t59 = t197 + t143 * (t145 - t200);
t58 = (t96 - t122) * t143;
t54 = t81 - t80;
t50 = -t111 - t80;
t48 = t89 * qJD(6) - t161;
t43 = -t80 - t81;
t30 = (-qJD(6) - t114) * t89 + t161;
t27 = -t139 * t184 - t207;
t26 = t142 * t184 - t211;
t24 = t142 * t50 - t227;
t23 = t139 * t50 + t226;
t17 = t139 * t154 - t142 * t35;
t16 = t148 + t178 + t220;
t1 = [0, 0, 0, 0, 0, qJDD(1), t168, t169, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t135 * qJDD(1) - t134 * t146) + t173, pkin(1) * (-t134 * qJDD(1) - t135 * t146) - t215, 0, pkin(1) * (t134 * t215 + t135 * t173), t61, t52, t59, t58, t62, 0, -t143 * t46 + t182, t140 * t46 - t183, t21 + t181, -pkin(2) * t46 + pkin(7) * t21 + pkin(1) * (t134 * t21 - t135 * t46), t61, t59, -t52, 0, -t62, t58, t97 * t198 + t143 * ((t97 - t122) * pkin(3) + t155) + t182, t143 * (pkin(3) * t101 + t25) + (t28 + t205) * t140 + t181, pkin(3) * t206 + t140 * t147 + t183, t176 * (t140 * t28 + t143 * t25) + (t165 - t177) * t148, t58, t52, t62, t61, t59, 0, t140 * (t147 + (t108 - t199) * qJ(5) - t220) + t143 * (-qJ(5) * t105 + t219 * t94) + t183, t143 * (-t217 + (-t109 + t199) * qJ(5) + t152) + t176 * t163 + t160 * t97 + (-qJ(4) * t97 - qJ(5) * t104 + t143 * t162) * t140, t143 * (-0.2e1 * t187 + (t96 + t185) * qJ(5) + t149) - t176 * t100 + t160 * t101 + (-qJ(5) * t174 - t158 - t205 + (qJ(5) * qJDD(1) - t172 + t179) * t140) * t140, (t160 - t198) * t16 + t222 * (t140 * t20 + t143 * t19), t180 + t143 * (-t142 * t49 - t89 * t202), t140 * t54 + t143 * (t139 * t159 + t142 * t30), t140 * t35 + t143 * (t139 * t69 - t226), -t180 + t143 * (t139 * t48 + t88 * t201), t140 * t154 + t143 * (-t142 * t68 + t211), t140 * t82 + t143 * (t139 * t89 - t142 * t88) * t114, t140 * (-qJ(5) * t24 - t5) + t143 * (-qJ(5) * t30 - t212) + t176 * (t140 * t24 + t143 * t30) + t150 * t23, t140 * (-qJ(5) * t27 - t6) + t143 * (-qJ(5) * t159 - t208) + t176 * (t140 * t27 + t143 * t159) + t150 * t26, -qJ(5) * t210 + t143 * (-qJ(5) * t43 + t2) + t176 * (t143 * t43 + t210) + t150 * t17, t150 * t2 + t222 * (t143 * t14 + t140 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, 0, 0, 0, 0, 0, 0, t63, -t60, 0, t140 * t38 - t143 * t37, 0, 0, 0, 0, 0, 0, t63, 0, t60, t140 * t25 - t143 * t28, 0, 0, 0, 0, 0, 0, t60, -t63, 0, t140 * t19 - t143 * t20, 0, 0, 0, 0, 0, 0, t140 * t30 - t143 * t24, t140 * t159 - t143 * t27, t140 * t43 - t143 * t18, t140 * t14 - t143 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t102, t186, t113, t185, qJDD(3), -t37, -t38, 0, 0, -t113, t186, -t102, qJDD(3), -t185, t113, pkin(3) * t104 - t204 - t28, (-pkin(3) * t140 + t203) * qJDD(1), t25 + t224, -pkin(3) * t28 + qJ(4) * t25, t113, t102, t185, -t113, t186, qJDD(3), pkin(4) * t108 + t19 + t224, -t219 * t104 + t20 + t204, (t219 * t140 - t203) * qJDD(1), qJ(4) * t19 - t219 * t20, -t139 * t49 + t89 * t201, t139 * t30 - t142 * t159, -t142 * t69 - t227, -t142 * t48 + t88 * t202, -t139 * t68 - t207, (-t139 * t88 - t142 * t89) * t114, -t189 * t24 + t216 * t30 + t208, t159 * t216 - t189 * t27 - t212, -t189 * t18 + t216 * t43 - t3, t216 * t14 - t189 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t186, -t108, t28, 0, 0, 0, 0, 0, 0, -t108, t104, -t186, t20, 0, 0, 0, 0, 0, 0, t24, t27, t18, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, -t97, -t101, -t16, 0, 0, 0, 0, 0, 0, t23, t26, t17, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t54, t35, -t55, t154, t82, -t5, -t6, 0, 0;];
tauJ_reg  = t1;
