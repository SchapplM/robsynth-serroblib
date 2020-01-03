% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:45:37
% EndTime: 2019-12-31 18:45:44
% DurationCPUTime: 2.13s
% Computational Cost: add. (3754->250), mult. (7264->305), div. (0->0), fcn. (4517->8), ass. (0->168)
t143 = cos(qJ(4));
t140 = sin(qJ(4));
t141 = sin(qJ(3));
t182 = qJD(1) * qJD(3);
t126 = t141 * t182;
t144 = cos(qJ(3));
t180 = t144 * qJDD(1);
t110 = -t126 + t180;
t101 = -qJDD(4) + t110;
t184 = qJD(1) * t141;
t102 = -t143 * qJD(3) + t140 * t184;
t104 = qJD(3) * t140 + t143 * t184;
t190 = t104 * t102;
t151 = t101 - t190;
t202 = t140 * t151;
t100 = t104 ^ 2;
t122 = qJD(1) * t144 - qJD(4);
t210 = t122 ^ 2;
t219 = -t100 - t210;
t30 = -t143 * t219 - t202;
t256 = pkin(2) * t30;
t255 = pkin(3) * t30;
t254 = pkin(7) * t30;
t196 = t143 * t151;
t37 = -t140 * t219 + t196;
t253 = pkin(7) * t37;
t137 = cos(pkin(8));
t252 = t137 * t30;
t251 = t141 * t37;
t250 = t144 * t37;
t175 = t144 * t182;
t181 = t141 * qJDD(1);
t109 = t175 + t181;
t168 = t143 * qJDD(3) - t109 * t140;
t155 = qJD(4) * t104 - t168;
t89 = t104 * t122;
t46 = t155 + t89;
t211 = t102 ^ 2;
t82 = t211 - t210;
t249 = t141 * (t143 * t82 + t202) + t144 * t46;
t157 = -qJDD(3) * t140 - t109 * t143;
t150 = -qJD(4) * t102 - t157;
t191 = t102 * t122;
t215 = t150 + t191;
t204 = t140 * t215;
t218 = t100 - t211;
t220 = t155 - t89;
t248 = t141 * (t143 * t220 + t204) + t144 * t218;
t136 = sin(pkin(8));
t60 = t101 + t190;
t195 = t143 * t60;
t214 = -t210 - t211;
t223 = t140 * t214 - t195;
t201 = t140 * t60;
t222 = t143 * t214 + t201;
t235 = t141 * t220 + t144 * t222;
t247 = pkin(1) * (t136 * t235 - t137 * t223) + pkin(6) * t235 - pkin(2) * t223;
t245 = pkin(3) * t223;
t244 = pkin(7) * t222;
t243 = pkin(7) * t223;
t83 = -t100 + t210;
t242 = t143 * t83 - t201;
t241 = qJ(5) * t215;
t238 = t140 * t82 - t196;
t236 = t141 * t222 - t144 * t220;
t216 = t150 - t191;
t234 = t141 * (-t140 * t83 - t195) - t144 * t216;
t217 = t100 + t211;
t233 = pkin(3) * t217;
t185 = -g(3) + qJDD(2);
t125 = t144 * t185;
t146 = qJD(1) ^ 2;
t164 = -pkin(3) * t144 - pkin(7) * t141;
t142 = sin(qJ(1));
t145 = cos(qJ(1));
t173 = t142 * g(1) - g(2) * t145;
t105 = qJDD(1) * pkin(1) + t173;
t161 = g(1) * t145 + g(2) * t142;
t106 = -pkin(1) * t146 - t161;
t207 = t136 * t105 + t137 * t106;
t66 = -pkin(2) * t146 + qJDD(1) * pkin(6) + t207;
t169 = t146 * t164 + t66;
t209 = qJD(3) ^ 2;
t40 = -qJDD(3) * pkin(3) - t209 * pkin(7) + t141 * t169 - t125;
t232 = pkin(4) * t155 - t241 + t40;
t230 = t141 * t217;
t226 = t144 * t217;
t221 = -t140 * t220 + t143 * t215;
t159 = -t110 + t126;
t160 = t109 + t175;
t171 = t137 * t105 - t136 * t106;
t65 = -qJDD(1) * pkin(2) - t146 * pkin(6) - t171;
t34 = pkin(3) * t159 - pkin(7) * t160 + t65;
t170 = t141 * t185;
t41 = -t209 * pkin(3) + qJDD(3) * pkin(7) + t169 * t144 + t170;
t18 = t140 * t34 + t143 * t41;
t72 = pkin(4) * t102 - qJ(5) * t104;
t172 = t101 * qJ(5) + t102 * t72 - t18;
t212 = -(t210 + t219) * pkin(4) - qJ(5) * t151 - t172;
t208 = pkin(4) * t143;
t206 = t140 * t40;
t203 = t140 * t216;
t199 = t143 * t40;
t197 = t143 * t216;
t192 = qJ(5) * t143;
t189 = t122 * t140;
t188 = t122 * t143;
t121 = t144 * t146 * t141;
t114 = qJDD(3) + t121;
t187 = t141 * t114;
t115 = qJDD(3) - t121;
t186 = t144 * t115;
t183 = qJD(5) * t122;
t179 = t102 * t188;
t178 = t144 * t190;
t176 = pkin(1) * t136 + pkin(6);
t174 = -qJ(5) * t140 - pkin(3);
t17 = t140 * t41 - t143 * t34;
t6 = t140 * t17 + t143 * t18;
t54 = t141 * t66 - t125;
t55 = t144 * t66 + t170;
t27 = t141 * t54 + t144 * t55;
t116 = -0.2e1 * t183;
t167 = t116 - t172;
t81 = t104 * t189;
t166 = t141 * (t143 * t150 + t81) - t178;
t165 = -t102 * t189 - t143 * t155;
t11 = -pkin(4) * t210 + t167;
t12 = t101 * pkin(4) - qJ(5) * t210 + t104 * t72 + qJDD(5) + t17;
t163 = -pkin(4) * t12 + qJ(5) * t11;
t162 = -pkin(4) * t216 - qJ(5) * t46;
t158 = t140 * t18 - t143 * t17;
t154 = (t102 * t140 + t104 * t143) * t122;
t153 = t141 * (-t81 + t179) + t144 * t101;
t152 = -pkin(1) * t137 - pkin(2) + t164;
t149 = t141 * (t140 * t155 - t179) + t178;
t148 = 0.2e1 * qJD(5) * t104 - t232;
t147 = -pkin(4) * t60 + qJ(5) * t214 - t12;
t133 = t144 ^ 2;
t132 = t141 ^ 2;
t130 = t133 * t146;
t128 = t132 * t146;
t119 = -t130 - t209;
t118 = -t128 - t209;
t113 = t128 + t130;
t112 = (t132 + t133) * qJDD(1);
t111 = -0.2e1 * t126 + t180;
t108 = 0.2e1 * t175 + t181;
t80 = -t118 * t141 - t186;
t79 = t119 * t144 - t187;
t52 = (qJD(4) - t122) * t102 + t157;
t47 = (-qJD(4) - t122) * t104 + t168;
t43 = -t104 * t188 + t140 * t150;
t26 = t143 * t47 + t203;
t25 = -t143 * t46 + t203;
t23 = -t140 * t46 - t197;
t21 = -t141 * t52 + t250;
t19 = -t141 * t215 - t250;
t14 = t144 * t25 - t230;
t13 = (-pkin(4) * t122 - 0.2e1 * qJD(5)) * t104 + t232;
t10 = (-t220 + t89) * pkin(4) + t148;
t9 = pkin(4) * t89 + t148 + t241;
t8 = qJ(5) * t217 + t12;
t7 = (-t210 + t217) * pkin(4) + t167;
t3 = t11 * t143 + t12 * t140;
t2 = t11 * t140 - t12 * t143;
t1 = t13 * t141 + t144 * t3;
t4 = [0, 0, 0, 0, 0, qJDD(1), t173, t161, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t137 - t136 * t146) + t171, pkin(1) * (-qJDD(1) * t136 - t137 * t146) - t207, 0, pkin(1) * (t136 * t207 + t137 * t171), t160 * t141, t108 * t144 + t111 * t141, t187 + t144 * (-t128 + t209), -t159 * t144, t141 * (t130 - t209) + t186, 0, -t144 * t65 + pkin(2) * t111 + pkin(6) * t79 + pkin(1) * (t111 * t137 + t136 * t79), t141 * t65 - pkin(2) * t108 + pkin(6) * t80 + pkin(1) * (-t108 * t137 + t136 * t80), pkin(2) * t113 + pkin(6) * t112 + pkin(1) * (t112 * t136 + t113 * t137) + t27, -pkin(2) * t65 + pkin(6) * t27 + pkin(1) * (t136 * t27 - t137 * t65), t166, -t248, t234, t149, t249, t153, t141 * (t206 - t243) + t144 * (t17 - t245) + t247, t141 * (t199 + t254) + t144 * (t18 + t255) + t256 + pkin(6) * t21 + pkin(1) * (t136 * t21 + t252), -t141 * t158 + t176 * (t144 * t26 - t230) + t152 * (t140 * t47 - t197), t176 * (t141 * t40 + t144 * t6) + t152 * t158, t166, t234, t248, t153, -t249, t149, t141 * (-t10 * t140 - t192 * t220 - t243) + t144 * (-t147 - t245) + t247, t141 * (-pkin(7) * t23 - t140 * t7 + t143 * t8) + t144 * (-pkin(3) * t23 - t162) - pkin(2) * t23 + pkin(6) * t14 + pkin(1) * (t136 * t14 - t137 * t23), t141 * (-pkin(4) * t204 + t143 * t9 - t254) + t144 * (0.2e1 * t183 - t212 - t255) - t256 + pkin(6) * t19 + pkin(1) * (t136 * t19 - t252), t141 * (-pkin(7) * t2 + (pkin(4) * t140 - t192) * t13) + t144 * (-pkin(3) * t2 - t163) - pkin(2) * t2 + pkin(6) * t1 + pkin(1) * (t1 * t136 - t137 * t2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, 0, 0, 0, 0, 0, 0, t114 * t144 + t119 * t141, -t115 * t141 + t118 * t144, 0, t141 * t55 - t144 * t54, 0, 0, 0, 0, 0, 0, t236, t144 * t52 + t251, t141 * t26 + t226, t141 * t6 - t144 * t40, 0, 0, 0, 0, 0, 0, t236, t141 * t25 + t226, t144 * t215 - t251, -t13 * t144 + t141 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t128 - t130, t181, t121, t180, qJDD(3), -t54, -t55, 0, 0, t43, t221, t242, t165, t238, t154, -pkin(3) * t220 - t199 + t244, pkin(3) * t52 + t206 + t253, pkin(7) * t26 + t233 + t6, -pkin(3) * t40 + pkin(7) * t6, t43, t242, -t221, t154, -t238, t165, t10 * t143 + t174 * t220 + t244, pkin(7) * t25 + t140 * t8 + t143 * t7 + t233, -t253 + t140 * t9 + (pkin(3) + t208) * t215, pkin(7) * t3 + (t174 - t208) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, t218, t216, -t190, -t46, -t101, -t17, -t18, 0, 0, t190, t216, -t218, -t101, t46, -t190, t147, t162, t116 + t212, t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t216, t219, t12;];
tauJ_reg = t4;
