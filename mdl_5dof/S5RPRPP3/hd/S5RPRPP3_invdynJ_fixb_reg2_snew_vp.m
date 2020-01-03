% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:54
% EndTime: 2019-12-31 18:13:01
% DurationCPUTime: 2.65s
% Computational Cost: add. (3112->275), mult. (7834->310), div. (0->0), fcn. (5231->6), ass. (0->153)
t148 = sin(pkin(7));
t149 = cos(pkin(7));
t153 = qJD(3) ^ 2;
t152 = cos(qJ(3));
t188 = qJD(1) * t152;
t150 = sin(qJ(3));
t189 = qJD(1) * t150;
t125 = t148 * t189 - t149 * t188;
t222 = t125 ^ 2;
t106 = t222 + t153;
t127 = t148 * t188 + t149 * t189;
t193 = t127 * t125;
t228 = qJDD(3) - t193;
t63 = t152 * t228;
t249 = -t106 * t150 + t63;
t205 = t150 * t228;
t33 = t106 * t152 + t205;
t265 = qJ(2) * (t148 * t249 + t149 * t33);
t117 = t127 * qJD(3);
t184 = t149 * qJDD(1);
t185 = t148 * qJDD(1);
t123 = t150 * t185 - t152 * t184;
t83 = t123 + 0.2e1 * t117;
t268 = -pkin(1) * t83 - t265;
t264 = pkin(6) * t33;
t267 = -pkin(2) * t83 - t264;
t221 = t127 ^ 2;
t230 = t221 + t153;
t227 = qJDD(3) + t193;
t240 = t150 * t227;
t235 = t152 * t230 + t240;
t239 = t152 * t227;
t45 = t150 * t230 - t239;
t266 = qJ(2) * (t148 * t235 + t149 * t45);
t255 = pkin(6) * t45;
t256 = pkin(6) * t235;
t109 = t221 - t153;
t248 = t148 * (t109 * t150 + t63) + t149 * (-t109 * t152 + t205);
t261 = pkin(6) * t249;
t103 = t222 - t153;
t257 = t148 * (-t103 * t152 + t240) + t149 * (-t103 * t150 - t239);
t252 = qJ(4) * t106;
t116 = t125 * qJD(3);
t191 = t148 * t152;
t124 = (t149 * t150 + t191) * qJDD(1);
t86 = t124 - t116;
t174 = -t86 + t116;
t233 = t174 * qJ(4);
t154 = qJD(1) ^ 2;
t151 = sin(qJ(1));
t217 = cos(qJ(1));
t168 = g(1) * t217 + g(2) * t151;
t238 = -pkin(1) * t154 + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t168;
t246 = 2 * qJD(4);
t177 = g(1) * t151 - g(2) * t217;
t170 = -qJDD(2) + t177;
t142 = t148 ^ 2;
t143 = t149 ^ 2;
t190 = t142 + t143;
t213 = t149 * pkin(2);
t75 = (pkin(6) * t190 + qJ(2)) * t154 + (pkin(1) + t213) * qJDD(1) + t170;
t159 = -pkin(3) * t117 + t127 * t246 + t75;
t158 = t159 - t233;
t247 = -pkin(1) * t174 - t266;
t231 = -t221 - t222;
t245 = pkin(1) * t231;
t244 = pkin(2) * t231;
t241 = qJ(4) * t231;
t237 = (-t153 - t231) * pkin(3);
t234 = qJ(4) * t227 + (-t153 + t230) * pkin(3);
t229 = t222 - t221;
t169 = pkin(4) * t127 - qJD(3) * qJ(5);
t220 = -2 * qJD(5);
t84 = t123 + t117;
t226 = -qJ(5) * t84 + t125 * t220 + t127 * t169;
t196 = qJDD(1) * pkin(1);
t197 = qJ(2) * t154;
t118 = t170 + t196 + t197;
t192 = t143 * t154;
t225 = qJ(2) * t192 + t142 * t197 - t118 - t196;
t99 = t150 * t117;
t25 = t148 * (t152 * t86 - t99) + t149 * (t117 * t152 + t150 * t86);
t218 = pkin(3) * t84;
t216 = pkin(3) * t152;
t215 = pkin(3) * t153;
t214 = g(3) * t149;
t212 = pkin(3) + qJ(5);
t162 = -t214 + (-pkin(6) * qJDD(1) + t154 * t213 - t238) * t148;
t173 = -g(3) * t148 + t149 * t238;
t58 = -pkin(2) * t192 + pkin(6) * t184 + t173;
t28 = t150 * t162 + t152 * t58;
t27 = t150 * t58 - t152 * t162;
t18 = t150 * t28 - t152 * t27;
t210 = t148 * t18;
t209 = t150 * t174;
t207 = t150 * t75;
t204 = t150 * t83;
t203 = t152 * t174;
t201 = t152 * t75;
t199 = t152 * t83;
t198 = t222 - t230;
t194 = t124 * t150;
t178 = qJ(4) * t150 + pkin(2);
t176 = t148 * (t148 * t238 + t214) + t149 * t173;
t19 = t150 * t27 + t152 * t28;
t70 = pkin(3) * t125 - qJ(4) * t127;
t172 = qJDD(3) * qJ(4) + qJD(3) * t246 - t125 * t70 + t28;
t166 = t149 * (-t125 * t150 - t127 * t152);
t165 = qJD(3) * t169 + qJDD(5) + t172;
t164 = -pkin(4) * t84 + t165;
t22 = -qJDD(3) * pkin(3) - qJ(4) * t153 + t127 * t70 + qJDD(4) + t27;
t163 = t86 * pkin(4) - qJ(5) * t228 + t22;
t161 = qJD(3) * t220 + t163;
t24 = t148 * (t116 * t152 + t150 * t84) + t149 * (t116 * t150 - t152 * t84);
t160 = t148 * t99 + (-t125 * t191 + t166) * qJD(3);
t14 = t159 - t218 - 0.2e1 * t233;
t157 = t158 + t226;
t138 = t143 * qJDD(1);
t137 = t142 * qJDD(1);
t130 = t190 * t154;
t115 = qJ(4) * t123;
t102 = t152 * t124;
t101 = t152 * t123;
t100 = t150 * t123;
t85 = t124 - 0.2e1 * t116;
t60 = -t101 + t194;
t59 = -t100 - t102;
t50 = t116 + t86;
t49 = t117 + t84;
t48 = -t117 + t84;
t39 = -pkin(4) * t228 - qJ(4) * t83;
t38 = t150 * t50 - t101;
t37 = -t152 * t48 + t194;
t36 = -t152 * t50 - t100;
t35 = -t150 * t48 - t102;
t23 = pkin(4) * t227 - t174 * t212;
t21 = t172 - t215;
t20 = t158 - t218;
t17 = t22 - t241;
t16 = t237 + t172;
t15 = (t49 + t84) * pkin(3) - t158;
t13 = -qJ(5) * t222 + t164 - t215;
t12 = (pkin(4) * t125 + t220) * qJD(3) + t163;
t11 = pkin(4) * t222 + t157 - t218;
t10 = -t241 + (t124 + t116) * pkin(4) + t161;
t7 = (-t222 - t231) * qJ(5) + (-t123 - t84) * pkin(4) + t237 + t165;
t6 = pkin(4) * t198 + t14 + t226;
t5 = t157 + (-t106 + t222) * pkin(4) - qJ(5) * t83 + (-t83 - t84) * pkin(3);
t4 = t12 * t150 + t13 * t152;
t3 = -t12 * t152 + t13 * t150;
t2 = pkin(4) * t12 + qJ(4) * t11;
t1 = pkin(4) * t13 + t11 * t212;
t8 = [0, 0, 0, 0, 0, qJDD(1), t177, t168, 0, 0, t137, 0.2e1 * t148 * t184, 0, t138, 0, 0, -t225 * t149, t225 * t148, pkin(1) * t130 + qJ(2) * (t138 + t137) + t176, pkin(1) * t118 + qJ(2) * t176, t25, t148 * (-t150 * t85 - t199) + t149 * (t152 * t85 - t204), t248, t24, -t257, t160, t148 * (-t207 - t261) + t149 * (t201 + t267) + t268, t148 * (-t201 + t256) + t149 * (-pkin(2) * t85 - t207 + t255) - pkin(1) * t85 + t266, t148 * (-pkin(6) * t35 - t18) + t149 * (pkin(6) * t37 + t19 - t244) - t245 + qJ(2) * (-t148 * t35 + t149 * t37), -pkin(6) * t210 + t149 * (pkin(2) * t75 + pkin(6) * t19) + pkin(1) * t75 + qJ(2) * (t149 * t19 - t210), (t148 * (-t125 * t152 + t127 * t150) + t166) * qJD(3), -t248, t257, t25, t148 * (-t152 * t49 + t209) + t149 * (-t150 * t49 - t203), t24, t148 * (-pkin(6) * t36 - t150 * t16 + t152 * t17) + t149 * (pkin(6) * t38 + t150 * t17 + t152 * t16 - t244) - t245 + qJ(2) * (-t148 * t36 + t149 * t38), t148 * (-t15 * t150 + t261) + t149 * (t15 * t152 + t264) + t265 + (qJ(4) * t191 + t149 * t178 + pkin(1)) * t49, t148 * (pkin(3) * t209 + t14 * t152 - t256) + t149 * (-t255 + t14 * t150 - (pkin(2) + t216) * t174) + t247, (t148 * (-pkin(3) * t150 + qJ(4) * t152) + t149 * (t178 + t216) + pkin(1)) * t20 + (qJ(2) + pkin(6)) * (-t148 * (t150 * t21 - t152 * t22) + t149 * (t150 * t22 + t152 * t21)), t160, t257, t248, t24, t148 * (t199 - t209) + t149 * (t203 + t204), t25, t148 * (-pkin(6) * t59 + t10 * t152 - t150 * t7) + t149 * (pkin(6) * t60 + t10 * t150 + t152 * t7 - t244) - t245 + qJ(2) * (-t148 * t59 + t149 * t60), t148 * (-t150 * t23 + t152 * t6 - t256) + t149 * (-pkin(2) * t174 + t150 * t6 + t152 * t23 - t255) + t247, t148 * (-t150 * t5 + t152 * t39 - t261) + t149 * (t150 * t39 + t152 * t5 + t267) + t268, t148 * (-pkin(6) * t3 - t1 * t150 + t152 * t2) + t149 * (pkin(2) * t11 + pkin(6) * t4 + t1 * t152 + t150 * t2) + pkin(1) * t11 + qJ(2) * (-t148 * t3 + t149 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t184, t185, -t130, -t118, 0, 0, 0, 0, 0, 0, t83, t85, t231, -t75, 0, 0, 0, 0, 0, 0, t231, -t49, t174, -t20, 0, 0, 0, 0, 0, 0, t231, t174, t83, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t229, t124, -t193, -t123, qJDD(3), -t27, -t28, 0, 0, qJDD(3), -t50, t48, t193, -t229, -t193, -pkin(3) * t50 - t115, -pkin(3) * t228 + t22 + t252, t172 + t234, -pkin(3) * t22 + qJ(4) * t21, qJDD(3), t123, t50, -t193, t229, t193, -t124 * t212 - t115, -qJ(5) * t198 + t164 + t234, -pkin(4) * t116 + t212 * t228 - t161 - t252, qJ(4) * t13 - t12 * t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t228, -t230, t22, 0, 0, 0, 0, 0, 0, t124, -t230, -t228, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t227, -t106, t13;];
tauJ_reg = t8;
