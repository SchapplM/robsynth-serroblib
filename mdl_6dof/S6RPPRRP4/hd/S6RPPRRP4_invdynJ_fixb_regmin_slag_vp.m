% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:06:22
% EndTime: 2019-03-09 02:06:29
% DurationCPUTime: 2.80s
% Computational Cost: add. (3460->401), mult. (6166->502), div. (0->0), fcn. (3822->8), ass. (0->186)
t130 = cos(pkin(9));
t205 = qJD(1) * qJD(2);
t106 = t130 * t205;
t129 = sin(pkin(9));
t201 = qJDD(1) * t130;
t136 = -pkin(1) - pkin(2);
t98 = t136 * qJDD(1) + qJDD(2);
t232 = qJ(2) * t201 + t129 * t98;
t53 = t106 + t232;
t48 = -qJDD(1) * pkin(7) + t53;
t267 = -qJD(3) * qJD(4) - t48;
t133 = sin(qJ(4));
t135 = cos(qJ(4));
t217 = qJ(2) * qJD(1);
t99 = t136 * qJD(1) + qJD(2);
t70 = t129 * t99 + t130 * t217;
t59 = -qJD(1) * pkin(7) + t70;
t44 = qJD(3) * t135 - t133 * t59;
t266 = qJD(4) * t44;
t204 = qJD(1) * qJD(4);
t182 = t135 * t204;
t198 = t133 * qJDD(1);
t265 = t182 + t198;
t132 = sin(qJ(5));
t207 = t132 * qJD(4);
t134 = cos(qJ(5));
t208 = qJD(5) * t134;
t151 = t133 * t208 + t135 * t207;
t250 = sin(qJ(1));
t251 = cos(qJ(1));
t82 = -t250 * t129 - t251 * t130;
t83 = t251 * t129 - t250 * t130;
t173 = g(1) * t82 + g(2) * t83;
t264 = g(3) * t135 - t133 * t173;
t101 = qJD(1) * t135 + qJD(5);
t117 = t135 * qJDD(1);
t183 = t133 * t204;
t262 = -t183 + t117;
t39 = -qJD(4) * pkin(4) - t44;
t211 = qJD(4) * t134;
t216 = qJD(1) * t133;
t84 = t132 * t216 + t211;
t85 = t134 * t216 - t207;
t15 = -pkin(5) * t84 + qJ(6) * t85 + t39;
t81 = -qJDD(5) - t262;
t254 = pkin(8) * t81;
t261 = t101 * t15 + t254;
t210 = qJD(4) * t135;
t186 = t134 * t210;
t209 = qJD(5) * t132;
t31 = qJD(1) * t151 - qJD(5) * t207 + qJDD(4) * t134 + t132 * t198;
t260 = -t133 * (-t134 * t31 + t84 * t209) + t84 * t186;
t239 = pkin(8) * qJD(5);
t259 = t101 * t239 - t264;
t233 = t134 * t81;
t155 = t101 * t209 + t233;
t188 = t101 * t211;
t230 = qJD(4) * t85;
t229 = qJD(5) * t84;
t30 = qJDD(4) * t132 - t134 * t265 + t229;
t258 = t133 * (t155 - t230) - t135 * (t30 + t188);
t190 = t130 * t216;
t212 = qJD(4) * t133;
t187 = t129 * t212;
t224 = t134 * t135;
t226 = t132 * t135;
t78 = t129 * t226 + t130 * t134;
t243 = -qJD(5) * t78 - t134 * t187 - (t129 * t132 + t130 * t224) * qJD(1);
t79 = t129 * t224 - t130 * t132;
t257 = t101 * t243 + t129 * (-t133 * t30 + t85 * t210) - t85 * t190 - t79 * t81;
t256 = t85 ^ 2;
t255 = pkin(5) * t81;
t249 = g(3) * t133;
t45 = t133 * qJD(3) + t135 * t59;
t40 = qJD(4) * pkin(8) + t45;
t171 = pkin(4) * t135 + pkin(8) * t133;
t69 = -t129 * t217 + t130 * t99;
t58 = qJD(1) * pkin(3) - t69;
t41 = t171 * qJD(1) + t58;
t13 = t132 * t41 + t134 * t40;
t7 = qJ(6) * t101 + t13;
t247 = t101 * t7;
t193 = t133 * t267 - t59 * t210;
t11 = -qJDD(4) * pkin(4) - qJDD(3) * t135 - t193;
t3 = -pkin(5) * t31 - qJ(6) * t30 + qJD(6) * t85 + t11;
t246 = t134 * t3;
t245 = t85 * t84;
t170 = -pkin(4) * t133 + pkin(8) * t135;
t89 = t170 * qJD(1);
t244 = t132 * t89 + t134 * t44;
t242 = qJD(5) * t79 - t132 * t187 - (-t129 * t134 + t130 * t226) * qJD(1);
t163 = pkin(3) + t171;
t90 = -t129 * qJ(2) + t130 * t136;
t60 = t163 - t90;
t91 = t130 * qJ(2) + t129 * t136;
t87 = -pkin(7) + t91;
t241 = t132 * t60 + t87 * t224;
t165 = pkin(5) * t132 - qJ(6) * t134;
t240 = -qJD(6) * t132 + t101 * t165 - t45;
t238 = qJ(6) * t81;
t237 = t101 * t13;
t236 = t101 * t84;
t235 = t101 * t85;
t234 = t132 * t87;
t231 = pkin(1) * qJDD(1);
t228 = t101 * t134;
t227 = t132 * t133;
t225 = t133 * t134;
t12 = -t132 * t40 + t134 * t41;
t223 = qJD(6) - t12;
t222 = qJDD(3) + g(3);
t221 = t251 * pkin(1) + t250 * qJ(2);
t220 = g(1) * t250 - g(2) * t251;
t126 = t133 ^ 2;
t219 = -t135 ^ 2 + t126;
t137 = qJD(4) ^ 2;
t138 = qJD(1) ^ 2;
t218 = t137 + t138;
t214 = qJD(2) * t130;
t206 = qJ(2) * qJDD(1);
t202 = qJDD(1) * t129;
t200 = qJDD(4) * t133;
t199 = qJDD(4) * t135;
t10 = qJDD(4) * pkin(8) + qJDD(3) * t133 + t135 * t48 + t266;
t104 = t129 * t205;
t180 = -qJ(2) * t202 + t130 * t98;
t52 = -t104 + t180;
t47 = qJDD(1) * pkin(3) - t52;
t21 = t262 * pkin(4) + pkin(8) * t265 + t47;
t197 = -t134 * t10 - t132 * t21 - t41 * t208;
t196 = -t151 * t85 + t30 * t227;
t189 = t135 * t214;
t74 = qJD(2) * t129 + t170 * qJD(4);
t195 = t132 * t74 + t134 * t189 + t60 * t208;
t194 = g(3) * t224 - t173 * t225;
t192 = 0.2e1 * t205;
t191 = t251 * pkin(2) + t221;
t181 = t132 * t10 - t134 * t21 + t40 * t208 + t41 * t209;
t179 = 0.2e1 * t182;
t178 = qJDD(2) - t231;
t32 = -t82 * t134 + t83 * t226;
t36 = -t83 * t134 - t82 * t226;
t176 = g(1) * t32 + g(2) * t36;
t33 = t132 * t82 + t83 * t224;
t37 = t132 * t83 - t82 * t224;
t175 = -g(1) * t33 - g(2) * t37;
t174 = -g(1) * t83 + g(2) * t82;
t172 = -t250 * pkin(1) + t251 * qJ(2);
t6 = -pkin(5) * t101 + t223;
t168 = t132 * t7 - t134 * t6;
t167 = t132 * t6 + t134 * t7;
t166 = -pkin(5) * t134 - qJ(6) * t132;
t164 = t129 * t69 - t130 * t70;
t162 = pkin(4) - t166;
t160 = -t165 + t87;
t157 = -t40 * t209 - t197;
t156 = t101 * t208 - t132 * t81;
t154 = g(1) * t251 + g(2) * t250;
t152 = qJD(1) * t58 - t173;
t150 = -t250 * pkin(2) + t172;
t149 = t101 * t39 + t254;
t148 = t173 * t135 + t249;
t146 = g(1) * t36 - g(2) * t32 - g(3) * t227 - t181;
t86 = pkin(3) - t90;
t145 = -qJDD(4) * t87 + (-qJD(1) * t86 - t214 - t58) * qJD(4);
t1 = qJD(6) * t101 + t157 - t238;
t2 = qJDD(6) + t181 + t255;
t144 = -t168 * qJD(5) + t1 * t134 + t2 * t132;
t143 = -t242 * t101 + t84 * t190 + t78 * t81 + (-t133 * t31 - t210 * t84) * t129;
t142 = -g(1) * t37 + g(2) * t33 + g(3) * t225 + t157;
t141 = qJDD(1) * t86 - t137 * t87 + t104 + t174 + t47;
t140 = -t15 * t85 + qJDD(6) - t146;
t139 = -t101 * t151 + t135 * t31 - t84 * t212 + t81 * t227;
t95 = -t133 * t137 + t199;
t94 = -t135 * t137 - t200;
t46 = -pkin(5) * t85 - qJ(6) * t84;
t43 = t160 * t133;
t23 = -t134 * t60 + (-pkin(5) + t234) * t135;
t22 = qJ(6) * t135 + t241;
t20 = pkin(5) * t216 + t132 * t44 - t134 * t89;
t19 = -qJ(6) * t216 + t244;
t18 = t30 - t236;
t14 = t160 * t210 + (t166 * qJD(5) + qJD(6) * t134 + t214) * t133;
t5 = pkin(5) * t212 + (qJD(5) * t135 * t87 - t74) * t134 + (qJD(5) * t60 - t87 * t212 + t189) * t132;
t4 = (-t87 * t209 + qJD(6)) * t135 + (-t134 * t87 - qJ(6)) * t212 + t195;
t8 = [qJDD(1), t220, t154, -qJDD(2) + t220 + 0.2e1 * t231, -t154 + t192 + 0.2e1 * t206, -t178 * pkin(1) - g(1) * t172 - g(2) * t221 + (t192 + t206) * qJ(2), -qJDD(1) * t90 + 0.2e1 * t104 + t174 - t180, qJDD(1) * t91 + 0.2e1 * t106 + t173 + t232, -g(1) * t150 - g(2) * t191 - qJD(2) * t164 + t52 * t90 + t53 * t91, qJDD(1) * t126 + t133 * t179, 0.2e1 * t117 * t133 - 0.2e1 * t204 * t219, t94, -t95, 0, t133 * t145 + t135 * t141, -t133 * t141 + t135 * t145, t85 * t186 + (-t134 * t30 - t209 * t85) * t133, t196 - t260 (t30 - t188) * t135 + (t155 + t230) * t133 (t101 * t207 + t31) * t135 + (-qJD(4) * t84 + t156) * t133, -t101 * t212 - t135 * t81 (t134 * t74 - t209 * t60) * t101 - t60 * t233 + ((-t132 * t214 - t208 * t87) * t101 + t81 * t234 + (-t39 * t132 - t87 * t84) * qJD(4) - t181) * t135 + (-t84 * t214 - t39 * t208 - t11 * t132 - t31 * t87 + (t101 * t234 - t12) * qJD(4)) * t133 + t175, -t195 * t101 + t241 * t81 + ((t101 * t87 + t40) * t209 + (-t39 * t134 - t85 * t87) * qJD(4) + t197) * t135 + (-t85 * t214 + t39 * t209 - t11 * t134 + t87 * t30 + (t228 * t87 + t13) * qJD(4)) * t133 + t176, -t101 * t5 - t14 * t84 + t23 * t81 - t31 * t43 + (-t15 * t207 - t2) * t135 + (qJD(4) * t6 - t132 * t3 - t15 * t208) * t133 + t175, t22 * t31 + t23 * t30 + t4 * t84 - t5 * t85 + t168 * t210 + (qJD(5) * t167 + t1 * t132 - t134 * t2 + t174) * t133, t101 * t4 + t14 * t85 - t22 * t81 - t30 * t43 + (t15 * t211 + t1) * t135 + (-qJD(4) * t7 - t15 * t209 + t246) * t133 - t176, t1 * t22 + t7 * t4 + t3 * t43 + t15 * t14 + t2 * t23 + t6 * t5 - g(1) * (t33 * pkin(5) + t32 * qJ(6) + t150) - g(2) * (pkin(5) * t37 + qJ(6) * t36 + t191) + (-g(2) * pkin(7) - g(1) * t163) * t83 + (-g(1) * pkin(7) + g(2) * t163) * t82; 0, 0, 0, -qJDD(1), -t138, -qJ(2) * t138 + t178 - t220, -t129 * t138 - t201, -t130 * t138 + t202, qJD(1) * t164 + t129 * t53 + t130 * t52 - t220, 0, 0, 0, 0, 0 (0.2e1 * t183 - t117) * t130 + (-t135 * t218 - t200) * t129 (t179 + t198) * t130 + (t133 * t218 - t199) * t129, 0, 0, 0, 0, 0, t143, -t257, t143, -t242 * t85 + t243 * t84 + t30 * t78 + t31 * t79, t257, -t15 * t190 + t1 * t79 + t2 * t78 + t243 * t7 + t242 * t6 + (t133 * t3 + t15 * t210) * t129 - t220; 0, 0, 0, 0, 0, 0, 0, 0, t222, 0, 0, 0, 0, 0, t95, t94, 0, 0, 0, 0, 0, t139, t258, t139, t196 + t260, -t258, g(3) + (qJD(4) * t167 - t3) * t135 + (qJD(4) * t15 + t144) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133 * t138 * t135, t219 * t138, -t198, -t117, qJDD(4), qJD(4) * t45 + t133 * t152 + t135 * t222 + t193, t266 + (qJD(4) * t59 - t222) * t133 + (t152 + t267) * t135, t132 * t30 - t228 * t85 (t30 + t236) * t134 + (t31 + t235) * t132 (t101 * t224 - t133 * t85) * qJD(1) + t156 (-t101 * t226 + t133 * t84) * qJD(1) - t155, t101 * t216, t12 * t216 + pkin(4) * t31 + t45 * t84 + (-t11 + (-t89 - t239) * t101) * t134 + (t44 * t101 + t149) * t132 + t194, -pkin(4) * t30 + t244 * t101 - t13 * t216 + t45 * t85 + t149 * t134 + (t11 + t259) * t132, -t6 * t216 - t246 + t31 * t162 - t240 * t84 + (-pkin(8) * t208 + t20) * t101 + t261 * t132 + t194, -t19 * t84 + t20 * t85 + (t1 + t101 * t6 + (-qJD(5) * t85 + t31) * pkin(8)) * t134 + (t2 - t247 + (t30 - t229) * pkin(8)) * t132 + t148, t7 * t216 - t101 * t19 + t30 * t162 + t240 * t85 - t261 * t134 + (-t259 - t3) * t132, -t7 * t19 - t6 * t20 + t240 * t15 + (t144 + t148) * pkin(8) + (-t3 + t264) * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, -t84 ^ 2 + t256, t18, t31 - t235, -t81, t39 * t85 + t146 + t237, t101 * t12 - t39 * t84 - t142, t46 * t84 - t140 + t237 - 0.2e1 * t255, -pkin(5) * t30 + qJ(6) * t31 + (t13 - t7) * t85 + (-t6 + t223) * t84, -0.2e1 * t238 + t15 * t84 - t46 * t85 + (0.2e1 * qJD(6) - t12) * t101 + t142, t1 * qJ(6) - t2 * pkin(5) - t15 * t46 - t6 * t13 - g(1) * (-pkin(5) * t36 + qJ(6) * t37) - g(2) * (pkin(5) * t32 - qJ(6) * t33) + t223 * t7 - t165 * t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 + t245, t18, -t101 ^ 2 - t256, t140 - t247 + t255;];
tau_reg  = t8;
