% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:48
% EndTime: 2019-03-09 07:05:05
% DurationCPUTime: 8.27s
% Computational Cost: add. (22749->526), mult. (61937->721), div. (0->0), fcn. (50021->10), ass. (0->242)
t217 = sin(pkin(11));
t222 = sin(qJ(3));
t218 = cos(pkin(11));
t226 = cos(qJ(3));
t286 = t218 * t226;
t200 = -t217 * t222 + t286;
t192 = t200 * qJD(1);
t201 = t217 * t226 + t218 * t222;
t193 = t201 * qJD(1);
t221 = sin(qJ(4));
t225 = cos(qJ(4));
t168 = t192 * t221 + t193 * t225;
t220 = sin(qJ(5));
t224 = cos(qJ(5));
t266 = t225 * t192 - t193 * t221;
t352 = -t168 * t220 + t224 * t266;
t114 = qJD(6) - t352;
t351 = t224 * t168 + t220 * t266;
t81 = pkin(5) * t351 - pkin(10) * t352;
t162 = Ifges(5,4) * t266;
t216 = qJD(3) + qJD(4);
t109 = t168 * Ifges(5,1) + t216 * Ifges(5,5) + t162;
t194 = t200 * qJD(3);
t184 = qJD(1) * t194;
t195 = t201 * qJD(3);
t185 = qJD(1) * t195;
t248 = -t184 * t225 + t185 * t221;
t115 = qJD(4) * t266 - t248;
t116 = -qJD(4) * t168 - t184 * t221 - t185 * t225;
t270 = -pkin(2) * t218 - pkin(1);
t204 = qJD(1) * t270 + qJD(2);
t176 = -pkin(3) * t192 + t204;
t213 = qJD(5) + t216;
t219 = sin(qJ(6));
t223 = cos(qJ(6));
t100 = t213 * t223 - t219 * t351;
t63 = qJD(5) * t352 + t115 * t224 + t116 * t220;
t35 = qJD(6) * t100 + t223 * t63;
t101 = t213 * t219 + t223 * t351;
t36 = -qJD(6) * t101 - t219 * t63;
t64 = qJD(5) * t351 + t115 * t220 - t224 * t116;
t12 = t35 * Ifges(7,4) + t36 * Ifges(7,2) + t64 * Ifges(7,6);
t13 = t35 * Ifges(7,1) + t36 * Ifges(7,4) + t64 * Ifges(7,5);
t255 = Ifges(7,5) * t219 + Ifges(7,6) * t223;
t305 = Ifges(7,4) * t219;
t257 = Ifges(7,2) * t223 + t305;
t304 = Ifges(7,4) * t223;
t259 = Ifges(7,1) * t219 + t304;
t262 = mrSges(7,1) * t223 - mrSges(7,2) * t219;
t95 = pkin(3) * t185 - pkin(4) * t116;
t17 = pkin(5) * t64 - pkin(10) * t63 + t95;
t346 = pkin(9) * t266;
t310 = pkin(7) + qJ(2);
t205 = t310 * t217;
t202 = qJD(1) * t205;
t206 = t310 * t218;
t203 = qJD(1) * t206;
t175 = -t202 * t222 + t203 * t226;
t148 = pkin(8) * t192 + t175;
t145 = t225 * t148;
t339 = -t226 * t202 - t203 * t222;
t147 = -pkin(8) * t193 + t339;
t146 = qJD(3) * pkin(3) + t147;
t94 = t146 * t221 + t145;
t86 = t94 + t346;
t292 = t224 * t86;
t354 = pkin(9) * t168;
t143 = t221 * t148;
t93 = t225 * t146 - t143;
t85 = t93 - t354;
t82 = pkin(4) * t216 + t85;
t43 = t220 * t82 + t292;
t41 = pkin(10) * t213 + t43;
t131 = -pkin(4) * t266 + t176;
t67 = -pkin(5) * t352 - pkin(10) * t351 + t131;
t20 = -t219 * t41 + t223 * t67;
t173 = t200 * t221 + t201 * t225;
t275 = qJD(1) * qJD(2);
t62 = -t94 * qJD(4) + t248 * pkin(8) + (-t175 * t225 - t221 * t339) * qJD(3) - t173 * t275;
t228 = -t115 * pkin(9) + t62;
t141 = qJD(2) * t192 + qJD(3) * t339;
t242 = t201 * qJD(2);
t142 = -qJD(1) * t242 - qJD(3) * t175;
t281 = qJD(4) * t225;
t282 = qJD(4) * t221;
t61 = -t148 * t282 + t221 * (-pkin(8) * t184 + t142) + t225 * (-pkin(8) * t185 + t141) + t146 * t281;
t38 = pkin(9) * t116 + t61;
t293 = t220 * t86;
t42 = t224 * t82 - t293;
t8 = qJD(5) * t42 + t220 * t228 + t224 * t38;
t2 = qJD(6) * t20 + t17 * t219 + t223 * t8;
t314 = t2 * t223;
t315 = t223 / 0.2e1;
t331 = t64 / 0.2e1;
t332 = t36 / 0.2e1;
t333 = t35 / 0.2e1;
t256 = Ifges(7,5) * t223 - Ifges(7,6) * t219;
t258 = -Ifges(7,2) * t219 + t304;
t260 = Ifges(7,1) * t223 - t305;
t261 = mrSges(7,1) * t219 + mrSges(7,2) * t223;
t316 = -t219 / 0.2e1;
t327 = t101 / 0.2e1;
t40 = -pkin(5) * t213 - t42;
t297 = t101 * Ifges(7,4);
t53 = t100 * Ifges(7,2) + t114 * Ifges(7,6) + t297;
t99 = Ifges(7,4) * t100;
t54 = t101 * Ifges(7,1) + t114 * Ifges(7,5) + t99;
t350 = t40 * t261 + t53 * t316 + t54 * t315 + t100 * t258 / 0.2e1 + t260 * t327 + t114 * t256 / 0.2e1;
t9 = qJD(5) * t43 + t220 * t38 - t224 * t228;
t236 = -t8 * mrSges(6,2) + mrSges(7,3) * t314 + t219 * t13 / 0.2e1 + t12 * t315 + t259 * t333 + t257 * t332 + t255 * t331 - Ifges(6,6) * t64 + Ifges(6,5) * t63 + (-mrSges(6,1) - t262) * t9 + t350 * qJD(6);
t307 = Ifges(5,4) * t168;
t323 = -t168 / 0.2e1;
t357 = t62 * mrSges(5,1) - t61 * mrSges(5,2) + Ifges(5,5) * t115 + Ifges(5,6) * t116 + t236 - (Ifges(5,5) * t266 - Ifges(5,6) * t168) * t216 / 0.2e1 + (t168 * t94 + t266 * t93) * mrSges(5,3) - (-Ifges(5,2) * t168 + t109 + t162) * t266 / 0.2e1 - t176 * (mrSges(5,1) * t168 + mrSges(5,2) * t266) + (Ifges(5,1) * t266 - t307) * t323;
t356 = Ifges(6,2) / 0.2e1;
t355 = pkin(4) * t168;
t16 = -mrSges(7,1) * t36 + mrSges(7,2) * t35;
t353 = m(7) * t9 + t16;
t108 = Ifges(5,2) * t266 + t216 * Ifges(5,6) + t307;
t348 = t108 / 0.2e1;
t347 = t352 * t356;
t21 = t219 * t67 + t223 * t41;
t295 = t20 * t223;
t253 = t21 * t219 + t295;
t243 = t253 * mrSges(7,3);
t212 = pkin(3) * t225 + pkin(4);
t96 = -t147 * t221 - t145;
t246 = t96 - t346;
t279 = qJD(5) * t224;
t280 = qJD(5) * t220;
t284 = t221 * t224;
t97 = t225 * t147 - t143;
t87 = t97 - t354;
t341 = t220 * t87 - t224 * t246 - t212 * t280 - (t221 * t279 + (t220 * t225 + t284) * qJD(4)) * pkin(3);
t290 = mrSges(6,1) * t213 + mrSges(7,1) * t100 - mrSges(7,2) * t101 - mrSges(6,3) * t351;
t198 = t226 * t205;
t177 = -t206 * t222 - t198;
t160 = -pkin(8) * t201 + t177;
t178 = -t222 * t205 + t226 * t206;
t161 = pkin(8) * t200 + t178;
t106 = t221 * t160 + t225 * t161;
t113 = Ifges(6,4) * t352;
t273 = -t113 / 0.2e1;
t309 = Ifges(6,1) * t351;
t338 = t273 - t309 / 0.2e1;
t337 = -t20 * t219 + t21 * t223;
t3 = -qJD(6) * t21 + t17 * t223 - t219 * t8;
t336 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t35 + Ifges(7,6) * t36;
t335 = (m(3) * qJ(2) + mrSges(3,3)) * (t217 ^ 2 + t218 ^ 2);
t105 = t225 * t160 - t161 * t221;
t89 = -pkin(9) * t173 + t105;
t172 = t200 * t225 - t201 * t221;
t90 = pkin(9) * t172 + t106;
t65 = t220 * t90 - t224 * t89;
t330 = t65 * t9;
t329 = -t100 / 0.2e1;
t328 = -t101 / 0.2e1;
t326 = -t114 / 0.2e1;
t324 = t266 / 0.2e1;
t322 = t168 / 0.2e1;
t320 = -t193 / 0.2e1;
t319 = t194 / 0.2e1;
t318 = -t195 / 0.2e1;
t313 = t3 * t219;
t312 = t42 * mrSges(6,3);
t311 = t43 * mrSges(6,3);
t308 = Ifges(4,4) * t193;
t285 = t220 * t221;
t190 = pkin(3) * t284 + t220 * t212;
t278 = qJD(6) * t219;
t277 = qJD(6) * t223;
t269 = t64 * mrSges(6,1) + t63 * mrSges(6,2);
t268 = -t116 * mrSges(5,1) + t115 * mrSges(5,2);
t128 = -qJD(4) * t173 - t194 * t221 - t195 * t225;
t102 = pkin(3) * t195 - pkin(4) * t128;
t137 = pkin(3) * t193 + t355;
t263 = -t2 * t219 - t3 * t223;
t18 = mrSges(7,1) * t64 - mrSges(7,3) * t35;
t19 = -mrSges(7,2) * t64 + mrSges(7,3) * t36;
t254 = -t219 * t18 + t223 * t19;
t66 = t220 * t89 + t224 * t90;
t130 = t172 * t220 + t173 * t224;
t181 = -pkin(3) * t200 + t270;
t140 = -pkin(4) * t172 + t181;
t249 = t224 * t172 - t173 * t220;
t74 = -pkin(5) * t249 - pkin(10) * t130 + t140;
t30 = t219 * t74 + t223 * t66;
t29 = -t219 * t66 + t223 * t74;
t76 = -mrSges(7,2) * t114 + mrSges(7,3) * t100;
t77 = mrSges(7,1) * t114 - mrSges(7,3) * t101;
t251 = -t219 * t77 + t223 * t76;
t189 = -pkin(3) * t285 + t212 * t224;
t103 = -mrSges(6,2) * t213 + mrSges(6,3) * t352;
t245 = -t103 - t251;
t152 = -qJD(3) * t198 + qJD(2) * t286 + (-qJD(2) * t217 - qJD(3) * t206) * t222;
t135 = -pkin(8) * t195 + t152;
t153 = -qJD(3) * t178 - t242;
t136 = -pkin(8) * t194 + t153;
t71 = t225 * t135 + t221 * t136 + t160 * t281 - t161 * t282;
t238 = -qJD(6) * t253 - t313;
t72 = -qJD(4) * t106 - t135 * t221 + t225 * t136;
t127 = qJD(4) * t172 + t194 * t225 - t195 * t221;
t237 = -pkin(9) * t127 + t72;
t234 = -t131 * mrSges(6,1) - t20 * mrSges(7,1) + t21 * mrSges(7,2) + Ifges(6,4) * t351 - Ifges(7,5) * t101 + Ifges(6,6) * t213 - Ifges(7,6) * t100 - Ifges(7,3) * t114 + t347;
t233 = -t77 * t277 - t76 * t278 + m(7) * (-t20 * t277 - t21 * t278 - t313 + t314) + t254;
t232 = t347 + t234;
t231 = t131 * mrSges(6,2) + t113 / 0.2e1 + Ifges(6,5) * t213 + t309 / 0.2e1 + t350;
t230 = -t231 + t338;
t229 = -t231 + t243;
t188 = pkin(10) + t190;
t187 = -pkin(5) - t189;
t186 = Ifges(4,4) * t192;
t182 = t184 * mrSges(4,2);
t180 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t193;
t179 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t192;
t170 = t212 * t279 + (-t221 * t280 + (t224 * t225 - t285) * qJD(4)) * pkin(3);
t164 = t193 * Ifges(4,1) + Ifges(4,5) * qJD(3) + t186;
t163 = t192 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t308;
t150 = mrSges(5,1) * t216 - mrSges(5,3) * t168;
t149 = -mrSges(5,2) * t216 + mrSges(5,3) * t266;
t123 = -mrSges(5,1) * t266 + mrSges(5,2) * t168;
t80 = -mrSges(6,1) * t352 + mrSges(6,2) * t351;
t73 = t81 + t355;
t70 = t137 + t81;
t69 = qJD(5) * t130 + t127 * t220 - t224 * t128;
t68 = qJD(5) * t249 + t127 * t224 + t128 * t220;
t58 = Ifges(7,3) * t64;
t49 = pkin(9) * t128 + t71;
t48 = t220 * t246 + t224 * t87;
t45 = t224 * t85 - t293;
t44 = t220 * t85 + t292;
t28 = t219 * t81 + t223 * t42;
t27 = -t219 * t42 + t223 * t81;
t26 = t219 * t73 + t223 * t45;
t25 = -t219 * t45 + t223 * t73;
t24 = t219 * t70 + t223 * t48;
t23 = -t219 * t48 + t223 * t70;
t22 = pkin(5) * t69 - pkin(10) * t68 + t102;
t15 = qJD(5) * t66 + t220 * t49 - t224 * t237;
t14 = -qJD(5) * t65 + t220 * t237 + t224 * t49;
t5 = -qJD(6) * t30 - t14 * t219 + t22 * t223;
t4 = qJD(6) * t29 + t14 * t223 + t219 * t22;
t1 = [(-t42 * t68 - t43 * t69 + t63 * t65 - t64 * t66) * mrSges(6,3) - t232 * t69 + (t141 * t200 - t142 * t201 - t175 * t195 - t177 * t184 - t178 * t185 - t194 * t339) * mrSges(4,3) + (t185 * (-mrSges(5,1) * t172 + mrSges(5,2) * t173) + t195 * t123) * pkin(3) + m(5) * (t105 * t62 + t106 * t61 + t71 * t94 + t72 * t93 + (t176 * t195 + t181 * t185) * pkin(3)) + t270 * (t185 * mrSges(4,1) + t182) + (-t200 * t185 + t192 * t318) * Ifges(4,2) + (t200 * t184 - t201 * t185 + t192 * t319 + t193 * t318) * Ifges(4,4) - (t95 * mrSges(6,1) - Ifges(6,4) * t63 + t58 / 0.2e1 - t8 * mrSges(6,3) + (Ifges(6,2) + Ifges(7,3) / 0.2e1) * t64 + t336) * t249 + m(4) * (t141 * t178 + t142 * t177 + t152 * t175 + t153 * t339) + t71 * t149 + t72 * t150 + (-t105 * t115 + t106 * t116 - t127 * t93 + t128 * t94 + t172 * t61 - t173 * t62) * mrSges(5,3) + 0.2e1 * t335 * t275 + t128 * t348 + (-t243 - t230) * t68 + qJD(3) * (Ifges(4,5) * t194 - Ifges(4,6) * t195) / 0.2e1 + t204 * (mrSges(4,1) * t195 + mrSges(4,2) * t194) + t127 * t109 / 0.2e1 + t102 * t80 + t14 * t103 + t4 * t76 + t5 * t77 + t65 * t16 + t176 * (-mrSges(5,1) * t128 + mrSges(5,2) * t127) + t152 * t179 + t153 * t180 + t163 * t318 + t164 * t319 + t216 * (Ifges(5,5) * t127 + Ifges(5,6) * t128) / 0.2e1 + t181 * t268 + t140 * t269 - t290 * t15 + (t95 * mrSges(6,2) + Ifges(6,1) * t63 - Ifges(6,4) * t64 + t260 * t333 + t258 * t332 + t256 * t331 + t12 * t316 + t13 * t315 + (mrSges(6,3) + t261) * t9 + t263 * mrSges(7,3) + (-t223 * t53 / 0.2e1 + t54 * t316 + t255 * t326 + t257 * t329 + t259 * t328 + t40 * t262 - t337 * mrSges(7,3)) * qJD(6)) * t130 + (t201 * t184 + t193 * t319) * Ifges(4,1) + (t173 * t115 + t127 * t322) * Ifges(5,1) + (t172 * t116 + t128 * t324) * Ifges(5,2) + (t172 * t115 + t173 * t116 + t127 * t324 + t128 * t322) * Ifges(5,4) + m(7) * (t15 * t40 + t2 * t30 + t20 * t5 + t21 * t4 + t29 * t3 + t330) + m(6) * (t102 * t131 + t14 * t43 + t140 * t95 - t15 * t42 + t66 * t8 + t330) + t29 * t18 + t30 * t19; t269 + t268 + t182 + t168 * t150 - t266 * t149 - (-m(5) * pkin(3) - mrSges(4,1)) * t185 + t290 * t351 - m(5) * (-t168 * t93 + t266 * t94) - m(4) * (t175 * t192 - t193 * t339) - t192 * t179 + t193 * t180 + t245 * t352 + t219 * t19 + t223 * t18 + t251 * qJD(6) - t335 * qJD(1) ^ 2 + (t114 * t337 - t351 * t40 - t263) * m(7) + (t351 * t42 - t352 * t43 + t95) * m(6); -m(5) * (t93 * t96 + t94 * t97) + t168 * t348 + (t175 * t193 + t192 * t339) * mrSges(4,3) - t339 * t179 - t97 * t149 - t96 * t150 + t357 + ((-t219 * t76 - t223 * t77) * t188 - t243) * qJD(6) - (-Ifges(4,2) * t193 + t164 + t186) * t192 / 0.2e1 + (-t193 * t123 + (t149 * t225 - t150 * t221) * qJD(4) + (-t115 * t225 + t116 * t221) * mrSges(5,3) + (0.2e1 * t176 * t320 + t221 * t61 + t225 * t62 + t281 * t94 - t282 * t93) * m(5)) * pkin(3) - t137 * t80 - t141 * mrSges(4,2) + t142 * mrSges(4,1) - t48 * t103 - t24 * t76 - t23 * t77 + t175 * t180 + Ifges(4,5) * t184 - Ifges(4,6) * t185 + t187 * t16 + (Ifges(4,1) * t192 - t308) * t320 + t193 * t163 / 0.2e1 - qJD(3) * (Ifges(4,5) * t192 - Ifges(4,6) * t193) / 0.2e1 - t204 * (mrSges(4,1) * t193 + mrSges(4,2) * t192) + (-t131 * t137 - t189 * t9 + t190 * t8 + (t170 - t48) * t43 + t341 * t42) * m(6) + t290 * t341 + ((t238 + t314) * t188 - t20 * t23 - t21 * t24 + t187 * t9 - t341 * t40 + t337 * t170) * m(7) - t245 * t170 + t254 * t188 + t232 * t351 + (t229 + t338) * t352 + (-t189 * t63 - t190 * t64 + t351 * t43 + t352 * t42) * mrSges(6,3) - mrSges(7,3) * t313; t233 * (pkin(4) * t220 + pkin(10)) - m(7) * (t20 * t25 + t21 * t26 + t40 * t44) - m(6) * (-t42 * t44 + t43 * t45) + (-t114 * t295 + (-t114 * t21 - t3) * t219) * mrSges(7,3) + t353 * (-pkin(4) * t224 - pkin(5)) - t93 * t149 + t94 * t150 + (t232 + t311) * t351 + (t230 + t312) * t352 - t45 * t103 - t26 * t76 - t25 * t77 + t108 * t322 + (-t168 * t80 + (-t220 * t64 - t224 * t63) * mrSges(6,3) + (-t290 * t220 - t245 * t224 + m(7) * (t220 * t40 + t224 * t337)) * qJD(5) + (0.2e1 * t131 * t323 + t220 * t8 - t224 * t9 + (-t220 * t42 + t224 * t43) * qJD(5)) * m(6)) * pkin(4) + t290 * t44 + t357; -m(7) * (t20 * t27 + t21 * t28 + t40 * t43) + t233 * pkin(10) + (t234 + t311) * t351 + (t273 + t312 + (-Ifges(6,1) / 0.2e1 + t356) * t351 + t229) * t352 + t236 - t42 * t103 - t28 * t76 - t27 * t77 + t238 * mrSges(7,3) + t290 * t43 - t353 * pkin(5); t58 - t40 * (mrSges(7,1) * t101 + mrSges(7,2) * t100) + (Ifges(7,1) * t100 - t297) * t328 + t53 * t327 + (Ifges(7,5) * t100 - Ifges(7,6) * t101) * t326 - t20 * t76 + t21 * t77 + (t100 * t20 + t101 * t21) * mrSges(7,3) + (-Ifges(7,2) * t101 + t54 + t99) * t329 + t336;];
tauc  = t1(:);
