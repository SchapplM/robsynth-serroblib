% Calculate time derivative of joint inertia matrix for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14V3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:33
% EndTime: 2019-04-12 15:03:54
% DurationCPUTime: 7.40s
% Computational Cost: add. (2418->572), mult. (7951->887), div. (0->0), fcn. (6714->8), ass. (0->273)
t149 = sin(qJ(6));
t154 = cos(qJ(5));
t200 = qJD(6) * t154 - qJD(4);
t150 = sin(qJ(5));
t151 = sin(qJ(4));
t243 = qJD(5) * t151;
t214 = t150 * t243;
t153 = cos(qJ(6));
t266 = t151 * t153;
t155 = cos(qJ(4));
t247 = qJD(4) * t154;
t199 = -qJD(6) + t247;
t331 = t199 * t155;
t341 = t149 * (t331 - t214) + t200 * t266;
t318 = 2 * qJ(3);
t152 = sin(qJ(2));
t145 = t152 ^ 2;
t157 = qJ(3) ^ 2;
t245 = qJD(4) * t157;
t179 = t145 * t151 * t155 * t245;
t144 = t151 ^ 2;
t254 = qJ(3) * qJD(3);
t137 = t144 * t254;
t197 = t145 * t137;
t340 = t179 + t197;
t156 = cos(qJ(2));
t251 = qJD(2) * t156;
t206 = t151 * t251;
t246 = qJD(4) * t155;
t167 = t152 * t246 + t206;
t242 = qJD(5) * t154;
t165 = t150 * t246 + t151 * t242;
t143 = t150 ^ 2;
t147 = t154 ^ 2;
t255 = t143 + t147;
t317 = m(6) / 0.2e1;
t338 = t317 * (0.1e1 - t255) - m(7) * t143 / 0.2e1;
t244 = qJD(5) * t150;
t337 = t151 * (t149 * t200 + t153 * t244) - t153 * t331;
t336 = 2 * mrSges(4,3);
t335 = -0.2e1 * t155;
t334 = 0.2e1 * t156;
t296 = t154 / 0.2e1;
t332 = m(7) * (t143 - t147);
t275 = mrSges(7,1) * t153 - mrSges(7,2) * t149 + mrSges(6,1);
t329 = t275 * t150;
t268 = t150 * t151;
t104 = mrSges(6,2) * t155 - mrSges(6,3) * t268;
t263 = t154 * t104;
t265 = t151 * t154;
t171 = t149 * t265 + t153 * t155;
t87 = -t149 * t155 + t153 * t265;
t272 = mrSges(6,1) * t155 + mrSges(7,1) * t171 + mrSges(7,2) * t87 + mrSges(6,3) * t265;
t169 = -t272 * t150 - t263;
t238 = qJD(6) * t153;
t328 = -t149 * t242 - t150 * t238;
t239 = qJD(6) * t150;
t327 = t149 * t239 - t153 * t242;
t249 = qJD(4) * t151;
t326 = mrSges(6,3) * t255 - mrSges(5,2);
t267 = t151 * t152;
t105 = mrSges(5,2) * t156 - mrSges(5,3) * t267;
t261 = t154 * t156;
t264 = t152 * t155;
t85 = t150 * t264 + t261;
t63 = -mrSges(6,2) * t267 - mrSges(6,3) * t85;
t278 = t154 * t63;
t262 = t154 * t155;
t86 = -t150 * t156 + t152 * t262;
t60 = -t149 * t86 + t152 * t266;
t61 = t149 * t267 + t153 * t86;
t293 = -mrSges(6,1) * t267 - mrSges(7,1) * t60 + mrSges(7,2) * t61 + mrSges(6,3) * t86;
t324 = -t293 * t150 - t105 - t278;
t241 = qJD(5) * t155;
t323 = t150 * t241 + t151 * t199;
t276 = mrSges(6,1) * t154 - mrSges(6,2) * t150 + mrSges(5,1);
t188 = mrSges(7,1) * t149 + mrSges(7,2) * t153;
t88 = t188 * t150;
t322 = -qJD(4) * t276 + t88 * t242;
t321 = 2 * m(5);
t320 = 0.2e1 * m(6);
t319 = 0.2e1 * m(7);
t286 = Ifges(6,4) * t154;
t122 = Ifges(6,1) * t150 + t286;
t287 = Ifges(6,4) * t150;
t186 = Ifges(6,1) * t154 - t287;
t54 = -t122 * t243 + (Ifges(6,5) * t151 + t155 * t186) * qJD(4);
t316 = t54 / 0.2e1;
t315 = t60 / 0.2e1;
t314 = t61 / 0.2e1;
t284 = Ifges(7,4) * t153;
t182 = -Ifges(7,2) * t149 + t284;
t76 = -Ifges(7,6) * t154 + t150 * t182;
t313 = t76 / 0.2e1;
t285 = Ifges(7,4) * t149;
t185 = Ifges(7,1) * t153 - t285;
t79 = -Ifges(7,5) * t154 + t150 * t185;
t312 = t79 / 0.2e1;
t80 = -Ifges(6,5) * t155 + t151 * t186;
t311 = t80 / 0.2e1;
t310 = -t171 / 0.2e1;
t309 = t87 / 0.2e1;
t181 = Ifges(7,5) * t153 - Ifges(7,6) * t149;
t95 = t181 * qJD(6);
t308 = t95 / 0.2e1;
t97 = t182 * qJD(6);
t307 = t97 / 0.2e1;
t100 = t185 * qJD(6);
t306 = t100 / 0.2e1;
t101 = t186 * qJD(5);
t305 = t101 / 0.2e1;
t116 = Ifges(7,5) * t149 + Ifges(7,6) * t153;
t304 = t116 / 0.2e1;
t118 = Ifges(7,2) * t153 + t285;
t303 = t118 / 0.2e1;
t121 = Ifges(7,1) * t149 + t284;
t302 = t121 / 0.2e1;
t301 = t122 / 0.2e1;
t300 = -t149 / 0.2e1;
t299 = t149 / 0.2e1;
t298 = -t153 / 0.2e1;
t297 = t153 / 0.2e1;
t159 = -qJD(6) * t86 + t167;
t218 = t151 * t247;
t41 = (qJD(2) * t155 - qJD(5)) * t261 + (-t218 + (qJD(2) - t241) * t150) * t152;
t193 = qJD(6) * t267 + t41;
t17 = -t149 * t193 + t153 * t159;
t18 = t149 * t159 + t153 * t193;
t295 = mrSges(6,1) * t167 + mrSges(7,1) * t17 - mrSges(7,2) * t18 - mrSges(6,3) * t41;
t205 = t154 * t246;
t166 = t205 - t214;
t294 = -mrSges(6,1) * t249 + mrSges(7,1) * t341 - mrSges(7,2) * t337 + mrSges(6,3) * t166;
t42 = Ifges(7,5) * t87 - Ifges(7,6) * t171 + Ifges(7,3) * t268;
t183 = -Ifges(6,2) * t150 + t286;
t77 = -Ifges(6,6) * t155 + t151 * t183;
t292 = t42 - t77;
t291 = mrSges(6,2) * t154;
t290 = mrSges(7,3) * t150;
t289 = Ifges(5,4) * t151;
t288 = Ifges(5,4) * t155;
t283 = Ifges(5,5) * t156;
t282 = Ifges(5,6) * t152;
t281 = Ifges(5,6) * t156;
t280 = t150 * t63;
t248 = qJD(4) * t152;
t219 = t151 * t248;
t252 = qJD(2) * t152;
t40 = -t150 * t219 - t156 * t244 - t154 * t252 + (t150 * t251 + t152 * t242) * t155;
t28 = -mrSges(6,2) * t167 - mrSges(6,3) * t40;
t279 = t154 * t28;
t69 = -mrSges(6,2) * t249 - mrSges(6,3) * t165;
t277 = t154 * t69;
t274 = -mrSges(5,1) * t156 - mrSges(6,1) * t85 - mrSges(6,2) * t86 - mrSges(5,3) * t264;
t217 = t152 * t245;
t221 = qJD(3) * t267;
t269 = t143 * t144;
t271 = qJ(3) * t155;
t273 = -0.2e1 * t143 * t221 * t271 + t217 * t269;
t270 = t104 * t150;
t148 = t155 ^ 2;
t196 = t152 * t157 * t251;
t223 = t148 * t254;
t260 = t145 * t223 + t148 * t196;
t259 = Ifges(6,5) * t205 + Ifges(6,3) * t249;
t222 = t155 * t251;
t258 = -Ifges(5,5) * t222 - Ifges(5,3) * t252;
t257 = t149 ^ 2 + t153 ^ 2;
t253 = qJ(3) * qJD(4);
t250 = qJD(3) * t155;
t240 = qJD(6) * t149;
t237 = m(4) * t318;
t236 = t304 - Ifges(6,6);
t235 = Ifges(6,5) * t150 / 0.2e1 + Ifges(6,6) * t296 - Ifges(5,6);
t234 = 0.2e1 * t150;
t1 = Ifges(7,5) * t18 + Ifges(7,6) * t17 + Ifges(7,3) * t40;
t12 = Ifges(6,4) * t41 - Ifges(6,2) * t40 + Ifges(6,6) * t167;
t233 = t1 / 0.2e1 - t12 / 0.2e1;
t119 = Ifges(6,2) * t154 + t287;
t51 = -t119 * t243 + (Ifges(6,6) * t151 + t155 * t183) * qJD(4);
t9 = -Ifges(7,5) * t337 - Ifges(7,6) * t341 + Ifges(7,3) * t165;
t232 = t9 / 0.2e1 - t51 / 0.2e1;
t231 = t143 * t340 + t196 * t269;
t21 = Ifges(7,5) * t61 + Ifges(7,6) * t60 + Ifges(7,3) * t85;
t45 = Ifges(6,4) * t86 - Ifges(6,2) * t85 + Ifges(6,6) * t267;
t229 = t21 / 0.2e1 - t45 / 0.2e1;
t228 = t42 / 0.2e1 - t77 / 0.2e1;
t48 = -t116 * t239 + (Ifges(7,3) * t150 + t154 * t181) * qJD(5);
t98 = t183 * qJD(5);
t227 = t48 / 0.2e1 - t98 / 0.2e1;
t74 = -Ifges(7,3) * t154 + t150 * t181;
t225 = -t119 / 0.2e1 + t74 / 0.2e1;
t120 = Ifges(5,2) * t155 + t289;
t75 = -Ifges(6,3) * t155 + (Ifges(6,5) * t154 - Ifges(6,6) * t150) * t151;
t224 = -t120 / 0.2e1 + t75 / 0.2e1;
t212 = t150 * t242;
t93 = t188 * qJD(6);
t203 = mrSges(6,2) * qJD(5) + t93;
t10 = Ifges(6,5) * t41 - Ifges(6,6) * t40 + Ifges(6,3) * t167;
t201 = 0.2e1 * t293;
t198 = 0.2e1 * t272;
t195 = t157 * t212;
t43 = Ifges(6,5) * t86 - Ifges(6,6) * t85 + Ifges(6,3) * t267;
t184 = -Ifges(5,2) * t151 + t288;
t78 = t152 * t184 - t281;
t194 = t43 - t78 + t281;
t164 = t152 * t87;
t7 = -qJD(3) * t164 + (t152 * t337 - t87 * t251) * qJ(3);
t163 = t171 * t152;
t8 = qJD(3) * t163 + (t152 * t341 + t171 * t251) * qJ(3);
t192 = -t149 * t8 + t153 * t7;
t190 = mrSges(5,1) * t151 + mrSges(5,2) * t155;
t189 = mrSges(6,1) * t150 + t291;
t187 = Ifges(5,1) * t155 - t289;
t123 = Ifges(5,1) * t151 + t288;
t173 = t149 * t151 + t153 * t262;
t178 = t200 * t155;
t24 = t173 * qJD(3) + (-t149 * t178 - t153 * t323) * qJ(3);
t172 = -t149 * t262 + t266;
t25 = t172 * qJD(3) + (t149 * t323 - t153 * t178) * qJ(3);
t180 = -t149 * t25 + t153 * t24;
t22 = Ifges(7,4) * t61 + Ifges(7,2) * t60 + Ifges(7,6) * t85;
t23 = Ifges(7,1) * t61 + Ifges(7,4) * t60 + Ifges(7,5) * t85;
t175 = t22 * t300 + t23 * t297;
t44 = Ifges(7,4) * t87 - Ifges(7,2) * t171 + Ifges(7,6) * t268;
t46 = Ifges(7,1) * t87 - Ifges(7,4) * t171 + Ifges(7,5) * t268;
t174 = t297 * t46 + t300 * t44;
t170 = qJ(3) * qJD(5) * t189 - qJD(3) * t276;
t168 = -t150 * t88 - t326;
t162 = t152 * (t291 + t329);
t161 = -t150 * t201 - 0.2e1 * t105 - 0.2e1 * t278;
t160 = t157 * t338;
t158 = mrSges(4,1) * t334 + t151 * t161 + t152 * t336 + t274 * t335;
t141 = Ifges(5,5) * t246;
t140 = Ifges(6,5) * t242;
t124 = t143 * t223;
t106 = -mrSges(7,1) * t154 - t153 * t290;
t103 = mrSges(7,2) * t154 - t149 * t290;
t102 = t187 * qJD(4);
t99 = t184 * qJD(4);
t96 = -Ifges(6,6) * t244 + t140;
t89 = t189 * t151;
t83 = t173 * qJ(3);
t82 = t172 * qJ(3);
t81 = t152 * t187 - t283;
t73 = qJ(3) * t164;
t72 = qJ(3) * t163;
t71 = -mrSges(5,2) * t252 - mrSges(5,3) * t167;
t70 = mrSges(5,1) * t252 + (t219 - t222) * mrSges(5,3);
t67 = -mrSges(7,2) * t244 + mrSges(7,3) * t328;
t66 = mrSges(7,1) * t244 + mrSges(7,3) * t327;
t65 = mrSges(7,1) * t268 - mrSges(7,3) * t87;
t62 = -mrSges(7,2) * t268 - mrSges(7,3) * t171;
t59 = mrSges(6,1) * t165 + mrSges(6,2) * t166;
t58 = mrSges(7,1) * t328 + mrSges(7,2) * t327;
t55 = -t123 * t248 + (Ifges(5,5) * t152 + t156 * t187) * qJD(2);
t53 = -t121 * t239 + (Ifges(7,5) * t150 + t154 * t185) * qJD(5);
t52 = -t120 * t248 + (t156 * t184 + t282) * qJD(2);
t50 = -t118 * t239 + (Ifges(7,6) * t150 + t154 * t182) * qJD(5);
t49 = -Ifges(6,5) * t214 - Ifges(6,6) * t165 + t259;
t47 = Ifges(6,1) * t86 - Ifges(6,4) * t85 + Ifges(6,5) * t267;
t32 = mrSges(7,1) * t85 - mrSges(7,3) * t61;
t31 = -mrSges(7,2) * t85 + mrSges(7,3) * t60;
t27 = -mrSges(7,2) * t165 - mrSges(7,3) * t341;
t26 = mrSges(7,1) * t165 + mrSges(7,3) * t337;
t20 = mrSges(6,1) * t40 + mrSges(6,2) * t41;
t14 = Ifges(6,1) * t41 - Ifges(6,4) * t40 + Ifges(6,5) * t167;
t13 = -Ifges(7,1) * t337 - Ifges(7,4) * t341 + Ifges(7,5) * t165;
t11 = -Ifges(7,4) * t337 - Ifges(7,2) * t341 + Ifges(7,6) * t165;
t6 = mrSges(7,1) * t40 - mrSges(7,3) * t18;
t5 = -mrSges(7,2) * t40 + mrSges(7,3) * t17;
t3 = Ifges(7,1) * t18 + Ifges(7,4) * t17 + Ifges(7,5) * t40;
t2 = Ifges(7,4) * t18 + Ifges(7,2) * t17 + Ifges(7,6) * t40;
t4 = [t145 * qJD(3) * t237 + t86 * t14 + t17 * t22 + t18 * t23 + t60 * t2 + t61 * t3 + 0.2e1 * t7 * t31 + 0.2e1 * t8 * t32 + t41 * t47 - 0.2e1 * t73 * t5 + 0.2e1 * t72 * t6 + (t1 - t12) * t85 + (t21 - t45) * t40 + (t197 + t260) * t321 + (t147 * t340 - t179 + t231 + t260) * t320 + (t144 * t145 * t195 - t7 * t73 + t72 * t8 + t231) * t319 + ((t155 * t81 + (Ifges(3,4) - Ifges(4,5)) * t334 + t194 * t151 + t158 * qJ(3)) * qJD(2) + t258) * t156 + (t155 * t55 + (t10 - t52) * t151 + ((-t81 + t283) * t151 + t194 * t155) * qJD(4) + t158 * qJD(3) + ((Ifges(5,5) * t155 - Ifges(5,6) * t151 - 0.2e1 * Ifges(3,4) + 0.2e1 * Ifges(4,5)) * t152 + (-(2 * Ifges(4,3)) - (2 * Ifges(3,2)) - Ifges(5,3) + (2 * Ifges(3,1)) + (2 * Ifges(4,1)) + ((2 * m(4)) + 0.4e1 * (m(5) / 0.2e1 + t147 * t317) * t144) * t157) * t156) * qJD(2) + ((qJD(4) * t161 + 0.2e1 * t20 - 0.2e1 * t70) * t155 + (-0.2e1 * t279 - 0.2e1 * t71 + t295 * t234 + (-t154 * t201 + 0.2e1 * t280) * qJD(5)) * t151 + 0.2e1 * (-mrSges(4,1) * t152 + mrSges(4,3) * t156) * qJD(2) + 0.2e1 * t274 * t249) * qJ(3)) * t152; 0.2e1 * ((t147 - 0.1e1) * t317 * t144 + t338 * t148) * t217 - t337 * t23 / 0.2e1 + t2 * t310 + t41 * t311 + t13 * t314 + t11 * t315 + t86 * t316 + t3 * t309 + ((Ifges(4,4) + Ifges(3,5)) * t156 + (-mrSges(4,2) * qJ(3) - Ifges(3,6) + Ifges(4,6)) * t152) * qJD(2) - t341 * t22 / 0.2e1 + (-t141 / 0.2e1 + (qJD(3) * mrSges(4,2))) * t156 + m(7) * (-t24 * t73 + t25 * t72 + t7 * t83 + t8 * t82 + t273) + m(6) * t273 + t232 * t85 + t228 * t40 + t24 * t31 + t25 * t32 + t17 * t44 / 0.2e1 + t18 * t46 / 0.2e1 + t7 * t62 + t8 * t65 + t72 * t26 - t73 * t27 + t82 * t6 + t83 * t5 + (-t10 / 0.2e1 + t152 * t102 / 0.2e1 + t52 / 0.2e1 + (t282 / 0.2e1 + t156 * t123 / 0.2e1) * qJD(2) + (t152 * t89 - t324) * qJD(3) + (t81 / 0.2e1 + t47 * t296 + t224 * t152 + t229 * t150) * qJD(4) + (t89 * t251 + t152 * t59 + t279 + t71 - t295 * t150 + (t154 * t293 - t280) * qJD(5) + (t152 * t169 - t274) * qJD(4)) * qJ(3)) * t155 + (t14 * t296 + t55 / 0.2e1 + t233 * t150 + (-t78 / 0.2e1 + t43 / 0.2e1) * qJD(4) - t274 * qJD(3) + (-t150 * t47 / 0.2e1 + t229 * t154) * qJD(5) + (Ifges(5,6) * qJD(4) / 0.2e1 + (0.2e1 * t155 * t160 + t224) * qJD(2)) * t156 + (qJD(2) * Ifges(5,5) / 0.2e1 + t49 / 0.2e1 - t99 / 0.2e1 - qJD(4) * t123 / 0.2e1 + m(7) * t195 * t335 + t169 * qJD(3)) * t152 + (t20 - t70 + t169 * t251 + t324 * qJD(4) + (-qJD(4) * t89 - t277 - t294 * t150 + m(6) * (-0.2e1 * t147 + 0.2e1) * t250 + (-t154 * t272 + t270) * qJD(5)) * t152) * qJ(3)) * t151; -t171 * t11 + t87 * t13 + 0.2e1 * t24 * t62 + 0.2e1 * t25 * t65 + 0.2e1 * t82 * t26 + 0.2e1 * t83 * t27 - t337 * t46 - t341 * t44 + (t237 + t336 + 0.2e1 * (t144 + t148) * mrSges(5,3)) * qJD(3) + (t147 * t223 + t124 + t137) * t320 + (t148 * t195 + t24 * t83 + t25 * t82 + t124) * t319 + (t137 + t223) * t321 + (t59 * t318 + 0.2e1 * qJD(3) * t89 + t154 * t54 + t102 + (-t51 + t9) * t150 + (-t150 * t80 + t154 * t292) * qJD(5) + (t169 * t318 - t120 + t75) * qJD(4)) * t151 + (-t49 + t99 + (t150 * t198 + 0.2e1 * t263) * qJD(3) + (0.2e1 * t277 + t294 * t234 + (t154 * t198 - 0.2e1 * t270) * qJD(5)) * qJ(3) + (t150 * t292 + 0.4e1 * t151 * t160 + t154 * t80 + t318 * t89 + t123) * qJD(4)) * t155; t295 * t154 + (t155 * mrSges(5,1) - t151 * mrSges(5,2)) * t248 + (mrSges(4,2) + t190) * t251 + (-qJ(3) * t267 * t332 + (m(7) * (-t149 * t72 - t153 * t73) + t153 * t31 - t149 * t32 + t63) * t154) * qJD(5) + (m(7) * (t154 * t221 - t72 * t238 + t73 * t240 + t192 + (t152 * t205 + t154 * t206) * qJ(3)) - t31 * t240 + t153 * t5 - t32 * t238 - t149 * t6 + t28 + t293 * qJD(5)) * t150; -t294 * t154 + t190 * qJD(4) + (t271 * t332 + (m(7) * (-t149 * t82 + t153 * t83) + t153 * t62 - t149 * t65 + t104) * t154) * qJD(5) + (m(7) * (qJ(3) * t218 - t154 * t250 - t238 * t82 - t240 * t83 + t180) - t62 * t240 + t153 * t27 - t65 * t238 - t149 * t26 + t69 + t272 * qJD(5)) * t150; (-0.1e1 + t257) * t212 * t319; t50 * t315 + t53 * t314 + t72 * t66 - t73 * t67 + t17 * t313 + t18 * t312 + t86 * t305 + t7 * t103 + t8 * t106 + t41 * t301 + t227 * t85 + t225 * t40 + ((t47 / 0.2e1 + t175) * qJD(5) - t233) * t154 + (t2 * t300 + t14 / 0.2e1 + t3 * t297 + (t22 * t298 + t23 * t300) * qJD(6) + t229 * qJD(5)) * t150 + (t235 * t151 + (t151 * t168 - t155 * t276) * qJ(3)) * t251 + (((qJ(3) * t168 + t235) * qJD(4) + t170) * t155 + (-Ifges(5,5) * qJD(4) + t96 / 0.2e1 + t168 * qJD(3) + (t150 * t58 - t322) * qJ(3)) * t151) * t152 - t258; t141 - t341 * t313 - t337 * t312 + t82 * t66 + t83 * t67 + t50 * t310 + t53 * t309 + t24 * t103 + t25 * t106 + ((t311 + t174) * qJD(5) - t232) * t154 + (-t96 / 0.2e1 + t247 * t301 + t326 * qJD(3) + t322 * qJ(3)) * t155 + ((qJD(5) * t225 + t305) * t154 + (-qJ(3) * t326 + t235) * qJD(4) + t170) * t151 + (t11 * t300 + t316 + t13 * t297 + (t298 * t44 + t300 * t46) * qJD(6) + t228 * qJD(5) + (-qJ(3) * t58 + qJD(3) * t88 + qJD(4) * t225) * t155 + (-qJD(5) * t122 / 0.2e1 - t88 * t253 + t227) * t151) * t150; (t58 + (t103 * t153 - t106 * t149) * qJD(5)) * t154 + (qJD(5) * t88 - t149 * t66 + t153 * t67 + (-t103 * t149 - t106 * t153) * qJD(6)) * t150; (-t48 + t98 + (-t149 * t76 + t153 * t79 + t122) * qJD(5)) * t154 + (-t149 * t50 + t153 * t53 + t101 + (-t149 * t79 - t153 * t76) * qJD(6) + (-t119 + t74) * qJD(5)) * t150; t85 * t308 + t60 * t307 + t61 * t306 + t40 * t304 + t17 * t303 + t18 * t302 + t3 * t299 + t2 * t297 + t175 * qJD(6) + t151 * qJD(3) * t162 + ((t149 * t73 - t153 * t72) * qJD(6) + t192) * mrSges(7,3) + (t162 * t246 + ((qJD(5) * t152 * t275 + mrSges(6,2) * t251) * t154 + (-t152 * t203 + t251 * t275) * t150) * t151) * qJ(3) + t10; -t171 * t307 + t87 * t306 - t341 * t303 - t337 * t302 + t13 * t299 + t11 * t297 + t174 * qJD(6) + ((-t149 * t83 - t153 * t82) * qJD(6) + t180) * mrSges(7,3) + ((qJ(3) * t249 - t250) * mrSges(6,2) + (t151 * t236 - t271 * t275) * qJD(5)) * t154 + ((-Ifges(6,5) * qJD(5) + t253 * t275 + t308) * t151 + (qJ(3) * t203 - qJD(3) * t275 + qJD(4) * t236) * t155) * t150 + t259; -t154 * t93 + (-t329 + (mrSges(7,3) * t257 - mrSges(6,2)) * t154) * qJD(5); -t154 * t95 / 0.2e1 + t140 + (qJD(6) * t312 + t50 / 0.2e1 + t242 * t302) * t153 + (t53 / 0.2e1 - qJD(6) * t76 / 0.2e1 - t118 * t242 / 0.2e1) * t149 + (t100 * t297 + t97 * t300 + (t118 * t298 + t121 * t300) * qJD(6) + t236 * qJD(5)) * t150; t100 * t149 + t153 * t97 + (-t118 * t149 + t121 * t153) * qJD(6); mrSges(7,1) * t8 - mrSges(7,2) * t7 + t1; mrSges(7,1) * t25 - t24 * mrSges(7,2) + t9; t58; t48; t95; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t4(1) t4(2) t4(4) t4(7) t4(11) t4(16); t4(2) t4(3) t4(5) t4(8) t4(12) t4(17); t4(4) t4(5) t4(6) t4(9) t4(13) t4(18); t4(7) t4(8) t4(9) t4(10) t4(14) t4(19); t4(11) t4(12) t4(13) t4(14) t4(15) t4(20); t4(16) t4(17) t4(18) t4(19) t4(20) t4(21);];
Mq  = res;
