% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:09
% EndTime: 2019-12-31 20:32:29
% DurationCPUTime: 9.21s
% Computational Cost: add. (7059->533), mult. (17970->773), div. (0->0), fcn. (12882->8), ass. (0->243)
t223 = sin(pkin(9));
t224 = cos(pkin(9));
t226 = sin(qJ(4));
t229 = cos(qJ(4));
t190 = t223 * t229 + t224 * t226;
t230 = cos(qJ(2));
t237 = t190 * t230;
t155 = qJD(1) * t237;
t172 = t190 * qJD(4);
t270 = -t155 + t172;
t242 = t223 * t226 - t224 * t229;
t236 = t242 * t230;
t156 = qJD(1) * t236;
t171 = t242 * qJD(4);
t269 = -t156 + t171;
t227 = sin(qJ(2));
t244 = pkin(2) * t227 - qJ(3) * t230;
t192 = t244 * qJD(1);
t268 = qJD(1) * t227;
t255 = t223 * t268;
t147 = pkin(6) * t255 + t224 * t192;
t271 = t224 * t230;
t241 = pkin(3) * t227 - pkin(7) * t271;
t114 = qJD(1) * t241 + t147;
t173 = t223 * t192;
t272 = t224 * t227;
t273 = t223 * t230;
t238 = -pkin(6) * t272 - pkin(7) * t273;
t132 = qJD(1) * t238 + t173;
t298 = pkin(7) + qJ(3);
t197 = t298 * t223;
t198 = t298 * t224;
t262 = qJD(4) * t229;
t264 = qJD(3) * t224;
t265 = qJD(3) * t223;
t346 = -t197 * t262 + (-t132 + t264) * t229 + (-qJD(4) * t198 - t114 - t265) * t226;
t145 = -t226 * t197 + t229 * t198;
t345 = -t190 * qJD(3) - qJD(4) * t145 - t229 * t114 + t132 * t226;
t352 = -pkin(4) * t268 + t269 * pkin(8) + t345;
t351 = t270 * pkin(8) - t346;
t225 = sin(qJ(5));
t228 = cos(qJ(5));
t184 = t224 * qJD(2) - t255;
t254 = t224 * t268;
t261 = t223 * qJD(2);
t185 = t254 + t261;
t249 = t229 * t184 - t185 * t226;
t196 = -pkin(2) * t230 - t227 * qJ(3) - pkin(1);
t177 = t196 * qJD(1);
t267 = qJD(1) * t230;
t220 = pkin(6) * t267;
t203 = qJD(2) * qJ(3) + t220;
t134 = t224 * t177 - t223 * t203;
t257 = pkin(3) * t267;
t91 = -t185 * pkin(7) + t134 - t257;
t135 = t223 * t177 + t224 * t203;
t95 = pkin(7) * t184 + t135;
t48 = t226 * t91 + t229 * t95;
t40 = pkin(8) * t249 + t48;
t280 = t228 * t40;
t216 = qJD(4) - t267;
t129 = t184 * t226 + t185 * t229;
t47 = -t226 * t95 + t229 * t91;
t39 = -pkin(8) * t129 + t47;
t38 = pkin(4) * t216 + t39;
t10 = t225 * t38 + t280;
t260 = qJD(1) * qJD(2);
t252 = t227 * t260;
t343 = -t129 * t225 + t228 * t249;
t233 = qJD(2) * t236;
t84 = -qJD(1) * t233 + qJD(4) * t249;
t234 = qJD(2) * t237;
t85 = -qJD(1) * t234 - qJD(4) * t129;
t27 = qJD(5) * t343 + t225 * t85 + t228 * t84;
t68 = t129 * t228 + t225 * t249;
t28 = -qJD(5) * t68 - t225 * t84 + t228 * t85;
t259 = Ifges(6,5) * t27 + Ifges(6,6) * t28 + Ifges(6,3) * t252;
t305 = Ifges(6,4) * t68;
t209 = qJD(5) + t216;
t312 = -t209 / 0.2e1;
t327 = -t68 / 0.2e1;
t329 = -t343 / 0.2e1;
t169 = qJD(2) * t244 - t227 * qJD(3);
t157 = t169 * qJD(1);
t219 = pkin(6) * t268;
t194 = (qJD(3) - t219) * qJD(2);
t112 = t224 * t157 - t223 * t194;
t235 = t241 * qJD(2);
t88 = qJD(1) * t235 + t112;
t113 = t223 * t157 + t224 * t194;
t251 = t230 * t260;
t248 = t223 * t251;
t92 = -pkin(7) * t248 + t113;
t18 = -qJD(4) * t48 - t226 * t92 + t229 * t88;
t13 = pkin(4) * t252 - pkin(8) * t84 + t18;
t263 = qJD(4) * t226;
t17 = t226 * t88 + t229 * t92 + t91 * t262 - t263 * t95;
t14 = pkin(8) * t85 + t17;
t281 = t225 * t40;
t9 = t228 * t38 - t281;
t2 = qJD(5) * t9 + t13 * t225 + t14 * t228;
t3 = -qJD(5) * t10 + t13 * t228 - t14 * t225;
t338 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t62 = Ifges(6,4) * t343;
t34 = Ifges(6,1) * t68 + Ifges(6,5) * t209 + t62;
t195 = -qJD(2) * pkin(2) + qJD(3) + t219;
t146 = -pkin(3) * t184 + t195;
t81 = -pkin(4) * t249 + t146;
t350 = t259 + t338 + (Ifges(6,5) * t343 - Ifges(6,6) * t68) * t312 + (t10 * t68 + t343 * t9) * mrSges(6,3) + (-Ifges(6,2) * t68 + t34 + t62) * t329 - t81 * (mrSges(6,1) * t68 + mrSges(6,2) * t343) + (Ifges(6,1) * t343 - t305) * t327;
t349 = Ifges(3,5) / 0.2e1;
t144 = -t229 * t197 - t198 * t226;
t107 = -pkin(8) * t190 + t144;
t108 = -pkin(8) * t242 + t145;
t52 = t107 * t225 + t108 * t228;
t348 = -qJD(5) * t52 + t351 * t225 + t352 * t228;
t51 = t107 * t228 - t108 * t225;
t347 = qJD(5) * t51 + t352 * t225 - t351 * t228;
t178 = t223 * t257 + t220;
t344 = t270 * pkin(4) - t178;
t33 = Ifges(6,2) * t343 + Ifges(6,6) * t209 + t305;
t341 = t33 / 0.2e1;
t253 = -Ifges(3,6) * qJD(2) / 0.2e1;
t340 = qJD(2) * t349;
t183 = t224 * t196;
t133 = -pkin(7) * t272 + t183 + (-pkin(6) * t223 - pkin(3)) * t230;
t152 = pkin(6) * t271 + t223 * t196;
t274 = t223 * t227;
t141 = -pkin(7) * t274 + t152;
t75 = t226 * t133 + t229 * t141;
t333 = -t18 * mrSges(5,1) + t17 * mrSges(5,2);
t266 = qJD(2) * t227;
t218 = Ifges(3,4) * t267;
t288 = Ifges(4,2) * t223;
t292 = Ifges(4,4) * t224;
t245 = -t288 + t292;
t293 = Ifges(4,4) * t223;
t246 = Ifges(4,1) * t224 - t293;
t295 = mrSges(4,2) * t224;
t307 = t224 / 0.2e1;
t308 = -t223 / 0.2e1;
t332 = -(t134 * t224 + t135 * t223) * mrSges(4,3) + t195 * (mrSges(4,1) * t223 + t295) + Ifges(3,1) * t268 / 0.2e1 + t218 / 0.2e1 + t340 + t184 * t245 / 0.2e1 + t185 * t246 / 0.2e1 + (Ifges(4,4) * t185 + Ifges(4,2) * t184 - Ifges(4,6) * t267) * t308 + (Ifges(4,1) * t185 + Ifges(4,4) * t184 - Ifges(4,5) * t267) * t307;
t331 = t27 / 0.2e1;
t330 = t28 / 0.2e1;
t328 = t343 / 0.2e1;
t326 = t68 / 0.2e1;
t325 = t84 / 0.2e1;
t324 = t85 / 0.2e1;
t323 = pkin(1) * mrSges(3,1);
t322 = pkin(1) * mrSges(3,2);
t164 = t190 * t227;
t165 = t242 * t227;
t102 = -t164 * t228 + t165 * t225;
t320 = t102 / 0.2e1;
t103 = -t164 * t225 - t165 * t228;
t319 = t103 / 0.2e1;
t318 = -t249 / 0.2e1;
t317 = t249 / 0.2e1;
t316 = -t129 / 0.2e1;
t315 = t129 / 0.2e1;
t314 = -t164 / 0.2e1;
t313 = -t165 / 0.2e1;
t311 = t209 / 0.2e1;
t310 = -t216 / 0.2e1;
t309 = t216 / 0.2e1;
t130 = -t190 * t225 - t228 * t242;
t69 = qJD(5) * t130 - t171 * t228 - t172 * t225;
t94 = -t155 * t225 - t156 * t228;
t297 = t69 - t94;
t131 = t190 * t228 - t225 * t242;
t70 = -qJD(5) * t131 + t171 * t225 - t172 * t228;
t93 = -t155 * t228 + t156 * t225;
t296 = t70 - t93;
t294 = Ifges(3,4) * t227;
t291 = Ifges(5,4) * t129;
t289 = Ifges(4,5) * t224;
t286 = Ifges(4,6) * t223;
t277 = qJD(2) * mrSges(3,2);
t159 = mrSges(4,1) * t248 + t251 * t295;
t256 = pkin(6) * t266;
t139 = t224 * t169 + t223 * t256;
t215 = pkin(6) * t251;
t168 = pkin(3) * t248 + t215;
t179 = (pkin(3) * t261 + pkin(6) * qJD(2)) * t230;
t193 = pkin(3) * t274 + t227 * pkin(6);
t258 = Ifges(5,5) * t84 + Ifges(5,6) * t85 + Ifges(5,3) * t252;
t217 = -pkin(3) * t224 - pkin(2);
t45 = -t85 * mrSges(5,1) + t84 * mrSges(5,2);
t8 = -t28 * mrSges(6,1) + t27 * mrSges(6,2);
t74 = t229 * t133 - t226 * t141;
t247 = m(4) * t195 - qJD(2) * mrSges(3,1) - mrSges(4,1) * t184 + mrSges(4,2) * t185 + mrSges(3,3) * t268;
t50 = -pkin(4) * t230 + t165 * pkin(8) + t74;
t53 = -pkin(8) * t164 + t75;
t30 = -t225 * t53 + t228 * t50;
t31 = t225 * t50 + t228 * t53;
t239 = t289 / 0.2e1 - t286 / 0.2e1;
t104 = t235 + t139;
t160 = t223 * t169;
t115 = qJD(2) * t238 + t160;
t35 = t226 * t104 + t229 * t115 + t133 * t262 - t141 * t263;
t36 = -qJD(4) * t75 + t229 * t104 - t115 * t226;
t231 = t134 * mrSges(4,1) + t47 * mrSges(5,1) + t9 * mrSges(6,1) - Ifges(4,3) * t267 / 0.2e1 + Ifges(4,6) * t184 + Ifges(4,5) * t185 + t253 - (t230 * Ifges(3,2) + t294) * qJD(1) / 0.2e1 + t209 * Ifges(6,3) + t68 * Ifges(6,5) + t343 * Ifges(6,6) + t216 * Ifges(5,3) + t129 * Ifges(5,5) + t249 * Ifges(5,6) - t10 * mrSges(6,2) - t135 * mrSges(4,2) - t48 * mrSges(5,2);
t204 = mrSges(3,3) * t267 - t277;
t167 = (mrSges(4,1) * t227 - mrSges(4,3) * t271) * t260;
t166 = (-mrSges(4,2) * t227 - mrSges(4,3) * t273) * t260;
t158 = pkin(4) * t242 + t217;
t154 = -mrSges(4,1) * t267 - t185 * mrSges(4,3);
t153 = mrSges(4,2) * t267 + t184 * mrSges(4,3);
t151 = -pkin(6) * t273 + t183;
t148 = -pkin(6) * t254 + t173;
t143 = (Ifges(4,5) * t227 + t230 * t246) * t260;
t142 = (Ifges(4,6) * t227 + t230 * t245) * t260;
t140 = -t224 * t256 + t160;
t137 = pkin(4) * t164 + t193;
t125 = Ifges(5,4) * t249;
t110 = t171 * t227 - t234;
t109 = -t172 * t227 - t233;
t101 = mrSges(5,1) * t216 - mrSges(5,3) * t129;
t100 = -mrSges(5,2) * t216 + mrSges(5,3) * t249;
t86 = -pkin(4) * t110 + t179;
t73 = -mrSges(5,2) * t252 + mrSges(5,3) * t85;
t72 = mrSges(5,1) * t252 - mrSges(5,3) * t84;
t71 = -mrSges(5,1) * t249 + mrSges(5,2) * t129;
t61 = -pkin(4) * t85 + t168;
t58 = t129 * Ifges(5,1) + t216 * Ifges(5,5) + t125;
t57 = Ifges(5,2) * t249 + t216 * Ifges(5,6) + t291;
t55 = mrSges(6,1) * t209 - mrSges(6,3) * t68;
t54 = -mrSges(6,2) * t209 + mrSges(6,3) * t343;
t44 = t84 * Ifges(5,1) + t85 * Ifges(5,4) + Ifges(5,5) * t252;
t43 = t84 * Ifges(5,4) + t85 * Ifges(5,2) + Ifges(5,6) * t252;
t42 = -qJD(5) * t103 - t109 * t225 + t110 * t228;
t41 = qJD(5) * t102 + t109 * t228 + t110 * t225;
t37 = -mrSges(6,1) * t343 + mrSges(6,2) * t68;
t29 = pkin(8) * t110 + t35;
t23 = pkin(4) * t266 - pkin(8) * t109 + t36;
t22 = -mrSges(6,2) * t252 + mrSges(6,3) * t28;
t21 = mrSges(6,1) * t252 - mrSges(6,3) * t27;
t12 = t228 * t39 - t281;
t11 = -t225 * t39 - t280;
t7 = t27 * Ifges(6,1) + t28 * Ifges(6,4) + Ifges(6,5) * t252;
t6 = t27 * Ifges(6,4) + t28 * Ifges(6,2) + Ifges(6,6) * t252;
t5 = -qJD(5) * t31 - t225 * t29 + t228 * t23;
t4 = qJD(5) * t30 + t225 * t23 + t228 * t29;
t1 = [t30 * t21 + t31 * t22 + t193 * t45 + (Ifges(6,4) * t103 + Ifges(6,2) * t102) * t330 + ((-pkin(6) * t204 + t231 + t253) * qJD(2) + t142 * t308 + t143 * t307 + pkin(6) * t159 + (-t112 * t224 - t113 * t223) * mrSges(4,3) + (-0.2e1 * t323 + Ifges(5,5) * t313 + Ifges(5,6) * t314 + Ifges(6,5) * t319 + Ifges(6,6) * t320 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t239) * t227) * t260) * t227 + (t10 * t42 + t102 * t2 - t103 * t3 - t41 * t9) * mrSges(6,3) + t42 * t341 + (Ifges(6,1) * t103 + Ifges(6,4) * t102) * t331 + m(4) * (t112 * t151 + t113 * t152 + t134 * t139 + t135 * t140) + t179 * t71 + t152 * t166 + t151 * t167 + t140 * t153 + t139 * t154 + t137 * t8 + t146 * (-mrSges(5,1) * t110 + mrSges(5,2) * t109) + t109 * t58 / 0.2e1 + t110 * t57 / 0.2e1 + t35 * t100 + t36 * t101 + t61 * (-mrSges(6,1) * t102 + mrSges(6,2) * t103) + t81 * (-mrSges(6,1) * t42 + mrSges(6,2) * t41) + t86 * t37 + t74 * t72 + t75 * t73 + t4 * t54 + t5 * t55 + m(6) * (t10 * t4 + t137 * t61 + t2 * t31 + t3 * t30 + t5 * t9 + t81 * t86) + m(5) * (t146 * t179 + t168 * t193 + t17 * t75 + t18 * t74 + t35 * t48 + t36 * t47) + (-t112 * mrSges(4,1) + t113 * mrSges(4,2) - Ifges(5,6) * t324 - Ifges(5,5) * t325 - Ifges(6,6) * t330 - Ifges(6,5) * t331 + (t247 * pkin(6) + t332 + t340) * qJD(2) + (-0.2e1 * t322 + (-0.3e1 / 0.2e1 * t289 + 0.3e1 / 0.2e1 * t286 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t230 + (-0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + Ifges(4,1) * t224 ^ 2 / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(6,3) / 0.2e1 + (m(4) * pkin(6) + t295) * pkin(6) + (pkin(6) * mrSges(4,1) - t292 + t288 / 0.2e1) * t223) * t227) * t260 + t333 - t338) * t230 - (t259 + t258) * t230 / 0.2e1 + (-Ifges(5,4) * t165 - Ifges(5,2) * t164) * t324 + t168 * (mrSges(5,1) * t164 - mrSges(5,2) * t165) + (-t109 * t47 + t110 * t48 - t164 * t17 + t165 * t18) * mrSges(5,3) + (-Ifges(5,1) * t165 - Ifges(5,4) * t164) * t325 + t41 * t34 / 0.2e1 + (Ifges(5,5) * t109 + Ifges(5,6) * t110) * t309 + (Ifges(6,5) * t41 + Ifges(6,6) * t42) * t311 + t44 * t313 + t43 * t314 + (Ifges(5,1) * t109 + Ifges(5,4) * t110) * t315 + (Ifges(5,4) * t109 + Ifges(5,2) * t110) * t317 + t7 * t319 + t6 * t320 + (Ifges(6,1) * t41 + Ifges(6,4) * t42) * t326 + (Ifges(6,4) * t41 + Ifges(6,2) * t42) * t328; (-t93 / 0.2e1 + t70 / 0.2e1) * t33 + t346 * t100 + (t144 * t18 + t145 * t17 - t146 * t178 + t168 * t217 + t345 * t47 + t346 * t48) * m(5) + t347 * t54 + (t347 * t10 + t158 * t61 + t2 * t52 + t3 * t51 + t344 * t81 + t348 * t9) * m(6) + t348 * t55 + m(4) * (-t134 * t265 + t135 * t264 + (-t112 * t223 + t113 * t224) * qJ(3)) + t217 * t45 + t190 * t44 / 0.2e1 + ((-t218 / 0.2e1 + (t230 * t239 + t322) * qJD(1) + (t349 + (Ifges(4,1) * t223 + t292) * t307 + (Ifges(4,2) * t224 + t293) * t308) * qJD(2) + ((-m(4) * pkin(2) - mrSges(4,1) * t224 + mrSges(4,2) * t223 - mrSges(3,1)) * qJD(2) - t247) * pkin(6) - t332) * t230 + ((t323 + t294 / 0.2e1 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t230) * qJD(1) - t231 + (t204 + t277) * pkin(6) + t253) * t227 + (Ifges(4,5) * t223 + Ifges(5,5) * t190 + Ifges(6,5) * t131 + Ifges(4,6) * t224 - Ifges(5,6) * t242 + Ifges(6,6) * t130) * t266 / 0.2e1) * qJD(1) - m(4) * (t134 * t147 + t135 * t148) + t345 * t101 + (-t94 / 0.2e1 + t69 / 0.2e1) * t34 - t178 * t71 + t158 * t8 - pkin(2) * t159 - t148 * t153 - t147 * t154 + t144 * t72 + t145 * t73 + t130 * t6 / 0.2e1 + t61 * (-mrSges(6,1) * t130 + mrSges(6,2) * t131) + t131 * t7 / 0.2e1 + t51 * t21 + t52 * t22 + (t142 / 0.2e1 + t113 * mrSges(4,3) + qJD(3) * t153 + qJ(3) * t166) * t224 + (t143 / 0.2e1 - qJD(3) * t154 - qJ(3) * t167 - t112 * mrSges(4,3)) * t223 + (-Ifges(5,1) * t156 - Ifges(5,4) * t155) * t316 + (-Ifges(5,4) * t156 - Ifges(5,2) * t155) * t318 + (-t171 / 0.2e1 + t156 / 0.2e1) * t58 - t242 * t43 / 0.2e1 + t168 * (mrSges(5,1) * t242 + mrSges(5,2) * t190) + (Ifges(5,4) * t190 - Ifges(5,2) * t242) * t324 + (Ifges(5,1) * t190 - Ifges(5,4) * t242) * t325 + (-t17 * t242 - t18 * t190 + t269 * t47 - t270 * t48) * mrSges(5,3) + (-t172 / 0.2e1 + t155 / 0.2e1) * t57 + (-Ifges(5,5) * t156 - Ifges(5,6) * t155) * t310 + t344 * t37 + (-Ifges(5,5) * t171 - Ifges(5,6) * t172) * t309 + (-Ifges(5,1) * t171 - Ifges(5,4) * t172) * t315 + (-Ifges(5,4) * t171 - Ifges(5,2) * t172) * t317 + (Ifges(6,5) * t69 + Ifges(6,6) * t70) * t311 + (Ifges(6,5) * t94 + Ifges(6,6) * t93) * t312 + (Ifges(6,1) * t69 + Ifges(6,4) * t70) * t326 + (Ifges(6,1) * t94 + Ifges(6,4) * t93) * t327 + (Ifges(6,4) * t69 + Ifges(6,2) * t70) * t328 + (Ifges(6,4) * t94 + Ifges(6,2) * t93) * t329 + (Ifges(6,4) * t131 + Ifges(6,2) * t130) * t330 + (Ifges(6,1) * t131 + Ifges(6,4) * t130) * t331 + (mrSges(5,1) * t270 - mrSges(5,2) * t269) * t146 + (-mrSges(6,1) * t296 + mrSges(6,2) * t297) * t81 + (t10 * t296 + t130 * t2 - t131 * t3 - t297 * t9) * mrSges(6,3); -t249 * t100 + t129 * t101 - t184 * t153 + t185 * t154 - t343 * t54 + t68 * t55 + t159 + t45 + t8 + (-t10 * t343 + t68 * t9 + t61) * m(6) + (t129 * t47 - t249 * t48 + t168) * m(5) + (t134 * t185 - t135 * t184 + t215) * m(4); -t333 + (-t129 * t37 + t228 * t21 + t225 * t22 + (-t225 * t55 + t228 * t54) * qJD(5) + (-t129 * t81 + t2 * t225 + t228 * t3 + (t10 * t228 - t225 * t9) * qJD(5)) * m(6)) * pkin(4) - m(6) * (t10 * t12 + t11 * t9) + t258 + (-Ifges(5,2) * t129 + t125 + t58) * t318 - t47 * t100 + t48 * t101 - t12 * t54 - t11 * t55 + t68 * t341 + (Ifges(5,5) * t249 - Ifges(5,6) * t129) * t310 - t146 * (mrSges(5,1) * t129 + mrSges(5,2) * t249) + (Ifges(5,1) * t249 - t291) * t316 + (t129 * t48 + t249 * t47) * mrSges(5,3) + t57 * t315 + t350; t10 * t55 + t33 * t326 - t9 * t54 + t350;];
tauc = t1(:);
