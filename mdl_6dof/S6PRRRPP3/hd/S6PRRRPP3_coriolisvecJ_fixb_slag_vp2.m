% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:38
% EndTime: 2019-03-08 22:54:54
% DurationCPUTime: 8.73s
% Computational Cost: add. (3648->566), mult. (9387->703), div. (0->0), fcn. (5937->8), ass. (0->262)
t354 = -Ifges(5,4) + Ifges(7,6);
t352 = Ifges(5,1) + Ifges(7,3);
t348 = Ifges(7,4) + Ifges(6,5);
t347 = Ifges(5,5) + Ifges(7,5);
t351 = Ifges(7,2) + Ifges(6,3);
t181 = cos(qJ(3));
t259 = qJD(2) * t181;
t168 = Ifges(4,4) * t259;
t353 = -t168 / 0.2e1;
t177 = sin(qJ(4));
t350 = t354 * t177;
t178 = sin(qJ(3));
t261 = qJD(2) * t178;
t349 = t261 / 0.2e1;
t237 = Ifges(4,5) * qJD(3) / 0.2e1;
t253 = qJD(2) * qJD(3);
t235 = t178 * t253;
t180 = cos(qJ(4));
t234 = t181 * t253;
t256 = qJD(4) * t177;
t239 = t178 * t256;
t252 = qJD(3) * qJD(4);
t96 = qJD(2) * t239 + (-t234 - t252) * t180;
t255 = qJD(4) * t180;
t257 = qJD(3) * t181;
t336 = t177 * t257 + t178 * t255;
t97 = t336 * qJD(2) + t177 * t252;
t346 = t351 * t97 + (Ifges(6,6) - Ifges(7,6)) * t96 + t348 * t235;
t345 = t347 * t235 - t352 * t96 + t354 * t97;
t276 = qJ(5) * t180;
t198 = qJ(6) * t177 - t276;
t195 = t198 * t181;
t231 = pkin(4) * t256 - qJD(5) * t177;
t179 = sin(qJ(2));
t174 = sin(pkin(6));
t263 = qJD(1) * t174;
t245 = t179 * t263;
t137 = qJD(2) * pkin(8) + t245;
t175 = cos(pkin(6));
t262 = qJD(1) * t175;
t108 = t181 * t137 + t178 * t262;
t247 = t177 * pkin(4) * t259 + t108;
t344 = -qJD(2) * t195 + qJD(4) * t198 - qJD(6) * t180 + t231 - t247;
t144 = -pkin(3) * t181 - pkin(9) * t178 - pkin(2);
t194 = t144 * qJD(2);
t182 = cos(qJ(2));
t244 = t182 * t263;
t190 = t194 - t244;
t274 = t174 * t182;
t242 = qJD(2) * t274;
t192 = qJD(1) * (qJD(3) * t175 + t242);
t258 = qJD(3) * t178;
t58 = -t137 * t258 + t181 * t192;
t343 = qJD(4) * t190 + t58;
t316 = pkin(5) + pkin(9);
t148 = t316 * t180;
t270 = t180 * t181;
t250 = pkin(5) * t270;
t299 = pkin(4) + qJ(6);
t107 = -t178 * t137 + t181 * t262;
t225 = pkin(3) * t178 - pkin(9) * t181;
t133 = t225 * qJD(2);
t44 = -t177 * t107 + t133 * t180;
t342 = -(-t178 * t299 + t250) * qJD(2) + t44 + qJD(4) * t148;
t272 = t177 * t181;
t251 = pkin(5) * t272;
t45 = t180 * t107 + t177 * t133;
t341 = -(qJ(5) * t178 - t251) * qJD(2) - t45 - t316 * t256;
t286 = qJD(2) * pkin(2);
t138 = -t244 - t286;
t162 = -qJD(4) + t259;
t229 = t177 * t244;
t85 = qJD(3) * pkin(9) + t108;
t31 = t177 * t194 + t180 * t85 - t229;
t20 = t162 * qJ(5) - t31;
t254 = t180 * qJD(3);
t130 = t177 * t261 - t254;
t334 = -t130 * pkin(5) + qJD(6);
t10 = -t20 + t334;
t131 = qJD(3) * t177 + t180 * t261;
t84 = -qJD(3) * pkin(3) - t107;
t191 = -qJ(5) * t131 + t84;
t18 = t130 * t299 + t191;
t30 = t177 * t85 - t180 * t190;
t199 = t31 * t177 - t30 * t180;
t328 = -qJD(5) - t30;
t19 = pkin(4) * t162 - t328;
t200 = t177 * t20 + t180 * t19;
t291 = Ifges(6,6) * t177;
t208 = -Ifges(6,2) * t180 + t291;
t292 = Ifges(5,4) * t180;
t213 = -Ifges(5,2) * t177 + t292;
t217 = -mrSges(7,2) * t180 + mrSges(7,3) * t177;
t218 = -mrSges(6,2) * t177 - mrSges(6,3) * t180;
t220 = mrSges(5,1) * t177 + mrSges(5,2) * t180;
t306 = t180 / 0.2e1;
t307 = -t180 / 0.2e1;
t308 = t177 / 0.2e1;
t309 = -t177 / 0.2e1;
t311 = -t162 / 0.2e1;
t312 = t131 / 0.2e1;
t313 = -t131 / 0.2e1;
t314 = t130 / 0.2e1;
t315 = -t130 / 0.2e1;
t32 = pkin(4) * t130 + t191;
t323 = -Ifges(5,6) + t348;
t324 = Ifges(6,4) - t347;
t126 = Ifges(7,6) * t131;
t282 = t131 * Ifges(6,6);
t330 = t351 * t130 - t348 * t162 + t126 - t282;
t128 = Ifges(5,4) * t130;
t289 = Ifges(7,6) * t130;
t331 = t352 * t131 - t347 * t162 - t128 + t289;
t287 = Ifges(7,6) * t180;
t290 = Ifges(6,6) * t180;
t338 = t351 * t177 + t287 - t290;
t339 = t352 * t180 + t350;
t127 = Ifges(6,6) * t130;
t53 = -Ifges(6,4) * t162 - Ifges(6,2) * t131 + t127;
t283 = t131 * Ifges(5,4);
t54 = -Ifges(5,2) * t130 - Ifges(5,6) * t162 + t283;
t197 = pkin(5) * t131 + t30;
t335 = qJD(5) + t197;
t9 = t162 * t299 + t335;
t184 = t200 * mrSges(6,1) - (t10 * t177 - t180 * t9) * mrSges(7,1) - t199 * mrSges(5,3) + t18 * t217 + t208 * t313 + t213 * t315 + t32 * t218 + t84 * t220 + t53 * t307 + t54 * t309 + t338 * t314 + t339 * t312 + t330 * t308 + t331 * t306 + (t177 * t323 - t180 * t324) * t311;
t340 = -t138 * mrSges(4,2) + t107 * mrSges(4,3) - Ifges(4,1) * t349 - t184 - t237 + t353;
t337 = t290 + t292 + (Ifges(5,1) + Ifges(6,2)) * t177;
t333 = t291 - t350 + (Ifges(5,2) + t351) * t180;
t236 = -Ifges(4,6) * qJD(3) / 0.2e1;
t332 = pkin(8) * (t178 * t254 + t181 * t256);
t275 = t174 * t179;
t114 = -t175 * t181 + t178 * t275;
t74 = -qJD(3) * t114 + t181 * t242;
t329 = -qJD(4) * t274 + t74;
t327 = -t31 - t334;
t136 = t225 * qJD(3);
t100 = (t136 + t245) * qJD(2);
t6 = t177 * t100 + t343 * t180 - t256 * t85;
t7 = t100 * t180 - t343 * t177 - t85 * t255;
t326 = -t177 * t7 + t180 * t6;
t4 = -qJ(5) * t235 + qJD(5) * t162 - t6;
t5 = -pkin(4) * t235 - t7;
t325 = t177 * t5 - t180 * t4;
t226 = -Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t227 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t228 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,1) / 0.2e1;
t321 = -t227 * t130 + t226 * t131 + t228 * t162 - t10 * mrSges(7,2) - t138 * mrSges(4,1) - t19 * mrSges(6,2) - t236 + (Ifges(4,4) * t178 + Ifges(4,2) * t181) * qJD(2) / 0.2e1 - Ifges(5,6) * t315 - Ifges(6,4) * t313 + t108 * mrSges(4,3) + t20 * mrSges(6,3) + t30 * mrSges(5,1) + t31 * mrSges(5,2) + t9 * mrSges(7,3) - t348 * t314 - t347 * t312 - (Ifges(5,3) + Ifges(7,1) + Ifges(6,1)) * t311;
t320 = -t96 / 0.2e1;
t319 = t96 / 0.2e1;
t318 = -t97 / 0.2e1;
t310 = t162 / 0.2e1;
t305 = pkin(8) * t177;
t62 = mrSges(6,1) * t97 - mrSges(6,3) * t235;
t66 = -mrSges(5,2) * t235 - mrSges(5,3) * t97;
t298 = -t62 + t66;
t64 = -t96 * mrSges(6,1) + mrSges(6,2) * t235;
t65 = mrSges(5,1) * t235 + mrSges(5,3) * t96;
t297 = t64 - t65;
t69 = -mrSges(7,2) * t131 + mrSges(7,3) * t130;
t72 = -mrSges(6,2) * t130 - mrSges(6,3) * t131;
t296 = t69 + t72;
t295 = mrSges(5,3) * t130;
t294 = mrSges(5,3) * t131;
t59 = t137 * t257 + t178 * t192;
t284 = t114 * t59;
t280 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t130 + mrSges(5,2) * t131 + mrSges(4,3) * t261;
t277 = qJ(5) * t130;
t273 = t177 * t178;
t271 = t178 * t180;
t269 = t180 * t182;
t101 = mrSges(5,2) * t162 - t295;
t104 = mrSges(6,1) * t130 + mrSges(6,3) * t162;
t268 = t101 - t104;
t102 = -mrSges(5,1) * t162 - t294;
t106 = mrSges(6,1) * t131 - mrSges(6,2) * t162;
t267 = t102 - t106;
t105 = -mrSges(7,1) * t130 - mrSges(7,2) * t162;
t266 = -t104 + t105;
t265 = t177 * t136 + t144 * t255;
t264 = pkin(4) * t273 + t178 * pkin(8);
t165 = pkin(8) * t270;
t110 = t177 * t144 + t165;
t260 = qJD(2) * t179;
t164 = pkin(8) * t272;
t249 = t101 + t266;
t103 = mrSges(7,1) * t131 + mrSges(7,3) * t162;
t248 = -t103 + t267;
t246 = -pkin(4) - t305;
t243 = t174 * t260;
t146 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t259;
t233 = -m(4) * t108 - t146;
t63 = -t97 * mrSges(7,1) + mrSges(7,2) * t235;
t232 = -qJ(5) * t177 - pkin(3);
t109 = t144 * t180 - t164;
t230 = t336 * pkin(4) + pkin(8) * t257 + qJ(5) * t239;
t94 = qJ(5) * t181 - t110;
t224 = qJD(4) * t165 - t136 * t180 + t144 * t256;
t222 = -qJD(5) * t181 + t265;
t221 = mrSges(5,1) * t180 - mrSges(5,2) * t177;
t219 = mrSges(6,2) * t180 - mrSges(6,3) * t177;
t216 = mrSges(7,2) * t177 + mrSges(7,3) * t180;
t211 = Ifges(6,4) * t177 + Ifges(6,5) * t180;
t210 = Ifges(7,4) * t180 - Ifges(7,5) * t177;
t209 = Ifges(5,5) * t177 + Ifges(5,6) * t180;
t202 = -Ifges(7,3) * t177 + t287;
t61 = -t96 * mrSges(7,1) - mrSges(7,3) * t235;
t115 = t175 * t178 + t181 * t275;
t196 = -m(4) * t107 + m(5) * t84 + t280;
t1 = -pkin(5) * t96 + qJD(6) * t162 - t235 * t299 - t7;
t2 = -pkin(5) * t97 - t4;
t189 = -t7 * mrSges(5,1) + t6 * mrSges(5,2) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3);
t187 = qJ(5) * t96 - qJD(5) * t131 + t59;
t183 = qJD(2) ^ 2;
t173 = t181 * pkin(4);
t161 = Ifges(6,1) * t235;
t160 = Ifges(7,1) * t235;
t159 = Ifges(5,3) * t235;
t147 = t316 * t177;
t139 = -pkin(4) * t180 + t232;
t132 = (-mrSges(4,1) * t181 + mrSges(4,2) * t178) * qJD(2);
t121 = -t180 * t299 + t232;
t120 = (mrSges(4,1) * t178 + mrSges(4,2) * t181) * t253;
t113 = -qJ(5) * t255 + t231;
t112 = -qJ(5) * t271 + t264;
t99 = (t177 * t179 + t181 * t269) * t263;
t98 = -t180 * t245 + t181 * t229;
t95 = -t109 + t173;
t92 = Ifges(6,4) * t96;
t91 = Ifges(7,4) * t97;
t90 = Ifges(5,5) * t96;
t89 = Ifges(6,5) * t97;
t88 = Ifges(7,5) * t96;
t87 = Ifges(5,6) * t97;
t83 = t178 * t198 + t264;
t76 = t115 * t180 - t177 * t274;
t75 = t115 * t177 + t174 * t269;
t73 = qJD(3) * t115 + t178 * t242;
t70 = pkin(4) * t131 + t277;
t67 = -pkin(5) * t273 - t94;
t60 = qJ(6) * t181 + t164 + t173 + (pkin(5) * t178 - t144) * t180;
t48 = -t259 * t276 + t247;
t43 = t131 * t299 + t277;
t42 = t258 * t305 - t224;
t41 = t265 - t332;
t40 = (-qJ(5) * t257 - qJD(5) * t178) * t180 + t230;
t39 = -pkin(4) * t261 - t44;
t38 = -qJ(5) * t261 - t45;
t36 = t246 * t258 + t224;
t35 = mrSges(7,2) * t96 + mrSges(7,3) * t97;
t34 = mrSges(5,1) * t97 - mrSges(5,2) * t96;
t33 = -mrSges(6,2) * t97 + mrSges(6,3) * t96;
t28 = -qJ(5) * t258 - t222 + t332;
t26 = -t96 * Ifges(5,4) - t97 * Ifges(5,2) + Ifges(5,6) * t235;
t25 = Ifges(6,4) * t235 + t96 * Ifges(6,2) + t97 * Ifges(6,6);
t15 = qJD(3) * t195 + (qJD(6) * t177 + (qJ(6) * qJD(4) - qJD(5)) * t180) * t178 + t230;
t14 = (-pkin(5) * t271 - t164) * qJD(4) + (-t251 + (-pkin(8) * t180 + qJ(5)) * t178) * qJD(3) + t222;
t13 = -t115 * t256 + t177 * t243 + t180 * t329;
t12 = t115 * t255 + t177 * t329 - t180 * t243;
t11 = -pkin(5) * t239 + qJD(6) * t181 + (t250 + (-qJ(6) + t246) * t178) * qJD(3) + t224;
t8 = pkin(4) * t97 + t187;
t3 = qJD(6) * t130 + t299 * t97 + t187;
t16 = [-t115 * mrSges(4,3) * t235 + t74 * t146 + (t63 + t298) * t76 + (t61 + t297) * t75 + t249 * t13 - t248 * t12 + (t280 + t296) * t73 + ((-mrSges(3,2) * t183 - t120) * t182 + (-mrSges(3,1) * t183 + qJD(2) * t132) * t179) * t174 + (mrSges(4,3) * t234 + t33 + t34 + t35) * t114 + m(7) * (t1 * t75 + t10 * t13 + t114 * t3 + t12 * t9 + t18 * t73 + t2 * t76) + m(6) * (t114 * t8 + t12 * t19 - t13 * t20 + t32 * t73 - t4 * t76 + t5 * t75) + m(5) * (t12 * t30 + t13 * t31 + t6 * t76 - t7 * t75 + t73 * t84 + t284) + m(4) * (-t107 * t73 + t108 * t74 + t284 + t115 * t58 + (t138 - t244) * t243); m(6) * (t112 * t8 + t19 * t36 + t20 * t28 + t32 * t40 + t4 * t94 + t5 * t95) + m(7) * (t1 * t60 + t10 * t14 + t11 * t9 + t15 * t18 + t2 * t67 + t3 * t83) + t95 * t64 + t83 * t35 + t94 * t62 + t40 * t72 + t60 * t61 + t67 * t63 + t15 * t69 + (t8 * t218 + t3 * t217 + t208 * t319 + t213 * t318 + t26 * t309 + t25 * t307 + (-t177 * t6 - t180 * t7) * mrSges(5,3) + (t1 * t180 - t177 * t2) * mrSges(7,1) + (t177 * t4 + t180 * t5) * mrSges(6,1) + (mrSges(4,3) + t220) * t59 + (mrSges(4,2) * t260 + (-m(6) * t32 - m(7) * t18 - t196 - t296) * t182) * t263 + (((-0.3e1 / 0.2e1 * Ifges(4,4) - t226 * t180 + t227 * t177) * t178 + (-0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(4,1) - t228) * t181) * qJD(2) + t236 - t321) * qJD(3) + (t54 * t307 + t202 * t312 + t209 * t310 + t84 * t221 - t32 * t219 + t18 * t216 + (-t177 * t30 - t180 * t31) * mrSges(5,3) + (-t10 * t180 - t177 * t9) * mrSges(7,1) + (-t177 * t19 + t180 * t20) * mrSges(6,1) + t337 * t313 + (t211 + t210) * t311 + t331 * t309 + t333 * t314) * qJD(4) + t339 * t320 + t338 * t97 / 0.2e1 + (qJD(4) * t53 + t346) * t308 + (qJD(4) * t330 + t345) * t306 + (t34 + (m(5) + m(4)) * t59 + t233 * qJD(3)) * pkin(8)) * t178 + m(5) * (t109 * t7 + t110 * t6 - t30 * t42 + t31 * t41) + (0.2e1 * (-t138 / 0.2e1 - t286 / 0.2e1) * m(4) - t132) * t245 + (-t159 / 0.2e1 - t160 / 0.2e1 - t161 / 0.2e1 + (-mrSges(4,1) * t260 + t182 * t233) * t263 + (t196 * pkin(8) + t237 + 0.3e1 / 0.2e1 * t168 - t340) * qJD(3) - t227 * t97 - t226 * t96 + (m(4) * pkin(8) + mrSges(4,3)) * t58 + t189 + t87 / 0.2e1 + t88 / 0.2e1 - t89 / 0.2e1 + t90 / 0.2e1 - t91 / 0.2e1 - t92 / 0.2e1) * t181 - m(7) * (t10 * t99 + t9 * t98) - m(5) * (t30 * t98 + t31 * t99) - m(6) * (t19 * t98 - t20 * t99) + t248 * t98 - t249 * t99 + t41 * t101 + t42 * t102 + t11 * t103 + t28 * t104 + t14 * t105 + t36 * t106 + t109 * t65 + t110 * t66 + t112 * t33 - pkin(2) * t120; -t58 * mrSges(4,2) - pkin(3) * t34 + t341 * t105 + t342 * t103 + t344 * t69 + (t1 * t147 + t341 * t10 + t121 * t3 + t148 * t2 + t344 * t18 + t342 * t9) * m(7) + t345 * t308 + t346 * t307 - m(5) * (t108 * t84 - t30 * t44 + t31 * t45) - m(6) * (t19 * t39 + t20 * t38 + t32 * t48) + t202 * t319 + t8 * t219 + (-mrSges(4,1) - t221) * t59 - t3 * t216 + t333 * t318 + m(6) * (pkin(9) * t325 + t113 * t32 + t139 * t8) + t325 * mrSges(6,1) + m(5) * (-pkin(3) * t59 + pkin(9) * t326) + t326 * mrSges(5,3) + t337 * t320 + t26 * t306 + t25 * t309 + (t177 * t297 + t180 * t298) * pkin(9) - t280 * t108 + ((-m(5) * t199 + m(6) * t200 - t177 * t268 - t180 * t267) * pkin(9) + t184) * qJD(4) + (-t48 + t113) * t72 - t45 * t101 - t44 * t102 - t38 * t104 - t39 * t106 + t121 * t35 + t139 * t33 - t107 * t146 + t147 * t61 + t148 * t63 + ((Ifges(4,4) * t349 + t236 + t321) * t178 + (t353 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t261 + t237 + t340) * t181 + (-t211 / 0.2e1 + t209 / 0.2e1 - t210 / 0.2e1) * t258) * qJD(2) + (t1 * t177 + t180 * t2) * mrSges(7,1); t159 + t160 + t161 - t70 * t72 - pkin(4) * t64 - t43 * t69 + (-t62 + t63) * qJ(5) + t327 * t103 + (-pkin(4) * t5 - qJ(5) * t4 - t19 * t31 + t20 * t328 - t32 * t70) * m(6) + (Ifges(6,2) * t130 + t282 + t54) * t312 + t197 * t105 - t299 * t61 + (-Ifges(5,2) * t131 - t128 + t331) * t314 + (t10 * t131 + t130 * t9) * mrSges(7,1) + (t130 * t19 - t131 * t20) * mrSges(6,1) + (qJ(5) * t2 - t1 * t299 + t335 * t10 - t18 * t43 + t327 * t9) * m(7) - t189 + (t130 * t324 + t131 * t323) * t310 + (t267 + t294) * t31 + (t268 + t295) * t30 + t266 * qJD(5) - t32 * (-mrSges(6,2) * t131 + mrSges(6,3) * t130) - t18 * (mrSges(7,2) * t130 + mrSges(7,3) * t131) - t87 - t88 + t89 - t90 + t91 + t92 - t84 * (mrSges(5,1) * t131 - mrSges(5,2) * t130) + (t351 * t131 + t127 - t289 + t53) * t315 + (-t352 * t130 + t126 - t283 + t330) * t313; t296 * t131 + t266 * t162 + t61 + t64 + (t10 * t162 + t131 * t18 + t1) * m(7) + (t131 * t32 - t162 * t20 + t5) * m(6); -t162 * t103 - t130 * t69 + 0.2e1 * (t2 / 0.2e1 + t18 * t315 + t9 * t311) * m(7) + t63;];
tauc  = t16(:);
