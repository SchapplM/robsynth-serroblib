% Calculate vector of centrifugal and coriolis load on the joints for
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:21:48
% EndTime: 2018-11-23 15:21:56
% DurationCPUTime: 8.42s
% Computational Cost: add. (3648->563), mult. (9387->697), div. (0->0), fcn. (5937->8), ass. (0->258)
t358 = -Ifges(5,4) + Ifges(7,6);
t356 = Ifges(5,1) + Ifges(7,3);
t351 = Ifges(7,4) + Ifges(6,5);
t350 = Ifges(5,5) + Ifges(7,5);
t355 = Ifges(7,2) + Ifges(6,3);
t180 = cos(qJ(3));
t257 = qJD(2) * t180;
t167 = Ifges(4,4) * t257;
t357 = -t167 / 0.2e1;
t176 = sin(qJ(4));
t354 = t358 * t176;
t179 = cos(qJ(4));
t252 = t179 * qJD(3);
t177 = sin(qJ(3));
t259 = qJD(2) * t177;
t129 = t176 * t259 - t252;
t130 = qJD(3) * t176 + t179 * t259;
t178 = sin(qJ(2));
t173 = sin(pkin(6));
t261 = qJD(1) * t173;
t244 = t178 * t261;
t136 = qJD(2) * pkin(8) + t244;
t174 = cos(pkin(6));
t260 = qJD(1) * t174;
t107 = -t177 * t136 + t180 * t260;
t84 = -qJD(3) * pkin(3) - t107;
t191 = -qJ(5) * t130 + t84;
t296 = pkin(4) + qJ(6);
t18 = t129 * t296 + t191;
t279 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t129 - mrSges(5,2) * t130 - mrSges(4,3) * t259;
t196 = -m(4) * t107 + m(5) * t84 - t279;
t70 = -mrSges(7,2) * t130 + mrSges(7,3) * t129;
t73 = -mrSges(6,2) * t129 - mrSges(6,3) * t130;
t295 = t70 + t73;
t32 = pkin(4) * t129 + t191;
t353 = m(6) * t32 + m(7) * t18 + t196 + t295;
t352 = t259 / 0.2e1;
t237 = Ifges(4,5) * qJD(3) / 0.2e1;
t251 = qJD(2) * qJD(3);
t235 = t177 * t251;
t234 = t180 * t251;
t254 = qJD(4) * t176;
t239 = t177 * t254;
t250 = qJD(3) * qJD(4);
t96 = qJD(2) * t239 + (-t234 - t250) * t179;
t253 = qJD(4) * t179;
t255 = qJD(3) * t180;
t338 = t176 * t255 + t177 * t253;
t97 = qJD(2) * t338 + t176 * t250;
t349 = t355 * t97 + (Ifges(6,6) - Ifges(7,6)) * t96 + t351 * t235;
t348 = t350 * t235 - t356 * t96 + t358 * t97;
t275 = qJ(5) * t179;
t198 = qJ(6) * t176 - t275;
t195 = t198 * t180;
t231 = pkin(4) * t254 - qJD(5) * t176;
t108 = t180 * t136 + t177 * t260;
t246 = t176 * pkin(4) * t257 + t108;
t347 = -qJD(2) * t195 + qJD(4) * t198 - qJD(6) * t179 + t231 - t246;
t143 = -pkin(3) * t180 - pkin(9) * t177 - pkin(2);
t194 = t143 * qJD(2);
t181 = cos(qJ(2));
t243 = t181 * t261;
t190 = t194 - t243;
t273 = t173 * t181;
t241 = qJD(2) * t273;
t192 = qJD(1) * (qJD(3) * t174 + t241);
t256 = qJD(3) * t177;
t58 = -t136 * t256 + t180 * t192;
t346 = qJD(4) * t190 + t58;
t313 = pkin(5) + pkin(9);
t147 = t313 * t179;
t268 = t179 * t180;
t248 = pkin(5) * t268;
t225 = pkin(3) * t177 - pkin(9) * t180;
t132 = t225 * qJD(2);
t44 = -t176 * t107 + t132 * t179;
t345 = -(-t177 * t296 + t248) * qJD(2) + t44 + qJD(4) * t147;
t270 = t176 * t180;
t249 = pkin(5) * t270;
t45 = t179 * t107 + t176 * t132;
t344 = -(qJ(5) * t177 - t249) * qJD(2) - t45 - t313 * t254;
t161 = -qJD(4) + t257;
t103 = mrSges(7,1) * t130 + mrSges(7,3) * t161;
t294 = mrSges(5,3) * t130;
t102 = -mrSges(5,1) * t161 - t294;
t106 = mrSges(6,1) * t130 - mrSges(6,2) * t161;
t265 = t102 - t106;
t343 = -t103 + t265;
t286 = qJD(2) * pkin(2);
t137 = -t243 - t286;
t229 = t176 * t243;
t85 = qJD(3) * pkin(9) + t108;
t31 = t176 * t194 + t179 * t85 - t229;
t20 = t161 * qJ(5) - t31;
t335 = -t129 * pkin(5) + qJD(6);
t10 = -t20 + t335;
t30 = t176 * t85 - t179 * t190;
t199 = t31 * t176 - t30 * t179;
t328 = -qJD(5) - t30;
t19 = pkin(4) * t161 - t328;
t200 = t20 * t176 + t19 * t179;
t291 = Ifges(6,6) * t176;
t208 = -Ifges(6,2) * t179 + t291;
t292 = Ifges(5,4) * t179;
t213 = -Ifges(5,2) * t176 + t292;
t217 = -mrSges(7,2) * t179 + mrSges(7,3) * t176;
t218 = -mrSges(6,2) * t176 - mrSges(6,3) * t179;
t220 = mrSges(5,1) * t176 + mrSges(5,2) * t179;
t303 = t179 / 0.2e1;
t304 = -t179 / 0.2e1;
t305 = t176 / 0.2e1;
t306 = -t176 / 0.2e1;
t308 = -t161 / 0.2e1;
t309 = t130 / 0.2e1;
t310 = -t130 / 0.2e1;
t311 = t129 / 0.2e1;
t312 = -t129 / 0.2e1;
t322 = -Ifges(5,6) + t351;
t323 = Ifges(6,4) - t350;
t125 = Ifges(7,6) * t130;
t282 = t130 * Ifges(6,6);
t331 = t129 * t355 - t161 * t351 + t125 - t282;
t127 = Ifges(5,4) * t129;
t289 = Ifges(7,6) * t129;
t332 = t130 * t356 - t161 * t350 - t127 + t289;
t287 = Ifges(7,6) * t179;
t290 = Ifges(6,6) * t179;
t340 = t176 * t355 + t287 - t290;
t341 = t179 * t356 + t354;
t126 = Ifges(6,6) * t129;
t53 = -Ifges(6,4) * t161 - Ifges(6,2) * t130 + t126;
t283 = t130 * Ifges(5,4);
t54 = -Ifges(5,2) * t129 - Ifges(5,6) * t161 + t283;
t197 = pkin(5) * t130 + t30;
t336 = qJD(5) + t197;
t9 = t161 * t296 + t336;
t183 = t200 * mrSges(6,1) - (t10 * t176 - t9 * t179) * mrSges(7,1) - t199 * mrSges(5,3) + t18 * t217 + t208 * t310 + t213 * t312 + t32 * t218 + t84 * t220 + t53 * t304 + t54 * t306 + t340 * t311 + t341 * t309 + t331 * t305 + t332 * t303 + (t176 * t322 - t179 * t323) * t308;
t342 = -t137 * mrSges(4,2) + t107 * mrSges(4,3) - Ifges(4,1) * t352 - t183 - t237 + t357;
t339 = t290 + t292 + (Ifges(5,1) + Ifges(6,2)) * t176;
t334 = t291 - t354 + (Ifges(5,2) + t355) * t179;
t236 = -Ifges(4,6) * qJD(3) / 0.2e1;
t333 = pkin(8) * (t177 * t252 + t180 * t254);
t62 = mrSges(6,1) * t97 - mrSges(6,3) * t235;
t66 = -mrSges(5,2) * t235 - mrSges(5,3) * t97;
t330 = -t62 + t66;
t64 = -t96 * mrSges(6,1) + mrSges(6,2) * t235;
t65 = mrSges(5,1) * t235 + mrSges(5,3) * t96;
t329 = t64 - t65;
t284 = t129 * mrSges(5,3);
t101 = mrSges(5,2) * t161 - t284;
t104 = mrSges(6,1) * t129 + mrSges(6,3) * t161;
t105 = -mrSges(7,1) * t129 - mrSges(7,2) * t161;
t264 = -t104 + t105;
t327 = -t101 - t264;
t326 = -t31 - t335;
t135 = t225 * qJD(3);
t100 = (t135 + t244) * qJD(2);
t6 = t176 * t100 + t179 * t346 - t254 * t85;
t7 = t100 * t179 - t176 * t346 - t85 * t253;
t325 = -t176 * t7 + t179 * t6;
t4 = -qJ(5) * t235 + qJD(5) * t161 - t6;
t5 = -pkin(4) * t235 - t7;
t324 = t176 * t5 - t179 * t4;
t226 = -Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1;
t227 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1;
t228 = Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,1) / 0.2e1;
t318 = -t227 * t129 + t226 * t130 + t228 * t161 - t10 * mrSges(7,2) - t137 * mrSges(4,1) - t19 * mrSges(6,2) - t236 + (Ifges(4,4) * t177 + Ifges(4,2) * t180) * qJD(2) / 0.2e1 - Ifges(5,6) * t312 - Ifges(6,4) * t310 + t108 * mrSges(4,3) + t20 * mrSges(6,3) + t30 * mrSges(5,1) + t31 * mrSges(5,2) + t9 * mrSges(7,3) - t351 * t311 - t350 * t309 - (Ifges(5,3) + Ifges(7,1) + Ifges(6,1)) * t308;
t317 = -t96 / 0.2e1;
t316 = t96 / 0.2e1;
t315 = -t97 / 0.2e1;
t307 = t161 / 0.2e1;
t302 = pkin(8) * t176;
t274 = t173 * t178;
t113 = -t174 * t180 + t177 * t274;
t59 = t136 * t255 + t177 * t192;
t280 = t59 * t113;
t276 = qJ(5) * t129;
t271 = t176 * t177;
t269 = t177 * t179;
t267 = t180 * t181;
t266 = -t101 + t104;
t263 = t176 * t135 + t143 * t253;
t262 = pkin(4) * t271 + t177 * pkin(8);
t164 = pkin(8) * t268;
t110 = t176 * t143 + t164;
t258 = qJD(2) * t178;
t163 = pkin(8) * t270;
t247 = t176 * t273;
t245 = -pkin(4) - t302;
t242 = t173 * t258;
t145 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t257;
t233 = -m(4) * t108 - t145;
t63 = -t97 * mrSges(7,1) + mrSges(7,2) * t235;
t232 = -qJ(5) * t176 - pkin(3);
t109 = t143 * t179 - t163;
t230 = pkin(4) * t338 + pkin(8) * t255 + qJ(5) * t239;
t94 = qJ(5) * t180 - t110;
t224 = qJD(4) * t164 - t135 * t179 + t143 * t254;
t222 = -qJD(5) * t180 + t263;
t221 = mrSges(5,1) * t179 - mrSges(5,2) * t176;
t219 = mrSges(6,2) * t179 - mrSges(6,3) * t176;
t216 = mrSges(7,2) * t176 + mrSges(7,3) * t179;
t211 = Ifges(6,4) * t176 + Ifges(6,5) * t179;
t210 = Ifges(7,4) * t179 - Ifges(7,5) * t176;
t209 = Ifges(5,5) * t176 + Ifges(5,6) * t179;
t202 = -Ifges(7,3) * t176 + t287;
t61 = -t96 * mrSges(7,1) - mrSges(7,3) * t235;
t114 = t174 * t177 + t180 * t274;
t75 = t114 * t176 + t179 * t273;
t1 = -pkin(5) * t96 + qJD(6) * t161 - t235 * t296 - t7;
t2 = -pkin(5) * t97 - t4;
t189 = -t7 * mrSges(5,1) + t6 * mrSges(5,2) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3);
t187 = -qJD(3) * t113 + t180 * t241;
t186 = qJ(5) * t96 - qJD(5) * t130 + t59;
t172 = t180 * pkin(4);
t160 = Ifges(6,1) * t235;
t159 = Ifges(7,1) * t235;
t158 = Ifges(5,3) * t235;
t146 = t313 * t176;
t138 = -pkin(4) * t179 + t232;
t131 = (-mrSges(4,1) * t180 + mrSges(4,2) * t177) * qJD(2);
t120 = -t179 * t296 + t232;
t119 = (mrSges(4,1) * t177 + mrSges(4,2) * t180) * t251;
t112 = -qJ(5) * t253 + t231;
t111 = -qJ(5) * t269 + t262;
t99 = (t176 * t178 + t179 * t267) * t261;
t98 = -t179 * t244 + t180 * t229;
t95 = -t109 + t172;
t92 = Ifges(6,4) * t96;
t91 = Ifges(7,4) * t97;
t90 = Ifges(5,5) * t96;
t89 = Ifges(6,5) * t97;
t88 = Ifges(7,5) * t96;
t87 = Ifges(5,6) * t97;
t83 = t177 * t198 + t262;
t71 = pkin(4) * t130 + t276;
t68 = -pkin(5) * t271 - t94;
t60 = qJ(6) * t180 + t163 + t172 + (pkin(5) * t177 - t143) * t179;
t48 = -t257 * t275 + t246;
t43 = t130 * t296 + t276;
t42 = t256 * t302 - t224;
t41 = t263 - t333;
t40 = (-qJ(5) * t255 - qJD(5) * t177) * t179 + t230;
t39 = -pkin(4) * t259 - t44;
t38 = -qJ(5) * t259 - t45;
t36 = t245 * t256 + t224;
t35 = mrSges(7,2) * t96 + mrSges(7,3) * t97;
t34 = mrSges(5,1) * t97 - mrSges(5,2) * t96;
t33 = -mrSges(6,2) * t97 + mrSges(6,3) * t96;
t28 = -qJ(5) * t256 - t222 + t333;
t26 = -t96 * Ifges(5,4) - t97 * Ifges(5,2) + Ifges(5,6) * t235;
t25 = Ifges(6,4) * t235 + t96 * Ifges(6,2) + t97 * Ifges(6,6);
t15 = qJD(3) * t195 + (qJD(6) * t176 + (qJ(6) * qJD(4) - qJD(5)) * t179) * t177 + t230;
t14 = (-pkin(5) * t269 - t163) * qJD(4) + (-t249 + (-pkin(8) * t179 + qJ(5)) * t177) * qJD(3) + t222;
t11 = -pkin(5) * t239 + qJD(6) * t180 + (t248 + (-qJ(6) + t245) * t177) * qJD(3) + t224;
t8 = pkin(4) * t97 + t186;
t3 = qJD(6) * t129 + t296 * t97 + t186;
t12 = [t187 * t145 + t131 * t242 + m(5) * t280 + m(4) * (t108 * t174 * t255 + t280 + t58 * t114 + (-t108 * t178 * t256 + (t108 * t267 + (t137 - t243) * t178) * qJD(2)) * t173) - t119 * t273 - t114 * mrSges(4,3) * t235 + (-mrSges(3,1) * t178 - mrSges(3,2) * t181) * t173 * qJD(2) ^ 2 + (m(5) * t6 - m(6) * t4 + m(7) * t2 + t330 + t63) * (t114 * t179 - t247) + (-m(5) * t7 + m(6) * t5 + m(7) * t1 + t329 + t61) * t75 + (m(5) * t30 + m(6) * t19 + m(7) * t9 - t343) * (-qJD(4) * t247 + t114 * t253 + t176 * t187 - t179 * t242) + (-m(5) * t31 + m(6) * t20 - m(7) * t10 + t327) * (qJD(4) * t75 - t176 * t242 - t179 * t187) + t353 * (qJD(3) * t114 + t177 * t241) + (m(6) * t8 + m(7) * t3 + mrSges(4,3) * t234 + t33 + t34 + t35) * t113; (-t158 / 0.2e1 - t159 / 0.2e1 - t160 / 0.2e1 + t189 - t227 * t97 - t226 * t96 + (m(4) * pkin(8) + mrSges(4,3)) * t58 + (-mrSges(4,1) * t258 + t181 * t233) * t261 + t87 / 0.2e1 + t88 / 0.2e1 - t89 / 0.2e1 + t90 / 0.2e1 - t91 / 0.2e1 - t92 / 0.2e1 + (t196 * pkin(8) + t237 + 0.3e1 / 0.2e1 * t167 - t342) * qJD(3)) * t180 + t343 * t98 + (t213 * t315 + t25 * t304 + t8 * t218 + t3 * t217 + t208 * t316 + t26 * t306 + (-t176 * t6 - t179 * t7) * mrSges(5,3) + (t1 * t179 - t176 * t2) * mrSges(7,1) + (t176 * t4 + t179 * t5) * mrSges(6,1) + (mrSges(4,3) + t220) * t59 + (mrSges(4,2) * t258 - t181 * t353) * t261 + (t236 + ((-0.3e1 / 0.2e1 * Ifges(4,4) - t226 * t179 + t227 * t176) * t177 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t228) * t180) * qJD(2) - t318) * qJD(3) + (t54 * t304 + t202 * t309 + t209 * t307 + t84 * t221 - t32 * t219 + t18 * t216 + (-t176 * t30 - t179 * t31) * mrSges(5,3) + (-t10 * t179 - t176 * t9) * mrSges(7,1) + (-t176 * t19 + t179 * t20) * mrSges(6,1) + t339 * t310 + (t211 + t210) * t308 + t332 * t306 + t334 * t311) * qJD(4) + t341 * t317 + t340 * t97 / 0.2e1 + (qJD(4) * t53 + t349) * t305 + (qJD(4) * t331 + t348) * t303 + (t34 + (m(4) + m(5)) * t59 + t233 * qJD(3)) * pkin(8)) * t177 + t83 * t35 + t94 * t62 + t95 * t64 + t40 * t73 + t60 * t61 + t68 * t63 + t15 * t70 + (-t131 + 0.2e1 * (-t286 / 0.2e1 - t137 / 0.2e1) * m(4)) * t244 - m(5) * (t30 * t98 + t31 * t99) + (t1 * t60 + t15 * t18 + t2 * t68 + t3 * t83 + (-t98 + t11) * t9 + (t14 - t99) * t10) * m(7) + m(5) * (t109 * t7 + t110 * t6 - t30 * t42 + t31 * t41) + (t111 * t8 + t32 * t40 + t4 * t94 + t5 * t95 + (t28 + t99) * t20 + (t36 - t98) * t19) * m(6) + t327 * t99 + t41 * t101 + t42 * t102 + t11 * t103 + t28 * t104 + t14 * t105 + t36 * t106 + t109 * t65 + t110 * t66 + t111 * t33 - pkin(2) * t119; (t1 * t146 + t10 * t344 + t120 * t3 + t147 * t2 + t18 * t347 + t345 * t9) * m(7) + t347 * t70 + t348 * t305 + t349 * t304 + t345 * t103 + t344 * t105 + t339 * t317 + t334 * t315 - t58 * mrSges(4,2) - pkin(3) * t34 - m(5) * (t108 * t84 - t30 * t44 + t31 * t45) - m(6) * (t19 * t39 + t20 * t38 + t32 * t48) + (t1 * t176 + t179 * t2) * mrSges(7,1) + (-t48 + t112) * t73 + ((Ifges(4,4) * t352 + t236 + t318) * t177 + (t357 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t259 + t237 + t342) * t180 + (-t210 / 0.2e1 - t211 / 0.2e1 + t209 / 0.2e1) * t256) * qJD(2) + t202 * t316 + t25 * t306 + t26 * t303 + t279 * t108 + ((-m(5) * t199 + m(6) * t200 + t176 * t266 - t179 * t265) * pkin(9) + t183) * qJD(4) + (t176 * t329 + t179 * t330) * pkin(9) + m(6) * (pkin(9) * t324 + t112 * t32 + t138 * t8) + t324 * mrSges(6,1) + m(5) * (-pkin(3) * t59 + pkin(9) * t325) + t325 * mrSges(5,3) - t45 * t101 - t44 * t102 - t38 * t104 - t39 * t106 + t120 * t35 + t138 * t33 - t107 * t145 + t146 * t61 + t147 * t63 - t3 * t216 + t8 * t219 + (-mrSges(4,1) - t221) * t59; (qJ(5) * t2 - t1 * t296 + t10 * t336 - t18 * t43 + t326 * t9) * m(7) + t158 + t159 + t160 + (t130 * t355 + t126 - t289 + t53) * t312 - t71 * t73 - pkin(4) * t64 - t43 * t70 + (-t129 * t356 + t125 - t283 + t331) * t310 - t189 + (Ifges(6,2) * t129 + t282 + t54) * t309 + (t10 * t130 + t129 * t9) * mrSges(7,1) + (t129 * t19 - t130 * t20) * mrSges(6,1) + (-t62 + t63) * qJ(5) + t197 * t105 - t296 * t61 + (t265 + t294) * t31 + (-t266 + t284) * t30 + t264 * qJD(5) + (-Ifges(5,2) * t130 - t127 + t332) * t311 + (-pkin(4) * t5 - qJ(5) * t4 - t19 * t31 + t20 * t328 - t32 * t71) * m(6) + t326 * t103 + (t129 * t323 + t130 * t322) * t307 - t87 - t88 + t89 - t90 + t91 + t92 - t32 * (-mrSges(6,2) * t130 + mrSges(6,3) * t129) - t18 * (mrSges(7,2) * t129 + mrSges(7,3) * t130) - t84 * (mrSges(5,1) * t130 - mrSges(5,2) * t129); t295 * t130 + t264 * t161 + t61 + t64 + (t10 * t161 + t130 * t18 + t1) * m(7) + (t130 * t32 - t161 * t20 + t5) * m(6); -t161 * t103 - t129 * t70 + 0.2e1 * (t2 / 0.2e1 + t18 * t312 + t9 * t308) * m(7) + t63;];
tauc  = t12(:);
