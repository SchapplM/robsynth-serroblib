% Calculate vector of inverse dynamics joint torques for
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:08
% EndTime: 2019-12-05 16:30:42
% DurationCPUTime: 13.03s
% Computational Cost: add. (3740->520), mult. (8907->747), div. (0->0), fcn. (6660->14), ass. (0->253)
t192 = sin(qJ(3));
t195 = cos(qJ(3));
t226 = mrSges(4,1) * t195 - mrSges(4,2) * t192;
t290 = pkin(8) + qJ(4);
t326 = -m(5) * qJ(4) - m(6) * t290 - mrSges(5,3) - mrSges(6,3);
t188 = cos(pkin(10));
t180 = pkin(4) * t188 + pkin(3);
t185 = pkin(10) + qJ(5);
t182 = sin(t185);
t183 = cos(t185);
t186 = sin(pkin(10));
t225 = -mrSges(5,1) * t188 + mrSges(5,2) * t186;
t351 = m(5) * pkin(3) + m(6) * t180 + mrSges(6,1) * t183 - mrSges(6,2) * t182 - t225;
t360 = t192 * t326 - t351 * t195 - mrSges(3,1) - t226;
t217 = pkin(3) * t192 - qJ(4) * t195;
t261 = qJD(4) * t192;
t127 = qJD(3) * t217 - t261;
t193 = sin(qJ(2));
t263 = qJD(3) * t192;
t256 = pkin(7) * t263;
t187 = sin(pkin(5));
t269 = qJD(1) * t187;
t196 = cos(qJ(2));
t271 = t195 * t196;
t335 = t188 * t127 + t186 * t256 - (-t186 * t271 + t188 * t193) * t269;
t359 = (t186 * t193 + t188 * t271) * t269 - t186 * t127;
t264 = qJD(2) * t195;
t178 = qJD(5) - t264;
t297 = t178 / 0.2e1;
t266 = qJD(2) * t192;
t147 = qJD(3) * t188 - t186 * t266;
t148 = qJD(3) * t186 + t188 * t266;
t191 = sin(qJ(5));
t194 = cos(qJ(5));
t75 = t147 * t191 + t148 * t194;
t303 = t75 / 0.2e1;
t228 = t194 * t147 - t148 * t191;
t305 = t228 / 0.2e1;
t357 = Ifges(6,5) * t303 + Ifges(6,6) * t305 + Ifges(6,3) * t297;
t239 = -t264 / 0.2e1;
t272 = t188 * t195;
t215 = pkin(4) * t192 - pkin(8) * t272;
t356 = qJD(3) * t215 + t335;
t273 = t188 * t192;
t276 = t186 * t195;
t355 = -(-pkin(7) * t273 - pkin(8) * t276) * qJD(3) + t359;
t151 = t186 * t194 + t188 * t191;
t211 = t151 * t195;
t111 = qJD(2) * t211;
t134 = t151 * qJD(5);
t353 = t111 - t134;
t216 = t186 * t191 - t188 * t194;
t210 = t216 * t195;
t112 = qJD(2) * t210;
t133 = t216 * qJD(5);
t352 = t112 - t133;
t295 = Ifges(6,4) * t75;
t27 = Ifges(6,2) * t228 + Ifges(6,6) * t178 + t295;
t309 = t27 / 0.2e1;
t71 = Ifges(6,4) * t228;
t28 = Ifges(6,1) * t75 + Ifges(6,5) * t178 + t71;
t308 = t28 / 0.2e1;
t267 = qJD(2) * t187;
t242 = qJD(1) * t267;
t169 = t196 * t242;
t259 = qJDD(1) * t187;
t124 = t193 * t259 + t169;
t189 = cos(pkin(5));
t268 = qJD(1) * t189;
t350 = qJDD(2) * pkin(7) + qJD(3) * t268 + t124;
t323 = m(6) + m(5) + m(4);
t349 = pkin(2) * t323 - t360;
t159 = -pkin(3) * t195 - qJ(4) * t192 - pkin(2);
t249 = t196 * t269;
t106 = qJD(2) * t159 - t249;
t250 = t193 * t269;
t157 = qJD(2) * pkin(7) + t250;
t248 = t192 * t268;
t104 = t157 * t195 + t248;
t98 = qJD(3) * qJ(4) + t104;
t46 = t188 * t106 - t186 * t98;
t29 = -pkin(4) * t264 - pkin(8) * t148 + t46;
t47 = t186 * t106 + t188 * t98;
t31 = pkin(8) * t147 + t47;
t11 = -t191 * t31 + t194 * t29;
t12 = t191 * t29 + t194 * t31;
t289 = Ifges(4,4) * t192;
t221 = t195 * Ifges(4,2) + t289;
t348 = t11 * mrSges(6,1) - Ifges(4,6) * qJD(3) / 0.2e1 - qJD(2) * t221 / 0.2e1 - t12 * mrSges(6,2) + t148 * Ifges(5,5) / 0.2e1 + t147 * Ifges(5,6) / 0.2e1 + Ifges(5,3) * t239 + t357;
t224 = t186 * mrSges(5,1) + t188 * mrSges(5,2);
t347 = -t182 * mrSges(6,1) - t183 * mrSges(6,2) + mrSges(3,2) - mrSges(4,3) - t224;
t260 = qJD(2) * qJD(3);
t155 = -t195 * qJDD(2) + t192 * t260;
t299 = t155 / 0.2e1;
t156 = qJDD(2) * t192 + t195 * t260;
t114 = qJDD(3) * t186 + t156 * t188;
t301 = t114 / 0.2e1;
t346 = Ifges(5,1) * t301 + Ifges(5,5) * t299;
t345 = m(4) * pkin(7);
t113 = qJDD(3) * t188 - t156 * t186;
t23 = qJD(5) * t228 + t113 * t191 + t114 * t194;
t311 = t23 / 0.2e1;
t24 = -qJD(5) * t75 + t113 * t194 - t114 * t191;
t310 = t24 / 0.2e1;
t302 = t113 / 0.2e1;
t149 = qJDD(5) + t155;
t300 = t149 / 0.2e1;
t344 = -t155 / 0.2e1;
t343 = t156 / 0.2e1;
t146 = t188 * t159;
t76 = -pkin(8) * t273 + t146 + (-pkin(7) * t186 - pkin(4)) * t195;
t108 = pkin(7) * t272 + t186 * t159;
t277 = t186 * t192;
t90 = -pkin(8) * t277 + t108;
t32 = -t191 * t90 + t194 * t76;
t342 = qJD(5) * t32 + t356 * t191 - t355 * t194;
t341 = qJD(3) / 0.2e1;
t33 = t191 * t76 + t194 * t90;
t339 = -qJD(5) * t33 + t355 * t191 + t356 * t194;
t144 = t192 * t157;
t103 = t195 * t268 - t144;
t153 = t217 * qJD(2);
t60 = -t103 * t186 + t188 * t153;
t43 = qJD(2) * t215 + t60;
t247 = t186 * t264;
t61 = t188 * t103 + t186 * t153;
t51 = -pkin(8) * t247 + t61;
t160 = t290 * t186;
t161 = t290 * t188;
t91 = -t160 * t194 - t161 * t191;
t338 = -qJD(4) * t216 + qJD(5) * t91 - t191 * t43 - t194 * t51;
t92 = -t160 * t191 + t161 * t194;
t337 = -qJD(4) * t151 - qJD(5) * t92 + t191 * t51 - t194 * t43;
t52 = -t113 * mrSges(5,1) + t114 * mrSges(5,2);
t333 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t156 + t52;
t254 = mrSges(4,3) * t266;
t332 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t147 + mrSges(5,2) * t148 + t254;
t331 = -t188 * t256 - t359;
t286 = Ifges(5,4) * t188;
t220 = -Ifges(5,2) * t186 + t286;
t287 = Ifges(5,4) * t186;
t222 = Ifges(5,1) * t188 - t287;
t328 = t147 * (Ifges(5,6) * t192 + t195 * t220) + t148 * (Ifges(5,5) * t192 + t195 * t222);
t258 = qJDD(1) * t189;
t252 = t192 * t258 + t195 * t350;
t41 = -t157 * t263 + t252;
t262 = qJD(3) * t195;
t42 = -t157 * t262 - t192 * t350 + t195 * t258;
t325 = -t192 * t42 + t195 * t41;
t37 = qJDD(3) * qJ(4) + (qJD(4) - t144) * qJD(3) + t252;
t168 = t193 * t242;
t123 = t196 * t259 - t168;
t115 = -qJDD(2) * pkin(2) - t123;
t50 = pkin(3) * t155 - qJ(4) * t156 - qJD(2) * t261 + t115;
t14 = -t186 * t37 + t188 * t50;
t15 = t186 * t50 + t188 * t37;
t324 = -t14 * t186 + t15 * t188;
t322 = mrSges(4,1) + t351;
t321 = mrSges(4,2) + t326;
t158 = -qJD(2) * pkin(2) - t249;
t94 = -qJD(3) * pkin(3) + qJD(4) - t103;
t318 = -t195 * t94 * t224 - t47 * (-mrSges(5,2) * t192 - mrSges(5,3) * t276) - t46 * (mrSges(5,1) * t192 - mrSges(5,3) * t272) - t158 * (mrSges(4,1) * t192 + mrSges(4,2) * t195);
t317 = -m(5) * t94 - t332;
t13 = pkin(8) * t113 + t15;
t8 = pkin(4) * t155 - pkin(8) * t114 + t14;
t1 = qJD(5) * t11 + t13 * t194 + t191 * t8;
t2 = -qJD(5) * t12 - t13 * t191 + t194 * t8;
t316 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t294 = pkin(4) * t186;
t251 = pkin(7) + t294;
t315 = -m(5) * pkin(7) - m(6) * t251 - t345 + t347;
t197 = qJD(2) ^ 2;
t313 = Ifges(6,4) * t311 + Ifges(6,2) * t310 + Ifges(6,6) * t300;
t312 = Ifges(6,1) * t311 + Ifges(6,4) * t310 + Ifges(6,5) * t300;
t307 = Ifges(5,4) * t302 + t346;
t306 = -t228 / 0.2e1;
t304 = -t75 / 0.2e1;
t298 = -t178 / 0.2e1;
t288 = Ifges(4,4) * t195;
t38 = -qJDD(3) * pkin(3) + qJDD(4) - t42;
t283 = t192 * t38;
t280 = cos(pkin(9));
t279 = sin(pkin(9));
t275 = t187 * t193;
t274 = t187 * t196;
t265 = qJD(2) * t193;
t257 = Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t149;
t253 = mrSges(4,3) * t264;
t246 = t187 * t265;
t245 = t196 * t267;
t5 = -t24 * mrSges(6,1) + t23 * mrSges(6,2);
t236 = t187 * t280;
t235 = t187 * t279;
t234 = t280 * t193;
t233 = t280 * t196;
t232 = t279 * t193;
t231 = t279 * t196;
t230 = -t260 / 0.2e1;
t229 = t260 / 0.2e1;
t219 = Ifges(4,5) * t195 - Ifges(4,6) * t192;
t218 = Ifges(5,5) * t188 - Ifges(5,6) * t186;
t137 = t189 * t192 + t195 * t275;
t82 = -t137 * t186 - t188 * t274;
t83 = t137 * t188 - t186 * t274;
t35 = -t191 * t83 + t194 * t82;
t36 = t191 * t82 + t194 * t83;
t136 = -t189 * t195 + t192 * t275;
t213 = t192 * (Ifges(4,1) * t195 - t289);
t130 = t189 * t234 + t231;
t84 = t130 * t192 + t195 * t236;
t132 = -t189 * t232 + t233;
t86 = t132 * t192 - t195 * t235;
t212 = -g(1) * t86 - g(2) * t84 - g(3) * t136;
t200 = t195 * (Ifges(5,3) * t192 + t195 * t218);
t181 = Ifges(4,4) * t264;
t163 = -qJD(3) * mrSges(4,2) + t253;
t154 = t251 * t192;
t152 = t226 * qJD(2);
t143 = t251 * t262;
t141 = Ifges(4,1) * t266 + Ifges(4,5) * qJD(3) + t181;
t131 = t189 * t231 + t234;
t129 = -t189 * t233 + t232;
t125 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t155;
t120 = t216 * t192;
t119 = t151 * t192;
t110 = -mrSges(5,1) * t264 - mrSges(5,3) * t148;
t109 = mrSges(5,2) * t264 + mrSges(5,3) * t147;
t107 = -pkin(7) * t276 + t146;
t93 = mrSges(4,1) * t155 + mrSges(4,2) * t156;
t89 = -qJD(3) * t136 + t195 * t245;
t88 = qJD(3) * t137 + t192 * t245;
t87 = t132 * t195 + t192 * t235;
t85 = t130 * t195 - t192 * t236;
t77 = t248 + (qJD(2) * t294 + t157) * t195;
t70 = t148 * Ifges(5,1) + t147 * Ifges(5,4) - Ifges(5,5) * t264;
t69 = t148 * Ifges(5,4) + t147 * Ifges(5,2) - Ifges(5,6) * t264;
t67 = mrSges(5,1) * t155 - mrSges(5,3) * t114;
t66 = -mrSges(5,2) * t155 + mrSges(5,3) * t113;
t64 = -qJD(3) * t211 + t133 * t192;
t63 = -qJD(3) * t210 - t134 * t192;
t62 = -pkin(4) * t147 + t94;
t59 = t186 * t246 + t188 * t89;
t58 = -t186 * t89 + t188 * t246;
t56 = mrSges(6,1) * t178 - mrSges(6,3) * t75;
t55 = -mrSges(6,2) * t178 + mrSges(6,3) * t228;
t39 = t114 * Ifges(5,4) + t113 * Ifges(5,2) + t155 * Ifges(5,6);
t30 = -mrSges(6,1) * t228 + mrSges(6,2) * t75;
t25 = -pkin(4) * t113 + t38;
t19 = -mrSges(6,2) * t149 + mrSges(6,3) * t24;
t18 = mrSges(6,1) * t149 - mrSges(6,3) * t23;
t7 = -qJD(5) * t36 - t191 * t59 + t194 * t58;
t6 = qJD(5) * t35 + t191 * t58 + t194 * t59;
t3 = [m(2) * qJDD(1) + t59 * t109 + t58 * t110 + t137 * t125 + t89 * t163 + t35 * t18 + t36 * t19 + t6 * t55 + t7 * t56 + t83 * t66 + t82 * t67 + (t30 + t332) * t88 + (t5 + t333) * t136 + (-m(2) - m(3) - t323) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t197 - t93) * t196 + (-mrSges(3,1) * t197 - mrSges(3,2) * qJDD(2) - qJD(2) * t152) * t193) * t187 + m(6) * (t1 * t36 + t11 * t7 + t12 * t6 + t136 * t25 + t2 * t35 + t62 * t88) + m(5) * (t136 * t38 + t14 * t82 + t15 * t83 + t46 * t58 + t47 * t59 + t88 * t94) + m(4) * (-t103 * t88 + t104 * t89 - t136 * t42 + t137 * t41 + (-t115 * t196 + t158 * t265) * t187) + m(3) * (qJDD(1) * t189 ^ 2 + (t123 * t196 + t124 * t193) * t187); t273 * t307 + t63 * t308 + t64 * t309 - t120 * t312 - t119 * t313 + t224 * t283 + t62 * (-mrSges(6,1) * t64 + mrSges(6,2) * t63) + t32 * t18 + t33 * t19 - t39 * t277 / 0.2e1 + (Ifges(6,5) * t63 + Ifges(6,6) * t64) * t297 + (t188 * t70 + t141) * t262 / 0.2e1 - (Ifges(5,5) * t114 + Ifges(5,6) * t113 + Ifges(5,3) * t155 + t257) * t195 / 0.2e1 + (-Ifges(6,5) * t120 - Ifges(6,6) * t119) * t300 + (-t1 * t119 - t11 * t63 + t12 * t64 + t120 * t2) * mrSges(6,3) + (-Ifges(6,1) * t120 - Ifges(6,4) * t119) * t311 + t25 * (mrSges(6,1) * t119 - mrSges(6,2) * t120) + (-Ifges(6,4) * t120 - Ifges(6,2) * t119) * t310 + t213 * t229 + t200 * t230 + (-t14 * t273 - t15 * t277) * mrSges(5,3) - t163 * t256 + (t123 + t168) * mrSges(3,1) + (t129 * t349 + t130 * t315) * g(2) + (t131 * t349 + t132 * t315) * g(1) + (-t103 * t262 + t325) * mrSges(4,3) + (-t124 + t169) * mrSges(3,2) + (-m(6) * t62 - t30 + t317) * t192 * t249 + t152 * t250 + t332 * pkin(7) * t262 + t333 * pkin(7) * t192 + (-t104 * mrSges(4,3) + t348 + t357) * t263 + (-pkin(2) * t115 + t325 * pkin(7) - (t158 * t193 + (-t103 * t192 + t104 * t195) * t196) * t269) * m(4) - pkin(2) * t93 + t107 * t67 + t108 * t66 + (Ifges(6,1) * t63 + Ifges(6,4) * t64) * t303 + (-t323 * (pkin(2) * t274 + pkin(7) * t275) + (t360 * t196 + (-m(6) * t294 + t347) * t193) * t187) * g(3) + t331 * t109 - t186 * t69 * t262 / 0.2e1 + t335 * t110 + (t107 * t14 + t108 * t15 + (t262 * t94 + t283) * pkin(7) + t331 * t47 + t335 * t46) * m(5) + t339 * t56 + t143 * t30 + (t1 * t33 + t11 * t339 + t12 * t342 + t143 * t62 + t154 * t25 + t2 * t32) * m(6) + t342 * t55 + (-Ifges(6,6) * t310 - Ifges(6,5) * t311 - Ifges(5,3) * t299 - Ifges(6,3) * t300 - Ifges(5,5) * t301 - Ifges(5,6) * t302 + t15 * mrSges(5,2) - t14 * mrSges(5,1) + pkin(7) * t125 - t163 * t249 + (-Ifges(4,2) * t192 + t288) * t229 + Ifges(4,4) * t343 + Ifges(4,2) * t344 - t316) * t195 + (Ifges(4,1) * t156 + Ifges(4,4) * t344 + t218 * t299 + t220 * t302 + t222 * t301) * t192 + t154 * t5 + ((-t103 * t195 - t104 * t192) * t345 + t219 * t341 - t318) * qJD(3) + qJDD(3) * (Ifges(4,5) * t192 + Ifges(4,6) * t195) - t115 * t226 + Ifges(3,3) * qJDD(2) + (Ifges(6,4) * t63 + Ifges(6,2) * t64) * t305 + t328 * t341 + t288 * t343 + t221 * t344; t151 * t312 - t77 * t30 - pkin(3) * t52 - t41 * mrSges(4,2) + t42 * mrSges(4,1) + (-t1 * t216 - t11 * t352 + t12 * t353 - t151 * t2) * mrSges(6,3) + (-mrSges(6,1) * t353 + mrSges(6,2) * t352) * t62 + (-Ifges(6,4) * t133 - Ifges(6,2) * t134) * t305 + (-Ifges(6,5) * t133 - Ifges(6,6) * t134) * t297 + (-Ifges(6,1) * t133 - Ifges(6,4) * t134) * t303 + (t200 / 0.2e1 - t213 / 0.2e1) * t197 + (-t328 / 0.2e1 + t318) * qJD(2) + (Ifges(5,6) * t299 + Ifges(5,2) * t302 + t70 * t239 + qJ(4) * t66 + t39 / 0.2e1 + (m(5) * t47 + t109) * qJD(4)) * t188 + t219 * t230 + (Ifges(6,5) * t304 - Ifges(4,2) * t239 + Ifges(6,6) * t306 + Ifges(6,3) * t298 - t348) * t266 + (t141 + t181) * t239 + t286 * t301 + (t321 * t85 + t322 * t84) * g(2) + (t321 * t87 + t322 * t86) * g(1) + (t136 * t322 + t137 * t321) * g(3) + (Ifges(6,4) * t151 - Ifges(6,2) * t216) * t310 + (Ifges(6,1) * t151 - Ifges(6,4) * t216) * t311 + (Ifges(6,5) * t151 - Ifges(6,6) * t216) * t300 + t25 * (mrSges(6,1) * t216 + mrSges(6,2) * t151) - t216 * t313 + (-Ifges(6,1) * t112 - Ifges(6,4) * t111) * t304 + (-Ifges(6,4) * t112 - Ifges(6,2) * t111) * t306 + (-Ifges(6,5) * t112 - Ifges(6,6) * t111) * t298 + (t254 + t317) * t104 + t287 * t302 + (-qJ(4) * t67 + t307 + (-m(5) * t46 - t110) * qJD(4) + t346) * t186 + (-pkin(3) * t38 + qJ(4) * t324 - t46 * t60 - t47 * t61) * m(5) + t324 * mrSges(5,3) + t91 * t18 + t92 * t19 + t352 * t308 + t353 * t309 - t61 * t109 - t60 * t110 + t337 * t56 + t338 * t55 + (t1 * t92 + t11 * t337 + t12 * t338 - t180 * t25 + t2 * t91 - t62 * t77) * m(6) - Ifges(4,6) * t155 + Ifges(4,5) * t156 - t180 * t5 + t38 * t225 + (t253 - t163) * t103 + Ifges(4,3) * qJDD(3) + t69 * t247 / 0.2e1; -t147 * t109 + t148 * t110 - t228 * t55 + t75 * t56 + t5 + t52 + (t11 * t75 - t12 * t228 + t212 + t25) * m(6) + (-t147 * t47 + t148 * t46 + t212 + t38) * m(5); -t62 * (mrSges(6,1) * t75 + mrSges(6,2) * t228) + (Ifges(6,1) * t228 - t295) * t304 + t27 * t303 + (Ifges(6,5) * t228 - Ifges(6,6) * t75) * t298 + t12 * t56 - t11 * t55 - g(1) * ((t131 * t183 - t182 * t87) * mrSges(6,1) + (-t131 * t182 - t183 * t87) * mrSges(6,2)) - g(2) * ((t129 * t183 - t182 * t85) * mrSges(6,1) + (-t129 * t182 - t183 * t85) * mrSges(6,2)) - g(3) * ((-t137 * t182 - t183 * t274) * mrSges(6,1) + (-t137 * t183 + t182 * t274) * mrSges(6,2)) + (t11 * t228 + t12 * t75) * mrSges(6,3) + t257 + (-Ifges(6,2) * t75 + t28 + t71) * t306 + t316;];
tau = t3;
