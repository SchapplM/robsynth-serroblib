% Calculate vector of inverse dynamics joint torques for
% S5RRRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:40
% EndTime: 2019-12-31 21:50:51
% DurationCPUTime: 5.12s
% Computational Cost: add. (4750->410), mult. (7129->531), div. (0->0), fcn. (4236->12), ass. (0->202)
t376 = mrSges(5,2) - mrSges(6,3);
t375 = mrSges(5,1) + mrSges(6,1);
t366 = Ifges(5,1) + Ifges(6,1);
t364 = Ifges(6,4) + Ifges(5,5);
t215 = qJD(3) + qJD(4);
t221 = sin(qJ(4));
t222 = sin(qJ(3));
t225 = cos(qJ(4));
t226 = cos(qJ(3));
t150 = t221 * t226 + t222 * t225;
t216 = qJD(1) + qJD(2);
t132 = t150 * t216;
t308 = t132 * mrSges(5,3);
t317 = mrSges(6,2) * t132;
t347 = t375 * t215 - t308 - t317;
t223 = sin(qJ(2));
t311 = pkin(1) * qJD(1);
t280 = t223 * t311;
t156 = pkin(7) * t216 + t280;
t265 = pkin(8) * t216 + t156;
t118 = t265 * t222;
t111 = qJD(3) * pkin(3) - t118;
t119 = t265 * t226;
t302 = t119 * t221;
t70 = t111 * t225 - t302;
t287 = -t70 + qJD(5);
t55 = -pkin(4) * t215 + t287;
t374 = -m(6) * t55 + t347;
t365 = -Ifges(5,4) + Ifges(6,5);
t358 = t226 * mrSges(4,1) - t222 * mrSges(4,2);
t372 = -mrSges(3,1) - t358;
t219 = qJ(3) + qJ(4);
t207 = sin(t219);
t209 = cos(t219);
t371 = t209 * mrSges(6,1) - t376 * t207;
t370 = m(5) * t70 + t374;
t369 = -mrSges(6,2) - mrSges(5,3) - mrSges(4,3) + mrSges(3,2);
t214 = qJDD(1) + qJDD(2);
t285 = qJD(3) * t222;
t139 = t214 * t226 - t216 * t285;
t367 = t139 / 0.2e1;
t333 = t150 / 0.2e1;
t363 = -Ifges(5,6) + Ifges(6,6);
t149 = t221 * t222 - t225 * t226;
t131 = t149 * t216;
t126 = Ifges(5,4) * t131;
t312 = Ifges(6,5) * t131;
t362 = t366 * t132 + t364 * t215 - t126 + t312;
t361 = -m(6) * qJ(5) - mrSges(6,3);
t282 = qJD(4) * t225;
t73 = -t118 * t225 - t302;
t360 = pkin(3) * t282 - t73;
t227 = cos(qJ(2));
t310 = pkin(1) * qJD(2);
t275 = qJD(1) * t310;
t303 = pkin(1) * qJDD(1);
t147 = t223 * t303 + t227 * t275;
t136 = pkin(7) * t214 + t147;
t89 = t226 * t136 - t156 * t285;
t284 = qJD(3) * t226;
t90 = -t136 * t222 - t156 * t284;
t359 = -t222 * t90 + t226 * t89;
t316 = mrSges(5,3) * t131;
t106 = -mrSges(5,2) * t215 - t316;
t309 = t131 * mrSges(6,2);
t109 = mrSges(6,3) * t215 - t309;
t293 = t225 * t119;
t71 = t111 * t221 + t293;
t61 = qJ(5) * t215 + t71;
t355 = m(5) * t71 + m(6) * t61 + t106 + t109;
t213 = qJDD(3) + qJDD(4);
t140 = t214 * t222 + t216 * t284;
t238 = t149 * qJD(4);
t63 = t139 * t221 + t140 * t225 - t216 * t238;
t45 = mrSges(5,1) * t213 - mrSges(5,3) * t63;
t46 = -t213 * mrSges(6,1) + t63 * mrSges(6,2);
t353 = t46 - t45;
t239 = t150 * qJD(4);
t64 = -t225 * t139 + t140 * t221 + t216 * t239;
t47 = -mrSges(5,2) * t213 - mrSges(5,3) * t64;
t48 = -mrSges(6,2) * t64 + mrSges(6,3) * t213;
t352 = t47 + t48;
t349 = qJD(5) + t360;
t289 = t209 * pkin(4) + t207 * qJ(5);
t346 = -t209 * mrSges(5,1) - t371;
t220 = qJ(1) + qJ(2);
t208 = sin(t220);
t210 = cos(t220);
t345 = g(1) * t210 + g(2) * t208;
t146 = -t223 * t275 + t227 * t303;
t135 = -pkin(2) * t214 - t146;
t344 = m(4) * t135 - mrSges(4,1) * t139 + mrSges(4,2) * t140;
t66 = qJDD(3) * pkin(3) - pkin(8) * t140 + t90;
t69 = pkin(8) * t139 + t89;
t10 = -qJD(4) * t71 - t221 * t69 + t225 * t66;
t324 = pkin(3) * t226;
t199 = pkin(2) + t324;
t300 = t208 * t209;
t343 = mrSges(5,1) * t300 + t369 * t210 + (-m(6) * (-t199 - t289) + t371 - t372) * t208;
t299 = t209 * t210;
t301 = t207 * t210;
t342 = t369 * t208 + t372 * t210 - t299 * t375 + t376 * t301;
t279 = t227 * t311;
t157 = -pkin(2) * t216 - t279;
t328 = pkin(1) * t227;
t329 = pkin(1) * t223;
t341 = m(4) * pkin(1) * (t157 * t223 + (t222 ^ 2 + t226 ^ 2) * t227 * t156) + (-mrSges(3,1) * t329 - mrSges(3,2) * t328) * t216;
t297 = t216 * t222;
t154 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t297;
t296 = t216 * t226;
t155 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t296;
t340 = m(4) * t359 - t154 * t284 - t155 * t285 + t226 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t139) - t222 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t140);
t229 = -pkin(8) - pkin(7);
t337 = -t131 / 0.2e1;
t336 = t131 / 0.2e1;
t334 = t132 / 0.2e1;
t330 = t215 / 0.2e1;
t327 = pkin(3) * t221;
t326 = pkin(3) * t222;
t325 = pkin(3) * t225;
t224 = sin(qJ(1));
t321 = t224 * pkin(1);
t228 = cos(qJ(1));
t212 = t228 * pkin(1);
t197 = pkin(7) + t329;
t320 = -pkin(8) - t197;
t318 = mrSges(5,2) * t209;
t315 = Ifges(4,4) * t222;
t314 = Ifges(4,4) * t226;
t313 = Ifges(5,4) * t132;
t298 = t210 * t229;
t288 = t210 * pkin(2) + t208 * pkin(7);
t286 = qJD(3) * t216;
t283 = qJD(4) * t221;
t281 = pkin(3) * t297;
t278 = t227 * t310;
t203 = t223 * t310;
t202 = pkin(3) * t285;
t273 = qJD(3) * t229;
t268 = t361 * t300;
t267 = t361 * t299;
t266 = -pkin(2) * t208 + t210 * pkin(7);
t264 = qJD(3) * t320;
t259 = t210 * t199 - t208 * t229;
t258 = t222 * t278;
t257 = t226 * t278;
t256 = t226 * t273;
t254 = mrSges(4,1) * t222 + mrSges(4,2) * t226;
t253 = t226 * Ifges(4,2) + t315;
t252 = Ifges(4,5) * t226 - Ifges(4,6) * t222;
t83 = pkin(4) * t132 + qJ(5) * t131;
t144 = t320 * t222;
t211 = t226 * pkin(8);
t145 = t197 * t226 + t211;
t249 = t225 * t144 - t145 * t221;
t95 = t144 * t221 + t145 * t225;
t172 = t229 * t222;
t173 = pkin(7) * t226 + t211;
t248 = t225 * t172 - t173 * t221;
t115 = t172 * t221 + t173 * t225;
t247 = -t208 * t199 - t298;
t242 = pkin(4) * t299 + qJ(5) * t301 + t259;
t9 = t111 * t282 - t119 * t283 + t221 * t66 + t225 * t69;
t241 = t157 * t254;
t240 = t222 * (Ifges(4,1) * t226 - t315);
t98 = pkin(4) * t149 - qJ(5) * t150 - t199;
t236 = t318 + (m(6) * pkin(4) + t375) * t207;
t134 = -t199 * t216 - t279;
t100 = -qJD(3) * t149 - t238;
t101 = qJD(3) * t150 + t239;
t31 = pkin(4) * t101 - qJ(5) * t100 - qJD(5) * t150 + t202;
t235 = m(6) * (-pkin(4) * t207 - t326) - t207 * mrSges(6,1);
t92 = -pkin(3) * t139 + t135;
t233 = t226 * t264 - t258;
t125 = Ifges(6,5) * t132;
t5 = qJ(5) * t213 + qJD(5) * t215 + t9;
t6 = -pkin(4) * t213 + qJDD(5) - t10;
t62 = pkin(4) * t131 - qJ(5) * t132 + t134;
t77 = Ifges(6,6) * t215 + Ifges(6,3) * t131 + t125;
t78 = -Ifges(5,2) * t131 + Ifges(5,6) * t215 + t313;
t231 = t61 * t317 + t78 * t334 + (Ifges(6,3) * t132 - t312) * t337 - t6 * mrSges(6,1) - t9 * mrSges(5,2) + t10 * mrSges(5,1) + t5 * mrSges(6,3) - t70 * t316 + t309 * t55 - t62 * (mrSges(6,1) * t132 + mrSges(6,3) * t131) - t134 * (mrSges(5,1) * t132 - mrSges(5,2) * t131) + t363 * t64 + t364 * t63 - (-t131 * t364 + t132 * t363) * t215 / 0.2e1 + (Ifges(6,2) + Ifges(5,3)) * t213 + (-Ifges(5,2) * t132 - t126 + t362) * t336 - (-t366 * t131 + t125 - t313 + t77) * t132 / 0.2e1;
t129 = Ifges(4,6) * qJD(3) + t216 * t253;
t176 = Ifges(4,4) * t296;
t130 = Ifges(4,1) * t297 + Ifges(4,5) * qJD(3) + t176;
t7 = pkin(4) * t64 - qJ(5) * t63 - qJD(5) * t132 + t92;
t230 = qJDD(3) * (Ifges(4,5) * t222 + Ifges(4,6) * t226) + Ifges(3,3) * t214 + t92 * (mrSges(5,1) * t149 + mrSges(5,2) * t150) + t7 * (mrSges(6,1) * t149 - mrSges(6,3) * t150) + t146 * mrSges(3,1) - t147 * mrSges(3,2) + (Ifges(4,1) * t140 + Ifges(4,4) * t367) * t222 + (t363 * t149 + t364 * t150) * t213 / 0.2e1 + (Ifges(6,3) * t336 - Ifges(5,2) * t337 + t134 * mrSges(5,1) + t77 / 0.2e1 + t62 * mrSges(6,1) - t78 / 0.2e1 + t365 * t334 + t363 * t330) * t101 + (t364 * t213 + t366 * t63) * t333 + (t365 * t149 + t366 * t150) * t63 / 0.2e1 + ((Ifges(5,2) + Ifges(6,3)) * t149 + 0.2e1 * t365 * t333) * t64 - t129 * t285 / 0.2e1 + t240 * t286 / 0.2e1 + (-t61 * t101 - t5 * t149 + t150 * t6) * mrSges(6,2) + (-t10 * t150 - t101 * t71 - t149 * t9) * mrSges(5,3) + (t55 * mrSges(6,2) + t134 * mrSges(5,2) - t62 * mrSges(6,3) + Ifges(5,4) * t337 + Ifges(6,5) * t336 + t364 * t330 + t366 * t334 + t362 / 0.2e1 - t70 * mrSges(5,3)) * t100 - t149 * (Ifges(5,4) * t63 + Ifges(5,6) * t213) / 0.2e1 + t359 * mrSges(4,3) + t149 * (Ifges(6,5) * t63 + Ifges(6,6) * t213) / 0.2e1 + (t216 * (-Ifges(4,2) * t222 + t314) + t130) * t284 / 0.2e1 + t226 * (Ifges(4,4) * t140 + Ifges(4,2) * t139) / 0.2e1 + (t241 + t252 * qJD(3) / 0.2e1) * qJD(3) + t253 * t367 + t140 * t314 / 0.2e1 - t135 * t358;
t198 = -pkin(4) - t325;
t191 = qJ(5) + t327;
t169 = -t199 - t328;
t153 = t203 + t202;
t152 = t222 * t273;
t138 = t358 * t216;
t116 = t222 * t264 + t257;
t91 = t98 - t328;
t87 = mrSges(5,1) * t131 + mrSges(5,2) * t132;
t86 = mrSges(6,1) * t131 - mrSges(6,3) * t132;
t74 = t83 + t281;
t72 = -t118 * t221 + t293;
t25 = t203 + t31;
t17 = mrSges(5,1) * t64 + mrSges(5,2) * t63;
t16 = mrSges(6,1) * t64 - mrSges(6,3) * t63;
t1 = [m(3) * (t146 * t227 + t147 * t223) * pkin(1) + t169 * t17 + t153 * t87 + t91 * t16 + t230 + t25 * t86 + m(5) * (t134 * t153 + t169 * t92) + m(6) * (t25 * t62 + t7 * t91) - t138 * t203 + t155 * t257 - t154 * t258 + Ifges(2,3) * qJDD(1) + (m(5) * t9 + m(6) * t5 + t352) * t95 - (-m(5) * t10 + m(6) * t6 + t353) * t249 + (mrSges(3,1) * t328 - mrSges(3,2) * t329) * t214 + t344 * (-pkin(2) - t328) - t370 * (qJD(4) * t95 + t116 * t221 - t225 * t233) + t355 * (qJD(4) * t249 + t225 * t116 + t221 * t233) + t341 * qJD(2) + t340 * t197 + (-mrSges(2,1) * t228 + mrSges(2,2) * t224 - m(3) * t212 - m(4) * (t212 + t288) - m(6) * (t212 + t242) - m(5) * (t212 + t259) + t342) * g(2) + (mrSges(2,1) * t224 + mrSges(2,2) * t228 + m(3) * t321 - m(6) * (-t298 - t321) - m(5) * (t247 - t321) - m(4) * (t266 - t321) + t343) * g(1); t98 * t16 - t199 * t17 + t87 * t202 + t31 * t86 + t230 + (t154 * t222 - t155 * t226) * t279 + t352 * t115 - t353 * t248 - t344 * pkin(2) + (t115 * t5 - t248 * t6 + t31 * t62 + t7 * t98) * m(6) + (t10 * t248 + t115 * t9 + t134 * t202 - t199 * t92) * m(5) + (-m(5) * t134 - m(6) * t62 + t138 - t86 - t87) * t280 - t341 * qJD(1) + (-m(4) * t288 - m(5) * t259 - m(6) * t242 + t342) * g(2) + (-m(4) * t266 - m(5) * t247 + m(6) * t298 + t343) * g(1) + t340 * pkin(7) + t370 * (-qJD(4) * t115 + t150 * t279 - t152 * t221 + t225 * t256) + t355 * (qJD(4) * t248 + t149 * t279 + t225 * t152 + t221 * t256); t45 * t325 + t47 * t327 + t191 * t48 + t198 * t46 + Ifges(4,5) * t140 + Ifges(4,6) * t139 - t89 * mrSges(4,2) + t90 * mrSges(4,1) + (t191 * t5 + t198 * t6 + t349 * t61 - t62 * t74) * m(6) + ((t10 * t225 + t221 * t9 + (-t221 * t70 + t225 * t71) * qJD(4)) * pkin(3) - t134 * t281 + t70 * t72 - t71 * t73) * m(5) + (m(5) * t326 + mrSges(5,1) * t207 + t254 + t318) * t345 + t129 * t297 / 0.2e1 - t252 * t286 / 0.2e1 - t87 * t281 + t231 + t374 * (-pkin(3) * t283 + t72) + t360 * t106 - t74 * t86 + t71 * t308 - g(1) * (t210 * t235 - t267) - g(2) * (t208 * t235 - t268) + t349 * t109 + (-t241 - t240 * t216 / 0.2e1) * t216 + (t154 * t226 + t155 * t222) * t156 + (-t358 - m(6) * (t289 + t324) - m(5) * t324 + t346) * g(3) + Ifges(4,3) * qJDD(3) - (-Ifges(4,2) * t297 + t130 + t176) * t296 / 0.2e1; (t308 + t347) * t71 - t70 * t106 + t231 - t83 * t86 - pkin(4) * t46 + qJ(5) * t48 + (t210 * t236 + t267) * g(1) + (t208 * t236 + t268) * g(2) + t346 * g(3) + t287 * t109 + (-pkin(4) * t6 - t289 * g(3) + qJ(5) * t5 + t287 * t61 - t55 * t71 - t62 * t83) * m(6); -t215 * t109 + t132 * t86 + (g(3) * t209 + t132 * t62 - t207 * t345 - t215 * t61 + t6) * m(6) + t46;];
tau = t1;
