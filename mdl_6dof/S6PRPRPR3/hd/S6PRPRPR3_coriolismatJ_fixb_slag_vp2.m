% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRPR3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:37
% EndTime: 2019-03-08 19:34:43
% DurationCPUTime: 3.00s
% Computational Cost: add. (4672->335), mult. (11564->493), div. (0->0), fcn. (11468->10), ass. (0->202)
t196 = sin(qJ(6));
t191 = t196 ^ 2;
t199 = cos(qJ(6));
t193 = t199 ^ 2;
t257 = t191 + t193;
t327 = m(6) / 0.2e1;
t195 = sin(pkin(11));
t198 = sin(qJ(2));
t279 = sin(pkin(6));
t280 = cos(pkin(11));
t239 = t280 * t279;
t311 = cos(qJ(2));
t240 = t311 * t279;
t125 = t195 * t240 + t198 * t239;
t197 = sin(qJ(4));
t200 = cos(qJ(4));
t281 = cos(pkin(6));
t98 = t125 * t200 + t281 * t197;
t319 = m(7) * t98;
t246 = t198 * t279;
t124 = t195 * t246 - t311 * t239;
t97 = t125 * t197 - t281 * t200;
t60 = -t124 * t196 + t199 * t97;
t288 = t199 * t60;
t61 = t124 * t199 + t196 * t97;
t295 = t196 * t61;
t339 = t288 + t295 - t97;
t308 = m(7) * t197;
t337 = pkin(2) * t195;
t336 = t280 * pkin(2);
t298 = Ifges(7,6) * t199;
t301 = Ifges(7,5) * t196;
t222 = t301 / 0.2e1 + t298 / 0.2e1;
t335 = Ifges(5,4) + Ifges(6,6) - t222;
t289 = t199 * mrSges(7,2);
t297 = t196 * mrSges(7,1);
t166 = t289 + t297;
t334 = t166 + mrSges(6,3);
t304 = mrSges(7,3) * t200;
t333 = t257 * t304 / 0.2e1;
t262 = t200 * qJ(5);
t332 = t197 * pkin(4) - t262;
t185 = pkin(8) + t337;
t306 = pkin(5) + t185;
t150 = t306 * t200;
t157 = pkin(9) * t197 + t332;
t83 = t199 * t150 - t157 * t196;
t286 = t199 * t83;
t84 = t196 * t150 + t157 * t199;
t293 = t196 * t84;
t229 = t286 + t293;
t264 = t199 * t197;
t68 = -t124 * t264 - t125 * t196;
t287 = t199 * t68;
t270 = t196 * t197;
t69 = -t124 * t270 + t125 * t199;
t294 = t196 * t69;
t231 = t287 + t294;
t186 = -pkin(3) - t336;
t268 = t197 * qJ(5);
t211 = t186 - t268;
t307 = t200 * pkin(4);
t145 = t211 - t307;
t165 = t200 * mrSges(6,2) - t197 * mrSges(6,3);
t330 = -m(6) * t145 - t165;
t326 = m(7) / 0.2e1;
t329 = -t262 * t326 + t327 * t332;
t302 = Ifges(7,4) * t199;
t234 = -Ifges(7,2) * t196 + t302;
t303 = Ifges(7,4) * t196;
t236 = Ifges(7,1) * t199 - t303;
t312 = t199 / 0.2e1;
t315 = t196 / 0.2e1;
t328 = t234 * t312 + t236 * t315;
t325 = mrSges(7,1) / 0.2e1;
t324 = -mrSges(7,2) / 0.2e1;
t322 = Ifges(7,3) / 0.2e1;
t321 = t60 / 0.2e1;
t201 = -pkin(4) - pkin(9);
t275 = t124 * t197;
t320 = m(7) * (-t231 - t275) * t200;
t148 = t200 * t166;
t318 = -t148 / 0.2e1;
t317 = t150 / 0.2e1;
t316 = -t196 / 0.2e1;
t314 = t197 / 0.2e1;
t313 = -t199 / 0.2e1;
t310 = m(6) * t332;
t242 = (0.1e1 - t257) * t200;
t309 = t242 * t308;
t290 = t199 * mrSges(7,1);
t296 = t196 * mrSges(7,2);
t237 = t290 - t296;
t147 = t237 * t200;
t146 = t197 * t237;
t149 = t306 * t197;
t269 = t196 * t200;
t292 = t197 * mrSges(7,1);
t159 = mrSges(7,3) * t269 + t292;
t263 = t199 * t200;
t291 = t197 * mrSges(7,2);
t161 = -mrSges(7,3) * t263 - t291;
t129 = t201 * t200 + t211;
t78 = -t129 * t196 + t149 * t199;
t79 = t129 * t199 + t149 * t196;
t230 = t196 * t79 + t199 * t78;
t204 = -t146 / 0.2e1 + t159 * t312 + (-t149 + t230) * t326 + t161 * t315;
t285 = t200 * mrSges(7,1);
t158 = -mrSges(7,3) * t270 + t285;
t266 = t199 * t158;
t284 = t200 * mrSges(7,2);
t160 = mrSges(7,3) * t264 - t284;
t271 = t196 * t160;
t216 = -t271 / 0.2e1 - t266 / 0.2e1;
t18 = (t147 / 0.2e1 + (t150 - t229) * t326 + t216) * t200 + t204 * t197;
t194 = t200 ^ 2;
t265 = t199 * t161;
t272 = t196 * t159;
t217 = t272 / 0.2e1 - t265 / 0.2e1;
t248 = -t193 / 0.2e1 - t191 / 0.2e1;
t241 = mrSges(7,3) * t248;
t25 = t194 * t241 + t197 * t318 + t217 * t200;
t305 = t18 * qJD(4) + t25 * qJD(6);
t300 = Ifges(7,5) * t197;
t299 = Ifges(7,6) * t197;
t283 = t97 * qJ(5);
t278 = qJ(5) * t148;
t274 = t124 * t200;
t235 = Ifges(7,1) * t196 + t302;
t214 = t235 * t200;
t133 = -t214 + t300;
t273 = t196 * t133;
t233 = Ifges(7,2) * t199 + t303;
t213 = t233 * t200;
t131 = -t213 + t299;
t267 = t199 * t131;
t261 = t201 * t159;
t260 = t201 * t161;
t259 = t25 * qJD(2);
t28 = (t291 / 0.2e1 - t161 / 0.2e1) * t199 + (t292 / 0.2e1 + t159 / 0.2e1) * t196 + t200 * t241;
t258 = t28 * qJD(2);
t256 = qJD(4) * t197;
t255 = qJD(4) * t200;
t254 = -t320 / 0.2e1;
t253 = t320 / 0.2e1;
t252 = mrSges(5,2) - t334;
t250 = -t275 / 0.2e1;
t247 = m(7) * t257;
t245 = t257 * t201;
t238 = -t200 * mrSges(5,1) + t197 * mrSges(5,2);
t54 = t193 * Ifges(7,4) - qJ(5) * t237 + (-t303 + (Ifges(7,1) - Ifges(7,2)) * t199) * t196;
t224 = t84 * t324 + t83 * t325;
t8 = t278 / 0.2e1 + (t201 * t241 + t322) * t200 + (0.3e1 / 0.4e1 * t299 - t150 * mrSges(7,1) / 0.2e1 + t131 / 0.4e1 - t260 / 0.2e1 + (Ifges(7,1) / 0.2e1 - Ifges(7,2) / 0.4e1) * t263) * t199 + (0.3e1 / 0.4e1 * t300 + mrSges(7,2) * t317 + t133 / 0.4e1 + t261 / 0.2e1 + (-0.3e1 / 0.2e1 * t302 + (Ifges(7,2) / 0.2e1 - Ifges(7,1) / 0.4e1) * t196) * t200) * t196 + t224;
t228 = t8 * qJD(2) + t54 * qJD(4);
t218 = m(7) * t231;
t209 = -t218 / 0.2e1;
t210 = (t196 * t60 - t199 * t61) * t326;
t20 = t197 * t210 + t209;
t26 = (t272 + m(7) * (t196 * t78 - t199 * t79) - t265 + t330) * t197;
t227 = qJD(1) * t20 + qJD(2) * t26;
t11 = t339 * t319;
t13 = (t197 * t339 + t98 * t242) * t326;
t226 = t11 * qJD(1) + t13 * qJD(3);
t225 = t69 * t324 + t68 * t325;
t223 = -t257 * mrSges(7,3) - mrSges(5,1) + mrSges(6,2);
t221 = -t296 / 0.2e1 + t290 / 0.2e1;
t168 = -t197 * mrSges(6,2) - t200 * mrSges(6,3);
t169 = t197 * mrSges(5,1) + t200 * mrSges(5,2);
t220 = t310 / 0.2e1 + t169 / 0.2e1 + t168 / 0.2e1;
t65 = t98 * t274;
t80 = t124 * t125;
t7 = m(6) * (t80 + (-t197 * t97 - t200 * t98) * t124) + m(7) * (t60 * t68 + t61 * t69 - t65) + m(5) * (-t97 * t275 - t65 + t80);
t215 = -t7 * qJD(1) + qJD(3) * t254;
t202 = t204 * t98 + (-t150 * t97 + t83 * t60 + t84 * t61) * t326 + t158 * t321 + t61 * t160 / 0.2e1 - t97 * t147 / 0.2e1;
t1 = (t294 / 0.2e1 + t287 / 0.2e1) * mrSges(7,3) + t201 * t209 + ((-mrSges(5,1) / 0.2e1 + mrSges(6,2) / 0.2e1) * t197 + (t166 / 0.2e1 - mrSges(5,2) / 0.2e1 + mrSges(6,3) / 0.2e1) * t200 + t220 - t329) * t124 + t202;
t130 = Ifges(7,6) * t200 + t197 * t233;
t132 = Ifges(7,5) * t200 + t197 * t235;
t5 = t186 * t169 + t332 * t165 + t78 * t158 + t83 * t159 + t79 * t160 + t84 * t161 - t149 * t147 - t150 * t146 + m(7) * (-t149 * t150 + t78 * t83 + t79 * t84) + (t168 + t310) * t145 + (t267 / 0.2e1 + t273 / 0.2e1 - t335 * t197) * t197 + (t130 * t313 + t132 * t316 + (Ifges(6,2) - Ifges(6,3) + Ifges(5,1) - Ifges(5,2) + Ifges(7,3)) * t197 + t335 * t200) * t200;
t208 = t1 * qJD(1) + t5 * qJD(2) + t18 * qJD(3);
t183 = Ifges(7,6) * t269;
t10 = -t150 * t148 + t183 * t314 - t79 * t159 + t78 * t161 + (-Ifges(7,5) * t264 / 0.2e1 + t133 * t313 + t131 * t315 + t328 * t200 + t230 * mrSges(7,3)) * t200;
t203 = (t295 / 0.2e1 + t288 / 0.2e1) * t304 + t161 * t321 - t61 * t159 / 0.2e1 + t98 * t318;
t3 = t203 - t225;
t207 = t3 * qJD(1) + t10 * qJD(2) + t25 * qJD(3);
t206 = t13 * qJD(1) + t18 * qJD(2) + qJD(3) * t309;
t112 = (-0.1e1 / 0.2e1 - t248) * t308;
t144 = (m(6) + m(7)) * qJ(5) + t334;
t30 = (t285 / 0.2e1 - t158 / 0.2e1) * t199 + (-t284 / 0.2e1 - t160 / 0.2e1) * t196 + 0.2e1 * (t150 / 0.4e1 - t293 / 0.4e1 - t286 / 0.4e1) * m(7);
t41 = 0.2e1 * (t191 / 0.4e1 + t193 / 0.4e1 - 0.1e1 / 0.4e1) * t319;
t205 = qJD(1) * t41 - qJD(2) * t30 + qJD(3) * t112 - qJD(4) * t144;
t101 = t308 / 0.2e1 + (m(6) + t247 / 0.2e1) * t197;
t96 = t221 * t197 + t237 * t314;
t31 = t319 / 0.2e1 + 0.2e1 * (t327 + t247 / 0.4e1) * t98;
t29 = (t297 / 0.2e1 + t289 / 0.2e1) * t197 - t217 + t333;
t27 = m(7) * t317 + t229 * t326 + (m(6) * t185 + mrSges(6,1) + t221) * t200 - t216;
t19 = m(6) * t250 + t218 / 0.2e1 + (-m(6) * t124 / 0.2e1 + t210) * t197;
t16 = (t237 / 0.2e1 + t221) * t98;
t9 = t237 * t317 - t273 / 0.4e1 - t236 * t263 / 0.2e1 - t278 / 0.2e1 + t197 * (-t298 - t301) / 0.4e1 - t267 / 0.4e1 + t234 * t269 / 0.2e1 + t196 * t214 / 0.4e1 + t199 * t213 / 0.4e1 + t261 * t316 + t260 * t312 + t200 * t322 + t222 * t197 + t224 + t333 * t201;
t6 = qJD(2) * t253 + t13 * qJD(4);
t4 = t203 + t225;
t2 = mrSges(5,1) * t275 / 0.2e1 + mrSges(6,2) * t250 + t202 + (t220 + t329) * t124 + (mrSges(5,2) / 0.2e1 - t334 / 0.2e1) * t274 + (t201 * t326 - mrSges(7,3) / 0.2e1) * t231;
t12 = [t7 * qJD(2) + t11 * qJD(4), t2 * qJD(4) + t19 * qJD(5) + t4 * qJD(6) - t215 + (-mrSges(3,1) * t246 - mrSges(3,2) * t240 + m(7) * (-t150 * t274 + t78 * t68 + t79 * t69) + t69 * t161 + t68 * t159 - t147 * t274 + (-m(4) * t336 + m(5) * t186 - mrSges(4,1) + t238 - t330) * t125 + (-m(4) * t337 + mrSges(4,2) + (mrSges(5,3) + mrSges(6,1) + (m(6) + m(5)) * t185) * (-t197 ^ 2 - t194)) * t124) * qJD(2), t6, t2 * qJD(2) + t31 * qJD(5) + t16 * qJD(6) + t226 + (t252 * t97 + t223 * t98 + 0.2e1 * (-pkin(4) * t98 - t283) * t327 + 0.2e1 * (t245 * t98 - t283) * t326) * qJD(4), qJD(2) * t19 + qJD(4) * t31, t4 * qJD(2) + t16 * qJD(4) + (-mrSges(7,1) * t61 - mrSges(7,2) * t60) * qJD(6); t1 * qJD(4) + t20 * qJD(5) + t3 * qJD(6) + t215, qJD(4) * t5 + qJD(5) * t26 + qJD(6) * t10, qJD(1) * t254 + t305, t27 * qJD(5) + t9 * qJD(6) + (-pkin(4) * mrSges(6,1) + Ifges(7,5) * t312 + Ifges(7,6) * t316 - Ifges(6,4) + Ifges(5,5)) * t255 + (-qJ(5) * mrSges(6,1) + Ifges(6,5) - Ifges(5,6) + t328) * t256 + t208 + (-t149 * t166 + t132 * t312 + t130 * t316 + (m(6) * (-t268 - t307) + t165 + t238) * t185 + (-m(7) * t149 - t146) * qJ(5) + (m(7) * t229 + t266 + t271) * t201 - t229 * mrSges(7,3)) * qJD(4), qJD(4) * t27 + qJD(6) * t29 + t227, t9 * qJD(4) + t29 * qJD(5) + (-mrSges(7,1) * t79 - mrSges(7,2) * t78 - Ifges(7,5) * t263 + t183) * qJD(6) + t207; t6, qJD(1) * t253 + t305, qJD(4) * t309, -qJD(4) * t310 + t101 * qJD(5) + t96 * qJD(6) + (m(7) * qJ(5) - t252) * t255 + (m(7) * t245 + t223) * t256 + t206, t101 * qJD(4), t96 * qJD(4) + t148 * qJD(6) + t259; -qJD(2) * t1 - qJD(5) * t41 - t226, qJD(5) * t30 - qJD(6) * t8 - t208, -t112 * qJD(5) - t206, qJD(5) * t144 - qJD(6) * t54, -t205 ((-mrSges(7,2) * t201 - Ifges(7,6)) * t199 + (-mrSges(7,1) * t201 - Ifges(7,5)) * t196) * qJD(6) - t228; -qJD(2) * t20 + qJD(4) * t41, -qJD(4) * t30 - qJD(6) * t28 - t227, t112 * qJD(4), t205, 0, -t166 * qJD(6) - t258; -t3 * qJD(2), qJD(4) * t8 + qJD(5) * t28 - t207, -t259, t228, t258, 0;];
Cq  = t12;
