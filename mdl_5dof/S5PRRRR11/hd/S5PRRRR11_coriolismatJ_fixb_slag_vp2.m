% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR11_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:16
% EndTime: 2024-09-27 21:45:20
% DurationCPUTime: 2.41s
% Computational Cost: add. (10587->304), mult. (28763->413), div. (0->0), fcn. (30031->8), ass. (0->182)
t191 = sin(pkin(5));
t194 = sin(qJ(4));
t195 = sin(qJ(3));
t197 = cos(qJ(4));
t198 = cos(qJ(3));
t166 = (-t194 * t195 + t197 * t198) * t191;
t249 = t191 * t198;
t250 = t191 * t195;
t167 = -t194 * t249 - t197 * t250;
t193 = sin(qJ(5));
t196 = cos(qJ(5));
t136 = t166 * t193 - t167 * t196;
t323 = t136 * mrSges(6,1);
t297 = -t323 / 0.2e1;
t228 = t196 * t166 + t167 * t193;
t315 = t228 * mrSges(6,2);
t224 = t297 - t315 / 0.2e1;
t179 = (-pkin(3) * t198 - pkin(2)) * t191;
t138 = -t166 * pkin(4) + t179;
t192 = cos(pkin(5));
t270 = Ifges(6,4) * t136;
t291 = t136 / 0.2e1;
t318 = t228 / 0.2e1;
t328 = -t136 / 0.2e1;
t331 = t323 + t315;
t284 = pkin(2) * t192;
t187 = t198 * t284;
t298 = pkin(7) + pkin(8);
t152 = -t250 * t298 + t187;
t143 = pkin(3) * t192 + t152;
t236 = t195 * t284;
t153 = t249 * t298 + t236;
t148 = t197 * t153;
t100 = t143 * t194 + t148;
t280 = pkin(9) * t166;
t76 = t100 + t280;
t259 = t193 * t76;
t282 = pkin(4) * t192;
t163 = t167 * pkin(9);
t146 = t194 * t153;
t99 = t197 * t143 - t146;
t75 = t163 + t99;
t73 = t282 + t75;
t38 = t196 * t73 - t259;
t258 = t196 * t76;
t39 = t193 * t73 + t258;
t332 = t138 * t331 + (-t136 * t39 - t228 * t38) * mrSges(6,3) + (Ifges(6,2) * t228 + Ifges(6,6) * t192 + t270) * t328 + (Ifges(6,1) * t228 - t270) * t291 + (Ifges(6,5) * t192 + 0.2e1 * Ifges(6,4) * t228 + (Ifges(6,1) - Ifges(6,2)) * t136) * t318;
t316 = Ifges(6,5) * t228;
t324 = Ifges(6,6) * t136;
t240 = t316 - t324;
t310 = t324 / 0.2e1 - t316 / 0.2e1;
t325 = mrSges(6,3) * t228;
t105 = -t192 * mrSges(6,2) + t325;
t263 = t136 * mrSges(6,3);
t106 = mrSges(6,1) * t192 - t263;
t330 = t105 * t318 + t106 * t328;
t283 = pkin(3) * t197;
t190 = pkin(4) + t283;
t247 = t193 * t194;
t171 = -pkin(3) * t247 + t190 * t196;
t313 = t171 * t228;
t327 = -t313 / 0.2e1;
t321 = qJD(3) + qJD(4);
t161 = Ifges(5,5) * t166;
t159 = t167 * mrSges(5,1);
t262 = t166 * mrSges(5,2);
t231 = -t159 + t262;
t312 = t331 + t231;
t160 = Ifges(5,6) * t167;
t226 = t160 + t161 + t240;
t294 = -t228 / 0.2e1;
t205 = t297 + t323 / 0.2e1 + (t294 + t318) * mrSges(6,2);
t311 = qJD(1) * t205;
t278 = t39 * mrSges(6,1);
t279 = t38 * mrSges(6,2);
t309 = -t278 / 0.2e1 - t279 / 0.2e1;
t245 = t194 * t196;
t176 = (-t193 * t197 - t245) * pkin(3);
t173 = t176 * mrSges(6,1);
t177 = (t196 * t197 - t247) * pkin(3);
t260 = t177 * mrSges(6,2);
t308 = (mrSges(5,1) * t194 + mrSges(5,2) * t197) * pkin(3) - t173 + t260;
t234 = -pkin(4) * t196 / 0.2e1;
t43 = t196 * t75 - t259;
t276 = t43 * mrSges(6,2);
t42 = -t193 * t75 - t258;
t277 = t42 * mrSges(6,1);
t307 = t234 * t325 + t277 / 0.2e1 - t276 / 0.2e1;
t306 = qJD(5) * t205;
t212 = 0.2e1 * t224;
t305 = t212 * qJD(5);
t304 = m(5) / 0.2e1;
t303 = m(6) / 0.2e1;
t302 = pkin(4) / 0.2e1;
t301 = m(6) * pkin(4);
t300 = t42 / 0.2e1;
t299 = -t43 / 0.2e1;
t261 = t166 * mrSges(5,3);
t144 = -t192 * mrSges(5,2) + t261;
t290 = t144 / 0.2e1;
t289 = -t171 / 0.2e1;
t172 = pkin(3) * t245 + t190 * t193;
t288 = t172 / 0.2e1;
t287 = t177 / 0.2e1;
t286 = t192 / 0.2e1;
t285 = t193 / 0.2e1;
t281 = pkin(4) * t193;
t103 = -t152 * t194 - t148;
t78 = t103 - t280;
t104 = t197 * t152 - t146;
t79 = t163 + t104;
t46 = -t193 * t79 + t196 * t78;
t275 = t46 * mrSges(6,1);
t47 = t193 * t78 + t196 * t79;
t274 = t47 * mrSges(6,2);
t273 = t99 * mrSges(5,2);
t272 = t275 / 0.2e1 - t274 / 0.2e1;
t271 = -t136 * t171 + t172 * t228;
t269 = t100 * mrSges(5,1);
t268 = t103 * mrSges(5,1);
t267 = t104 * mrSges(5,2);
t237 = pkin(3) * t250;
t142 = -pkin(4) * t167 + t237;
t145 = mrSges(5,1) * t192 + t167 * mrSges(5,3);
t201 = (-t167 ^ 2 / 0.2e1 - t166 ^ 2 / 0.2e1) * mrSges(5,3) + (-t136 * t291 + t228 * t294) * mrSges(6,3) + t166 * t290 + t167 * t145 / 0.2e1 + t330;
t219 = -t136 * t38 + t228 * t39;
t7 = (t136 * t47 + t228 * t46 + t219) * t303 + ((-t104 + t99) * t167 + (t100 + t103) * t166) * t304 + (t262 / 0.2e1 - t159 / 0.2e1 + t142 * t303 + t237 * t304 - t224) * t192 + t201;
t257 = t7 * qJD(2);
t8 = (t136 * t43 - t167 * t282 + t228 * t42 + t219) * t303 + t201 + t312 * t286;
t256 = t8 * qJD(2);
t168 = t172 * mrSges(6,1);
t56 = t171 * mrSges(6,2) + t168;
t255 = qJD(5) * t56;
t13 = t331 * t286 + (-t228 ^ 2 / 0.2e1 - t136 ^ 2 / 0.2e1) * mrSges(6,3) + t330;
t254 = t13 * qJD(2);
t251 = t172 * t136;
t248 = t193 * t136;
t246 = t194 * t167;
t244 = t196 * t105;
t243 = t196 * t228;
t242 = t197 * t166;
t235 = -t281 / 0.2e1;
t233 = -t263 / 0.2e1;
t225 = t289 + t234;
t223 = -(pkin(7) * t249 + t236) * mrSges(4,1) - (-pkin(7) * t250 + t187) * mrSges(4,2);
t162 = Ifges(5,4) * t166;
t185 = Ifges(4,5) * t249;
t200 = -t99 * t261 + t332;
t203 = Ifges(5,6) * t286 + t100 * mrSges(5,3) - Ifges(5,4) * t167 + (-Ifges(5,1) + Ifges(5,2)) * t166;
t69 = -mrSges(6,1) * t228 + mrSges(6,2) * t136;
t2 = ((-Ifges(4,6) * t192 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t195) * t191 + (m(5) * t179 - mrSges(5,1) * t166 - mrSges(5,2) * t167) * pkin(3)) * t195 + (Ifges(4,4) * t249 + Ifges(4,5) * t286 + (-pkin(2) * mrSges(4,2) + (Ifges(4,1) - Ifges(4,2)) * t195) * t191) * t198) * t191 + t203 * t167 + (t185 / 0.2e1 + t161 + t160 / 0.2e1 + t223 - t310) * t192 + m(6) * (t138 * t142 + t38 * t46 + t39 * t47) + t200 + t179 * t231 + t142 * t69 + t104 * t144 + t103 * t145 + m(5) * (t100 * t104 + t103 * t99) + t47 * t105 + t46 * t106 + t166 * t162;
t221 = t7 * qJD(1) + t2 * qJD(2);
t3 = m(6) * (t38 * t42 + t39 * t43) + ((-m(6) * t138 - t69) * pkin(4) + t203) * t167 + (t179 * mrSges(5,2) + t162) * t166 + t200 - t179 * t159 + t99 * t144 - t100 * t145 + t43 * t105 + t42 * t106 + (t226 + t161) * t286;
t220 = t8 * qJD(1) + t3 * qJD(2);
t6 = t38 * t105 - t39 * t106 + t240 * t286 + t332;
t218 = t13 * qJD(1) + t6 * qJD(2);
t217 = -t168 / 0.2e1 + mrSges(6,1) * t235;
t215 = t105 * t289 + t106 * t288;
t214 = -t251 / 0.2e1 + t327;
t211 = (t193 * t47 + t196 * t46) * t301;
t210 = -0.2e1 * t310 + t309;
t208 = (-t136 * t196 + t193 * t228) * t301;
t9 = (t136 * t288 + t313 / 0.2e1) * mrSges(6,3) + t215 + t272 - t309;
t207 = -t9 * qJD(2) - t56 * qJD(3) - t311;
t202 = (t177 * t136 + t176 * t228 + t271) * t303;
t22 = -t208 / 0.2e1 + t202;
t199 = (t197 * t290 - t194 * t145 / 0.2e1 + (t246 / 0.2e1 - t242 / 0.2e1) * mrSges(5,3)) * pkin(3) + (t171 * t42 + t172 * t43 + t176 * t38 + t177 * t39) * t303 + t176 * t106 / 0.2e1 + t105 * t287;
t5 = (t47 / 0.2e1 + t299) * mrSges(6,2) + (t104 / 0.2e1 - t99 / 0.2e1) * mrSges(5,2) + (-t46 / 0.2e1 + t300) * mrSges(6,1) + (-t103 / 0.2e1 - t100 / 0.2e1) * mrSges(5,1) - t211 / 0.2e1 + ((t248 / 0.2e1 + t243 / 0.2e1) * pkin(4) + t214) * mrSges(6,3) + t199;
t55 = -m(6) * (t171 * t176 + t172 * t177) + t308;
t206 = t22 * qJD(1) + t5 * qJD(2) - t55 * qJD(3);
t11 = (t299 + t38 / 0.2e1) * mrSges(6,2) + (t300 + t39 / 0.2e1) * mrSges(6,1) + (-t244 / 0.2e1 + t106 * t285 + (t136 * t285 + t196 * t318) * mrSges(6,3)) * pkin(4);
t180 = (t193 * mrSges(6,1) + t196 * mrSges(6,2)) * pkin(4);
t52 = -t173 / 0.2e1 + (t287 + t225) * mrSges(6,2) + t217;
t204 = qJD(2) * t11 - qJD(3) * t52 + qJD(4) * t180 + t311;
t178 = t180 * qJD(5);
t53 = -t260 / 0.2e1 + t173 / 0.2e1 + t225 * mrSges(6,2) + t217;
t19 = t208 / 0.2e1 + t202 - t312;
t12 = t106 * t235 + t233 * t281 + t244 * t302 + t210 + t307;
t10 = mrSges(6,3) * t327 + t172 * t233 + t210 - t215 + t272;
t4 = t211 / 0.2e1 + t199 - t267 / 0.2e1 - t269 / 0.2e1 + t268 / 0.2e1 - t273 / 0.2e1 + t226 + t272 + (-t248 * t302 + t214) * mrSges(6,3) + t307;
t1 = qJD(3) * t7 + qJD(4) * t8 + qJD(5) * t13;
t14 = [0, t1, t257 + (-mrSges(4,1) * t250 - mrSges(4,2) * t249 - t312) * qJD(3) + t19 * qJD(4) + 0.2e1 * (t271 * t303 + (pkin(3) * t166 * t194 + t167 * t283) * t304) * qJD(3) + t305, t256 + t19 * qJD(3) + (t208 - t312) * qJD(4) + t305, -qJD(5) * t331 + t321 * t212 + t254; t1, qJD(3) * t2 + qJD(4) * t3 + qJD(5) * t6, t4 * qJD(4) + t10 * qJD(5) + t221 + (t185 - Ifges(4,6) * t250 - t274 + m(6) * (t171 * t46 + t172 * t47) + t275 + t268 - t267 + t223 + t226 + (m(5) * (t103 * t197 + t104 * t194) + (-t242 + t246) * mrSges(5,3)) * pkin(3) + (-t251 - t313) * mrSges(6,3)) * qJD(3), t4 * qJD(3) + t12 * qJD(5) + t220 + (t226 - t269 - t273 - t276 + t277 + (m(6) * (t193 * t43 + t196 * t42) + (-t243 - t248) * mrSges(6,3)) * pkin(4)) * qJD(4), t10 * qJD(3) + t12 * qJD(4) + (t240 - t278 - t279) * qJD(5) + t218; qJD(4) * t22 - t257 - t306, qJD(4) * t5 - qJD(5) * t9 - t221, -qJD(4) * t55 - t255, ((t176 * t196 + t177 * t193) * t301 - t308) * qJD(4) + t53 * qJD(5) + t206, t53 * qJD(4) + t207 - t255; -qJD(3) * t22 - t256 - t306, -qJD(3) * t5 - qJD(5) * t11 - t220, qJD(5) * t52 - t206, -t178, -t178 - t204; t321 * t205 - t254, qJD(3) * t9 + qJD(4) * t11 - t218, -qJD(4) * t52 - t207, t204, 0;];
Cq = t14;
