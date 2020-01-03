% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPPR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:31:39
% EndTime: 2019-12-31 19:31:49
% DurationCPUTime: 4.39s
% Computational Cost: add. (9214->324), mult. (18477->467), div. (0->0), fcn. (20730->8), ass. (0->172)
t274 = sin(pkin(8));
t275 = cos(pkin(8));
t301 = sin(qJ(2));
t303 = cos(qJ(2));
t196 = -t274 * t303 - t275 * t301;
t215 = sin(pkin(9));
t216 = cos(pkin(9));
t300 = sin(qJ(5));
t302 = cos(qJ(5));
t327 = t300 * t215 - t302 * t216;
t334 = t327 * t196;
t336 = -t334 / 0.2e1;
t337 = mrSges(6,1) * t336;
t257 = t301 * pkin(2);
t335 = m(4) * t257;
t195 = t274 * t301 - t275 * t303;
t143 = t327 * t195;
t314 = -t143 / 0.2e1;
t307 = -t327 / 0.2e1;
t198 = -t215 * t302 - t216 * t300;
t223 = t198 * t196;
t332 = t223 / 0.2e1;
t333 = mrSges(6,2) * t332;
t214 = t216 ^ 2;
t261 = t215 ^ 2 + t214;
t331 = t261 * mrSges(5,3);
t256 = t301 * pkin(6);
t203 = -qJ(3) * t301 - t256;
t258 = t303 * pkin(6);
t204 = qJ(3) * t303 + t258;
t235 = -t275 * t203 + t204 * t274;
t330 = t196 * t235;
t328 = t274 * t203 + t275 * t204;
t158 = -t196 * pkin(3) + t195 * qJ(4) + t257;
t96 = t216 * t158 + t215 * t235;
t97 = t215 * t158 - t216 * t235;
t326 = -t215 * t96 + t216 * t97;
t103 = -mrSges(6,2) * t195 - mrSges(6,3) * t223;
t104 = mrSges(6,1) * t195 - mrSges(6,3) * t334;
t298 = t103 * t307 + t198 * t104 / 0.2e1;
t325 = t196 ^ 2;
t324 = -m(5) / 0.2e1;
t323 = m(5) / 0.2e1;
t322 = m(6) / 0.2e1;
t321 = m(4) * pkin(2);
t320 = -mrSges(6,3) / 0.2e1;
t140 = t198 * t195;
t317 = t140 / 0.2e1;
t250 = t274 * pkin(2);
t205 = t250 + qJ(4);
t299 = pkin(7) + t205;
t183 = t299 * t215;
t184 = t299 * t216;
t155 = -t183 * t300 + t184 * t302;
t312 = -t155 / 0.2e1;
t294 = Ifges(6,4) * t198;
t166 = -Ifges(6,2) * t327 - t294;
t311 = t166 / 0.2e1;
t190 = Ifges(6,4) * t327;
t168 = -Ifges(6,1) * t198 - t190;
t310 = t168 / 0.2e1;
t309 = -t196 / 0.2e1;
t308 = t196 / 0.2e1;
t306 = -t198 / 0.2e1;
t305 = t215 / 0.2e1;
t304 = t216 / 0.2e1;
t297 = mrSges(4,3) * t196;
t296 = Ifges(5,4) * t215;
t295 = Ifges(5,4) * t216;
t293 = Ifges(5,2) * t215;
t154 = -t183 * t302 - t184 * t300;
t164 = mrSges(6,1) * t327 - t198 * mrSges(6,2);
t252 = t275 * pkin(2);
t212 = -t252 - pkin(3);
t201 = -t216 * pkin(4) + t212;
t202 = -t216 * mrSges(5,1) + t215 * mrSges(5,2);
t283 = t198 * mrSges(6,3);
t254 = t283 / 0.2e1;
t284 = t327 * mrSges(6,3);
t217 = (-t195 * t205 * t261 - t196 * t212) * t323 + (-t140 * t154 + t143 * t155 - t196 * t201) * t322 + (-t195 * t274 + t196 * t275) * t321 / 0.2e1 - t140 * t254 + t284 * t314 + (t164 + t202) * t309 - t195 * t331 / 0.2e1;
t269 = t215 * t195;
t230 = mrSges(5,2) * t196 + mrSges(5,3) * t269;
t267 = t216 * t195;
t231 = -t196 * mrSges(5,1) + mrSges(5,3) * t267;
t239 = t196 * mrSges(6,2) - t140 * mrSges(6,3);
t240 = -t196 * mrSges(6,1) - t143 * mrSges(6,3);
t61 = -t196 * pkin(4) + pkin(7) * t267 + t96;
t73 = pkin(7) * t269 + t97;
t34 = -t300 * t73 + t302 * t61;
t35 = t300 * t61 + t302 * t73;
t218 = (t215 * t97 + t216 * t96) * t323 + (-t198 * t35 - t327 * t34) * t322 + t240 * t307 + t239 * t306 + t230 * t305 + t231 * t304 + t335 / 0.2e1;
t246 = -t196 * mrSges(4,1) - t195 * mrSges(4,2);
t5 = t217 - t218 - t246;
t292 = qJD(1) * t5;
t123 = -pkin(4) * t269 + t328;
t268 = t215 * t196;
t124 = -pkin(4) * t268 + t235;
t160 = -t195 * mrSges(5,2) + mrSges(5,3) * t268;
t266 = t216 * t196;
t161 = t195 * mrSges(5,1) + mrSges(5,3) * t266;
t280 = t216 * mrSges(5,2);
t282 = t215 * mrSges(5,1);
t243 = t280 + t282;
t227 = t243 * t195;
t229 = Ifges(6,5) * t314 + Ifges(6,6) * t317;
t232 = t216 * Ifges(5,5) - t215 * Ifges(5,6) - Ifges(4,4);
t237 = Ifges(6,4) * t143 - Ifges(6,2) * t140;
t238 = Ifges(6,1) * t143 - Ifges(6,4) * t140;
t241 = mrSges(6,1) * t223 + mrSges(6,2) * t334;
t288 = t143 * mrSges(6,2);
t290 = t140 * mrSges(6,1);
t242 = t288 + t290;
t245 = -pkin(2) * t303 - pkin(1);
t157 = t195 * pkin(3) + t196 * qJ(4) + t245;
t90 = t216 * t157 - t215 * t328;
t60 = t195 * pkin(4) + pkin(7) * t266 + t90;
t91 = t215 * t157 + t216 * t328;
t72 = pkin(7) * t268 + t91;
t32 = -t300 * t72 + t302 * t60;
t33 = t300 * t60 + t302 * t72;
t286 = t334 * Ifges(6,4);
t64 = -Ifges(6,2) * t223 + t195 * Ifges(6,6) + t286;
t137 = Ifges(6,4) * t223;
t65 = Ifges(6,1) * t334 + t195 * Ifges(6,5) - t137;
t1 = (Ifges(3,2) - Ifges(3,1)) * t303 * t301 - m(5) * (t235 * t328 + t90 * t96 + t91 * t97) + (-mrSges(4,1) * t257 + t232 * t195 + t229) * t195 + (mrSges(4,2) * t257 + Ifges(6,5) * t334 - Ifges(6,6) * t223 + (-t214 * Ifges(5,1) - Ifges(4,1) + Ifges(4,2) + Ifges(5,3) + Ifges(6,3) + (-t293 + 0.2e1 * t295) * t215) * t195 - t232 * t196 + (t243 + mrSges(4,3)) * t328) * t196 - t328 * t297 + t237 * t332 + t238 * t336 + t65 * t314 + t64 * t317 - m(6) * (t123 * t124 + t32 * t34 + t33 * t35) - t91 * t230 - t90 * t231 - t33 * t239 - t32 * t240 - t123 * t241 - t124 * t242 - t35 * t103 - t34 * t104 - t97 * t160 - t96 * t161 + pkin(1) * (mrSges(3,1) * t301 + mrSges(3,2) * t303) + (t301 ^ 2 - t303 ^ 2) * Ifges(3,4) + (-t246 - t335) * t245 + t235 * t227;
t291 = t1 * qJD(1);
t289 = t223 * mrSges(6,2);
t287 = t334 * mrSges(6,1);
t263 = -Ifges(6,5) * t223 - Ifges(6,6) * t334;
t75 = t287 - t289;
t76 = -Ifges(6,2) * t334 - t137;
t77 = -Ifges(6,1) * t223 - t286;
t4 = t195 * t263 / 0.2e1 + t32 * t103 - t33 * t104 + t124 * t75 + (-t33 * mrSges(6,3) - t64 / 0.2e1 + t77 / 0.2e1) * t334 - (-t32 * mrSges(6,3) + t76 / 0.2e1 + t65 / 0.2e1) * t223;
t278 = t4 * qJD(1);
t236 = t215 * t90 - t91 * t216;
t6 = -t143 * t103 + t140 * t104 - m(6) * (-t124 * t196 - t140 * t32 + t143 * t33) + t196 * t241 + t160 * t267 - t161 * t269 - m(5) * (t195 * t236 - t330) - m(4) * (-t195 * t328 - t330) - t325 * t243 + (-t195 ^ 2 - t325) * mrSges(4,3);
t277 = t6 * qJD(1);
t12 = -t223 * t103 - t334 * t104 + m(6) * (-t223 * t33 - t32 * t334) + (t215 * t160 + t216 * t161 + m(5) * (t215 * t91 + t216 * t90)) * t196;
t273 = qJD(1) * t12;
t272 = t223 * t327;
t225 = (t198 * t223 + t327 * t334) * t322;
t248 = m(5) * t261;
t39 = t225 + 0.2e1 * (t248 / 0.4e1 + m(5) / 0.4e1 + m(6) / 0.4e1) * t196;
t265 = t39 * qJD(1);
t262 = -Ifges(6,5) * t327 + Ifges(6,6) * t198;
t163 = -t198 * mrSges(6,1) - mrSges(6,2) * t327;
t260 = t163 * qJD(5);
t249 = -t267 / 0.2e1;
t36 = (t198 ^ 2 + t327 ^ 2) * mrSges(6,3) + t331 + m(6) * (t154 * t198 - t155 * t327) + t205 * t248;
t220 = (t198 * t336 - t223 * t307) * mrSges(6,3) + t236 * t324 + (-t154 * t334 - t155 * t223 + t198 * t32 - t327 * t33) * t322 - t215 * t161 / 0.2e1 + t160 * t304 + t298;
t228 = -t290 / 0.2e1 - t288 / 0.2e1;
t222 = t328 * t324 - m(6) * t123 / 0.2e1 + t228;
t9 = (t280 / 0.2e1 + t282 / 0.2e1) * t195 + t220 + t222;
t234 = t9 * qJD(1) + t36 * qJD(2);
t40 = 0.2e1 * t333 + 0.2e1 * t337;
t233 = qJD(1) * t40 - qJD(2) * t163;
t10 = (t334 * t306 + t272 / 0.2e1) * mrSges(6,3) + t228 - t298;
t226 = t10 * qJD(1);
t165 = Ifges(6,2) * t198 - t190;
t167 = -Ifges(6,1) * t327 + t294;
t18 = t201 * t163 + (-t167 / 0.2e1 + t311) * t198 - (t310 + t165 / 0.2e1) * t327;
t219 = -(t65 / 0.4e1 + t76 / 0.4e1) * t327 + (-t77 / 0.4e1 + t64 / 0.4e1) * t198 - (t154 * t320 + t165 / 0.4e1 + t168 / 0.4e1) * t223 + (mrSges(6,3) * t312 - t166 / 0.4e1 + t167 / 0.4e1) * t334 + t124 * t163 / 0.2e1 + t154 * t103 / 0.2e1 + t104 * t312 + t195 * t262 / 0.4e1 + t201 * t75 / 0.2e1;
t221 = Ifges(6,3) * t308 - t34 * mrSges(6,1) / 0.2e1 + t35 * mrSges(6,2) / 0.2e1 + t229;
t3 = t219 + t221;
t224 = -t3 * qJD(1) - t18 * qJD(2);
t41 = -t289 / 0.2e1 + t287 / 0.2e1 + t337 + t333;
t38 = t248 * t308 + t225 + (m(5) + m(6)) * t309;
t11 = t254 * t334 + t272 * t320 + t228 + t298;
t8 = mrSges(5,2) * t249 - mrSges(5,1) * t269 / 0.2e1 + t220 - t222;
t7 = t217 + t218;
t2 = t219 - t221;
t13 = [-qJD(2) * t1 - qJD(3) * t6 + qJD(4) * t12 + qJD(5) * t4, -t291 + (t326 * mrSges(5,3) + (m(5) * t326 - t215 * t231 + t216 * t230) * t205 + (Ifges(5,5) * t215 - Ifges(6,5) * t198 + Ifges(5,6) * t216 - Ifges(6,6) * t327) * t309 + m(6) * (t123 * t201 + t154 * t34 + t155 * t35) + (m(5) * t212 - t275 * t321 - mrSges(4,1) + t202) * t328 - (t274 * t321 - mrSges(4,2)) * t235 + t143 * t310 - t140 * t311 + t250 * t297 + t34 * t283 - t35 * t284 + ((t293 - t295) * t304 + (-Ifges(5,1) * t216 + t296) * t305 + t252 * mrSges(4,3) - Ifges(4,5)) * t195 - t212 * t227 + mrSges(3,2) * t256 + t237 * t307 + t238 * t306 + (-Ifges(5,5) * t305 - Ifges(6,5) * t306 - Ifges(5,6) * t304 - Ifges(6,6) * t307 + Ifges(4,6)) * t196 + t155 * t239 + t154 * t240 + t201 * t242 + t123 * t164 - mrSges(3,1) * t258 + (Ifges(5,1) * t215 + t295) * t249 + (Ifges(5,2) * t216 + t296) * t269 / 0.2e1 - Ifges(3,6) * t301 + Ifges(3,5) * t303) * qJD(2) + t7 * qJD(3) + t8 * qJD(4) + t2 * qJD(5), -t277 + t7 * qJD(2) + m(6) * (t140 * t327 - t143 * t198) * qJD(3) + t38 * qJD(4) + t11 * qJD(5), qJD(2) * t8 + qJD(3) * t38 + qJD(5) * t41 + t273, t278 + t2 * qJD(2) + t11 * qJD(3) + t41 * qJD(4) + (-mrSges(6,1) * t33 - mrSges(6,2) * t32 + t263) * qJD(5); qJD(3) * t5 + qJD(4) * t9 + qJD(5) * t3 + t291, qJD(4) * t36 + qJD(5) * t18, t292, t234, (-mrSges(6,1) * t155 - mrSges(6,2) * t154 + t262) * qJD(5) - t224; -qJD(2) * t5 + qJD(4) * t39 - qJD(5) * t10 + t277, -t292, 0, t265, -t226 - t260; -qJD(2) * t9 - qJD(3) * t39 - qJD(5) * t40 - t273, -t234 + t260, -t265, 0, -t233; -qJD(2) * t3 + qJD(3) * t10 + qJD(4) * t40 - t278, -t163 * qJD(4) + t224, t226, t233, 0;];
Cq = t13;
