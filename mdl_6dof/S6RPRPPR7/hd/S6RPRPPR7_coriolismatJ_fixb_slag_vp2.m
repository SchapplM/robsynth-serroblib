% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPPR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:16
% EndTime: 2019-03-09 02:56:19
% DurationCPUTime: 2.62s
% Computational Cost: add. (6900->327), mult. (12448->430), div. (0->0), fcn. (12949->6), ass. (0->195)
t297 = sin(qJ(3));
t175 = t297 * pkin(3) + qJ(2);
t332 = m(5) * t175;
t184 = cos(pkin(9));
t187 = cos(qJ(3));
t271 = sin(pkin(9));
t158 = t184 * t297 + t187 * t271;
t160 = t184 * t187 - t271 * t297;
t226 = -qJ(5) * t160 + t175;
t103 = pkin(4) * t158 + t226;
t331 = m(6) * t103;
t186 = cos(qJ(6));
t280 = t186 * mrSges(7,1);
t185 = sin(qJ(6));
t283 = t185 * mrSges(7,2);
t209 = t280 / 0.2e1 - t283 / 0.2e1;
t247 = t209 * t158;
t188 = -pkin(1) - pkin(7);
t238 = t297 * t188;
t162 = -qJ(4) * t297 + t238;
t163 = (-qJ(4) + t188) * t187;
t325 = t184 * t162 + t271 * t163;
t330 = -t158 * pkin(5) + t325;
t262 = t158 * t186;
t108 = -t160 * mrSges(7,2) + mrSges(7,3) * t262;
t254 = t186 * t108;
t263 = t158 * t185;
t107 = t160 * mrSges(7,1) - mrSges(7,3) * t263;
t257 = t185 * t107;
t207 = t254 / 0.2e1 - t257 / 0.2e1;
t182 = t185 ^ 2;
t183 = t186 ^ 2;
t235 = t183 / 0.2e1 + t182 / 0.2e1;
t227 = mrSges(7,3) * t235;
t329 = -t158 * t227 + t207;
t301 = -t185 / 0.4e1;
t328 = mrSges(6,2) - mrSges(5,1);
t295 = mrSges(5,3) + mrSges(6,1);
t327 = -mrSges(6,3) + mrSges(5,2);
t164 = -t280 + t283;
t326 = t330 * t164;
t246 = t183 + t182;
t285 = t183 * mrSges(7,3);
t286 = t182 * mrSges(7,3);
t324 = -t285 / 0.2e1 - t286 / 0.2e1;
t323 = t158 * t184 - t271 * t160;
t174 = -pkin(3) * t184 - pkin(4);
t171 = pkin(3) * t271 + qJ(5);
t261 = t160 * t171;
t322 = -t158 * t174 - t261;
t181 = t187 * pkin(3);
t229 = qJ(5) * t158 + t181;
t306 = pkin(4) + pkin(8);
t78 = t160 * t306 + t229;
t43 = -t185 * t78 + t186 * t330;
t276 = t186 * t43;
t44 = t185 * t330 + t186 * t78;
t282 = t185 * t44;
t321 = -t282 - t276;
t292 = Ifges(7,6) * t186;
t293 = Ifges(7,5) * t185;
t319 = -Ifges(6,6) - Ifges(5,4) + t293 / 0.2e1 + t292 / 0.2e1;
t73 = t158 * t306 + t226;
t113 = t162 * t271 - t184 * t163;
t83 = t160 * pkin(5) + t113;
t41 = -t185 * t73 + t186 * t83;
t42 = t185 * t83 + t186 * t73;
t219 = t41 * t185 - t42 * t186;
t318 = -m(7) * t219 + t254 - t257 + t331;
t125 = t158 ^ 2;
t317 = t160 ^ 2;
t316 = 0.2e1 * t158;
t315 = 0.2e1 * t160;
t313 = m(6) / 0.2e1;
t312 = m(6) / 0.4e1;
t311 = m(7) / 0.2e1;
t310 = m(7) / 0.4e1;
t309 = m(5) * pkin(3);
t308 = -mrSges(7,1) / 0.2e1;
t307 = mrSges(7,2) / 0.2e1;
t304 = -t158 / 0.2e1;
t303 = -t160 / 0.2e1;
t302 = -t185 / 0.2e1;
t300 = t185 / 0.2e1;
t298 = t186 / 0.2e1;
t294 = Ifges(7,4) * t185;
t220 = t185 * t42 + t186 * t41;
t255 = t186 * t107;
t256 = t185 * t108;
t290 = t158 * t330;
t7 = (t160 * t295 + t255 + t256) * t160 + (-t164 + t295) * t125 + m(7) * (t160 * t220 - t290) + (m(6) + m(5)) * (t113 * t160 - t158 * t325);
t291 = qJD(1) * t7;
t289 = t160 * mrSges(6,3);
t288 = t160 * Ifges(7,5);
t287 = t160 * Ifges(7,6);
t284 = t185 * mrSges(7,1);
t149 = Ifges(7,4) * t262;
t76 = Ifges(7,1) * t263 + t149 + t288;
t281 = t185 * t76;
t279 = t186 * mrSges(7,2);
t278 = t186 * Ifges(7,4);
t277 = t186 * Ifges(7,2);
t224 = t277 + t294;
t74 = t158 * t224 + t287;
t275 = t186 * t74;
t274 = t187 * mrSges(4,2);
t165 = t279 + t284;
t169 = Ifges(7,1) * t186 - t294;
t223 = Ifges(7,5) * t186 - Ifges(7,6) * t185;
t4 = -t165 * t290 - t41 * t108 + t42 * t107 + (t74 * t300 + t223 * t303 + (t169 * t302 + t277 * t300) * t158 + t220 * mrSges(7,3) - (t149 + t76) * t186 / 0.2e1) * t158;
t273 = t4 * qJD(1);
t151 = t158 * mrSges(6,3);
t153 = t158 * mrSges(5,2);
t170 = -pkin(8) + t174;
t230 = t246 * t170;
t244 = t309 / 0.2e1;
t191 = t165 * t304 + (t174 * t313 - t184 * t244 + t230 * t311 + t324) * t160 + (-t271 * t244 + (-t313 - t311) * t171) * t158;
t104 = pkin(4) * t160 + t229;
t259 = t160 * t186;
t212 = t158 * mrSges(7,2) + mrSges(7,3) * t259;
t260 = t160 * t185;
t213 = -t158 * mrSges(7,1) - mrSges(7,3) * t260;
t193 = t104 * t313 + (-t185 * t43 + t186 * t44) * t311 + t213 * t302 + t212 * t298 + t187 * t244;
t8 = t160 * t328 - t151 + t153 + t191 - t193;
t272 = t8 * qJD(1);
t152 = t158 * mrSges(6,2);
t112 = -t152 - t289;
t15 = (-t112 - t318) * t160;
t270 = qJD(1) * t15;
t154 = t158 * mrSges(5,1);
t195 = mrSges(4,1) * t297 - t152 + t154 + t274;
t17 = mrSges(3,3) + t327 * t160 + (m(4) + m(3)) * qJ(2) + t332 + t195 + t318;
t269 = qJD(1) * t17;
t205 = t125 * t165;
t210 = t284 / 0.2e1 + t279 / 0.2e1;
t11 = -t205 / 0.2e1 + t329 * t160 - t210;
t268 = t11 * qJD(1);
t126 = t158 * t160;
t267 = t158 * t164;
t199 = t210 * t160;
t18 = t199 - t329;
t258 = t18 * qJD(1);
t198 = t209 * t160;
t208 = t256 / 0.2e1 + t255 / 0.2e1;
t20 = t198 + t208;
t253 = t20 * qJD(1);
t194 = (-t246 * t317 - t125) * t311 + 0.2e1 * (t312 + m(5) / 0.4e1) * (-t125 - t317);
t234 = m(7) * t246;
t203 = -m(5) / 0.2e1 - m(6) / 0.2e1 - t234 / 0.2e1;
t25 = t194 + t203;
t252 = t25 * qJD(1);
t35 = (m(7) * t235 + t313) * t315;
t251 = t35 * qJD(1);
t248 = t262 * t308 + t263 * t307;
t245 = qJD(6) * t165;
t239 = -t326 / 0.2e1;
t237 = t260 / 0.2e1;
t236 = t259 / 0.2e1;
t231 = t246 * t160;
t225 = t185 * Ifges(7,1) + t278;
t222 = -t292 - t293;
t75 = -Ifges(7,6) * t158 + t160 * t224;
t77 = -Ifges(7,5) * t158 + t160 * t225;
t1 = -qJ(2) * t297 * mrSges(4,2) + t103 * t151 + t43 * t107 + t44 * t108 - t175 * t153 + (qJ(2) * mrSges(4,1) + pkin(3) * t154 + (-Ifges(4,1) + Ifges(4,2)) * t297) * t187 + t181 * t332 + m(7) * (-t330 * t83 + t41 * t43 + t42 * t44) + (-t41 * mrSges(7,1) + t42 * mrSges(7,2) - t319 * t158 - t164 * t83 + t75 * t298 + t77 * t300) * t158 + (-t187 ^ 2 + t297 ^ 2) * Ifges(4,4) + (mrSges(5,2) * t181 + t275 / 0.2e1 + t281 / 0.2e1 - t103 * mrSges(6,2) + t175 * mrSges(5,1) + t326 + t319 * t160 - t219 * mrSges(7,3) + (Ifges(6,3) - Ifges(5,1) - Ifges(6,2) + Ifges(5,2) - Ifges(7,3)) * t158) * t160 + (t112 + t331) * t104;
t204 = t160 * t164;
t3 = -m(7) * (t321 + t330) * t315 / 0.4e1 + (-t204 / 0.2e1 - m(7) * (t220 - t83) / 0.2e1 - t208) * t158;
t221 = t1 * qJD(1) - t3 * qJD(2);
t30 = m(7) * (-t158 * t231 + t126);
t217 = -t3 * qJD(1) + t30 * qJD(2);
t214 = t307 * t44 + t308 * t43;
t206 = t158 * t223;
t202 = t185 * t212;
t201 = t186 * t213;
t167 = -Ifges(7,2) * t185 + t278;
t49 = -t171 * t164 + (-t225 / 0.2e1 - t167 / 0.2e1) * t186 + (-t169 / 0.2e1 + t224 / 0.2e1) * t185;
t52 = t267 / 0.2e1 + t247;
t192 = -t186 * t224 / 0.4e1 + t169 * t298 + t171 * t165 / 0.2e1 + t182 * Ifges(7,2) / 0.4e1 - t170 * t227 + (t167 + t225) * t301;
t6 = t239 + (t170 * t108 / 0.2e1 - 0.3e1 / 0.4e1 * t287 - t74 / 0.4e1) * t186 + (-t170 * t107 / 0.2e1 - 0.3e1 / 0.4e1 * t288 - t76 / 0.4e1 - t149 / 0.4e1) * t185 + (Ifges(7,3) / 0.2e1 + t192) * t158 + t214;
t197 = t6 * qJD(1) - t52 * qJD(2) + t49 * qJD(3);
t111 = mrSges(6,3) + 0.4e1 * (t312 + t310) * t171 + t165;
t14 = 0.2e1 * (t330 / 0.4e1 - t282 / 0.4e1 - t276 / 0.4e1) * m(7) + t247 + t248;
t57 = (0.1e1 / 0.4e1 - t183 / 0.4e1 - t182 / 0.4e1) * m(7) * t316;
t196 = qJD(1) * t14 + qJD(2) * t57 + qJD(3) * t111;
t53 = -t267 / 0.2e1 + t247;
t47 = t158 * t234 / 0.2e1 + (t313 + t310) * t316;
t36 = m(6) * t303 - t231 * t311 + (t312 + t234 / 0.4e1) * t315;
t24 = t194 - t203;
t21 = t198 - t208;
t19 = t158 * t324 + t199 + t207;
t13 = -t158 * mrSges(6,1) + t201 / 0.2e1 + t202 / 0.2e1 + t248 + 0.2e1 * t325 * t313 + (-t321 + t330) * t311;
t12 = t205 / 0.2e1 + t254 * t303 + t107 * t237 - t210 + t246 * mrSges(7,3) * t126 / 0.2e1;
t10 = t191 + t193;
t5 = t160 * t222 / 0.4e1 + t239 - t281 / 0.4e1 - t275 / 0.4e1 + t149 * t301 + Ifges(7,6) * t236 + Ifges(7,5) * t237 + Ifges(7,3) * t304 + t207 * t170 + t192 * t158 - t214;
t2 = t3 * qJD(3);
t9 = [qJD(2) * t17 + qJD(3) * t1 + qJD(4) * t7 + qJD(5) * t15 - qJD(6) * t4, qJD(4) * t24 + qJD(6) * t12 - t2 + t269 (-t206 / 0.2e1 + t77 * t298 + t75 * t302 + t167 * t236 + t169 * t237 - mrSges(4,1) * t238 - Ifges(4,5) * t297 - t188 * t274 - t83 * t165 - Ifges(4,6) * t187 + t323 * mrSges(5,3) * pkin(3) + (-m(7) * t83 + t204) * t171 + (-m(7) * t321 + t201 + t202) * t170 + (Ifges(6,5) - Ifges(5,6)) * t160 + (Ifges(6,4) - Ifges(5,5)) * t158 - (m(6) * t171 + t271 * t309 - t327) * t113 + (m(6) * t174 - t184 * t309 + t328) * t325 + t321 * mrSges(7,3) + t322 * mrSges(6,1)) * qJD(3) + t10 * qJD(4) + t13 * qJD(5) + t5 * qJD(6) + t221, qJD(2) * t24 + qJD(3) * t10 + qJD(5) * t36 + qJD(6) * t21 + t291, qJD(3) * t13 + qJD(4) * t36 + qJD(6) * t19 + t270, -t273 + t12 * qJD(2) + t5 * qJD(3) + t21 * qJD(4) + t19 * qJD(5) + (-mrSges(7,1) * t42 - mrSges(7,2) * t41 + t206) * qJD(6); qJD(4) * t25 - qJD(6) * t11 - t2 - t269, t30 * qJD(3) (-t323 * t309 - m(6) * t322 + m(7) * t261 + t289 - t195 + (t165 - mrSges(5,2)) * t160 + (m(7) * t230 - t285 - t286) * t158) * qJD(3) + t47 * qJD(5) + t53 * qJD(6) + t217, t252, t47 * qJD(3), t53 * qJD(3) + t160 * t245 - t268; qJD(4) * t8 + qJD(5) * t14 + qJD(6) * t6 - t221, qJD(5) * t57 - qJD(6) * t52 - t217, qJD(5) * t111 + qJD(6) * t49, t272, t196 (-t165 * t170 + t222) * qJD(6) + t197; -qJD(2) * t25 - qJD(3) * t8 - qJD(5) * t35 - qJD(6) * t20 - t291, -t252, -t272, 0, -t251, qJD(6) * t164 - t253; -qJD(3) * t14 + qJD(4) * t35 - qJD(6) * t18 - t270, -t57 * qJD(3), -t196, t251, 0, -t245 - t258; qJD(2) * t11 - qJD(3) * t6 + qJD(4) * t20 + qJD(5) * t18 + t273, qJD(3) * t52 + t268, -t197, t253, t258, 0;];
Cq  = t9;
