% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPPRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:49
% EndTime: 2019-03-08 19:17:53
% DurationCPUTime: 2.18s
% Computational Cost: add. (4365->299), mult. (11069->451), div. (0->0), fcn. (10906->10), ass. (0->197)
t187 = sin(pkin(11));
t191 = sin(qJ(2));
t271 = sin(pkin(6));
t272 = cos(pkin(11));
t235 = t272 * t271;
t295 = cos(qJ(2));
t236 = t295 * t271;
t127 = t187 * t236 + t191 * t235;
t323 = m(6) * t127;
t189 = sin(qJ(6));
t183 = t189 ^ 2;
t192 = cos(qJ(6));
t185 = t192 ^ 2;
t255 = t183 + t185;
t317 = -t255 * mrSges(7,3) + mrSges(6,2);
t278 = t192 * mrSges(7,1);
t282 = t189 * mrSges(7,2);
t158 = -t278 + t282;
t273 = -mrSges(6,1) + t158;
t322 = -m(7) * pkin(5) + t273;
t320 = 0.1e1 - t255;
t319 = t187 * pkin(2);
t318 = t272 * pkin(2);
t190 = sin(qJ(5));
t193 = cos(qJ(5));
t316 = t190 * mrSges(6,1) + t193 * mrSges(6,2) + mrSges(5,3);
t182 = Ifges(7,4) * t192;
t315 = -Ifges(7,2) * t189 + t182;
t161 = Ifges(7,1) * t189 + t182;
t277 = t192 * mrSges(7,2);
t283 = t189 * mrSges(7,1);
t159 = t277 + t283;
t217 = t277 / 0.2e1 + t283 / 0.2e1;
t314 = t217 + t159 / 0.2e1;
t186 = t193 ^ 2;
t103 = t186 * t127;
t184 = t190 ^ 2;
t266 = t184 * t127;
t313 = -t266 - t103;
t294 = pkin(5) * t193;
t165 = t190 * pkin(9) + t294;
t238 = -pkin(3) - t318;
t174 = -pkin(8) + t238;
t262 = t189 * t193;
t95 = t192 * t165 - t174 * t262;
t258 = t192 * t193;
t96 = t189 * t165 + t174 * t258;
t228 = -t189 * t95 + t192 * t96;
t141 = t190 * t159;
t249 = mrSges(7,3) * t262;
t153 = -t190 * mrSges(7,2) - t249;
t248 = mrSges(7,3) * t258;
t155 = t190 * mrSges(7,1) - t248;
t298 = -t192 / 0.2e1;
t300 = t189 / 0.2e1;
t205 = t153 * t298 + t155 * t300 - t141 / 0.2e1;
t312 = -m(7) / 0.2e1;
t311 = m(7) / 0.2e1;
t244 = t191 * t271;
t126 = t187 * t244 - t295 * t235;
t188 = cos(pkin(6));
t91 = t126 * t190 + t188 * t193;
t59 = t127 * t189 + t192 * t91;
t310 = -t59 / 0.2e1;
t254 = t184 - t186;
t84 = t320 * t254;
t309 = -t84 / 0.2e1;
t261 = t190 * t192;
t66 = -t126 * t189 + t127 * t261;
t275 = t192 * t66;
t263 = t189 * t190;
t65 = -t126 * t192 - t127 * t263;
t280 = t189 * t65;
t230 = t275 - t280;
t270 = t127 * t190;
t308 = m(7) * (t230 - t270) * t193;
t140 = -mrSges(7,1) * t258 + mrSges(7,2) * t262;
t307 = t140 / 0.2e1;
t305 = -t153 / 0.2e1;
t304 = -t159 / 0.2e1;
t302 = -t189 / 0.2e1;
t299 = -t190 / 0.2e1;
t297 = t192 / 0.2e1;
t296 = t193 / 0.2e1;
t293 = t190 * pkin(5);
t276 = t192 * t59;
t58 = t127 * t192 - t189 * t91;
t281 = t189 * t58;
t222 = -t276 + t91 + t281;
t90 = -t126 * t193 + t188 * t190;
t239 = t320 * t90;
t19 = (t190 * t239 - t222 * t193) * t311;
t292 = t19 * qJD(5);
t142 = t193 * t159;
t152 = -mrSges(7,2) * t193 + mrSges(7,3) * t263;
t260 = t192 * t152;
t154 = mrSges(7,1) * t193 + mrSges(7,3) * t261;
t264 = t189 * t154;
t204 = t142 / 0.2e1 - t264 / 0.2e1 + t260 / 0.2e1;
t214 = m(7) * t228;
t176 = qJ(4) + t319;
t139 = -pkin(9) * t193 + t176 + t293;
t82 = t139 * t192 - t174 * t263;
t83 = t139 * t189 + t174 * t261;
t229 = t189 * t82 - t192 * t83;
t22 = t254 * t174 * t311 + (t214 / 0.2e1 + t204) * t193 + (t229 * t311 + t205) * t190;
t259 = t192 * t155;
t265 = t189 * t153;
t213 = -t265 / 0.2e1 - t259 / 0.2e1;
t237 = mrSges(7,3) * (-t185 / 0.2e1 - t183 / 0.2e1);
t33 = t140 * t299 + t186 * t237 + t213 * t193;
t291 = t22 * qJD(5) + t33 * qJD(6);
t290 = m(7) * qJD(5);
t288 = Ifges(7,4) * t189;
t287 = Ifges(7,5) * t190;
t181 = Ifges(7,5) * t192;
t285 = Ifges(7,6) * t189;
t284 = Ifges(7,6) * t190;
t269 = t127 * t193;
t268 = t174 * t193;
t267 = t176 * t126;
t196 = (t193 * t237 + t213) * t190 + t140 * t296;
t218 = t282 / 0.2e1 - t278 / 0.2e1;
t29 = t196 + t218;
t257 = t29 * qJD(2);
t256 = t33 * qJD(2);
t253 = qJD(5) * t190;
t252 = -t308 / 0.2e1;
t251 = t308 / 0.2e1;
t240 = t255 * t193;
t100 = (-t193 + t240) * t190;
t250 = t100 * t290;
t245 = pkin(9) * t255;
t241 = qJD(5) * t273;
t233 = Ifges(7,1) * t192 - t288;
t160 = Ifges(7,2) * t192 + t288;
t231 = Ifges(7,5) * t189 + Ifges(7,6) * t192;
t18 = (t222 * t190 + t193 * t239) * t311;
t227 = -t18 * qJD(1) - t22 * qJD(2);
t24 = (t229 * t312 - t205) * t193 + ((t228 - 0.2e1 * t268) * t311 + t204) * t190;
t226 = -t19 * qJD(1) - t24 * qJD(2);
t197 = -m(6) * t313 / 0.2e1 + (t230 * t190 + t103) * t311;
t199 = -t323 / 0.2e1 + (t189 * t59 + t192 * t58) * t312;
t20 = t197 + t199;
t34 = t265 + t259 + m(7) * (t189 * t83 + t192 * t82) + 0.4e1 * (m(6) / 0.4e1 + m(5) / 0.4e1) * t176 + t316;
t225 = -qJD(1) * t20 + qJD(2) * t34;
t224 = t65 * mrSges(7,1) / 0.2e1 - t66 * mrSges(7,2) / 0.2e1;
t223 = -t95 * mrSges(7,1) / 0.2e1 + t96 * mrSges(7,2) / 0.2e1;
t221 = t100 * qJD(3) + qJD(4) * t309;
t220 = qJD(3) * t309 - t100 * qJD(4);
t219 = pkin(5) * t307 + t190 * t181 / 0.4e1;
t129 = t193 * t315 + t284;
t144 = t193 * t161;
t216 = -t129 / 0.4e1 - t144 / 0.4e1 + pkin(9) * t305;
t131 = t233 * t193 + t287;
t143 = t193 * t160;
t215 = -pkin(9) * t155 / 0.2e1 + t131 / 0.4e1 - t143 / 0.4e1;
t212 = t160 * t300 + t161 * t298;
t6 = m(7) * (-t90 * t269 + t58 * t65 + t59 * t66) + (t190 * t91 - t193 * t90 - t126) * t323;
t211 = -t6 * qJD(1) + qJD(3) * t252;
t210 = t181 / 0.2e1 - t285 / 0.2e1 - Ifges(6,4);
t209 = t174 * t190 + t229;
t208 = t304 + t217;
t10 = t82 * t153 - t83 * t155 + (t229 * mrSges(7,3) + t129 * t298 + t174 * t140 - t144 * t297 + t231 * t299 + (t131 - t143) * t302) * t193;
t202 = t58 * t305 + t59 * t155 / 0.2e1 + t90 * t307;
t3 = (t276 / 0.2e1 - t281 / 0.2e1) * t193 * mrSges(7,3) + t202 + t224;
t207 = -t3 * qJD(1) + t10 * qJD(2) + t33 * qJD(3);
t12 = m(7) * t222 * t90;
t206 = t12 * qJD(1) + t18 * qJD(3) + t19 * qJD(4);
t203 = -t91 * t268 + t95 * t58 + t96 * t59;
t201 = -t58 * t154 / 0.2e1 + t152 * t310 - t91 * t142 / 0.2e1;
t195 = (t275 / 0.2e1 - t280 / 0.2e1) * mrSges(7,3) + (pkin(5) * t269 + t230 * pkin(9)) * t311 - t158 * t269 / 0.2e1;
t1 = t203 * t312 + (t209 * t312 - t205) * t90 + t195 + t201;
t128 = Ifges(7,6) * t193 - t190 * t315;
t130 = Ifges(7,5) * t193 - t233 * t190;
t7 = m(7) * (t82 * t95 + t83 * t96) + t96 * t153 + t83 * t152 + t95 * t155 + t82 * t154 + (t176 * mrSges(6,1) + t128 * t302 + t130 * t297 + t174 * t141 + t210 * t193) * t193 + (-t176 * mrSges(6,2) + t174 * t142 + t131 * t298 + t129 * t300 - t210 * t190 + (-m(7) * t174 ^ 2 - Ifges(6,1) + Ifges(6,2) + Ifges(7,3)) * t193) * t190;
t200 = -t1 * qJD(1) + t7 * qJD(2) + t22 * qJD(3) + t24 * qJD(4);
t14 = t208 * t90;
t53 = pkin(5) * t159 + t233 * t302 + t298 * t315 + t212;
t194 = pkin(9) * t237 + t174 * t304 - (t161 + t315) * t189 / 0.4e1 + (t233 / 0.4e1 - t160 / 0.4e1) * t192;
t8 = (t287 / 0.2e1 + t215) * t192 + (-0.3e1 / 0.4e1 * t284 + t216) * t189 + (-Ifges(7,3) / 0.2e1 + t194) * t193 + t219 + t223;
t86 = t208 * t190;
t87 = t208 * t193;
t198 = t14 * qJD(1) - t8 * qJD(2) + t86 * qJD(3) - t87 * qJD(4) + t53 * qJD(5);
t92 = t174 * t103;
t89 = (-t217 + t304) * t193;
t88 = t314 * t190;
t76 = t84 * t290 / 0.2e1;
t28 = t196 - t218;
t23 = t24 * qJD(5);
t16 = m(5) * t127 + t197 - t199;
t15 = t314 * t90;
t9 = -Ifges(7,5) * t261 / 0.2e1 + Ifges(7,6) * t263 / 0.2e1 + Ifges(7,3) * t296 + (-t284 / 0.4e1 + t216) * t189 + t215 * t192 + t194 * t193 + t219 - t223;
t5 = qJD(2) * t251 + t18 * qJD(5);
t4 = t248 * t310 + t58 * t249 / 0.2e1 - t202 + t224;
t2 = t203 * t311 + t127 * (t193 * mrSges(6,1) - t190 * mrSges(6,2)) / 0.2e1 - mrSges(6,2) * t270 / 0.2e1 + mrSges(6,1) * t269 / 0.2e1 + t195 - t201 + (t209 * t311 + t205) * t90;
t11 = [t6 * qJD(2) + t12 * qJD(5) (-mrSges(3,1) * t244 - mrSges(3,2) * t236 + t66 * t153 + t65 * t155 + m(7) * (t65 * t82 + t66 * t83 + t92) - t142 * t269 + m(6) * (t174 * t266 - t267 + t92) - m(5) * t267 + (-m(4) * t318 + m(5) * t238 - mrSges(4,1) + mrSges(5,2)) * t127 + (-m(4) * t319 + mrSges(4,2) - t316) * t126 + t313 * mrSges(6,3)) * qJD(2) + t16 * qJD(4) + t2 * qJD(5) + t4 * qJD(6) - t211, t5, qJD(2) * t16 + t292, t2 * qJD(2) + (t322 * t91 + (-m(7) * t245 + t317) * t90) * qJD(5) + t15 * qJD(6) + t206, t4 * qJD(2) + t15 * qJD(5) + (-mrSges(7,1) * t59 - mrSges(7,2) * t58) * qJD(6); -t20 * qJD(4) - t1 * qJD(5) - t3 * qJD(6) + t211, qJD(4) * t34 + qJD(5) * t7 + qJD(6) * t10, qJD(1) * t252 + t291, qJD(6) * t28 + t225 + t23, t9 * qJD(6) + (t322 * t174 - Ifges(6,5) + t212) * t253 + t200 + (-mrSges(6,2) * t268 - Ifges(6,6) * t193 + pkin(5) * t141 + t128 * t297 + t130 * t300 + t231 * t296 + (t214 + t260 - t264) * pkin(9) + t228 * mrSges(7,3)) * qJD(5), t28 * qJD(4) + t9 * qJD(5) + (-t83 * mrSges(7,1) - t82 * mrSges(7,2) - t193 * t231) * qJD(6) + t207; t5, qJD(1) * t251 + t291, -t250, t76, t88 * qJD(6) + t193 * t241 + t317 * t253 + ((-t190 * t245 - t294) * qJD(5) - t221) * m(7) - t227, t88 * qJD(5) + t140 * qJD(6) + t256; qJD(2) * t20 + t292, qJD(6) * t29 - t225 + t23, t76, t250, t89 * qJD(6) + t190 * t241 - t317 * qJD(5) * t193 + ((pkin(9) * t240 - t293) * qJD(5) - t220) * m(7) - t226, qJD(6) * t158 * t190 + t89 * qJD(5) + t257; qJD(2) * t1 - qJD(6) * t14 - t206, qJD(6) * t8 - t200, t221 * m(7) - t86 * qJD(6) + t227, t220 * m(7) + t87 * qJD(6) + t226, -t53 * qJD(6) (pkin(9) * t158 + t181 - t285) * qJD(6) - t198; t3 * qJD(2) + t14 * qJD(5), -qJD(4) * t29 - qJD(5) * t8 - t207, t86 * qJD(5) - t256, -t87 * qJD(5) - t257, t198, 0;];
Cq  = t11;
