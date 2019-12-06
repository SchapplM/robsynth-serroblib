% Calculate vector of inverse dynamics joint torques for
% S5RRRRP1
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
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:46
% EndTime: 2019-12-05 18:45:11
% DurationCPUTime: 11.39s
% Computational Cost: add. (6786->478), mult. (16064->614), div. (0->0), fcn. (11290->12), ass. (0->214)
t227 = sin(qJ(3));
t228 = sin(qJ(2));
t231 = cos(qJ(3));
t232 = cos(qJ(2));
t178 = -t227 * t228 + t231 * t232;
t165 = t178 * qJD(1);
t179 = t227 * t232 + t228 * t231;
t166 = t179 * qJD(1);
t226 = sin(qJ(4));
t230 = cos(qJ(4));
t126 = t165 * t226 + t166 * t230;
t221 = t232 * pkin(2);
t209 = t221 + pkin(1);
t194 = t209 * qJD(1);
t145 = -pkin(3) * t165 - t194;
t222 = qJDD(2) + qJDD(3);
t215 = qJDD(4) + t222;
t265 = qJD(1) * qJD(2);
t188 = qJDD(1) * t232 - t228 * t265;
t189 = qJDD(1) * t228 + t232 * t265;
t246 = t178 * qJD(3);
t103 = qJD(1) * t246 + t188 * t227 + t189 * t231;
t247 = t179 * qJD(3);
t104 = -qJD(1) * t247 + t188 * t231 - t189 * t227;
t253 = t230 * t165 - t166 * t226;
t33 = qJD(4) * t253 + t103 * t230 + t104 * t226;
t234 = -pkin(7) - pkin(6);
t196 = t234 * t232;
t185 = qJD(1) * t196;
t170 = t231 * t185;
t195 = t234 * t228;
t184 = qJD(1) * t195;
t173 = qJD(2) * pkin(2) + t184;
t133 = t173 * t227 - t170;
t177 = t189 * pkin(6);
t143 = qJDD(2) * pkin(2) - pkin(7) * t189 - t177;
t176 = t188 * pkin(6);
t144 = pkin(7) * t188 + t176;
t58 = -qJD(3) * t133 + t231 * t143 - t144 * t227;
t16 = pkin(3) * t222 - pkin(8) * t103 + t58;
t270 = qJD(3) * t231;
t271 = qJD(3) * t227;
t57 = t227 * t143 + t231 * t144 + t173 * t270 + t185 * t271;
t25 = pkin(8) * t104 + t57;
t223 = qJD(2) + qJD(3);
t167 = t227 * t185;
t132 = t231 * t173 + t167;
t160 = t166 * pkin(8);
t99 = t132 - t160;
t88 = pkin(3) * t223 + t99;
t302 = pkin(8) * t165;
t100 = t133 + t302;
t91 = t230 * t100;
t49 = t226 * t88 + t91;
t6 = -qJD(4) * t49 + t230 * t16 - t226 * t25;
t2 = pkin(4) * t215 - qJ(5) * t33 - qJD(5) * t126 + t6;
t216 = qJD(4) + t223;
t34 = -qJD(4) * t126 - t103 * t226 + t104 * t230;
t268 = qJD(4) * t230;
t269 = qJD(4) * t226;
t5 = -t100 * t269 + t226 * t16 + t230 * t25 + t88 * t268;
t3 = qJ(5) * t34 + qJD(5) * t253 + t5;
t315 = t126 / 0.2e1;
t337 = Ifges(5,6) + Ifges(6,6);
t338 = Ifges(5,2) + Ifges(6,2);
t339 = Ifges(5,5) + Ifges(6,5);
t341 = Ifges(5,1) + Ifges(6,1);
t340 = Ifges(5,4) + Ifges(6,4);
t357 = t340 * t253;
t353 = t341 * t126 + t339 * t216 + t357;
t356 = t340 * t126;
t354 = t337 * t216 + t338 * t253 + t356;
t76 = -pkin(4) * t253 + qJD(5) + t145;
t359 = t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) + t337 * t34 + t339 * t33 + (Ifges(5,3) + Ifges(6,3)) * t215 - (-t337 * t126 + t339 * t253) * t216 / 0.2e1 - t145 * (mrSges(5,1) * t126 + mrSges(5,2) * t253) - t76 * (mrSges(6,1) * t126 + mrSges(6,2) * t253) + t354 * t315 - (-t338 * t126 + t353 + t357) * t253 / 0.2e1 - (t341 * t253 - t356) * t126 / 0.2e1;
t295 = mrSges(5,1) + mrSges(6,1);
t193 = -mrSges(3,1) * t232 + mrSges(3,2) * t228;
t225 = qJ(2) + qJ(3);
t217 = sin(t225);
t218 = cos(t225);
t220 = qJ(4) + t225;
t205 = sin(t220);
t206 = cos(t220);
t294 = mrSges(5,2) + mrSges(6,2);
t250 = t205 * t294 - t295 * t206;
t241 = -t218 * mrSges(4,1) + mrSges(4,2) * t217 + t250;
t358 = t193 + t241;
t138 = -t184 * t227 + t170;
t105 = t138 - t302;
t139 = t231 * t184 + t167;
t106 = -t160 + t139;
t208 = pkin(2) * t231 + pkin(3);
t275 = t227 * t230;
t350 = -t230 * t105 + t106 * t226 - t208 * t269 + (-t227 * t268 + (-t226 * t231 - t275) * qJD(3)) * pkin(2);
t276 = t226 * t227;
t349 = -t226 * t105 - t230 * t106 + t208 * t268 + (-t227 * t269 + (t230 * t231 - t276) * qJD(3)) * pkin(2);
t281 = qJ(5) * t253;
t352 = t281 + t350;
t118 = qJ(5) * t126;
t351 = t118 + t349;
t89 = t226 * t100;
t48 = t230 * t88 - t89;
t17 = t48 - t118;
t13 = pkin(4) * t216 + t17;
t347 = t48 * mrSges(5,3) + t13 * mrSges(6,3);
t229 = sin(qJ(1));
t233 = cos(qJ(1));
t346 = g(1) * t233 + g(2) * t229;
t307 = t228 / 0.2e1;
t280 = qJDD(1) * pkin(1);
t146 = t231 * t195 + t196 * t227;
t121 = -pkin(8) * t179 + t146;
t147 = t227 * t195 - t231 * t196;
t122 = pkin(8) * t178 + t147;
t71 = t226 * t121 + t230 * t122;
t330 = t176 * t232 + t177 * t228;
t18 = t49 + t281;
t328 = -t49 * mrSges(5,3) - t18 * mrSges(6,3);
t200 = pkin(4) * t206;
t204 = pkin(3) * t218;
t274 = t204 + t221;
t262 = t200 + t274;
t327 = mrSges(2,1) + m(6) * (pkin(1) + t262) + m(5) * (pkin(1) + t274) + m(4) * t209 + m(3) * pkin(1) - t358;
t224 = -pkin(8) + t234;
t326 = mrSges(2,2) + m(6) * (-qJ(5) + t224) - mrSges(6,3) + m(5) * t224 - mrSges(5,3) + m(4) * t234 - mrSges(4,3) - m(3) * pkin(6) - mrSges(3,3);
t325 = m(4) * pkin(2);
t311 = t166 / 0.2e1;
t305 = pkin(2) * t228;
t304 = pkin(3) * t166;
t303 = pkin(3) * t217;
t22 = -mrSges(6,2) * t215 + mrSges(6,3) * t34;
t23 = -mrSges(5,2) * t215 + mrSges(5,3) * t34;
t293 = t22 + t23;
t51 = t230 * t99 - t89;
t292 = mrSges(3,2) * t232;
t290 = mrSges(4,3) * t166;
t289 = Ifges(3,4) * t228;
t288 = Ifges(3,4) * t232;
t285 = t132 * mrSges(4,3);
t284 = t133 * mrSges(4,3);
t283 = t166 * Ifges(4,4);
t282 = t232 * Ifges(3,2);
t273 = qJD(2) * t228;
t272 = qJD(2) * t232;
t267 = t228 * qJD(1);
t266 = t232 * qJD(1);
t213 = pkin(2) * t273;
t261 = qJD(2) * t234;
t257 = -t34 * mrSges(6,1) + t33 * mrSges(6,2);
t256 = t294 * t206;
t141 = -qJD(2) * t179 - t247;
t127 = -pkin(3) * t141 + t213;
t50 = -t226 * t99 - t91;
t70 = t230 * t121 - t122 * t226;
t151 = -pkin(3) * t178 - t209;
t161 = -pkin(2) * t276 + t230 * t208;
t85 = pkin(4) * t126 + t304;
t174 = -pkin(4) * t205 - t303;
t252 = t282 + t289;
t251 = Ifges(3,5) * t232 - Ifges(3,6) * t228;
t158 = -pkin(2) * t188 - t280;
t136 = t178 * t230 - t179 * t226;
t137 = t178 * t226 + t179 * t230;
t249 = pkin(1) * (mrSges(3,1) * t228 + t292);
t248 = t228 * (Ifges(3,1) * t232 - t289);
t186 = t228 * t261;
t187 = t232 * t261;
t81 = t231 * t186 + t227 * t187 + t195 * t270 + t196 * t271;
t62 = pkin(8) * t141 + t81;
t140 = qJD(2) * t178 + t246;
t82 = -qJD(3) * t147 - t186 * t227 + t231 * t187;
t63 = -pkin(8) * t140 + t82;
t9 = t121 * t268 - t122 * t269 + t226 * t63 + t230 * t62;
t77 = -pkin(3) * t104 + t158;
t10 = -qJD(4) * t71 - t226 * t62 + t230 * t63;
t238 = mrSges(4,1) * t217 + mrSges(4,2) * t218 + t205 * t295 + t256;
t116 = t165 * Ifges(4,2) + t223 * Ifges(4,6) + t283;
t159 = Ifges(4,4) * t165;
t117 = t166 * Ifges(4,1) + t223 * Ifges(4,5) + t159;
t235 = -(-Ifges(4,2) * t166 + t117 + t159) * t165 / 0.2e1 - t328 * t126 + t194 * (mrSges(4,1) * t166 + mrSges(4,2) * t165) + t347 * t253 + t116 * t311 + t165 * t285 - t166 * (Ifges(4,1) * t165 - t283) / 0.2e1 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + Ifges(4,5) * t103 + Ifges(4,6) * t104 + Ifges(4,3) * t222 - t223 * (Ifges(4,5) * t165 - Ifges(4,6) * t166) / 0.2e1 + t359;
t212 = pkin(2) * t267;
t211 = Ifges(3,4) * t266;
t192 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t266;
t191 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t267;
t164 = Ifges(3,1) * t267 + Ifges(3,5) * qJD(2) + t211;
t163 = Ifges(3,6) * qJD(2) + qJD(1) * t252;
t162 = pkin(2) * t275 + t208 * t226;
t157 = pkin(4) + t161;
t150 = mrSges(4,1) * t223 - t290;
t149 = -mrSges(4,2) * t223 + mrSges(4,3) * t165;
t148 = t212 + t304;
t131 = -mrSges(4,1) * t165 + mrSges(4,2) * t166;
t110 = mrSges(5,1) * t216 - mrSges(5,3) * t126;
t109 = mrSges(6,1) * t216 - mrSges(6,3) * t126;
t108 = -mrSges(5,2) * t216 + mrSges(5,3) * t253;
t107 = -mrSges(6,2) * t216 + mrSges(6,3) * t253;
t92 = -pkin(4) * t136 + t151;
t87 = -mrSges(4,2) * t222 + mrSges(4,3) * t104;
t86 = mrSges(4,1) * t222 - mrSges(4,3) * t103;
t78 = t212 + t85;
t75 = -mrSges(5,1) * t253 + mrSges(5,2) * t126;
t74 = -mrSges(6,1) * t253 + mrSges(6,2) * t126;
t55 = -qJD(4) * t137 - t140 * t226 + t141 * t230;
t54 = qJD(4) * t136 + t140 * t230 + t141 * t226;
t43 = qJ(5) * t136 + t71;
t42 = -qJ(5) * t137 + t70;
t37 = -pkin(4) * t55 + t127;
t27 = -t118 + t51;
t26 = t50 - t281;
t21 = mrSges(5,1) * t215 - mrSges(5,3) * t33;
t20 = mrSges(6,1) * t215 - mrSges(6,3) * t33;
t11 = -pkin(4) * t34 + qJDD(5) + t77;
t8 = -qJ(5) * t54 - qJD(5) * t137 + t10;
t7 = qJ(5) * t55 + qJD(5) * t136 + t9;
t1 = [t353 * t54 / 0.2e1 + t354 * t55 / 0.2e1 + m(4) * (t132 * t82 + t133 * t81 + t146 * t58 + t147 * t57 - t158 * t209 - t194 * t213) + (m(3) * t280 + mrSges(3,1) * t188 - mrSges(3,2) * t189) * pkin(1) + t330 * mrSges(3,3) + (-t191 * t272 - t192 * t273 - t228 * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t189) + m(3) * t330 + t232 * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t188)) * pkin(6) + (-t145 * mrSges(5,1) - t76 * mrSges(6,1) - t328) * t55 + (t229 * t327 + t233 * t326) * g(1) + (t229 * t326 - t233 * t327) * g(2) + (mrSges(5,2) * t77 + mrSges(6,2) * t11 - mrSges(5,3) * t6 - mrSges(6,3) * t2 + t215 * t339 + t33 * t341 + t34 * t340) * t137 + (-mrSges(5,1) * t77 - mrSges(6,1) * t11 + mrSges(5,3) * t5 + mrSges(6,3) * t3 + t215 * t337 + t33 * t340 + t338 * t34) * t136 + (t145 * mrSges(5,2) + t76 * mrSges(6,2) - t347) * t54 + (t248 + t232 * (-Ifges(3,2) * t228 + t288)) * t265 / 0.2e1 + t232 * (Ifges(3,4) * t189 + Ifges(3,2) * t188) / 0.2e1 + (-mrSges(4,1) * t158 + mrSges(4,3) * t57 + Ifges(4,4) * t103 + Ifges(4,2) * t104 + Ifges(4,6) * t222) * t178 + m(6) * (t11 * t92 + t13 * t8 + t18 * t7 + t2 * t42 + t3 * t43 + t37 * t76) + m(5) * (t10 * t48 + t127 * t145 + t151 * t77 + t49 * t9 + t5 * t71 + t6 * t70) + (Ifges(4,1) * t140 + Ifges(4,4) * t141) * t311 + t141 * t284 + (mrSges(4,2) * t158 - mrSges(4,3) * t58 + Ifges(4,1) * t103 + Ifges(4,4) * t104 + Ifges(4,5) * t222) * t179 + (Ifges(3,1) * t189 + Ifges(3,4) * t188) * t307 + t131 * t213 + (t337 * t55 + t339 * t54) * t216 / 0.2e1 + (t338 * t55 + t340 * t54) * t253 / 0.2e1 + qJD(2) ^ 2 * t251 / 0.2e1 + (t340 * t55 + t341 * t54) * t315 + t188 * t252 / 0.2e1 + (0.2e1 * Ifges(3,5) * t307 + Ifges(3,6) * t232) * qJDD(2) + t92 * t257 - t249 * t265 + t164 * t272 / 0.2e1 - t163 * t273 / 0.2e1 - t193 * t280 - t140 * t285 + t189 * (t228 * Ifges(3,1) + t288) / 0.2e1 + t42 * t20 + t43 * t22 + Ifges(2,3) * qJDD(1) + t70 * t21 + t71 * t23 + t37 * t74 + t7 * t107 + t9 * t108 + t8 * t109 + t10 * t110 + t127 * t75 + t140 * t117 / 0.2e1 + t141 * t116 / 0.2e1 + t146 * t86 + t147 * t87 + t81 * t149 + t82 * t150 + t151 * (-mrSges(5,1) * t34 + mrSges(5,2) * t33) + t165 * (Ifges(4,4) * t140 + Ifges(4,2) * t141) / 0.2e1 - t194 * (-mrSges(4,1) * t141 + mrSges(4,2) * t140) - t209 * (-mrSges(4,1) * t104 + mrSges(4,2) * t103) + t223 * (Ifges(4,5) * t140 + Ifges(4,6) * t141) / 0.2e1; (-qJD(2) * t251 / 0.2e1 + t163 * t307 + (-t248 / 0.2e1 + t249 + t282 * t307) * qJD(1) + (t232 * t191 + t228 * t192) * pkin(6) + (m(4) * t194 - t131) * t305 - (t164 + t211) * t232 / 0.2e1) * qJD(1) + t349 * t108 + (-t274 * g(3) - t145 * t148 + t161 * t6 + t162 * t5 + t349 * t49 + t350 * t48) * m(5) + t350 * t110 + t351 * t107 + (-t262 * g(3) + t352 * t13 + t157 * t2 + t162 * t3 + t351 * t18 - t76 * t78) * m(6) + t352 * t109 + t346 * (-m(5) * (-t303 - t305) - m(6) * (t174 - t305) + t292 + (mrSges(3,1) + t325) * t228 + t238) + (-g(3) * m(4) * t232 + t227 * t87 + t231 * t86 + (t149 * t231 - t150 * t227) * qJD(3)) * pkin(2) + (t227 * t57 + t231 * t58 + (-t132 * t227 + t133 * t231) * qJD(3)) * t325 + t166 * t284 + t358 * g(3) + t235 - m(4) * (t132 * t138 + t133 * t139) + Ifges(3,3) * qJDD(2) + t293 * t162 - t78 * t74 - t148 * t75 - t139 * t149 - t138 * t150 + t157 * t20 + t161 * t21 - t176 * mrSges(3,2) - t177 * mrSges(3,1) + Ifges(3,6) * t188 + Ifges(3,5) * t189; (-m(6) * (t200 + t204) + t241) * g(3) - m(6) * (t13 * t26 + t18 * t27 + t76 * t85) + (t150 + t290) * t133 + (-t166 * t75 + t230 * t21 + (m(6) * t3 + t293) * t226 + ((m(6) * t18 + t107 + t108) * t230 + (-m(6) * t13 - t109 - t110) * t226) * qJD(4) + (-g(3) * t218 - t145 * t166 + t346 * t217 + t226 * t5 + t230 * t6 + t49 * t268 - t48 * t269) * m(5)) * pkin(3) - m(5) * (t48 * t50 + t49 * t51) + t235 - t85 * t74 - t27 * t107 - t51 * t108 - t26 * t109 - t50 * t110 - t132 * t149 + t346 * (-m(6) * t174 + t238) + (m(6) * t2 + t20) * (pkin(3) * t230 + pkin(4)); (t126 * t18 + t13 * t253) * mrSges(6,3) + (-(-t13 + t17) * t18 - g(3) * t200 + (-t126 * t76 + t2) * pkin(4)) * m(6) + (t126 * t49 + t253 * t48) * mrSges(5,3) + t250 * g(3) - t17 * t107 - t48 * t108 + t18 * t109 + t49 * t110 + (-t126 * t74 + t20) * pkin(4) + t346 * (t256 + (m(6) * pkin(4) + t295) * t205) + t359; -t253 * t107 + t126 * t109 + (-g(1) * t229 + g(2) * t233 + t13 * t126 - t18 * t253 + t11) * m(6) + t257;];
tau = t1;
