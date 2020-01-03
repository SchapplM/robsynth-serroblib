% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:58
% EndTime: 2019-12-31 20:55:09
% DurationCPUTime: 6.90s
% Computational Cost: add. (83864->309), mult. (188312->380), div. (0->0), fcn. (128060->8), ass. (0->115)
t263 = sin(qJ(2));
t266 = cos(qJ(2));
t287 = qJD(1) * qJD(2);
t241 = t263 * qJDD(1) + t266 * t287;
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t248 = -t267 * g(1) - t264 * g(2);
t268 = qJD(1) ^ 2;
t236 = -t268 * pkin(1) + qJDD(1) * pkin(6) + t248;
t292 = t263 * t236;
t295 = pkin(2) * t268;
t198 = qJDD(2) * pkin(2) - t241 * pkin(7) - t292 + (pkin(7) * t287 + t263 * t295 - g(3)) * t266;
t223 = -t263 * g(3) + t266 * t236;
t242 = t266 * qJDD(1) - t263 * t287;
t289 = qJD(1) * t263;
t246 = qJD(2) * pkin(2) - pkin(7) * t289;
t260 = t266 ^ 2;
t199 = t242 * pkin(7) - qJD(2) * t246 - t260 * t295 + t223;
t262 = sin(qJ(3));
t265 = cos(qJ(3));
t166 = t265 * t198 - t262 * t199;
t233 = (-t262 * t263 + t265 * t266) * qJD(1);
t207 = t233 * qJD(3) + t265 * t241 + t262 * t242;
t234 = (t262 * t266 + t263 * t265) * qJD(1);
t257 = qJDD(2) + qJDD(3);
t258 = qJD(2) + qJD(3);
t158 = (t233 * t258 - t207) * qJ(4) + (t233 * t234 + t257) * pkin(3) + t166;
t167 = t262 * t198 + t265 * t199;
t206 = -t234 * qJD(3) - t262 * t241 + t265 * t242;
t225 = t258 * pkin(3) - t234 * qJ(4);
t229 = t233 ^ 2;
t160 = -t229 * pkin(3) + t206 * qJ(4) - t258 * t225 + t167;
t261 = sin(pkin(8));
t293 = cos(pkin(8));
t219 = -t293 * t233 + t261 * t234;
t296 = -2 * qJD(4);
t156 = t261 * t158 + t293 * t160 + t219 * t296;
t178 = -t293 * t206 + t261 * t207;
t220 = t261 * t233 + t293 * t234;
t210 = t258 * mrSges(5,1) - t220 * mrSges(5,3);
t191 = t219 * pkin(4) - t220 * qJ(5);
t256 = t258 ^ 2;
t149 = -t256 * pkin(4) + t257 * qJ(5) + 0.2e1 * qJD(5) * t258 - t219 * t191 + t156;
t211 = -t258 * mrSges(6,1) + t220 * mrSges(6,2);
t286 = m(6) * t149 + t257 * mrSges(6,3) + t258 * t211;
t192 = t219 * mrSges(6,1) - t220 * mrSges(6,3);
t290 = -t219 * mrSges(5,1) - t220 * mrSges(5,2) - t192;
t294 = -mrSges(5,3) - mrSges(6,2);
t139 = m(5) * t156 - t257 * mrSges(5,2) + t294 * t178 - t258 * t210 + t290 * t219 + t286;
t278 = t293 * t158 - t261 * t160;
t155 = t220 * t296 + t278;
t179 = t261 * t206 + t293 * t207;
t209 = -t258 * mrSges(5,2) - t219 * mrSges(5,3);
t151 = -t257 * pkin(4) - t256 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t191) * t220 - t278;
t212 = -t219 * mrSges(6,2) + t258 * mrSges(6,3);
t282 = -m(6) * t151 + t257 * mrSges(6,1) + t258 * t212;
t141 = m(5) * t155 + t257 * mrSges(5,1) + t294 * t179 + t258 * t209 + t290 * t220 + t282;
t134 = t261 * t139 + t293 * t141;
t221 = -t233 * mrSges(4,1) + t234 * mrSges(4,2);
t224 = -t258 * mrSges(4,2) + t233 * mrSges(4,3);
t129 = m(4) * t166 + t257 * mrSges(4,1) - t207 * mrSges(4,3) - t234 * t221 + t258 * t224 + t134;
t226 = t258 * mrSges(4,1) - t234 * mrSges(4,3);
t283 = t293 * t139 - t261 * t141;
t130 = m(4) * t167 - t257 * mrSges(4,2) + t206 * mrSges(4,3) + t233 * t221 - t258 * t226 + t283;
t125 = t265 * t129 + t262 * t130;
t222 = -t266 * g(3) - t292;
t231 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t263 + Ifges(3,2) * t266) * qJD(1);
t232 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t263 + Ifges(3,4) * t266) * qJD(1);
t215 = Ifges(4,4) * t234 + Ifges(4,2) * t233 + Ifges(4,6) * t258;
t216 = Ifges(4,1) * t234 + Ifges(4,4) * t233 + Ifges(4,5) * t258;
t184 = Ifges(5,4) * t220 - Ifges(5,2) * t219 + Ifges(5,6) * t258;
t186 = Ifges(5,1) * t220 - Ifges(5,4) * t219 + Ifges(5,5) * t258;
t181 = Ifges(6,5) * t220 + Ifges(6,6) * t258 + Ifges(6,3) * t219;
t185 = Ifges(6,1) * t220 + Ifges(6,4) * t258 + Ifges(6,5) * t219;
t274 = mrSges(6,1) * t151 - mrSges(6,3) * t149 - Ifges(6,4) * t179 - Ifges(6,2) * t257 - Ifges(6,6) * t178 + t220 * t181 - t219 * t185;
t272 = mrSges(5,2) * t156 - t219 * t186 - qJ(5) * (-t178 * mrSges(6,2) - t219 * t192 + t286) - pkin(4) * (-t179 * mrSges(6,2) - t220 * t192 + t282) - mrSges(5,1) * t155 - t220 * t184 + Ifges(5,6) * t178 - Ifges(5,5) * t179 - Ifges(5,3) * t257 + t274;
t270 = -mrSges(4,1) * t166 + mrSges(4,2) * t167 - Ifges(4,5) * t207 - Ifges(4,6) * t206 - Ifges(4,3) * t257 - pkin(3) * t134 - t234 * t215 + t233 * t216 + t272;
t297 = mrSges(3,1) * t222 - mrSges(3,2) * t223 + Ifges(3,5) * t241 + Ifges(3,6) * t242 + Ifges(3,3) * qJDD(2) + pkin(2) * t125 + (t263 * t231 - t266 * t232) * qJD(1) - t270;
t183 = Ifges(6,4) * t220 + Ifges(6,2) * t258 + Ifges(6,6) * t219;
t291 = -Ifges(5,5) * t220 + Ifges(5,6) * t219 - Ifges(5,3) * t258 - t183;
t288 = qJD(1) * t266;
t247 = t264 * g(1) - t267 * g(2);
t240 = (-mrSges(3,1) * t266 + mrSges(3,2) * t263) * qJD(1);
t245 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t288;
t123 = m(3) * t222 + qJDD(2) * mrSges(3,1) - t241 * mrSges(3,3) + qJD(2) * t245 - t240 * t289 + t125;
t244 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t289;
t284 = -t262 * t129 + t265 * t130;
t124 = m(3) * t223 - qJDD(2) * mrSges(3,2) + t242 * mrSges(3,3) - qJD(2) * t244 + t240 * t288 + t284;
t285 = -t263 * t123 + t266 * t124;
t279 = -qJDD(1) * pkin(1) - t247;
t208 = -t242 * pkin(2) + t246 * t289 + (-pkin(7) * t260 - pkin(6)) * t268 + t279;
t162 = -t206 * pkin(3) - t229 * qJ(4) + t234 * t225 + qJDD(4) + t208;
t153 = -0.2e1 * qJD(5) * t220 + (t219 * t258 - t179) * qJ(5) + (t220 * t258 + t178) * pkin(4) + t162;
t281 = -mrSges(6,1) * t153 + mrSges(6,2) * t149;
t146 = m(6) * t153 + t178 * mrSges(6,1) - t179 * mrSges(6,3) - t220 * t211 + t219 * t212;
t277 = mrSges(6,2) * t151 - mrSges(6,3) * t153 + Ifges(6,1) * t179 + Ifges(6,4) * t257 + Ifges(6,5) * t178 + t258 * t181;
t131 = -mrSges(5,1) * t162 + mrSges(5,3) * t156 - pkin(4) * t146 + (t185 + t186) * t258 + (Ifges(5,6) - Ifges(6,6)) * t257 + t291 * t220 + (Ifges(5,4) - Ifges(6,5)) * t179 + (-Ifges(5,2) - Ifges(6,3)) * t178 + t281;
t132 = mrSges(5,2) * t162 - mrSges(5,3) * t155 + Ifges(5,1) * t179 - Ifges(5,4) * t178 + Ifges(5,5) * t257 - qJ(5) * t146 - t258 * t184 + t291 * t219 + t277;
t214 = Ifges(4,5) * t234 + Ifges(4,6) * t233 + Ifges(4,3) * t258;
t275 = m(5) * t162 + t178 * mrSges(5,1) + t179 * mrSges(5,2) + t219 * t209 + t220 * t210 + t146;
t120 = -mrSges(4,1) * t208 + mrSges(4,3) * t167 + Ifges(4,4) * t207 + Ifges(4,2) * t206 + Ifges(4,6) * t257 - pkin(3) * t275 + qJ(4) * t283 + t293 * t131 + t261 * t132 - t234 * t214 + t258 * t216;
t121 = mrSges(4,2) * t208 - mrSges(4,3) * t166 + Ifges(4,1) * t207 + Ifges(4,4) * t206 + Ifges(4,5) * t257 - qJ(4) * t134 - t261 * t131 + t293 * t132 + t233 * t214 - t258 * t215;
t230 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t263 + Ifges(3,6) * t266) * qJD(1);
t235 = -t268 * pkin(6) + t279;
t273 = m(4) * t208 - t206 * mrSges(4,1) + t207 * mrSges(4,2) - t233 * t224 + t234 * t226 + t275;
t113 = -mrSges(3,1) * t235 + mrSges(3,3) * t223 + Ifges(3,4) * t241 + Ifges(3,2) * t242 + Ifges(3,6) * qJDD(2) - pkin(2) * t273 + pkin(7) * t284 + qJD(2) * t232 + t265 * t120 + t262 * t121 - t230 * t289;
t116 = mrSges(3,2) * t235 - mrSges(3,3) * t222 + Ifges(3,1) * t241 + Ifges(3,4) * t242 + Ifges(3,5) * qJDD(2) - pkin(7) * t125 - qJD(2) * t231 - t262 * t120 + t265 * t121 + t230 * t288;
t271 = -m(3) * t235 + t242 * mrSges(3,1) - t241 * mrSges(3,2) - t244 * t289 + t245 * t288 - t273;
t276 = mrSges(2,1) * t247 - mrSges(2,2) * t248 + Ifges(2,3) * qJDD(1) + pkin(1) * t271 + pkin(6) * t285 + t266 * t113 + t263 * t116;
t135 = m(2) * t247 + qJDD(1) * mrSges(2,1) - t268 * mrSges(2,2) + t271;
t119 = t266 * t123 + t263 * t124;
t117 = m(2) * t248 - t268 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t285;
t114 = mrSges(2,1) * g(3) + mrSges(2,3) * t248 + t268 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t119 - t297;
t111 = -mrSges(2,2) * g(3) - mrSges(2,3) * t247 + Ifges(2,5) * qJDD(1) - t268 * Ifges(2,6) - pkin(6) * t119 - t263 * t113 + t266 * t116;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t267 * t111 - t264 * t114 - pkin(5) * (t264 * t117 + t267 * t135), t111, t116, t121, t132, -t219 * t183 + t277; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t264 * t111 + t267 * t114 + pkin(5) * (t267 * t117 - t264 * t135), t114, t113, t120, t131, -t274; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t276, t276, t297, -t270, -t272, Ifges(6,5) * t179 + Ifges(6,6) * t257 + Ifges(6,3) * t178 + t220 * t183 - t258 * t185 - t281;];
m_new = t1;
