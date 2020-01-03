% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP6
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
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:00:05
% EndTime: 2019-12-31 21:00:18
% DurationCPUTime: 5.87s
% Computational Cost: add. (72154->306), mult. (147982->374), div. (0->0), fcn. (97474->8), ass. (0->114)
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t245 = -t257 * g(1) - t254 * g(2);
t259 = qJD(1) ^ 2;
t225 = -t259 * pkin(1) + qJDD(1) * pkin(6) + t245;
t253 = sin(qJ(2));
t256 = cos(qJ(2));
t215 = -t256 * g(3) - t253 * t225;
t238 = (-pkin(2) * t256 - pkin(7) * t253) * qJD(1);
t258 = qJD(2) ^ 2;
t279 = qJD(1) * t253;
t193 = -qJDD(2) * pkin(2) - t258 * pkin(7) + t238 * t279 - t215;
t252 = sin(qJ(3));
t255 = cos(qJ(3));
t236 = t252 * qJD(2) + t255 * t279;
t277 = qJD(1) * qJD(2);
t274 = t256 * t277;
t239 = t253 * qJDD(1) + t274;
t207 = -t236 * qJD(3) + t255 * qJDD(2) - t252 * t239;
t278 = t256 * qJD(1);
t247 = qJD(3) - t278;
t213 = t247 * pkin(3) - t236 * qJ(4);
t235 = t255 * qJD(2) - t252 * t279;
t233 = t235 ^ 2;
t158 = -t207 * pkin(3) - t233 * qJ(4) + t236 * t213 + qJDD(4) + t193;
t208 = t235 * qJD(3) + t252 * qJDD(2) + t255 * t239;
t251 = sin(pkin(8));
t282 = cos(pkin(8));
t179 = -t282 * t207 + t251 * t208;
t180 = t251 * t207 + t282 * t208;
t209 = -t282 * t235 + t251 * t236;
t210 = t251 * t235 + t282 * t236;
t153 = -0.2e1 * qJD(5) * t210 + (t209 * t247 - t180) * qJ(5) + (t210 * t247 + t179) * pkin(4) + t158;
t197 = -t247 * mrSges(6,1) + t210 * mrSges(6,2);
t198 = -t209 * mrSges(6,2) + t247 * mrSges(6,3);
t143 = m(6) * t153 + t179 * mrSges(6,1) - t180 * mrSges(6,3) - t210 * t197 + t209 * t198;
t244 = t254 * g(1) - t257 * g(2);
t224 = -qJDD(1) * pkin(1) - t259 * pkin(6) - t244;
t275 = t253 * t277;
t240 = t256 * qJDD(1) - t275;
t189 = (-t239 - t274) * pkin(7) + (-t240 + t275) * pkin(2) + t224;
t216 = -t253 * g(3) + t256 * t225;
t194 = -t258 * pkin(2) + qJDD(2) * pkin(7) + t238 * t278 + t216;
t161 = t255 * t189 - t252 * t194;
t234 = qJDD(3) - t240;
t155 = (t235 * t247 - t208) * qJ(4) + (t235 * t236 + t234) * pkin(3) + t161;
t162 = t252 * t189 + t255 * t194;
t157 = -t233 * pkin(3) + t207 * qJ(4) - t247 * t213 + t162;
t284 = -2 * qJD(4);
t151 = t251 * t155 + t282 * t157 + t209 * t284;
t177 = Ifges(6,1) * t210 + Ifges(6,4) * t247 + Ifges(6,5) * t209;
t178 = Ifges(5,1) * t210 - Ifges(5,4) * t209 + Ifges(5,5) * t247;
t184 = t209 * pkin(4) - t210 * qJ(5);
t246 = t247 ^ 2;
t146 = -t246 * pkin(4) + t234 * qJ(5) + 0.2e1 * qJD(5) * t247 - t209 * t184 + t151;
t270 = -mrSges(6,1) * t153 + mrSges(6,2) * t146;
t175 = Ifges(6,4) * t210 + Ifges(6,2) * t247 + Ifges(6,6) * t209;
t281 = -Ifges(5,5) * t210 + Ifges(5,6) * t209 - Ifges(5,3) * t247 - t175;
t129 = -mrSges(5,1) * t158 + mrSges(5,3) * t151 - pkin(4) * t143 + (t177 + t178) * t247 + (Ifges(5,6) - Ifges(6,6)) * t234 + t281 * t210 + (Ifges(5,4) - Ifges(6,5)) * t180 + (-Ifges(5,2) - Ifges(6,3)) * t179 + t270;
t268 = t282 * t155 - t251 * t157;
t150 = t210 * t284 + t268;
t176 = Ifges(5,4) * t210 - Ifges(5,2) * t209 + Ifges(5,6) * t247;
t148 = -t234 * pkin(4) - t246 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t184) * t210 - t268;
t173 = Ifges(6,5) * t210 + Ifges(6,6) * t247 + Ifges(6,3) * t209;
t267 = mrSges(6,2) * t148 - mrSges(6,3) * t153 + Ifges(6,1) * t180 + Ifges(6,4) * t234 + Ifges(6,5) * t179 + t247 * t173;
t130 = mrSges(5,2) * t158 - mrSges(5,3) * t150 + Ifges(5,1) * t180 - Ifges(5,4) * t179 + Ifges(5,5) * t234 - qJ(5) * t143 - t247 * t176 + t281 * t209 + t267;
t200 = Ifges(4,5) * t236 + Ifges(4,6) * t235 + Ifges(4,3) * t247;
t202 = Ifges(4,1) * t236 + Ifges(4,4) * t235 + Ifges(4,5) * t247;
t195 = -t247 * mrSges(5,2) - t209 * mrSges(5,3);
t196 = t247 * mrSges(5,1) - t210 * mrSges(5,3);
t264 = m(5) * t158 + t179 * mrSges(5,1) + t180 * mrSges(5,2) + t209 * t195 + t210 * t196 + t143;
t276 = m(6) * t146 + t234 * mrSges(6,3) + t247 * t197;
t185 = t209 * mrSges(6,1) - t210 * mrSges(6,3);
t280 = -t209 * mrSges(5,1) - t210 * mrSges(5,2) - t185;
t283 = -mrSges(5,3) - mrSges(6,2);
t134 = m(5) * t151 - t234 * mrSges(5,2) + t283 * t179 - t247 * t196 + t280 * t209 + t276;
t271 = -m(6) * t148 + t234 * mrSges(6,1) + t247 * t198;
t136 = m(5) * t150 + t234 * mrSges(5,1) + t283 * t180 + t247 * t195 + t280 * t210 + t271;
t272 = t282 * t134 - t251 * t136;
t115 = -mrSges(4,1) * t193 + mrSges(4,3) * t162 + Ifges(4,4) * t208 + Ifges(4,2) * t207 + Ifges(4,6) * t234 - pkin(3) * t264 + qJ(4) * t272 + t282 * t129 + t251 * t130 - t236 * t200 + t247 * t202;
t131 = t251 * t134 + t282 * t136;
t201 = Ifges(4,4) * t236 + Ifges(4,2) * t235 + Ifges(4,6) * t247;
t116 = mrSges(4,2) * t193 - mrSges(4,3) * t161 + Ifges(4,1) * t208 + Ifges(4,4) * t207 + Ifges(4,5) * t234 - qJ(4) * t131 - t251 * t129 + t282 * t130 + t235 * t200 - t247 * t201;
t211 = -t235 * mrSges(4,1) + t236 * mrSges(4,2);
t212 = -t247 * mrSges(4,2) + t235 * mrSges(4,3);
t127 = m(4) * t161 + t234 * mrSges(4,1) - t208 * mrSges(4,3) - t236 * t211 + t247 * t212 + t131;
t214 = t247 * mrSges(4,1) - t236 * mrSges(4,3);
t128 = m(4) * t162 - t234 * mrSges(4,2) + t207 * mrSges(4,3) + t235 * t211 - t247 * t214 + t272;
t125 = -t252 * t127 + t255 * t128;
t138 = -m(4) * t193 + t207 * mrSges(4,1) - t208 * mrSges(4,2) + t235 * t212 - t236 * t214 - t264;
t222 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t253 + Ifges(3,2) * t256) * qJD(1);
t223 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t253 + Ifges(3,4) * t256) * qJD(1);
t285 = mrSges(3,1) * t215 - mrSges(3,2) * t216 + Ifges(3,5) * t239 + Ifges(3,6) * t240 + Ifges(3,3) * qJDD(2) + pkin(2) * t138 + pkin(7) * t125 + t255 * t115 + t252 * t116 + (t253 * t222 - t256 * t223) * qJD(1);
t237 = (-mrSges(3,1) * t256 + mrSges(3,2) * t253) * qJD(1);
t242 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t279;
t123 = m(3) * t216 - qJDD(2) * mrSges(3,2) + t240 * mrSges(3,3) - qJD(2) * t242 + t237 * t278 + t125;
t243 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t278;
t137 = m(3) * t215 + qJDD(2) * mrSges(3,1) - t239 * mrSges(3,3) + qJD(2) * t243 - t237 * t279 + t138;
t273 = t256 * t123 - t253 * t137;
t124 = t255 * t127 + t252 * t128;
t221 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t253 + Ifges(3,6) * t256) * qJD(1);
t112 = mrSges(3,2) * t224 - mrSges(3,3) * t215 + Ifges(3,1) * t239 + Ifges(3,4) * t240 + Ifges(3,5) * qJDD(2) - pkin(7) * t124 - qJD(2) * t222 - t252 * t115 + t255 * t116 + t221 * t278;
t265 = mrSges(6,1) * t148 - mrSges(6,3) * t146 - Ifges(6,4) * t180 - Ifges(6,2) * t234 - Ifges(6,6) * t179 + t210 * t173 - t209 * t177;
t262 = mrSges(5,2) * t151 - t209 * t178 - qJ(5) * (-t179 * mrSges(6,2) - t209 * t185 + t276) - pkin(4) * (-t180 * mrSges(6,2) - t210 * t185 + t271) - mrSges(5,1) * t150 - t210 * t176 + Ifges(5,6) * t179 - Ifges(5,5) * t180 - Ifges(5,3) * t234 + t265;
t260 = mrSges(4,1) * t161 - mrSges(4,2) * t162 + Ifges(4,5) * t208 + Ifges(4,6) * t207 + Ifges(4,3) * t234 + pkin(3) * t131 + t236 * t201 - t235 * t202 - t262;
t114 = -mrSges(3,1) * t224 + mrSges(3,3) * t216 + Ifges(3,4) * t239 + Ifges(3,2) * t240 + Ifges(3,6) * qJDD(2) - pkin(2) * t124 + qJD(2) * t223 - t221 * t279 - t260;
t263 = -m(3) * t224 + t240 * mrSges(3,1) - t239 * mrSges(3,2) - t242 * t279 + t243 * t278 - t124;
t266 = mrSges(2,1) * t244 - mrSges(2,2) * t245 + Ifges(2,3) * qJDD(1) + pkin(1) * t263 + pkin(6) * t273 + t253 * t112 + t256 * t114;
t120 = m(2) * t244 + qJDD(1) * mrSges(2,1) - t259 * mrSges(2,2) + t263;
t119 = t253 * t123 + t256 * t137;
t117 = m(2) * t245 - t259 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t273;
t110 = mrSges(2,1) * g(3) + mrSges(2,3) * t245 + t259 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t119 - t285;
t109 = -mrSges(2,2) * g(3) - mrSges(2,3) * t244 + Ifges(2,5) * qJDD(1) - t259 * Ifges(2,6) - pkin(6) * t119 + t256 * t112 - t253 * t114;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t257 * t109 - t254 * t110 - pkin(5) * (t254 * t117 + t257 * t120), t109, t112, t116, t130, -t209 * t175 + t267; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t254 * t109 + t257 * t110 + pkin(5) * (t257 * t117 - t254 * t120), t110, t114, t115, t129, -t265; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t266, t266, t285, t260, -t262, Ifges(6,5) * t180 + Ifges(6,6) * t234 + Ifges(6,3) * t179 + t210 * t175 - t247 * t177 - t270;];
m_new = t1;
