% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP9
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:03:54
% EndTime: 2019-12-31 22:04:08
% DurationCPUTime: 5.97s
% Computational Cost: add. (76683->306), mult. (152198->372), div. (0->0), fcn. (101666->8), ass. (0->114)
t256 = sin(qJ(1));
t259 = cos(qJ(1));
t245 = -t259 * g(1) - t256 * g(2);
t261 = qJD(1) ^ 2;
t229 = -t261 * pkin(1) + qJDD(1) * pkin(6) + t245;
t255 = sin(qJ(2));
t258 = cos(qJ(2));
t213 = -t258 * g(3) - t255 * t229;
t238 = (-pkin(2) * t258 - pkin(7) * t255) * qJD(1);
t260 = qJD(2) ^ 2;
t279 = qJD(1) * t255;
t193 = -qJDD(2) * pkin(2) - t260 * pkin(7) + t238 * t279 - t213;
t254 = sin(qJ(3));
t257 = cos(qJ(3));
t236 = t254 * qJD(2) + t257 * t279;
t277 = qJD(1) * qJD(2);
t275 = t258 * t277;
t239 = t255 * qJDD(1) + t275;
t206 = -t236 * qJD(3) + t257 * qJDD(2) - t254 * t239;
t278 = t258 * qJD(1);
t248 = qJD(3) - t278;
t215 = t248 * pkin(3) - t236 * pkin(8);
t235 = t257 * qJD(2) - t254 * t279;
t233 = t235 ^ 2;
t158 = -t206 * pkin(3) - t233 * pkin(8) + t236 * t215 + t193;
t207 = t235 * qJD(3) + t254 * qJDD(2) + t257 * t239;
t253 = sin(qJ(4));
t283 = cos(qJ(4));
t209 = t253 * t235 + t283 * t236;
t170 = t209 * qJD(4) - t283 * t206 + t253 * t207;
t208 = -t283 * t235 + t253 * t236;
t171 = -t208 * qJD(4) + t253 * t206 + t283 * t207;
t247 = qJD(4) + t248;
t150 = -0.2e1 * qJD(5) * t209 + (t208 * t247 - t171) * qJ(5) + (t209 * t247 + t170) * pkin(4) + t158;
t195 = -t208 * mrSges(6,2) + t247 * mrSges(6,3);
t198 = -t247 * mrSges(6,1) + t209 * mrSges(6,2);
t139 = m(6) * t150 + t170 * mrSges(6,1) - t171 * mrSges(6,3) + t208 * t195 - t209 * t198;
t244 = t256 * g(1) - t259 * g(2);
t228 = -qJDD(1) * pkin(1) - t261 * pkin(6) - t244;
t249 = t255 * t277;
t240 = t258 * qJDD(1) - t249;
t189 = (-t239 - t275) * pkin(7) + (-t240 + t249) * pkin(2) + t228;
t214 = -t255 * g(3) + t258 * t229;
t194 = -t260 * pkin(2) + qJDD(2) * pkin(7) + t238 * t278 + t214;
t172 = t257 * t189 - t254 * t194;
t234 = qJDD(3) - t240;
t155 = (t235 * t248 - t207) * pkin(8) + (t235 * t236 + t234) * pkin(3) + t172;
t173 = t254 * t189 + t257 * t194;
t157 = -t233 * pkin(3) + t206 * pkin(8) - t248 * t215 + t173;
t153 = t253 * t155 + t283 * t157;
t179 = Ifges(6,1) * t209 + Ifges(6,4) * t247 + Ifges(6,5) * t208;
t180 = Ifges(5,1) * t209 - Ifges(5,4) * t208 + Ifges(5,5) * t247;
t230 = qJDD(4) + t234;
t184 = t208 * pkin(4) - t209 * qJ(5);
t246 = t247 ^ 2;
t146 = -t246 * pkin(4) + t230 * qJ(5) + 0.2e1 * qJD(5) * t247 - t208 * t184 + t153;
t271 = -mrSges(6,1) * t150 + mrSges(6,2) * t146;
t177 = Ifges(6,4) * t209 + Ifges(6,2) * t247 + Ifges(6,6) * t208;
t281 = -Ifges(5,5) * t209 + Ifges(5,6) * t208 - Ifges(5,3) * t247 - t177;
t126 = -mrSges(5,1) * t158 + mrSges(5,3) * t153 - pkin(4) * t139 + (t179 + t180) * t247 + (Ifges(5,6) - Ifges(6,6)) * t230 + t281 * t209 + (Ifges(5,4) - Ifges(6,5)) * t171 + (-Ifges(5,2) - Ifges(6,3)) * t170 + t271;
t152 = t283 * t155 - t253 * t157;
t178 = Ifges(5,4) * t209 - Ifges(5,2) * t208 + Ifges(5,6) * t247;
t148 = -t230 * pkin(4) - t246 * qJ(5) + t209 * t184 + qJDD(5) - t152;
t175 = Ifges(6,5) * t209 + Ifges(6,6) * t247 + Ifges(6,3) * t208;
t269 = mrSges(6,2) * t148 - mrSges(6,3) * t150 + Ifges(6,1) * t171 + Ifges(6,4) * t230 + Ifges(6,5) * t170 + t247 * t175;
t127 = mrSges(5,2) * t158 - mrSges(5,3) * t152 + Ifges(5,1) * t171 - Ifges(5,4) * t170 + Ifges(5,5) * t230 - qJ(5) * t139 - t247 * t178 + t281 * t208 + t269;
t200 = Ifges(4,5) * t236 + Ifges(4,6) * t235 + Ifges(4,3) * t248;
t202 = Ifges(4,1) * t236 + Ifges(4,4) * t235 + Ifges(4,5) * t248;
t196 = -t247 * mrSges(5,2) - t208 * mrSges(5,3);
t197 = t247 * mrSges(5,1) - t209 * mrSges(5,3);
t266 = m(5) * t158 + t170 * mrSges(5,1) + t171 * mrSges(5,2) + t208 * t196 + t209 * t197 + t139;
t276 = m(6) * t146 + t230 * mrSges(6,3) + t247 * t198;
t185 = t208 * mrSges(6,1) - t209 * mrSges(6,3);
t280 = -t208 * mrSges(5,1) - t209 * mrSges(5,2) - t185;
t282 = -mrSges(5,3) - mrSges(6,2);
t136 = m(5) * t153 - t230 * mrSges(5,2) + t282 * t170 - t247 * t197 + t280 * t208 + t276;
t272 = -m(6) * t148 + t230 * mrSges(6,1) + t247 * t195;
t138 = m(5) * t152 + t230 * mrSges(5,1) + t282 * t171 + t247 * t196 + t280 * t209 + t272;
t273 = t283 * t136 - t253 * t138;
t115 = -mrSges(4,1) * t193 + mrSges(4,3) * t173 + Ifges(4,4) * t207 + Ifges(4,2) * t206 + Ifges(4,6) * t234 - pkin(3) * t266 + pkin(8) * t273 + t283 * t126 + t253 * t127 - t236 * t200 + t248 * t202;
t131 = t253 * t136 + t283 * t138;
t201 = Ifges(4,4) * t236 + Ifges(4,2) * t235 + Ifges(4,6) * t248;
t116 = mrSges(4,2) * t193 - mrSges(4,3) * t172 + Ifges(4,1) * t207 + Ifges(4,4) * t206 + Ifges(4,5) * t234 - pkin(8) * t131 - t253 * t126 + t283 * t127 + t235 * t200 - t248 * t201;
t210 = -t235 * mrSges(4,1) + t236 * mrSges(4,2);
t211 = -t248 * mrSges(4,2) + t235 * mrSges(4,3);
t129 = m(4) * t172 + t234 * mrSges(4,1) - t207 * mrSges(4,3) - t236 * t210 + t248 * t211 + t131;
t212 = t248 * mrSges(4,1) - t236 * mrSges(4,3);
t130 = m(4) * t173 - t234 * mrSges(4,2) + t206 * mrSges(4,3) + t235 * t210 - t248 * t212 + t273;
t125 = -t254 * t129 + t257 * t130;
t133 = -m(4) * t193 + t206 * mrSges(4,1) - t207 * mrSges(4,2) + t235 * t211 - t236 * t212 - t266;
t226 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t255 + Ifges(3,2) * t258) * qJD(1);
t227 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t255 + Ifges(3,4) * t258) * qJD(1);
t284 = mrSges(3,1) * t213 - mrSges(3,2) * t214 + Ifges(3,5) * t239 + Ifges(3,6) * t240 + Ifges(3,3) * qJDD(2) + pkin(2) * t133 + pkin(7) * t125 + t257 * t115 + t254 * t116 + (t255 * t226 - t258 * t227) * qJD(1);
t237 = (-mrSges(3,1) * t258 + mrSges(3,2) * t255) * qJD(1);
t242 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t279;
t123 = m(3) * t214 - qJDD(2) * mrSges(3,2) + t240 * mrSges(3,3) - qJD(2) * t242 + t237 * t278 + t125;
t243 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t278;
t132 = m(3) * t213 + qJDD(2) * mrSges(3,1) - t239 * mrSges(3,3) + qJD(2) * t243 - t237 * t279 + t133;
t274 = t258 * t123 - t255 * t132;
t124 = t257 * t129 + t254 * t130;
t225 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t255 + Ifges(3,6) * t258) * qJD(1);
t112 = mrSges(3,2) * t228 - mrSges(3,3) * t213 + Ifges(3,1) * t239 + Ifges(3,4) * t240 + Ifges(3,5) * qJDD(2) - pkin(7) * t124 - qJD(2) * t226 - t254 * t115 + t257 * t116 + t225 * t278;
t267 = mrSges(6,1) * t148 - mrSges(6,3) * t146 - Ifges(6,4) * t171 - Ifges(6,2) * t230 - Ifges(6,6) * t170 + t209 * t175 - t208 * t179;
t264 = mrSges(5,2) * t153 - t208 * t180 - qJ(5) * (-t170 * mrSges(6,2) - t208 * t185 + t276) - pkin(4) * (-t171 * mrSges(6,2) - t209 * t185 + t272) - mrSges(5,1) * t152 + Ifges(5,6) * t170 - Ifges(5,5) * t171 - t209 * t178 - Ifges(5,3) * t230 + t267;
t262 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,5) * t207 + Ifges(4,6) * t206 + Ifges(4,3) * t234 + pkin(3) * t131 + t236 * t201 - t235 * t202 - t264;
t114 = -mrSges(3,1) * t228 + mrSges(3,3) * t214 + Ifges(3,4) * t239 + Ifges(3,2) * t240 + Ifges(3,6) * qJDD(2) - pkin(2) * t124 + qJD(2) * t227 - t225 * t279 - t262;
t265 = -m(3) * t228 + t240 * mrSges(3,1) - t239 * mrSges(3,2) - t242 * t279 + t243 * t278 - t124;
t268 = mrSges(2,1) * t244 - mrSges(2,2) * t245 + Ifges(2,3) * qJDD(1) + pkin(1) * t265 + pkin(6) * t274 + t255 * t112 + t258 * t114;
t120 = m(2) * t244 + qJDD(1) * mrSges(2,1) - t261 * mrSges(2,2) + t265;
t119 = t255 * t123 + t258 * t132;
t117 = m(2) * t245 - t261 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t274;
t110 = mrSges(2,1) * g(3) + mrSges(2,3) * t245 + t261 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t119 - t284;
t109 = -mrSges(2,2) * g(3) - mrSges(2,3) * t244 + Ifges(2,5) * qJDD(1) - t261 * Ifges(2,6) - pkin(6) * t119 + t258 * t112 - t255 * t114;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t259 * t109 - t256 * t110 - pkin(5) * (t256 * t117 + t259 * t120), t109, t112, t116, t127, -t208 * t177 + t269; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t256 * t109 + t259 * t110 + pkin(5) * (t259 * t117 - t256 * t120), t110, t114, t115, t126, -t267; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t268, t268, t284, t262, -t264, Ifges(6,5) * t171 + Ifges(6,6) * t230 + Ifges(6,3) * t170 + t209 * t177 - t247 * t179 - t271;];
m_new = t1;
