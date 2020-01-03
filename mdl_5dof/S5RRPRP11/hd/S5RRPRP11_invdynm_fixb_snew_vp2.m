% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:26
% EndTime: 2019-12-31 20:12:32
% DurationCPUTime: 2.55s
% Computational Cost: add. (23452->312), mult. (47396->362), div. (0->0), fcn. (24631->6), ass. (0->108)
t304 = -2 * qJD(3);
t262 = sin(qJ(1));
t265 = cos(qJ(1));
t246 = -g(1) * t265 - g(2) * t262;
t267 = qJD(1) ^ 2;
t215 = -pkin(1) * t267 + qJDD(1) * pkin(6) + t246;
t261 = sin(qJ(2));
t264 = cos(qJ(2));
t200 = -t264 * g(3) - t261 * t215;
t201 = -g(3) * t261 + t264 * t215;
t210 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t261 + Ifges(3,4) * t264) * qJD(1);
t233 = (mrSges(4,2) * t264 - mrSges(4,3) * t261) * qJD(1);
t291 = qJD(1) * qJD(2);
t288 = t264 * t291;
t235 = qJDD(1) * t261 + t288;
t287 = t261 * t291;
t236 = qJDD(1) * t264 - t287;
t292 = qJD(1) * t264;
t242 = -mrSges(4,1) * t292 - qJD(2) * mrSges(4,3);
t232 = (-pkin(2) * t264 - qJ(3) * t261) * qJD(1);
t266 = qJD(2) ^ 2;
t169 = pkin(2) * t266 - qJDD(2) * qJ(3) + qJD(2) * t304 - t232 * t292 - t201;
t293 = qJD(1) * t261;
t243 = mrSges(4,1) * t293 + qJD(2) * mrSges(4,2);
t244 = pkin(3) * t293 - qJD(2) * pkin(7);
t259 = t264 ^ 2;
t164 = -pkin(7) * t259 * t267 + pkin(3) * t236 + qJD(2) * t244 - t169;
t260 = sin(qJ(4));
t263 = cos(qJ(4));
t231 = qJD(2) * t263 - t260 * t292;
t189 = qJD(4) * t231 + qJDD(2) * t260 + t263 * t236;
t230 = qJD(2) * t260 + t263 * t292;
t190 = -qJD(4) * t230 + qJDD(2) * t263 - t236 * t260;
t249 = qJD(4) + t293;
t156 = -0.2e1 * qJD(5) * t231 + (t230 * t249 - t190) * qJ(5) + (t231 * t249 + t189) * pkin(4) + t164;
t197 = -mrSges(6,2) * t230 + mrSges(6,3) * t249;
t199 = -mrSges(6,1) * t249 + mrSges(6,2) * t231;
t149 = m(6) * t156 + t189 * mrSges(6,1) - mrSges(6,3) * t190 + t230 * t197 - t199 * t231;
t196 = -mrSges(5,2) * t249 - mrSges(5,3) * t230;
t198 = mrSges(5,1) * t249 - mrSges(5,3) * t231;
t273 = m(5) * t164 + mrSges(5,1) * t189 + t190 * mrSges(5,2) + t196 * t230 + t231 * t198 + t149;
t272 = -m(4) * t169 + qJDD(2) * mrSges(4,3) + qJD(2) * t243 + t233 * t292 + t273;
t245 = t262 * g(1) - t265 * g(2);
t283 = -qJDD(1) * pkin(1) - t245;
t274 = pkin(2) * t287 + t293 * t304 + (-t235 - t288) * qJ(3) + t283;
t161 = -t244 * t293 + (-pkin(3) * t259 - pkin(6)) * t267 + (-pkin(2) - pkin(7)) * t236 + t274;
t171 = -qJDD(2) * pkin(2) - t266 * qJ(3) + t232 * t293 + qJDD(3) - t200;
t165 = (-t261 * t264 * t267 - qJDD(2)) * pkin(7) + (t235 - t288) * pkin(3) + t171;
t159 = t263 * t161 + t260 * t165;
t177 = Ifges(6,1) * t231 + Ifges(6,4) * t249 + Ifges(6,5) * t230;
t178 = Ifges(5,1) * t231 - Ifges(5,4) * t230 + Ifges(5,5) * t249;
t229 = qJDD(4) + t235;
t193 = pkin(4) * t230 - qJ(5) * t231;
t247 = t249 ^ 2;
t152 = -pkin(4) * t247 + qJ(5) * t229 + 0.2e1 * qJD(5) * t249 - t193 * t230 + t159;
t284 = -mrSges(6,1) * t156 + mrSges(6,2) * t152;
t175 = Ifges(6,4) * t231 + Ifges(6,2) * t249 + Ifges(6,6) * t230;
t297 = -Ifges(5,5) * t231 + Ifges(5,6) * t230 - Ifges(5,3) * t249 - t175;
t132 = -mrSges(5,1) * t164 + mrSges(5,3) * t159 - pkin(4) * t149 + (t177 + t178) * t249 + t297 * t231 + (Ifges(5,6) - Ifges(6,6)) * t229 + (Ifges(5,4) - Ifges(6,5)) * t190 + (-Ifges(5,2) - Ifges(6,3)) * t189 + t284;
t158 = -t161 * t260 + t165 * t263;
t176 = Ifges(5,4) * t231 - Ifges(5,2) * t230 + Ifges(5,6) * t249;
t154 = -pkin(4) * t229 - qJ(5) * t247 + t193 * t231 + qJDD(5) - t158;
t173 = Ifges(6,5) * t231 + Ifges(6,6) * t249 + Ifges(6,3) * t230;
t281 = mrSges(6,2) * t154 - mrSges(6,3) * t156 + Ifges(6,1) * t190 + Ifges(6,4) * t229 + Ifges(6,5) * t189 + t249 * t173;
t134 = mrSges(5,2) * t164 - mrSges(5,3) * t158 + Ifges(5,1) * t190 - Ifges(5,4) * t189 + Ifges(5,5) * t229 - qJ(5) * t149 - t176 * t249 + t230 * t297 + t281;
t289 = m(6) * t152 + t229 * mrSges(6,3) + t249 * t199;
t194 = mrSges(6,1) * t230 - mrSges(6,3) * t231;
t296 = -mrSges(5,1) * t230 - mrSges(5,2) * t231 - t194;
t300 = -mrSges(5,3) - mrSges(6,2);
t140 = m(5) * t159 - t229 * mrSges(5,2) + t189 * t300 - t249 * t198 + t230 * t296 + t289;
t285 = -m(6) * t154 + t229 * mrSges(6,1) + t249 * t197;
t143 = m(5) * t158 + t229 * mrSges(5,1) + t190 * t300 + t249 * t196 + t231 * t296 + t285;
t135 = t260 * t140 + t263 * t143;
t212 = Ifges(4,4) * qJD(2) + (-Ifges(4,2) * t261 - Ifges(4,6) * t264) * qJD(1);
t276 = -mrSges(4,2) * t171 + mrSges(4,3) * t169 - Ifges(4,1) * qJDD(2) + Ifges(4,4) * t235 + Ifges(4,5) * t236 + pkin(7) * t135 + t260 * t132 - t263 * t134 - t212 * t292;
t279 = -m(4) * t171 - t235 * mrSges(4,1) - t135;
t211 = Ifges(4,5) * qJD(2) + (-Ifges(4,6) * t261 - Ifges(4,3) * t264) * qJD(1);
t294 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t261 + Ifges(3,2) * t264) * qJD(1) - t211;
t303 = (-t264 * t210 + t261 * t294) * qJD(1) + mrSges(3,1) * t200 - mrSges(3,2) * t201 + Ifges(3,5) * t235 + Ifges(3,6) * t236 + Ifges(3,3) * qJDD(2) + pkin(2) * (-qJDD(2) * mrSges(4,2) - qJD(2) * t242 - t233 * t293 + t279) + qJ(3) * (mrSges(4,1) * t236 + t272) - t276;
t301 = t267 * pkin(6);
t299 = Ifges(3,4) + Ifges(4,6);
t136 = t263 * t140 - t260 * t143;
t213 = Ifges(4,1) * qJD(2) + (-Ifges(4,4) * t261 - Ifges(4,5) * t264) * qJD(1);
t295 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t261 + Ifges(3,6) * t264) * qJD(1) + t213;
t234 = (-mrSges(3,1) * t264 + mrSges(3,2) * t261) * qJD(1);
t241 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t292;
t129 = m(3) * t200 - t235 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t241 - t242) * qJD(2) + (-t233 - t234) * t293 + t279;
t240 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t293;
t138 = t234 * t292 - qJDD(2) * mrSges(3,2) + t272 + m(3) * t201 - qJD(2) * t240 + (mrSges(3,3) + mrSges(4,1)) * t236;
t286 = -t129 * t261 + t264 * t138;
t166 = -t236 * pkin(2) + t274 - t301;
t282 = -m(4) * t166 - t236 * mrSges(4,2) + t243 * t293 - t136;
t130 = -t235 * mrSges(4,3) + t242 * t292 - t282;
t214 = t283 - t301;
t275 = -mrSges(4,1) * t169 + mrSges(4,2) * t166 + pkin(3) * t273 - pkin(7) * t136 - t263 * t132 - t260 * t134;
t121 = -mrSges(3,1) * t214 + mrSges(3,3) * t201 - pkin(2) * t130 + (Ifges(3,2) + Ifges(4,3)) * t236 + t299 * t235 + (Ifges(3,6) - Ifges(4,5)) * qJDD(2) + (t210 - t212) * qJD(2) - t295 * t293 + t275;
t277 = mrSges(6,1) * t154 - mrSges(6,3) * t152 - Ifges(6,4) * t190 - Ifges(6,2) * t229 - Ifges(6,6) * t189 + t231 * t173 - t230 * t177;
t270 = mrSges(5,2) * t159 - t230 * t178 - qJ(5) * (-t189 * mrSges(6,2) - t230 * t194 + t289) - pkin(4) * (-t190 * mrSges(6,2) - t231 * t194 + t285) - mrSges(5,1) * t158 - t231 * t176 + Ifges(5,6) * t189 - Ifges(5,5) * t190 - Ifges(5,3) * t229 + t277;
t269 = -mrSges(4,1) * t171 + mrSges(4,3) * t166 - pkin(3) * t135 + t270;
t123 = t299 * t236 + (Ifges(3,1) + Ifges(4,2)) * t235 + (Ifges(3,5) - Ifges(4,4)) * qJDD(2) - t294 * qJD(2) - t269 - mrSges(3,3) * t200 + mrSges(3,2) * t214 - qJ(3) * t130 + t295 * t292;
t271 = -m(3) * t214 + t241 * t292 + t236 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t235 + (-t240 * t261 - t242 * t264) * qJD(1) + t282;
t278 = mrSges(2,1) * t245 - mrSges(2,2) * t246 + Ifges(2,3) * qJDD(1) + pkin(1) * t271 + pkin(6) * t286 + t264 * t121 + t261 * t123;
t127 = m(2) * t245 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t267 + t271;
t126 = t129 * t264 + t138 * t261;
t124 = m(2) * t246 - mrSges(2,1) * t267 - qJDD(1) * mrSges(2,2) + t286;
t119 = mrSges(2,1) * g(3) + mrSges(2,3) * t246 + t267 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t126 - t303;
t118 = -mrSges(2,2) * g(3) - mrSges(2,3) * t245 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t267 - pkin(6) * t126 - t121 * t261 + t123 * t264;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t265 * t118 - t262 * t119 - pkin(5) * (t124 * t262 + t127 * t265), t118, t123, -t211 * t293 - t276, t134, -t175 * t230 + t281; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t262 * t118 + t265 * t119 + pkin(5) * (t124 * t265 - t127 * t262), t119, t121, Ifges(4,4) * qJDD(2) - Ifges(4,2) * t235 - Ifges(4,6) * t236 - qJD(2) * t211 - t213 * t292 + t269, t132, -t277; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t278, t278, t303, Ifges(4,5) * qJDD(2) - Ifges(4,6) * t235 - Ifges(4,3) * t236 + qJD(2) * t212 + t213 * t293 - t275, -t270, Ifges(6,5) * t190 + Ifges(6,6) * t229 + Ifges(6,3) * t189 + t175 * t231 - t177 * t249 - t284;];
m_new = t1;
