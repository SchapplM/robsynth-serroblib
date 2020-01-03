% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:23
% EndTime: 2019-12-31 17:53:27
% DurationCPUTime: 2.18s
% Computational Cost: add. (16729->258), mult. (39531->318), div. (0->0), fcn. (23376->6), ass. (0->102)
t263 = sin(qJ(1));
t265 = cos(qJ(1));
t230 = -g(1) * t265 - g(2) * t263;
t267 = qJD(1) ^ 2;
t312 = -pkin(1) * t267 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t230;
t260 = sin(pkin(7));
t250 = t260 ^ 2;
t261 = cos(pkin(7));
t251 = t261 ^ 2;
t300 = t251 * t267;
t311 = t250 * t267 + t300;
t201 = -t261 * g(3) - t312 * t260;
t219 = (-pkin(2) * t261 - qJ(3) * t260) * qJD(1);
t294 = qJD(1) * t260;
t171 = t219 * t294 + qJDD(3) - t201;
t160 = (-pkin(3) * t261 * t267 - pkin(6) * qJDD(1)) * t260 + t171;
t202 = -g(3) * t260 + t312 * t261;
t293 = qJD(1) * t261;
t173 = t219 * t293 + t202;
t288 = qJDD(1) * t261;
t162 = -pkin(3) * t300 - pkin(6) * t288 + t173;
t262 = sin(qJ(4));
t264 = cos(qJ(4));
t158 = t262 * t160 + t264 * t162;
t279 = t260 * t262 + t261 * t264;
t280 = t260 * t264 - t261 * t262;
t217 = t280 * qJD(1);
t291 = t217 * qJD(4);
t199 = t279 * qJDD(1) + t291;
t206 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t217;
t216 = t279 * qJD(1);
t185 = pkin(4) * t216 - qJ(5) * t217;
t266 = qJD(4) ^ 2;
t149 = -pkin(4) * t266 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t185 * t216 + t158;
t207 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t217;
t287 = m(6) * t149 + qJDD(4) * mrSges(6,3) + qJD(4) * t207;
t186 = mrSges(6,1) * t216 - mrSges(6,3) * t217;
t298 = -mrSges(5,1) * t216 - mrSges(5,2) * t217 - t186;
t307 = -mrSges(5,3) - mrSges(6,2);
t139 = m(5) * t158 - qJDD(4) * mrSges(5,2) - qJD(4) * t206 + t307 * t199 + t298 * t216 + t287;
t157 = t160 * t264 - t162 * t262;
t292 = t216 * qJD(4);
t200 = t280 * qJDD(1) - t292;
t205 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t216;
t152 = -qJDD(4) * pkin(4) - qJ(5) * t266 + t185 * t217 + qJDD(5) - t157;
t208 = -mrSges(6,2) * t216 + qJD(4) * mrSges(6,3);
t283 = -m(6) * t152 + qJDD(4) * mrSges(6,1) + qJD(4) * t208;
t140 = m(5) * t157 + qJDD(4) * mrSges(5,1) + qJD(4) * t205 + t307 * t200 + t298 * t217 + t283;
t134 = t264 * t139 - t262 * t140;
t220 = (-mrSges(4,1) * t261 - mrSges(4,3) * t260) * qJD(1);
t132 = m(4) * t173 + mrSges(4,2) * t288 + t220 * t293 + t134;
t281 = Ifges(4,5) * t260 - Ifges(4,3) * t261;
t222 = t281 * qJD(1);
t226 = (Ifges(4,1) * t260 - Ifges(4,5) * t261) * qJD(1);
t133 = t262 * t139 + t264 * t140;
t178 = Ifges(5,4) * t217 - Ifges(5,2) * t216 + Ifges(5,6) * qJD(4);
t180 = Ifges(5,1) * t217 - Ifges(5,4) * t216 + Ifges(5,5) * qJD(4);
t175 = Ifges(6,5) * t217 + Ifges(6,6) * qJD(4) + Ifges(6,3) * t216;
t179 = Ifges(6,1) * t217 + Ifges(6,4) * qJD(4) + Ifges(6,5) * t216;
t278 = -mrSges(6,1) * t152 + mrSges(6,3) * t149 + Ifges(6,4) * t200 + Ifges(6,2) * qJDD(4) + Ifges(6,6) * t199 - t217 * t175 + t216 * t179;
t271 = qJ(5) * (-mrSges(6,2) * t199 - t186 * t216 + t287) + pkin(4) * (-mrSges(6,2) * t200 - t186 * t217 + t283) - mrSges(5,2) * t158 + mrSges(5,1) * t157 + t216 * t180 + t217 * t178 - Ifges(5,6) * t199 + Ifges(5,5) * t200 + Ifges(5,3) * qJDD(4) + t278;
t289 = qJDD(1) * t260;
t269 = -mrSges(4,1) * t171 + mrSges(4,3) * t173 + Ifges(4,4) * t289 - pkin(3) * t133 - t271;
t276 = -m(4) * t171 - t133;
t303 = Ifges(3,1) * t260;
t310 = ((t222 - (Ifges(3,4) * t260 + Ifges(3,2) * t261) * qJD(1)) * t260 + (t226 + (Ifges(3,4) * t261 + t303) * qJD(1)) * t261) * qJD(1) - mrSges(3,1) * t201 + mrSges(3,2) * t202 - pkin(2) * ((-qJDD(1) * mrSges(4,2) - qJD(1) * t220) * t260 + t276) - qJ(3) * t132 - t269;
t229 = t263 * g(1) - t265 * g(2);
t215 = -qJDD(1) * pkin(1) - t267 * qJ(2) + qJDD(2) - t229;
t198 = -pkin(2) * t288 - qJ(3) * t289 - 0.2e1 * qJD(3) * t294 + t215;
t302 = Ifges(3,5) * t260;
t309 = (Ifges(3,6) - Ifges(4,6)) * t261 + t302;
t306 = Ifges(3,4) - Ifges(4,5);
t304 = mrSges(3,2) * t260;
t177 = Ifges(6,4) * t217 + Ifges(6,2) * qJD(4) + Ifges(6,6) * t216;
t299 = -Ifges(5,5) * t217 + Ifges(5,6) * t216 - Ifges(5,3) * qJD(4) - t177;
t221 = (-mrSges(3,1) * t261 + t304) * qJD(1);
t126 = m(3) * t201 + ((-mrSges(4,2) - mrSges(3,3)) * qJDD(1) + (-t220 - t221) * qJD(1)) * t260 + t276;
t127 = m(3) * t202 + (qJDD(1) * mrSges(3,3) + qJD(1) * t221) * t261 + t132;
t284 = -t126 * t260 + t261 * t127;
t168 = pkin(3) * t288 + (-t250 - t251) * t267 * pkin(6) - t198;
t154 = -0.2e1 * qJD(5) * t217 + (-t200 + t292) * qJ(5) + (t199 + t291) * pkin(4) + t168;
t282 = -mrSges(6,1) * t154 + mrSges(6,2) * t149;
t142 = m(6) * t154 + t199 * mrSges(6,1) - t200 * mrSges(6,3) - t217 * t207 + t216 * t208;
t275 = mrSges(6,2) * t152 - mrSges(6,3) * t154 + Ifges(6,1) * t200 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t199 + qJD(4) * t175;
t141 = -m(5) * t168 - t199 * mrSges(5,1) - t200 * mrSges(5,2) - t216 * t205 - t217 * t206 - t142;
t137 = m(4) * t198 - mrSges(4,1) * t288 - t311 * mrSges(4,2) - mrSges(4,3) * t289 + t141;
t223 = (Ifges(3,6) * t261 + t302) * qJD(1);
t224 = (Ifges(4,4) * t260 - Ifges(4,6) * t261) * qJD(1);
t129 = -mrSges(5,1) * t168 + mrSges(5,3) * t158 - pkin(4) * t142 + t299 * t217 + (Ifges(5,4) - Ifges(6,5)) * t200 + (-Ifges(5,2) - Ifges(6,3)) * t199 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t179 + t180) * qJD(4) + t282;
t130 = mrSges(5,2) * t168 - mrSges(5,3) * t157 + Ifges(5,1) * t200 - Ifges(5,4) * t199 + Ifges(5,5) * qJDD(4) - qJ(5) * t142 - qJD(4) * t178 + t299 * t216 + t275;
t272 = mrSges(4,1) * t198 - mrSges(4,2) * t173 + pkin(3) * t141 + pkin(6) * t134 + t129 * t264 + t130 * t262;
t119 = -mrSges(3,1) * t215 + mrSges(3,3) * t202 - pkin(2) * t137 + (-t223 - t224) * t294 + ((Ifges(3,2) + Ifges(4,3)) * t261 + t306 * t260) * qJDD(1) - t272;
t273 = mrSges(4,2) * t171 - mrSges(4,3) * t198 + Ifges(4,1) * t289 - pkin(6) * t133 - t129 * t262 + t264 * t130 + t224 * t293;
t121 = t223 * t293 + mrSges(3,2) * t215 - mrSges(3,3) * t201 - qJ(3) * t137 + (t306 * t261 + t303) * qJDD(1) + t273;
t270 = -m(3) * t215 + mrSges(3,1) * t288 + t311 * mrSges(3,3) - t137;
t274 = -mrSges(2,2) * t230 + qJ(2) * t284 + t261 * t119 + t260 * t121 + pkin(1) * (-mrSges(3,2) * t289 + t270) + mrSges(2,1) * t229 + Ifges(2,3) * qJDD(1);
t135 = (mrSges(2,1) - t304) * qJDD(1) + m(2) * t229 - mrSges(2,2) * t267 + t270;
t124 = t126 * t261 + t127 * t260;
t122 = m(2) * t230 - mrSges(2,1) * t267 - qJDD(1) * mrSges(2,2) + t284;
t117 = -pkin(1) * t124 + mrSges(2,1) * g(3) + mrSges(2,3) * t230 + (Ifges(2,6) - t309) * qJDD(1) + t267 * Ifges(2,5) + t310;
t116 = -mrSges(2,2) * g(3) - mrSges(2,3) * t229 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t267 - qJ(2) * t124 - t119 * t260 + t121 * t261;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t265 * t116 - t263 * t117 - pkin(5) * (t122 * t263 + t135 * t265), t116, t121, -Ifges(4,5) * t288 + t273, t130, -t177 * t216 + t275; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t263 * t116 + t265 * t117 + pkin(5) * (t122 * t265 - t135 * t263), t117, t119, -Ifges(4,6) * t288 + (-t260 * t222 - t261 * t226) * qJD(1) + t269, t129, t278; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t274, t274, t309 * qJDD(1) - t310, t281 * qJDD(1) + t224 * t294 + t272, t271, Ifges(6,5) * t200 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t199 - qJD(4) * t179 + t177 * t217 - t282;];
m_new = t1;
