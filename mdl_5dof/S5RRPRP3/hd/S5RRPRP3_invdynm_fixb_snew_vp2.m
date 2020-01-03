% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:59
% EndTime: 2019-12-31 19:51:04
% DurationCPUTime: 3.14s
% Computational Cost: add. (48775->244), mult. (67202->293), div. (0->0), fcn. (40544->8), ass. (0->102)
t219 = qJD(1) + qJD(2);
t215 = t219 ^ 2;
t226 = sin(pkin(8));
t227 = cos(pkin(8));
t228 = sin(qJ(4));
t270 = cos(qJ(4));
t272 = t226 * t228 - t227 * t270;
t230 = sin(qJ(1));
t232 = cos(qJ(1));
t208 = t230 * g(1) - t232 * g(2);
t201 = qJDD(1) * pkin(1) + t208;
t209 = -t232 * g(1) - t230 * g(2);
t234 = qJD(1) ^ 2;
t202 = -t234 * pkin(1) + t209;
t229 = sin(qJ(2));
t231 = cos(qJ(2));
t185 = t229 * t201 + t231 * t202;
t216 = qJDD(1) + qJDD(2);
t182 = -t215 * pkin(2) + t216 * qJ(3) + t185;
t260 = qJD(3) * t219;
t255 = -t227 * g(3) - 0.2e1 * t226 * t260;
t269 = pkin(3) * t227;
t149 = (-pkin(7) * t216 + t215 * t269 - t182) * t226 + t255;
t155 = -t226 * g(3) + (t182 + 0.2e1 * t260) * t227;
t264 = t216 * t227;
t218 = t227 ^ 2;
t265 = t215 * t218;
t150 = -pkin(3) * t265 + pkin(7) * t264 + t155;
t146 = t228 * t149 + t270 * t150;
t243 = t270 * t226 + t227 * t228;
t193 = t243 * t219;
t258 = t193 * qJD(4);
t180 = t272 * t216 + t258;
t189 = qJD(4) * mrSges(5,1) - t193 * mrSges(5,3);
t192 = t272 * t219;
t176 = t192 * pkin(4) - t193 * qJ(5);
t233 = qJD(4) ^ 2;
t139 = -t233 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t192 * t176 + t146;
t190 = -qJD(4) * mrSges(6,1) + t193 * mrSges(6,2);
t257 = m(6) * t139 + qJDD(4) * mrSges(6,3) + qJD(4) * t190;
t177 = t192 * mrSges(6,1) - t193 * mrSges(6,3);
t261 = -t192 * mrSges(5,1) - t193 * mrSges(5,2) - t177;
t268 = -mrSges(5,3) - mrSges(6,2);
t130 = m(5) * t146 - qJDD(4) * mrSges(5,2) - qJD(4) * t189 + t268 * t180 + t261 * t192 + t257;
t145 = t270 * t149 - t228 * t150;
t259 = t192 * qJD(4);
t181 = t216 * t243 - t259;
t188 = -qJD(4) * mrSges(5,2) - t192 * mrSges(5,3);
t141 = -qJDD(4) * pkin(4) - t233 * qJ(5) + t193 * t176 + qJDD(5) - t145;
t191 = -t192 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t251 = -m(6) * t141 + qJDD(4) * mrSges(6,1) + qJD(4) * t191;
t131 = m(5) * t145 + qJDD(4) * mrSges(5,1) + qJD(4) * t188 + t268 * t181 + t261 * t193 + t251;
t123 = t228 * t130 + t270 * t131;
t154 = -t226 * t182 + t255;
t160 = Ifges(5,4) * t193 - Ifges(5,2) * t192 + Ifges(5,6) * qJD(4);
t162 = Ifges(5,1) * t193 - Ifges(5,4) * t192 + Ifges(5,5) * qJD(4);
t157 = Ifges(6,5) * t193 + Ifges(6,6) * qJD(4) + Ifges(6,3) * t192;
t161 = Ifges(6,1) * t193 + Ifges(6,4) * qJD(4) + Ifges(6,5) * t192;
t240 = mrSges(6,1) * t141 - mrSges(6,3) * t139 - Ifges(6,4) * t181 - Ifges(6,2) * qJDD(4) - Ifges(6,6) * t180 + t193 * t157 - t192 * t161;
t236 = mrSges(5,2) * t146 - t192 * t162 - qJ(5) * (-t180 * mrSges(6,2) - t192 * t177 + t257) - pkin(4) * (-t181 * mrSges(6,2) - t193 * t177 + t251) - mrSges(5,1) * t145 - t193 * t160 + Ifges(5,6) * t180 - Ifges(5,5) * t181 - Ifges(5,3) * qJDD(4) + t240;
t248 = Ifges(4,4) * t226 + Ifges(4,2) * t227;
t249 = Ifges(4,1) * t226 + Ifges(4,4) * t227;
t271 = -mrSges(4,1) * t154 + mrSges(4,2) * t155 - pkin(3) * t123 - (t226 * t248 - t227 * t249) * t215 + t236;
t267 = mrSges(4,2) * t226;
t247 = Ifges(4,5) * t226 + Ifges(4,6) * t227;
t266 = t247 * t215;
t245 = mrSges(4,3) * t216 + (-mrSges(4,1) * t227 + t267) * t215;
t121 = m(4) * t154 - t226 * t245 + t123;
t252 = t270 * t130 - t228 * t131;
t122 = m(4) * t155 + t227 * t245 + t252;
t253 = -t226 * t121 + t227 * t122;
t113 = m(3) * t185 - t215 * mrSges(3,1) - t216 * mrSges(3,2) + t253;
t184 = t231 * t201 - t229 * t202;
t246 = qJDD(3) - t184;
t179 = -t216 * pkin(2) - t215 * qJ(3) + t246;
t217 = t226 ^ 2;
t153 = (-pkin(2) - t269) * t216 + (-qJ(3) + (-t217 - t218) * pkin(7)) * t215 + t246;
t143 = -0.2e1 * qJD(5) * t193 + (-t181 + t259) * qJ(5) + (t180 + t258) * pkin(4) + t153;
t132 = m(6) * t143 + t180 * mrSges(6,1) - t181 * mrSges(6,3) - t193 * t190 + t192 * t191;
t238 = m(5) * t153 + t180 * mrSges(5,1) + t181 * mrSges(5,2) + t192 * t188 + t193 * t189 + t132;
t237 = -m(4) * t179 + mrSges(4,1) * t264 - t238 + (t215 * t217 + t265) * mrSges(4,3);
t125 = -t215 * mrSges(3,2) + m(3) * t184 + t237 + (mrSges(3,1) - t267) * t216;
t110 = t229 * t113 + t231 * t125;
t115 = t227 * t121 + t226 * t122;
t159 = Ifges(6,4) * t193 + Ifges(6,2) * qJD(4) + Ifges(6,6) * t192;
t262 = -Ifges(5,5) * t193 + Ifges(5,6) * t192 - Ifges(5,3) * qJD(4) - t159;
t254 = t231 * t113 - t229 * t125;
t250 = -mrSges(6,1) * t143 + mrSges(6,2) * t139;
t242 = mrSges(6,2) * t141 - mrSges(6,3) * t143 + Ifges(6,1) * t181 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t180 + qJD(4) * t157;
t116 = -mrSges(5,1) * t153 + mrSges(5,3) * t146 - pkin(4) * t132 + t262 * t193 + (Ifges(5,4) - Ifges(6,5)) * t181 + (-Ifges(5,2) - Ifges(6,3)) * t180 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t161 + t162) * qJD(4) + t250;
t117 = mrSges(5,2) * t153 - mrSges(5,3) * t145 + Ifges(5,1) * t181 - Ifges(5,4) * t180 + Ifges(5,5) * qJDD(4) - qJ(5) * t132 - qJD(4) * t160 + t262 * t192 + t242;
t104 = -mrSges(4,1) * t179 + mrSges(4,3) * t155 - pkin(3) * t238 + pkin(7) * t252 + t270 * t116 + t228 * t117 + t248 * t216 - t226 * t266;
t106 = mrSges(4,2) * t179 - mrSges(4,3) * t154 - pkin(7) * t123 - t228 * t116 + t270 * t117 + t249 * t216 + t227 * t266;
t241 = -mrSges(3,2) * t185 + qJ(3) * t253 + t227 * t104 + t226 * t106 + pkin(2) * (-t216 * t267 + t237) + mrSges(3,1) * t184 + Ifges(3,3) * t216;
t239 = mrSges(2,1) * t208 - mrSges(2,2) * t209 + Ifges(2,3) * qJDD(1) + pkin(1) * t110 + t241;
t108 = m(2) * t209 - t234 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t254;
t107 = m(2) * t208 + qJDD(1) * mrSges(2,1) - t234 * mrSges(2,2) + t110;
t102 = (Ifges(3,6) - t247) * t216 + t215 * Ifges(3,5) + mrSges(3,3) * t185 + mrSges(3,1) * g(3) - pkin(2) * t115 + t271;
t101 = -mrSges(3,2) * g(3) - mrSges(3,3) * t184 + Ifges(3,5) * t216 - t215 * Ifges(3,6) - qJ(3) * t115 - t226 * t104 + t227 * t106;
t100 = -mrSges(2,2) * g(3) - mrSges(2,3) * t208 + Ifges(2,5) * qJDD(1) - t234 * Ifges(2,6) - pkin(6) * t110 + t231 * t101 - t229 * t102;
t99 = Ifges(2,6) * qJDD(1) + t234 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t209 + t229 * t101 + t231 * t102 - pkin(1) * (-m(3) * g(3) + t115) + pkin(6) * t254;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t232 * t100 - t230 * t99 - pkin(5) * (t232 * t107 + t230 * t108), t100, t101, t106, t117, -t192 * t159 + t242; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t230 * t100 + t232 * t99 + pkin(5) * (-t230 * t107 + t232 * t108), t99, t102, t104, t116, -t240; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t239, t239, t241, t247 * t216 - t271, -t236, Ifges(6,5) * t181 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t180 - qJD(4) * t161 + t193 * t159 - t250;];
m_new = t1;
