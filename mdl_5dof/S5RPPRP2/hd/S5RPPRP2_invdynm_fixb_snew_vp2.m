% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 17:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:49:08
% EndTime: 2019-12-31 17:49:12
% DurationCPUTime: 2.68s
% Computational Cost: add. (31085->244), mult. (67202->296), div. (0->0), fcn. (40544->8), ass. (0->103)
t234 = qJD(1) ^ 2;
t226 = sin(pkin(8));
t262 = qJD(1) * t226;
t228 = cos(pkin(8));
t272 = qJD(1) * t228;
t231 = sin(qJ(1));
t232 = cos(qJ(1));
t204 = t231 * g(1) - t232 * g(2);
t201 = qJDD(1) * pkin(1) + t204;
t205 = -t232 * g(1) - t231 * g(2);
t202 = -t234 * pkin(1) + t205;
t227 = sin(pkin(7));
t229 = cos(pkin(7));
t185 = t227 * t201 + t229 * t202;
t165 = -t234 * pkin(2) + qJDD(1) * qJ(3) + t185;
t225 = -g(3) + qJDD(2);
t259 = qJD(1) * qJD(3);
t263 = t228 * t225 - 0.2e1 * t226 * t259;
t269 = pkin(3) * t228;
t149 = (-pkin(6) * qJDD(1) + t234 * t269 - t165) * t226 + t263;
t153 = t226 * t225 + (t165 + 0.2e1 * t259) * t228;
t257 = qJDD(1) * t228;
t217 = t228 ^ 2;
t266 = t217 * t234;
t150 = -pkin(3) * t266 + pkin(6) * t257 + t153;
t230 = sin(qJ(4));
t270 = cos(qJ(4));
t146 = t230 * t149 + t270 * t150;
t255 = t228 * t270;
t258 = qJDD(1) * t226;
t243 = t270 * t226 + t228 * t230;
t193 = t243 * qJD(1);
t260 = t193 * qJD(4);
t181 = -qJDD(1) * t255 + t230 * t258 + t260;
t189 = qJD(4) * mrSges(5,1) - t193 * mrSges(5,3);
t192 = -qJD(1) * t255 + t230 * t262;
t169 = t192 * pkin(4) - t193 * qJ(5);
t233 = qJD(4) ^ 2;
t139 = -t233 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t192 * t169 + t146;
t190 = -qJD(4) * mrSges(6,1) + t193 * mrSges(6,2);
t256 = m(6) * t139 + qJDD(4) * mrSges(6,3) + qJD(4) * t190;
t170 = t192 * mrSges(6,1) - t193 * mrSges(6,3);
t264 = -t192 * mrSges(5,1) - t193 * mrSges(5,2) - t170;
t268 = -mrSges(5,3) - mrSges(6,2);
t130 = m(5) * t146 - qJDD(4) * mrSges(5,2) - qJD(4) * t189 + t268 * t181 + t264 * t192 + t256;
t145 = t270 * t149 - t230 * t150;
t261 = t192 * qJD(4);
t182 = t243 * qJDD(1) - t261;
t188 = -qJD(4) * mrSges(5,2) - t192 * mrSges(5,3);
t141 = -qJDD(4) * pkin(4) - t233 * qJ(5) + t193 * t169 + qJDD(5) - t145;
t191 = -t192 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t251 = -m(6) * t141 + qJDD(4) * mrSges(6,1) + qJD(4) * t191;
t131 = m(5) * t145 + qJDD(4) * mrSges(5,1) + qJD(4) * t188 + t268 * t182 + t264 * t193 + t251;
t123 = t230 * t130 + t270 * t131;
t152 = -t226 * t165 + t263;
t162 = Ifges(5,4) * t193 - Ifges(5,2) * t192 + Ifges(5,6) * qJD(4);
t164 = Ifges(5,1) * t193 - Ifges(5,4) * t192 + Ifges(5,5) * qJD(4);
t159 = Ifges(6,5) * t193 + Ifges(6,6) * qJD(4) + Ifges(6,3) * t192;
t163 = Ifges(6,1) * t193 + Ifges(6,4) * qJD(4) + Ifges(6,5) * t192;
t240 = mrSges(6,1) * t141 - mrSges(6,3) * t139 - Ifges(6,4) * t182 - Ifges(6,2) * qJDD(4) - Ifges(6,6) * t181 + t193 * t159 - t192 * t163;
t236 = mrSges(5,2) * t146 - t192 * t164 - qJ(5) * (-t181 * mrSges(6,2) - t192 * t170 + t256) - pkin(4) * (-t182 * mrSges(6,2) - t193 * t170 + t251) - mrSges(5,1) * t145 - t193 * t162 + Ifges(5,6) * t181 - Ifges(5,5) * t182 - Ifges(5,3) * qJDD(4) + t240;
t248 = Ifges(4,4) * t226 + Ifges(4,2) * t228;
t249 = Ifges(4,1) * t226 + Ifges(4,4) * t228;
t271 = -mrSges(4,1) * t152 + mrSges(4,2) * t153 - pkin(3) * t123 - (t248 * t262 - t249 * t272) * qJD(1) + t236;
t267 = mrSges(4,2) * t226;
t244 = mrSges(4,3) * qJDD(1) + t234 * (-mrSges(4,1) * t228 + t267);
t121 = m(4) * t152 - t244 * t226 + t123;
t252 = t270 * t130 - t230 * t131;
t122 = m(4) * t153 + t244 * t228 + t252;
t253 = -t226 * t121 + t228 * t122;
t113 = m(3) * t185 - t234 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t253;
t184 = t229 * t201 - t227 * t202;
t246 = qJDD(3) - t184;
t158 = -qJDD(1) * pkin(2) - t234 * qJ(3) + t246;
t216 = t226 ^ 2;
t151 = (-pkin(2) - t269) * qJDD(1) + (-qJ(3) + (-t216 - t217) * pkin(6)) * t234 + t246;
t143 = -0.2e1 * qJD(5) * t193 + (-t182 + t261) * qJ(5) + (t181 + t260) * pkin(4) + t151;
t132 = m(6) * t143 + t181 * mrSges(6,1) - t182 * mrSges(6,3) - t193 * t190 + t192 * t191;
t238 = m(5) * t151 + t181 * mrSges(5,1) + t182 * mrSges(5,2) + t192 * t188 + t193 * t189 + t132;
t237 = -m(4) * t158 + mrSges(4,1) * t257 - t238 + (t216 * t234 + t266) * mrSges(4,3);
t125 = -t234 * mrSges(3,2) + m(3) * t184 + t237 + (mrSges(3,1) - t267) * qJDD(1);
t110 = t227 * t113 + t229 * t125;
t115 = t228 * t121 + t226 * t122;
t161 = Ifges(6,4) * t193 + Ifges(6,2) * qJD(4) + Ifges(6,6) * t192;
t265 = -Ifges(5,5) * t193 + Ifges(5,6) * t192 - Ifges(5,3) * qJD(4) - t161;
t254 = t229 * t113 - t227 * t125;
t250 = -mrSges(6,1) * t143 + mrSges(6,2) * t139;
t247 = Ifges(4,5) * t226 + Ifges(4,6) * t228;
t242 = mrSges(6,2) * t141 - mrSges(6,3) * t143 + Ifges(6,1) * t182 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t181 + qJD(4) * t159;
t116 = -mrSges(5,1) * t151 + mrSges(5,3) * t146 - pkin(4) * t132 + t265 * t193 + (Ifges(5,4) - Ifges(6,5)) * t182 + (-Ifges(5,2) - Ifges(6,3)) * t181 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t163 + t164) * qJD(4) + t250;
t117 = mrSges(5,2) * t151 - mrSges(5,3) * t145 + Ifges(5,1) * t182 - Ifges(5,4) * t181 + Ifges(5,5) * qJDD(4) - qJ(5) * t132 - qJD(4) * t162 + t265 * t192 + t242;
t198 = t247 * qJD(1);
t104 = -mrSges(4,1) * t158 + mrSges(4,3) * t153 - pkin(3) * t238 + pkin(6) * t252 + t248 * qJDD(1) + t270 * t116 + t230 * t117 - t198 * t262;
t106 = mrSges(4,2) * t158 - mrSges(4,3) * t152 - pkin(6) * t123 + t249 * qJDD(1) - t230 * t116 + t270 * t117 + t198 * t272;
t241 = -mrSges(3,2) * t185 + qJ(3) * t253 + t228 * t104 + t226 * t106 + pkin(2) * (-mrSges(4,2) * t258 + t237) + mrSges(3,1) * t184 + Ifges(3,3) * qJDD(1);
t239 = mrSges(2,1) * t204 - mrSges(2,2) * t205 + Ifges(2,3) * qJDD(1) + pkin(1) * t110 + t241;
t108 = m(2) * t205 - t234 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t254;
t107 = m(2) * t204 + qJDD(1) * mrSges(2,1) - t234 * mrSges(2,2) + t110;
t102 = t234 * Ifges(3,5) - mrSges(3,1) * t225 + mrSges(3,3) * t185 - pkin(2) * t115 + (Ifges(3,6) - t247) * qJDD(1) + t271;
t101 = mrSges(3,2) * t225 - mrSges(3,3) * t184 + Ifges(3,5) * qJDD(1) - t234 * Ifges(3,6) - qJ(3) * t115 - t226 * t104 + t228 * t106;
t100 = -mrSges(2,2) * g(3) - mrSges(2,3) * t204 + Ifges(2,5) * qJDD(1) - t234 * Ifges(2,6) - qJ(2) * t110 + t229 * t101 - t227 * t102;
t99 = Ifges(2,6) * qJDD(1) + t234 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t205 + t227 * t101 + t229 * t102 - pkin(1) * (m(3) * t225 + t115) + qJ(2) * t254;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t232 * t100 - t231 * t99 - pkin(5) * (t232 * t107 + t231 * t108), t100, t101, t106, t117, -t192 * t161 + t242; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t231 * t100 + t232 * t99 + pkin(5) * (-t231 * t107 + t232 * t108), t99, t102, t104, t116, -t240; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t239, t239, t241, t247 * qJDD(1) - t271, -t236, Ifges(6,5) * t182 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t181 - qJD(4) * t163 + t193 * t161 - t250;];
m_new = t1;
