% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:36
% EndTime: 2019-12-05 15:30:42
% DurationCPUTime: 2.62s
% Computational Cost: add. (22810->230), mult. (48985->298), div. (0->0), fcn. (29382->8), ass. (0->102)
t217 = sin(pkin(7));
t219 = cos(pkin(7));
t198 = t217 * g(1) - t219 * g(2);
t199 = -t219 * g(1) - t217 * g(2);
t221 = sin(qJ(2));
t223 = cos(qJ(2));
t168 = t221 * t198 + t223 * t199;
t224 = qJD(2) ^ 2;
t266 = -t224 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t168;
t215 = -g(3) + qJDD(1);
t216 = sin(pkin(8));
t218 = cos(pkin(8));
t148 = t218 * t215 - t266 * t216;
t258 = qJD(2) * t216;
t256 = t218 * qJD(2);
t149 = t216 * t215 + t266 * t218;
t241 = -pkin(3) * t218 - pkin(6) * t216;
t196 = t241 * qJD(2);
t147 = t196 * t256 + t149;
t167 = t223 * t198 - t221 * t199;
t234 = -t224 * qJ(3) + qJDD(3) - t167;
t152 = (-pkin(2) + t241) * qJDD(2) + t234;
t220 = sin(qJ(4));
t222 = cos(qJ(4));
t144 = t222 * t147 + t220 * t152;
t146 = t196 * t258 - t148;
t205 = qJD(4) - t256;
t166 = Ifges(5,5) * t205 + (Ifges(5,1) * t222 - Ifges(5,4) * t220) * t258;
t250 = t220 * t258;
t179 = -t205 * mrSges(6,2) - mrSges(6,3) * t250;
t249 = t222 * t258;
t182 = t205 * mrSges(6,1) - mrSges(6,3) * t249;
t184 = (mrSges(6,1) * t220 + mrSges(6,2) * t222) * t258;
t254 = qJD(2) * qJD(4);
t186 = (-qJDD(2) * t220 - t222 * t254) * t216;
t187 = (qJDD(2) * t222 - t220 * t254) * t216;
t253 = t218 * qJDD(2);
t204 = qJDD(4) - t253;
t181 = t205 * pkin(4) - qJ(5) * t249;
t242 = -0.2e1 * qJD(5) * t258;
t262 = t216 ^ 2 * t224;
t252 = t220 ^ 2 * t262;
t140 = -pkin(4) * t252 + t186 * qJ(5) - t205 * t181 + t220 * t242 + t144;
t142 = -t186 * pkin(4) - qJ(5) * t252 + t181 * t249 + qJDD(5) + t146;
t165 = Ifges(6,5) * t205 + (Ifges(6,1) * t222 - Ifges(6,4) * t220) * t258;
t233 = -mrSges(6,1) * t142 + mrSges(6,3) * t140 + Ifges(6,4) * t187 + Ifges(6,2) * t186 + Ifges(6,6) * t204 + t205 * t165;
t246 = -m(6) * t142 + t186 * mrSges(6,1);
t161 = Ifges(6,3) * t205 + (Ifges(6,5) * t222 - Ifges(6,6) * t220) * t258;
t260 = -t161 - Ifges(5,3) * t205 - (Ifges(5,5) * t222 - Ifges(5,6) * t220) * t258;
t261 = m(6) * t140 + t186 * mrSges(6,3);
t118 = Ifges(5,4) * t187 + Ifges(5,2) * t186 + Ifges(5,6) * t204 + t205 * t166 - mrSges(5,1) * t146 + mrSges(5,3) * t144 - pkin(4) * (t187 * mrSges(6,2) - t246) + qJ(5) * (-t204 * mrSges(6,2) - t205 * t182 + t261) + ((-pkin(4) * t179 - qJ(5) * t184) * t220 + (-pkin(4) * t182 + t260) * t222) * t258 + t233;
t151 = t222 * t152;
t138 = t222 * t242 + t204 * pkin(4) - t187 * qJ(5) + t151 + (-pkin(4) * t222 * t262 - qJ(5) * t205 * t258 - t147) * t220;
t251 = m(6) * t138 + t204 * mrSges(6,1) + t205 * t179;
t134 = -t187 * mrSges(6,3) - t184 * t249 + t251;
t143 = -t220 * t147 + t151;
t163 = Ifges(6,6) * t205 + (Ifges(6,4) * t222 - Ifges(6,2) * t220) * t258;
t164 = Ifges(5,6) * t205 + (Ifges(5,4) * t222 - Ifges(5,2) * t220) * t258;
t235 = mrSges(6,2) * t142 - mrSges(6,3) * t138 + Ifges(6,1) * t187 + Ifges(6,4) * t186 + Ifges(6,5) * t204;
t125 = mrSges(5,2) * t146 - mrSges(5,3) * t143 + Ifges(5,1) * t187 + Ifges(5,4) * t186 + Ifges(5,5) * t204 - qJ(5) * t134 + (-t163 - t164) * t205 + t260 * t250 + t235;
t180 = -t205 * mrSges(5,2) - mrSges(5,3) * t250;
t240 = (-t184 - (mrSges(5,1) * t220 + mrSges(5,2) * t222) * t258) * t258;
t129 = m(5) * t143 + t204 * mrSges(5,1) + t205 * t180 + (-mrSges(5,3) - mrSges(6,3)) * t187 + t222 * t240 + t251;
t259 = -t205 * mrSges(5,1) + mrSges(5,3) * t249 - t182;
t264 = -mrSges(5,2) - mrSges(6,2);
t130 = m(5) * t144 + t186 * mrSges(5,3) + t264 * t204 + t259 * t205 + t220 * t240 + t261;
t127 = -t220 * t129 + t222 * t130;
t229 = -m(5) * t146 + t186 * mrSges(5,1) + t264 * t187 + t246;
t232 = t259 * t222 + (-t179 - t180) * t220;
t239 = Ifges(4,1) * t216 + Ifges(4,4) * t218;
t265 = -((Ifges(4,4) * t216 + Ifges(4,2) * t218) * t258 - t239 * t256) * qJD(2) - mrSges(4,1) * t148 + mrSges(4,2) * t149 - pkin(3) * (t232 * t258 + t229) - pkin(6) * t127 - t222 * t118 - t220 * t125;
t263 = mrSges(4,2) * t216;
t192 = (-mrSges(4,1) * t218 + t263) * qJD(2);
t257 = qJDD(2) * mrSges(4,3);
t123 = m(4) * t149 + (qJD(2) * t192 + t257) * t218 + t127;
t132 = m(4) * t148 + (-t257 + (-t192 + t232) * qJD(2)) * t216 + t229;
t244 = t218 * t123 - t216 * t132;
t115 = m(3) * t168 - t224 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t244;
t126 = t222 * t129 + t220 * t130;
t155 = -qJDD(2) * pkin(2) + t234;
t228 = -m(4) * t155 + mrSges(4,1) * t253 - t126 + (t218 ^ 2 * t224 + t262) * mrSges(4,3);
t120 = m(3) * t167 - t224 * mrSges(3,2) + (mrSges(3,1) - t263) * qJDD(2) + t228;
t110 = t221 * t115 + t223 * t120;
t117 = t216 * t123 + t218 * t132;
t248 = t161 * t258;
t245 = t223 * t115 - t221 * t120;
t238 = Ifges(4,5) * t216 + Ifges(4,6) * t218;
t237 = t164 * t222 + t166 * t220;
t193 = t238 * qJD(2);
t106 = mrSges(4,2) * t155 - mrSges(4,3) * t148 - pkin(6) * t126 + t239 * qJDD(2) - t220 * t118 + t222 * t125 + t193 * t256;
t230 = -mrSges(6,1) * t138 + mrSges(6,2) * t140 - Ifges(6,5) * t187 - Ifges(6,6) * t186 - Ifges(6,3) * t204 - t163 * t249 - t165 * t250;
t225 = mrSges(5,1) * t143 - mrSges(5,2) * t144 + Ifges(5,5) * t187 + Ifges(5,6) * t186 + Ifges(5,3) * t204 + pkin(4) * t134 - t230;
t112 = -t225 + (Ifges(4,4) * qJDD(2) + (-t193 - t237) * qJD(2)) * t216 + mrSges(4,3) * t149 - mrSges(4,1) * t155 - pkin(3) * t126 + Ifges(4,2) * t253;
t231 = -mrSges(3,2) * t168 + qJ(3) * t244 + t216 * t106 + t218 * t112 + pkin(2) * (-qJDD(2) * t263 + t228) + mrSges(3,1) * t167 + Ifges(3,3) * qJDD(2);
t227 = mrSges(2,1) * t198 - mrSges(2,2) * t199 + pkin(1) * t110 + t231;
t108 = m(2) * t199 + t245;
t107 = m(2) * t198 + t110;
t104 = -mrSges(3,1) * t215 + mrSges(3,3) * t168 + t224 * Ifges(3,5) - pkin(2) * t117 + (Ifges(3,6) - t238) * qJDD(2) + t265;
t103 = mrSges(3,2) * t215 - mrSges(3,3) * t167 + Ifges(3,5) * qJDD(2) - t224 * Ifges(3,6) - qJ(3) * t117 + t218 * t106 - t216 * t112;
t102 = mrSges(2,2) * t215 - mrSges(2,3) * t198 - pkin(5) * t110 + t223 * t103 - t221 * t104;
t101 = -mrSges(2,1) * t215 + mrSges(2,3) * t199 + t221 * t103 + t223 * t104 - pkin(1) * (m(3) * t215 + t117) + pkin(5) * t245;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t219 * t102 - t217 * t101 - qJ(1) * (t219 * t107 + t217 * t108), t102, t103, t106, t125, -t205 * t163 - t220 * t248 + t235; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t217 * t102 + t219 * t101 + qJ(1) * (-t217 * t107 + t219 * t108), t101, t104, t112, t118, -t222 * t248 + t233; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t227, t227, t231, t238 * qJDD(2) - t265, t237 * t258 + t225, -t230;];
m_new = t1;
