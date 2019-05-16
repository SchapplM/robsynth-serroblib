% Calculate vector of cutting torques with Newton-Euler for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPPRR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:45:50
% EndTime: 2019-05-05 13:45:56
% DurationCPUTime: 2.78s
% Computational Cost: add. (38727->256), mult. (65705->302), div. (0->0), fcn. (28493->8), ass. (0->107)
t231 = qJD(1) ^ 2;
t225 = sin(qJ(1));
t228 = cos(qJ(1));
t199 = -g(1) * t228 - g(2) * t225;
t245 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t199;
t266 = -pkin(1) - pkin(2);
t175 = t231 * t266 + t245;
t198 = g(1) * t225 - t228 * g(2);
t244 = -qJ(2) * t231 + qJDD(2) - t198;
t178 = qJDD(1) * t266 + t244;
t221 = sin(pkin(9));
t222 = cos(pkin(9));
t160 = -t221 * t175 + t178 * t222;
t158 = qJDD(1) * pkin(3) - qJ(4) * t231 + qJDD(4) - t160;
t155 = qJDD(1) * pkin(7) + t158;
t215 = g(3) + qJDD(3);
t224 = sin(qJ(5));
t227 = cos(qJ(5));
t150 = t224 * t155 + t227 * t215;
t191 = (-mrSges(6,1) * t224 - mrSges(6,2) * t227) * qJD(1);
t256 = qJD(1) * qJD(5);
t253 = t227 * t256;
t193 = qJDD(1) * t224 + t253;
t257 = qJD(1) * t227;
t197 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t257;
t161 = t222 * t175 + t221 * t178;
t269 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t161;
t153 = (-pkin(3) - pkin(7)) * t231 + t269;
t254 = t224 * t256;
t194 = -qJDD(1) * t227 + t254;
t145 = (-t194 - t254) * pkin(8) + (-t193 - t253) * pkin(5) + t153;
t192 = (-pkin(5) * t224 + pkin(8) * t227) * qJD(1);
t230 = qJD(5) ^ 2;
t258 = qJD(1) * t224;
t147 = -pkin(5) * t230 + qJDD(5) * pkin(8) + t192 * t258 + t150;
t223 = sin(qJ(6));
t226 = cos(qJ(6));
t143 = t145 * t226 - t147 * t223;
t189 = qJD(5) * t226 + t223 * t257;
t168 = qJD(6) * t189 + qJDD(5) * t223 + t194 * t226;
t190 = qJD(5) * t223 - t226 * t257;
t169 = -mrSges(7,1) * t189 + mrSges(7,2) * t190;
t200 = qJD(6) - t258;
t173 = -mrSges(7,2) * t200 + mrSges(7,3) * t189;
t188 = qJDD(6) - t193;
t139 = m(7) * t143 + mrSges(7,1) * t188 - mrSges(7,3) * t168 - t169 * t190 + t173 * t200;
t144 = t145 * t223 + t147 * t226;
t167 = -qJD(6) * t190 + qJDD(5) * t226 - t194 * t223;
t174 = mrSges(7,1) * t200 - mrSges(7,3) * t190;
t140 = m(7) * t144 - mrSges(7,2) * t188 + mrSges(7,3) * t167 + t169 * t189 - t174 * t200;
t252 = -t139 * t223 + t226 * t140;
t126 = m(6) * t150 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t193 - qJD(5) * t197 + t191 * t258 + t252;
t259 = t215 * t224;
t149 = t155 * t227 - t259;
t196 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t258;
t146 = -qJDD(5) * pkin(5) - pkin(8) * t230 + t259 + (-qJD(1) * t192 - t155) * t227;
t243 = -m(7) * t146 + t167 * mrSges(7,1) - mrSges(7,2) * t168 + t189 * t173 - t174 * t190;
t135 = m(6) * t149 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t194 + qJD(5) * t196 + t191 * t257 + t243;
t119 = t126 * t224 + t135 * t227;
t162 = Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t200;
t164 = Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t200;
t132 = -mrSges(7,1) * t146 + mrSges(7,3) * t144 + Ifges(7,4) * t168 + Ifges(7,2) * t167 + Ifges(7,6) * t188 - t162 * t190 + t164 * t200;
t163 = Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t200;
t133 = mrSges(7,2) * t146 - mrSges(7,3) * t143 + Ifges(7,1) * t168 + Ifges(7,4) * t167 + Ifges(7,5) * t188 + t162 * t189 - t163 * t200;
t183 = (Ifges(6,6) * qJD(5)) + (-Ifges(6,4) * t227 + Ifges(6,2) * t224) * qJD(1);
t184 = (Ifges(6,5) * qJD(5)) + (-Ifges(6,1) * t227 + Ifges(6,4) * t224) * qJD(1);
t234 = mrSges(6,1) * t149 - mrSges(6,2) * t150 + Ifges(6,5) * t194 + Ifges(6,6) * t193 + Ifges(6,3) * qJDD(5) + pkin(5) * t243 + pkin(8) * t252 - (t227 * t183 + t224 * t184) * qJD(1) + t226 * t132 + t223 * t133;
t270 = -mrSges(5,1) * t158 - pkin(4) * t119 - t234;
t128 = t226 * t139 + t223 * t140;
t267 = -m(6) * t153 + mrSges(6,1) * t193 - t194 * mrSges(6,2) + (t196 * t224 + t197 * t227) * qJD(1) - t128;
t264 = mrSges(2,1) + mrSges(3,1);
t263 = mrSges(4,2) - mrSges(5,3);
t262 = Ifges(5,4) - Ifges(4,5);
t261 = Ifges(5,5) - Ifges(4,6);
t120 = t227 * t126 - t224 * t135;
t114 = -m(5) * t158 + qJDD(1) * mrSges(5,2) + t231 * mrSges(5,3) - t119;
t113 = m(4) * t160 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t231 + t114;
t156 = pkin(3) * t231 - t269;
t236 = -m(5) * t156 + t231 * mrSges(5,2) - t267;
t122 = m(4) * t161 - mrSges(4,1) * t231 + qJDD(1) * t263 + t236;
t108 = -t113 * t221 + t222 * t122;
t107 = t113 * t222 + t122 * t221;
t116 = (-m(4) - m(5)) * t215 - t120;
t179 = -pkin(1) * t231 + t245;
t248 = m(3) * t179 + qJDD(1) * mrSges(3,3) + t108;
t182 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t227 + Ifges(6,6) * t224) * qJD(1);
t110 = mrSges(6,2) * t153 - mrSges(6,3) * t149 + Ifges(6,1) * t194 + Ifges(6,4) * t193 + Ifges(6,5) * qJDD(5) - pkin(8) * t128 - qJD(5) * t183 - t132 * t223 + t133 * t226 + t182 * t258;
t235 = mrSges(7,1) * t143 - mrSges(7,2) * t144 + Ifges(7,5) * t168 + Ifges(7,6) * t167 + Ifges(7,3) * t188 + t163 * t190 - t164 * t189;
t112 = -mrSges(6,1) * t153 + mrSges(6,3) * t150 + Ifges(6,4) * t194 + Ifges(6,2) * t193 + Ifges(6,6) * qJDD(5) - pkin(5) * t128 + qJD(5) * t184 + t182 * t257 - t235;
t247 = mrSges(5,2) * t158 - mrSges(5,3) * t156 - Ifges(5,1) * qJDD(1) - pkin(7) * t119 + t227 * t110 - t224 * t112;
t181 = -qJDD(1) * pkin(1) + t244;
t242 = -m(3) * t181 + qJDD(1) * mrSges(3,1) + t231 * mrSges(3,3) - t107;
t117 = m(5) * t215 + t120;
t239 = -mrSges(5,1) * t156 - pkin(4) * t267 - pkin(7) * t120 - t110 * t224 - t112 * t227;
t100 = mrSges(4,3) * t161 - pkin(3) * t117 - t262 * t231 + (-mrSges(4,1) + mrSges(5,2)) * t215 + t261 * qJDD(1) + t239;
t102 = -mrSges(4,3) * t160 - qJ(4) * t117 + qJDD(1) * t262 + t215 * t263 + t231 * t261 - t270;
t241 = mrSges(3,2) * t181 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t231 * Ifges(3,6) - qJ(3) * t107 - t100 * t221 + t222 * t102;
t238 = mrSges(3,2) * t179 - pkin(2) * t116 - qJ(3) * t108 - t100 * t222 - t102 * t221;
t237 = mrSges(4,1) * t160 + pkin(3) * t114 + qJ(4) * (-qJDD(1) * mrSges(5,3) + t236) - mrSges(4,2) * t161 - Ifges(4,3) * qJDD(1) + t247;
t233 = -mrSges(3,1) * t181 + mrSges(3,3) * t179 + Ifges(3,2) * qJDD(1) - pkin(2) * t107 - t237;
t232 = -mrSges(2,2) * t199 + qJ(2) * (-mrSges(3,1) * t231 + t248) + pkin(1) * t242 + mrSges(2,1) * t198 + Ifges(2,3) * qJDD(1) + t233;
t115 = -m(3) * g(3) + t116;
t104 = m(2) * t198 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t231 + t242;
t103 = m(2) * t199 - qJDD(1) * mrSges(2,2) - t231 * t264 + t248;
t99 = -mrSges(2,2) * g(3) - mrSges(2,3) * t198 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t231 - qJ(2) * t115 + t241;
t98 = mrSges(2,3) * t199 - pkin(1) * t115 + (Ifges(3,4) + Ifges(2,5)) * t231 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t264 * g(3) + t238;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t228 * t99 - t225 * t98 - pkin(6) * (t103 * t225 + t104 * t228), t99, t241, t102, t247, t110, t133; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t225 * t99 + t228 * t98 + pkin(6) * (t103 * t228 - t225 * t104), t98, t233, t100, mrSges(5,3) * t215 - Ifges(5,4) * qJDD(1) - t231 * Ifges(5,5) + t270, t112, t132; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t232, t232, -mrSges(3,1) * g(3) - Ifges(3,4) * t231 + Ifges(3,6) * qJDD(1) - t238, t237, -mrSges(5,2) * t215 + Ifges(5,4) * t231 - Ifges(5,5) * qJDD(1) - t239, t234, t235;];
m_new  = t1;
