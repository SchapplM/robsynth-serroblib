% Calculate vector of cutting torques with Newton-Euler for
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-05-05 13:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:50:10
% EndTime: 2019-05-05 13:50:15
% DurationCPUTime: 2.84s
% Computational Cost: add. (42650->261), mult. (71727->300), div. (0->0), fcn. (29517->8), ass. (0->109)
t267 = -2 * qJD(1);
t223 = sin(qJ(1));
t226 = cos(qJ(1));
t199 = -t226 * g(1) - t223 * g(2);
t266 = -qJDD(1) * qJ(2) + (qJD(2) * t267) - t199;
t229 = qJD(1) ^ 2;
t198 = t223 * g(1) - t226 * g(2);
t246 = qJDD(2) - t198;
t261 = -pkin(1) - qJ(3);
t238 = (qJD(3) * t267) + t261 * qJDD(1) + t246;
t171 = (-pkin(3) - qJ(2)) * t229 + t238;
t177 = t261 * t229 + qJDD(3) - t266;
t172 = qJDD(1) * pkin(3) + t177;
t219 = sin(pkin(9));
t220 = cos(pkin(9));
t158 = -t219 * t171 + t220 * t172;
t155 = -qJDD(1) * pkin(4) - (t229 * pkin(7)) - t158;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t254 = qJD(1) * qJD(5);
t249 = t225 * t254;
t193 = t222 * qJDD(1) + t249;
t250 = t222 * t254;
t194 = t225 * qJDD(1) - t250;
t149 = (-t193 - t249) * pkin(8) + (-t194 + t250) * pkin(5) + t155;
t159 = t220 * t171 + t219 * t172;
t156 = -t229 * pkin(4) + qJDD(1) * pkin(7) + t159;
t214 = -g(3) + qJDD(4);
t153 = t225 * t156 + t222 * t214;
t192 = (-pkin(5) * t225 - pkin(8) * t222) * qJD(1);
t228 = qJD(5) ^ 2;
t255 = t225 * qJD(1);
t151 = -t228 * pkin(5) + qJDD(5) * pkin(8) + t192 * t255 + t153;
t221 = sin(qJ(6));
t224 = cos(qJ(6));
t147 = t224 * t149 - t221 * t151;
t256 = qJD(1) * t222;
t189 = t224 * qJD(5) - t221 * t256;
t166 = t189 * qJD(6) + t221 * qJDD(5) + t224 * t193;
t190 = t221 * qJD(5) + t224 * t256;
t170 = -t189 * mrSges(7,1) + t190 * mrSges(7,2);
t200 = qJD(6) - t255;
t178 = -t200 * mrSges(7,2) + t189 * mrSges(7,3);
t188 = qJDD(6) - t194;
t144 = m(7) * t147 + t188 * mrSges(7,1) - t166 * mrSges(7,3) - t190 * t170 + t200 * t178;
t148 = t221 * t149 + t224 * t151;
t165 = -t190 * qJD(6) + t224 * qJDD(5) - t221 * t193;
t179 = t200 * mrSges(7,1) - t190 * mrSges(7,3);
t145 = m(7) * t148 - t188 * mrSges(7,2) + t165 * mrSges(7,3) + t189 * t170 - t200 * t179;
t138 = -t221 * t144 + t224 * t145;
t259 = t225 * t214;
t150 = -qJDD(5) * pkin(5) - t228 * pkin(8) - t259 + (qJD(1) * t192 + t156) * t222;
t160 = Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t200;
t162 = Ifges(7,1) * t190 + Ifges(7,4) * t189 + Ifges(7,5) * t200;
t139 = -mrSges(7,1) * t150 + mrSges(7,3) * t148 + Ifges(7,4) * t166 + Ifges(7,2) * t165 + Ifges(7,6) * t188 - t190 * t160 + t200 * t162;
t161 = Ifges(7,4) * t190 + Ifges(7,2) * t189 + Ifges(7,6) * t200;
t140 = mrSges(7,2) * t150 - mrSges(7,3) * t147 + Ifges(7,1) * t166 + Ifges(7,4) * t165 + Ifges(7,5) * t188 + t189 * t160 - t200 * t161;
t146 = -m(7) * t150 + t165 * mrSges(7,1) - t166 * mrSges(7,2) + t189 * t178 - t190 * t179;
t152 = -t222 * t156 + t259;
t184 = (Ifges(6,6) * qJD(5)) + (Ifges(6,4) * t222 + Ifges(6,2) * t225) * qJD(1);
t185 = (Ifges(6,5) * qJD(5)) + (Ifges(6,1) * t222 + Ifges(6,4) * t225) * qJD(1);
t265 = mrSges(6,1) * t152 - mrSges(6,2) * t153 + Ifges(6,5) * t193 + Ifges(6,6) * t194 + Ifges(6,3) * qJDD(5) + pkin(5) * t146 + pkin(8) * t138 + t224 * t139 + t221 * t140 + (t222 * t184 - t225 * t185) * qJD(1);
t264 = mrSges(3,2) - mrSges(4,3);
t263 = Ifges(3,4) + Ifges(4,6);
t262 = (Ifges(4,4) - Ifges(3,5));
t260 = t229 * mrSges(4,3);
t258 = t229 * qJ(2);
t191 = (-mrSges(6,1) * t225 + mrSges(6,2) * t222) * qJD(1);
t196 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t256;
t135 = m(6) * t153 - qJDD(5) * mrSges(6,2) + t194 * mrSges(6,3) - qJD(5) * t196 + t191 * t255 + t138;
t197 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t255;
t142 = m(6) * t152 + qJDD(5) * mrSges(6,1) - t193 * mrSges(6,3) + qJD(5) * t197 - t191 * t256 + t146;
t248 = t225 * t135 - t222 * t142;
t121 = m(5) * t159 - (t229 * mrSges(5,1)) - qJDD(1) * mrSges(5,2) + t248;
t137 = t224 * t144 + t221 * t145;
t234 = -m(6) * t155 + t194 * mrSges(6,1) - t193 * mrSges(6,2) - t196 * t256 + t197 * t255 - t137;
t131 = m(5) * t158 + qJDD(1) * mrSges(5,1) - t229 * mrSges(5,2) + t234;
t113 = t219 * t121 + t220 * t131;
t257 = t220 * t121 - t219 * t131;
t127 = t222 * t135 + t225 * t142;
t252 = (Ifges(2,6) + t262);
t251 = -m(5) * t214 - t127;
t247 = m(4) * t177 + qJDD(1) * mrSges(4,1) + t113;
t176 = t238 - t258;
t109 = m(4) * t176 - t229 * mrSges(4,1) - qJDD(1) * mrSges(4,3) + t257;
t183 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t222 + Ifges(6,6) * t225) * qJD(1);
t116 = mrSges(6,2) * t155 - mrSges(6,3) * t152 + Ifges(6,1) * t193 + Ifges(6,4) * t194 + Ifges(6,5) * qJDD(5) - pkin(8) * t137 - qJD(5) * t184 - t221 * t139 + t224 * t140 + t183 * t255;
t233 = mrSges(7,1) * t147 - mrSges(7,2) * t148 + Ifges(7,5) * t166 + Ifges(7,6) * t165 + Ifges(7,3) * t188 + t190 * t161 - t189 * t162;
t118 = -mrSges(6,1) * t155 + mrSges(6,3) * t153 + Ifges(6,4) * t193 + Ifges(6,2) * t194 + Ifges(6,6) * qJDD(5) - pkin(5) * t137 + qJD(5) * t185 - t183 * t256 - t233;
t102 = mrSges(5,2) * t214 - mrSges(5,3) * t158 + Ifges(5,5) * qJDD(1) - (t229 * Ifges(5,6)) - pkin(7) * t127 + t225 * t116 - t222 * t118;
t104 = -mrSges(5,1) * t214 + mrSges(5,3) * t159 + t229 * Ifges(5,5) + Ifges(5,6) * qJDD(1) - pkin(4) * t127 - t265;
t243 = -mrSges(4,1) * g(3) + mrSges(4,2) * t176 - pkin(3) * t251 - qJ(4) * t257 - t219 * t102 - t220 * t104;
t242 = -mrSges(4,2) * t177 - (t229 * Ifges(4,6)) + qJ(4) * t113 - t220 * t102 + t219 * t104;
t180 = (t229 * pkin(1)) + t266;
t241 = -m(3) * t180 + (t229 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t247;
t240 = mrSges(5,1) * t158 - mrSges(5,2) * t159 + Ifges(5,3) * qJDD(1) + pkin(4) * t234 + pkin(7) * t248 + t222 * t116 + t225 * t118;
t182 = -qJDD(1) * pkin(1) + t246 - t258;
t239 = -m(3) * t182 + t229 * mrSges(3,3) - t109;
t237 = -mrSges(3,1) * t182 - pkin(2) * t109 + t243;
t236 = -mrSges(4,1) * t177 + mrSges(4,3) * t176 - Ifges(4,2) * qJDD(1) - pkin(3) * t113 - t240;
t235 = mrSges(3,1) * t180 + pkin(2) * (-t247 + t260) + qJ(3) * (-m(4) * g(3) - t251) - t242;
t232 = mrSges(3,2) * t182 - mrSges(3,3) * t180 + Ifges(3,1) * qJDD(1) - qJ(3) * t109 - t236;
t230 = -mrSges(2,2) * t199 + qJ(2) * (t241 - t260) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t239) + mrSges(2,1) * t198 + Ifges(2,3) * qJDD(1) + t232;
t122 = (-m(3) - m(4)) * g(3) - t251;
t106 = m(2) * t198 - t229 * mrSges(2,2) + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) + t239;
t105 = m(2) * t199 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t229 + t241;
t99 = -mrSges(2,3) * t198 - qJ(2) * t122 + (-mrSges(2,2) + mrSges(3,3)) * g(3) - (t252 * t229) + (Ifges(2,5) - t263) * qJDD(1) - t237;
t98 = mrSges(2,3) * t199 - pkin(1) * t122 + (-Ifges(3,4) + Ifges(2,5)) * t229 + t252 * qJDD(1) + (mrSges(2,1) - t264) * g(3) - t235;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t226 * t99 - t223 * t98 - pkin(6) * (t223 * t105 + t226 * t106), t99, t232, -mrSges(4,3) * g(3) - Ifges(4,4) * qJDD(1) - t242, t102, t116, t140; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t223 * t99 + t226 * t98 + pkin(6) * (t226 * t105 - t223 * t106), t98, -mrSges(3,3) * g(3) + t263 * qJDD(1) + (t262 * t229) + t237, t236, t104, t118, t139; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t230, t230, t229 * Ifges(3,4) + t264 * g(3) - t262 * qJDD(1) + t235, -t229 * Ifges(4,4) - Ifges(4,6) * qJDD(1) - t243, t240, t265, t233;];
m_new  = t1;
