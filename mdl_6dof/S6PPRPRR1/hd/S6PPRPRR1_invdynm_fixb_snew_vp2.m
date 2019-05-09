% Calculate vector of cutting torques with Newton-Euler for
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-05-04 19:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:52:24
% EndTime: 2019-05-04 19:52:54
% DurationCPUTime: 21.63s
% Computational Cost: add. (441663->248), mult. (761245->329), div. (0->0), fcn. (625086->16), ass. (0->124)
t235 = sin(pkin(11));
t240 = cos(pkin(11));
t227 = -t240 * g(1) - t235 * g(2);
t234 = sin(pkin(12));
t239 = cos(pkin(12));
t226 = t235 * g(1) - t240 * g(2);
t232 = -g(3) + qJDD(1);
t237 = sin(pkin(6));
t242 = cos(pkin(6));
t257 = t226 * t242 + t232 * t237;
t197 = -t234 * t227 + t257 * t239;
t198 = t239 * t227 + t257 * t234;
t212 = -t237 * t226 + t242 * t232 + qJDD(2);
t245 = sin(qJ(3));
t241 = cos(pkin(7));
t248 = cos(qJ(3));
t268 = t241 * t248;
t236 = sin(pkin(7));
t270 = t236 * t248;
t188 = t197 * t268 - t245 * t198 + t212 * t270;
t186 = qJDD(3) * pkin(3) + t188;
t269 = t241 * t245;
t271 = t236 * t245;
t189 = t197 * t269 + t248 * t198 + t212 * t271;
t250 = qJD(3) ^ 2;
t187 = -t250 * pkin(3) + t189;
t233 = sin(pkin(13));
t238 = cos(pkin(13));
t182 = t233 * t186 + t238 * t187;
t179 = -t250 * pkin(4) + qJDD(3) * pkin(9) + t182;
t193 = -t236 * t197 + t241 * t212;
t192 = qJDD(4) + t193;
t244 = sin(qJ(5));
t247 = cos(qJ(5));
t175 = t247 * t179 + t244 * t192;
t222 = (-pkin(5) * t247 - pkin(10) * t244) * qJD(3);
t249 = qJD(5) ^ 2;
t265 = t247 * qJD(3);
t173 = -t249 * pkin(5) + qJDD(5) * pkin(10) + t222 * t265 + t175;
t181 = t238 * t186 - t233 * t187;
t178 = -qJDD(3) * pkin(4) - t250 * pkin(9) - t181;
t264 = qJD(3) * qJD(5);
t261 = t247 * t264;
t223 = t244 * qJDD(3) + t261;
t262 = t244 * t264;
t224 = t247 * qJDD(3) - t262;
t176 = (-t223 - t261) * pkin(10) + (-t224 + t262) * pkin(5) + t178;
t243 = sin(qJ(6));
t246 = cos(qJ(6));
t169 = -t243 * t173 + t246 * t176;
t266 = qJD(3) * t244;
t219 = t246 * qJD(5) - t243 * t266;
t205 = t219 * qJD(6) + t243 * qJDD(5) + t246 * t223;
t220 = t243 * qJD(5) + t246 * t266;
t206 = -t219 * mrSges(7,1) + t220 * mrSges(7,2);
t230 = qJD(6) - t265;
t210 = -t230 * mrSges(7,2) + t219 * mrSges(7,3);
t218 = qJDD(6) - t224;
t167 = m(7) * t169 + t218 * mrSges(7,1) - t205 * mrSges(7,3) - t220 * t206 + t230 * t210;
t170 = t246 * t173 + t243 * t176;
t204 = -t220 * qJD(6) + t246 * qJDD(5) - t243 * t223;
t211 = t230 * mrSges(7,1) - t220 * mrSges(7,3);
t168 = m(7) * t170 - t218 * mrSges(7,2) + t204 * mrSges(7,3) + t219 * t206 - t230 * t211;
t161 = -t243 * t167 + t246 * t168;
t221 = (-mrSges(6,1) * t247 + mrSges(6,2) * t244) * qJD(3);
t228 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t266;
t158 = m(6) * t175 - qJDD(5) * mrSges(6,2) + t224 * mrSges(6,3) - qJD(5) * t228 + t221 * t265 + t161;
t267 = t247 * t192;
t172 = -qJDD(5) * pkin(5) - t249 * pkin(10) - t267 + (qJD(3) * t222 + t179) * t244;
t171 = -m(7) * t172 + t204 * mrSges(7,1) - t205 * mrSges(7,2) + t219 * t210 - t220 * t211;
t174 = -t244 * t179 + t267;
t229 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t265;
t165 = m(6) * t174 + qJDD(5) * mrSges(6,1) - t223 * mrSges(6,3) + qJD(5) * t229 - t221 * t266 + t171;
t259 = t247 * t158 - t244 * t165;
t149 = m(5) * t182 - t250 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t259;
t160 = t246 * t167 + t243 * t168;
t253 = -m(6) * t178 + t224 * mrSges(6,1) - t223 * mrSges(6,2) - t228 * t266 + t229 * t265 - t160;
t155 = m(5) * t181 + qJDD(3) * mrSges(5,1) - t250 * mrSges(5,2) + t253;
t142 = t233 * t149 + t238 * t155;
t140 = m(4) * t188 + qJDD(3) * mrSges(4,1) - t250 * mrSges(4,2) + t142;
t260 = t238 * t149 - t233 * t155;
t141 = m(4) * t189 - t250 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t260;
t153 = t244 * t158 + t247 * t165;
t263 = m(5) * t192 + t153;
t151 = m(4) * t193 + t263;
t128 = t140 * t268 + t141 * t269 - t236 * t151;
t125 = m(3) * t197 + t128;
t133 = -t245 * t140 + t248 * t141;
t132 = m(3) * t198 + t133;
t279 = t125 * t239 + t132 * t234;
t199 = Ifges(7,5) * t220 + Ifges(7,6) * t219 + Ifges(7,3) * t230;
t201 = Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t230;
t162 = -mrSges(7,1) * t172 + mrSges(7,3) * t170 + Ifges(7,4) * t205 + Ifges(7,2) * t204 + Ifges(7,6) * t218 - t220 * t199 + t230 * t201;
t200 = Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t230;
t163 = mrSges(7,2) * t172 - mrSges(7,3) * t169 + Ifges(7,1) * t205 + Ifges(7,4) * t204 + Ifges(7,5) * t218 + t219 * t199 - t230 * t200;
t213 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t244 + Ifges(6,6) * t247) * qJD(3);
t214 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t244 + Ifges(6,2) * t247) * qJD(3);
t144 = mrSges(6,2) * t178 - mrSges(6,3) * t174 + Ifges(6,1) * t223 + Ifges(6,4) * t224 + Ifges(6,5) * qJDD(5) - pkin(10) * t160 - qJD(5) * t214 - t243 * t162 + t246 * t163 + t213 * t265;
t215 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t244 + Ifges(6,4) * t247) * qJD(3);
t252 = mrSges(7,1) * t169 - mrSges(7,2) * t170 + Ifges(7,5) * t205 + Ifges(7,6) * t204 + Ifges(7,3) * t218 + t220 * t200 - t219 * t201;
t146 = -mrSges(6,1) * t178 + mrSges(6,3) * t175 + Ifges(6,4) * t223 + Ifges(6,2) * t224 + Ifges(6,6) * qJDD(5) - pkin(5) * t160 + qJD(5) * t215 - t213 * t266 - t252;
t254 = mrSges(5,1) * t181 - mrSges(5,2) * t182 + Ifges(5,3) * qJDD(3) + pkin(4) * t253 + pkin(9) * t259 + t244 * t144 + t247 * t146;
t123 = mrSges(4,1) * t188 - mrSges(4,2) * t189 + Ifges(4,3) * qJDD(3) + pkin(3) * t142 + t254;
t127 = t140 * t270 + t141 * t271 + t241 * t151;
t129 = mrSges(5,2) * t192 - mrSges(5,3) * t181 + Ifges(5,5) * qJDD(3) - t250 * Ifges(5,6) - pkin(9) * t153 + t247 * t144 - t244 * t146;
t277 = mrSges(6,1) * t174 - mrSges(6,2) * t175 + Ifges(6,5) * t223 + Ifges(6,6) * t224 + Ifges(6,3) * qJDD(5) + pkin(5) * t171 + pkin(10) * t161 + t246 * t162 + t243 * t163 + (t244 * t214 - t247 * t215) * qJD(3);
t134 = -mrSges(5,1) * t192 + mrSges(5,3) * t182 + t250 * Ifges(5,5) + Ifges(5,6) * qJDD(3) - pkin(4) * t153 - t277;
t118 = -mrSges(4,1) * t193 + mrSges(4,3) * t189 + t250 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t263 + qJ(4) * t260 + t233 * t129 + t238 * t134;
t119 = mrSges(4,2) * t193 - mrSges(4,3) * t188 + Ifges(4,5) * qJDD(3) - t250 * Ifges(4,6) - qJ(4) * t142 + t238 * t129 - t233 * t134;
t256 = pkin(8) * t133 + t118 * t248 + t119 * t245;
t111 = -mrSges(3,1) * t212 + mrSges(3,3) * t198 - pkin(2) * t127 - t236 * t123 + t256 * t241;
t113 = mrSges(3,2) * t212 - mrSges(3,3) * t197 - t245 * t118 + t248 * t119 + (-t127 * t236 - t128 * t241) * pkin(8);
t122 = -t234 * t125 + t239 * t132;
t278 = qJ(2) * t122 + t111 * t239 + t113 * t234;
t126 = m(3) * t212 + t127;
t117 = -t237 * t126 + t242 * t279;
t109 = mrSges(3,1) * t197 - mrSges(3,2) * t198 + pkin(2) * t128 + t241 * t123 + t256 * t236;
t255 = mrSges(2,1) * t226 - mrSges(2,2) * t227 + pkin(1) * t117 + t242 * t109 + t237 * t278;
t120 = m(2) * t227 + t122;
t116 = t242 * t126 + t237 * t279;
t114 = m(2) * t226 + t117;
t107 = mrSges(2,2) * t232 - mrSges(2,3) * t226 - t234 * t111 + t239 * t113 + (-t116 * t237 - t117 * t242) * qJ(2);
t106 = -mrSges(2,1) * t232 + mrSges(2,3) * t227 - pkin(1) * t116 - t237 * t109 + t242 * t278;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t240 * t107 - t235 * t106 - qJ(1) * (t240 * t114 + t235 * t120), t107, t113, t119, t129, t144, t163; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t235 * t107 + t240 * t106 + qJ(1) * (-t235 * t114 + t240 * t120), t106, t111, t118, t134, t146, t162; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t255, t255, t109, t123, t254, t277, t252;];
m_new  = t1;
