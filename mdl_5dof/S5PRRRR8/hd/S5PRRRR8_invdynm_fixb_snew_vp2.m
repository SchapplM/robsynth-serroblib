% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:15:15
% EndTime: 2019-12-05 17:15:31
% DurationCPUTime: 7.55s
% Computational Cost: add. (127785->266), mult. (250262->349), div. (0->0), fcn. (174654->12), ass. (0->117)
t235 = sin(pkin(10));
t237 = cos(pkin(10));
t223 = t235 * g(1) - t237 * g(2);
t224 = -t237 * g(1) - t235 * g(2);
t234 = -g(3) + qJDD(1);
t246 = cos(qJ(2));
t238 = cos(pkin(5));
t242 = sin(qJ(2));
t266 = t238 * t242;
t236 = sin(pkin(5));
t267 = t236 * t242;
t196 = t223 * t266 + t246 * t224 + t234 * t267;
t247 = qJD(2) ^ 2;
t191 = -t247 * pkin(2) + qJDD(2) * pkin(7) + t196;
t206 = -t236 * t223 + t238 * t234;
t241 = sin(qJ(3));
t245 = cos(qJ(3));
t177 = -t241 * t191 + t245 * t206;
t263 = qJD(2) * qJD(3);
t262 = t245 * t263;
t220 = t241 * qJDD(2) + t262;
t167 = (-t220 + t262) * pkin(8) + (t241 * t245 * t247 + qJDD(3)) * pkin(3) + t177;
t178 = t245 * t191 + t241 * t206;
t221 = t245 * qJDD(2) - t241 * t263;
t265 = qJD(2) * t241;
t228 = qJD(3) * pkin(3) - pkin(8) * t265;
t233 = t245 ^ 2;
t168 = -t233 * t247 * pkin(3) + t221 * pkin(8) - qJD(3) * t228 + t178;
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t164 = t240 * t167 + t244 * t168;
t213 = (t240 * t245 + t241 * t244) * qJD(2);
t186 = -t213 * qJD(4) - t240 * t220 + t244 * t221;
t212 = (t240 * t241 - t244 * t245) * qJD(2);
t198 = t212 * mrSges(5,1) + t213 * mrSges(5,2);
t232 = qJD(3) + qJD(4);
t205 = t232 * mrSges(5,1) - t213 * mrSges(5,3);
t231 = qJDD(3) + qJDD(4);
t199 = t212 * pkin(4) - t213 * pkin(9);
t230 = t232 ^ 2;
t160 = -t230 * pkin(4) + t231 * pkin(9) - t212 * t199 + t164;
t195 = -t242 * t224 + (t223 * t238 + t234 * t236) * t246;
t253 = -qJDD(2) * pkin(2) - t195;
t172 = -t221 * pkin(3) + t228 * t265 + (-pkin(8) * t233 - pkin(7)) * t247 + t253;
t187 = -t212 * qJD(4) + t244 * t220 + t240 * t221;
t161 = (t212 * t232 - t187) * pkin(9) + (t213 * t232 - t186) * pkin(4) + t172;
t239 = sin(qJ(5));
t243 = cos(qJ(5));
t157 = -t239 * t160 + t243 * t161;
t200 = -t239 * t213 + t243 * t232;
t171 = t200 * qJD(5) + t243 * t187 + t239 * t231;
t201 = t243 * t213 + t239 * t232;
t179 = -t200 * mrSges(6,1) + t201 * mrSges(6,2);
t184 = qJDD(5) - t186;
t207 = qJD(5) + t212;
t188 = -t207 * mrSges(6,2) + t200 * mrSges(6,3);
t153 = m(6) * t157 + t184 * mrSges(6,1) - t171 * mrSges(6,3) - t201 * t179 + t207 * t188;
t158 = t243 * t160 + t239 * t161;
t170 = -t201 * qJD(5) - t239 * t187 + t243 * t231;
t189 = t207 * mrSges(6,1) - t201 * mrSges(6,3);
t154 = m(6) * t158 - t184 * mrSges(6,2) + t170 * mrSges(6,3) + t200 * t179 - t207 * t189;
t259 = -t239 * t153 + t243 * t154;
t140 = m(5) * t164 - t231 * mrSges(5,2) + t186 * mrSges(5,3) - t212 * t198 - t232 * t205 + t259;
t163 = t244 * t167 - t240 * t168;
t204 = -t232 * mrSges(5,2) - t212 * mrSges(5,3);
t159 = -t231 * pkin(4) - t230 * pkin(9) + t213 * t199 - t163;
t254 = -m(6) * t159 + t170 * mrSges(6,1) - t171 * mrSges(6,2) + t200 * t188 - t201 * t189;
t149 = m(5) * t163 + t231 * mrSges(5,1) - t187 * mrSges(5,3) - t213 * t198 + t232 * t204 + t254;
t135 = t240 * t140 + t244 * t149;
t219 = (-mrSges(4,1) * t245 + mrSges(4,2) * t241) * qJD(2);
t264 = qJD(2) * t245;
t226 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t264;
t133 = m(4) * t177 + qJDD(3) * mrSges(4,1) - t220 * mrSges(4,3) + qJD(3) * t226 - t219 * t265 + t135;
t225 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t265;
t260 = t244 * t140 - t240 * t149;
t134 = m(4) * t178 - qJDD(3) * mrSges(4,2) + t221 * mrSges(4,3) - qJD(3) * t225 + t219 * t264 + t260;
t127 = t245 * t133 + t241 * t134;
t210 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t241 + Ifges(4,2) * t245) * qJD(2);
t211 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t241 + Ifges(4,4) * t245) * qJD(2);
t173 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t207;
t175 = Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t207;
t146 = -mrSges(6,1) * t159 + mrSges(6,3) * t158 + Ifges(6,4) * t171 + Ifges(6,2) * t170 + Ifges(6,6) * t184 - t201 * t173 + t207 * t175;
t174 = Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t207;
t147 = mrSges(6,2) * t159 - mrSges(6,3) * t157 + Ifges(6,1) * t171 + Ifges(6,4) * t170 + Ifges(6,5) * t184 + t200 * t173 - t207 * t174;
t193 = Ifges(5,4) * t213 - Ifges(5,2) * t212 + Ifges(5,6) * t232;
t194 = Ifges(5,1) * t213 - Ifges(5,4) * t212 + Ifges(5,5) * t232;
t251 = -mrSges(5,1) * t163 + mrSges(5,2) * t164 - Ifges(5,5) * t187 - Ifges(5,6) * t186 - Ifges(5,3) * t231 - pkin(4) * t254 - pkin(9) * t259 - t243 * t146 - t239 * t147 - t213 * t193 - t212 * t194;
t271 = mrSges(4,1) * t177 - mrSges(4,2) * t178 + Ifges(4,5) * t220 + Ifges(4,6) * t221 + Ifges(4,3) * qJDD(3) + pkin(3) * t135 + (t241 * t210 - t245 * t211) * qJD(2) - t251;
t113 = -mrSges(3,1) * t206 + mrSges(3,3) * t196 + t247 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t127 - t271;
t261 = -t241 * t133 + t245 * t134;
t125 = m(3) * t196 - t247 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t261;
t190 = -t247 * pkin(7) + t253;
t142 = t243 * t153 + t239 * t154;
t252 = m(5) * t172 - t186 * mrSges(5,1) + t187 * mrSges(5,2) + t212 * t204 + t213 * t205 + t142;
t249 = -m(4) * t190 + t221 * mrSges(4,1) - t220 * mrSges(4,2) - t225 * t265 + t226 * t264 - t252;
t137 = m(3) * t195 + qJDD(2) * mrSges(3,1) - t247 * mrSges(3,2) + t249;
t122 = t246 * t125 - t242 * t137;
t272 = pkin(6) * t122 + t113 * t246;
t268 = t137 * t246;
t126 = m(3) * t206 + t127;
t117 = t125 * t266 - t236 * t126 + t238 * t268;
t192 = Ifges(5,5) * t213 - Ifges(5,6) * t212 + Ifges(5,3) * t232;
t128 = mrSges(5,2) * t172 - mrSges(5,3) * t163 + Ifges(5,1) * t187 + Ifges(5,4) * t186 + Ifges(5,5) * t231 - pkin(9) * t142 - t239 * t146 + t243 * t147 - t212 * t192 - t232 * t193;
t250 = mrSges(6,1) * t157 - mrSges(6,2) * t158 + Ifges(6,5) * t171 + Ifges(6,6) * t170 + Ifges(6,3) * t184 + t201 * t174 - t200 * t175;
t129 = -mrSges(5,1) * t172 + mrSges(5,3) * t164 + Ifges(5,4) * t187 + Ifges(5,2) * t186 + Ifges(5,6) * t231 - pkin(4) * t142 - t213 * t192 + t232 * t194 - t250;
t209 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t241 + Ifges(4,6) * t245) * qJD(2);
t118 = -mrSges(4,1) * t190 + mrSges(4,3) * t178 + Ifges(4,4) * t220 + Ifges(4,2) * t221 + Ifges(4,6) * qJDD(3) - pkin(3) * t252 + pkin(8) * t260 + qJD(3) * t211 + t240 * t128 + t244 * t129 - t209 * t265;
t119 = mrSges(4,2) * t190 - mrSges(4,3) * t177 + Ifges(4,1) * t220 + Ifges(4,4) * t221 + Ifges(4,5) * qJDD(3) - pkin(8) * t135 - qJD(3) * t210 + t244 * t128 - t240 * t129 + t209 * t264;
t109 = mrSges(3,1) * t195 - mrSges(3,2) * t196 + Ifges(3,3) * qJDD(2) + pkin(2) * t249 + pkin(7) * t261 + t245 * t118 + t241 * t119;
t111 = mrSges(3,2) * t206 - mrSges(3,3) * t195 + Ifges(3,5) * qJDD(2) - t247 * Ifges(3,6) - pkin(7) * t127 - t241 * t118 + t245 * t119;
t255 = mrSges(2,1) * t223 - mrSges(2,2) * t224 + pkin(1) * t117 + t238 * t109 + t111 * t267 + t272 * t236;
t120 = m(2) * t224 + t122;
t116 = t238 * t126 + (t125 * t242 + t268) * t236;
t114 = m(2) * t223 + t117;
t107 = mrSges(2,2) * t234 - mrSges(2,3) * t223 + t246 * t111 - t242 * t113 + (-t116 * t236 - t117 * t238) * pkin(6);
t106 = -mrSges(2,1) * t234 + mrSges(2,3) * t224 - pkin(1) * t116 - t236 * t109 + (t111 * t242 + t272) * t238;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t237 * t107 - t235 * t106 - qJ(1) * (t237 * t114 + t235 * t120), t107, t111, t119, t128, t147; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t235 * t107 + t237 * t106 + qJ(1) * (-t235 * t114 + t237 * t120), t106, t113, t118, t129, t146; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t255, t255, t109, t271, -t251, t250;];
m_new = t1;
