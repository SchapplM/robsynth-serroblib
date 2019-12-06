% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:41:45
% EndTime: 2019-12-05 17:41:54
% DurationCPUTime: 6.70s
% Computational Cost: add. (82225->247), mult. (182364->311), div. (0->0), fcn. (122991->10), ass. (0->110)
t242 = qJD(1) ^ 2;
t238 = sin(qJ(1));
t241 = cos(qJ(1));
t212 = t241 * g(2) + t238 * g(3);
t208 = qJDD(1) * pkin(1) + t212;
t211 = t238 * g(2) - t241 * g(3);
t209 = -t242 * pkin(1) + t211;
t233 = sin(pkin(8));
t235 = cos(pkin(8));
t194 = t233 * t208 + t235 * t209;
t184 = -t242 * pkin(2) + qJDD(1) * qJ(3) + t194;
t232 = sin(pkin(9));
t231 = -g(1) + qJDD(2);
t234 = cos(pkin(9));
t264 = qJD(1) * qJD(3);
t267 = t234 * t231 - 0.2e1 * t232 * t264;
t270 = pkin(3) * t234;
t169 = (-pkin(6) * qJDD(1) + t242 * t270 - t184) * t232 + t267;
t173 = t232 * t231 + (t184 + 0.2e1 * t264) * t234;
t263 = qJDD(1) * t234;
t226 = t234 ^ 2;
t268 = t226 * t242;
t170 = -pkin(3) * t268 + pkin(6) * t263 + t173;
t237 = sin(qJ(4));
t240 = cos(qJ(4));
t151 = t240 * t169 - t237 * t170;
t252 = t232 * t240 + t234 * t237;
t251 = -t232 * t237 + t234 * t240;
t199 = t251 * qJD(1);
t265 = t199 * qJD(4);
t191 = t252 * qJDD(1) + t265;
t200 = t252 * qJD(1);
t146 = (-t191 + t265) * pkin(7) + (t199 * t200 + qJDD(4)) * pkin(4) + t151;
t152 = t237 * t169 + t240 * t170;
t190 = -t200 * qJD(4) + t251 * qJDD(1);
t197 = qJD(4) * pkin(4) - t200 * pkin(7);
t198 = t199 ^ 2;
t147 = -t198 * pkin(4) + t190 * pkin(7) - qJD(4) * t197 + t152;
t236 = sin(qJ(5));
t239 = cos(qJ(5));
t144 = t239 * t146 - t236 * t147;
t182 = t239 * t199 - t236 * t200;
t159 = t182 * qJD(5) + t236 * t190 + t239 * t191;
t183 = t236 * t199 + t239 * t200;
t165 = -t182 * mrSges(6,1) + t183 * mrSges(6,2);
t227 = qJD(4) + qJD(5);
t175 = -t227 * mrSges(6,2) + t182 * mrSges(6,3);
t224 = qJDD(4) + qJDD(5);
t141 = m(6) * t144 + t224 * mrSges(6,1) - t159 * mrSges(6,3) - t183 * t165 + t227 * t175;
t145 = t236 * t146 + t239 * t147;
t158 = -t183 * qJD(5) + t239 * t190 - t236 * t191;
t176 = t227 * mrSges(6,1) - t183 * mrSges(6,3);
t142 = m(6) * t145 - t224 * mrSges(6,2) + t158 * mrSges(6,3) + t182 * t165 - t227 * t176;
t132 = t239 * t141 + t236 * t142;
t186 = -t199 * mrSges(5,1) + t200 * mrSges(5,2);
t195 = -qJD(4) * mrSges(5,2) + t199 * mrSges(5,3);
t129 = m(5) * t151 + qJDD(4) * mrSges(5,1) - t191 * mrSges(5,3) + qJD(4) * t195 - t200 * t186 + t132;
t196 = qJD(4) * mrSges(5,1) - t200 * mrSges(5,3);
t259 = -t236 * t141 + t239 * t142;
t130 = m(5) * t152 - qJDD(4) * mrSges(5,2) + t190 * mrSges(5,3) - qJD(4) * t196 + t199 * t186 + t259;
t125 = t240 * t129 + t237 * t130;
t172 = -t232 * t184 + t267;
t180 = Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * qJD(4);
t181 = Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * qJD(4);
t161 = Ifges(6,4) * t183 + Ifges(6,2) * t182 + Ifges(6,6) * t227;
t162 = Ifges(6,1) * t183 + Ifges(6,4) * t182 + Ifges(6,5) * t227;
t248 = -mrSges(6,1) * t144 + mrSges(6,2) * t145 - Ifges(6,5) * t159 - Ifges(6,6) * t158 - Ifges(6,3) * t224 - t183 * t161 + t182 * t162;
t244 = -mrSges(5,1) * t151 + mrSges(5,2) * t152 - Ifges(5,5) * t191 - Ifges(5,6) * t190 - Ifges(5,3) * qJDD(4) - pkin(4) * t132 - t200 * t180 + t199 * t181 + t248;
t257 = Ifges(4,4) * t232 + Ifges(4,2) * t234;
t258 = Ifges(4,1) * t232 + Ifges(4,4) * t234;
t271 = -mrSges(4,1) * t172 + mrSges(4,2) * t173 - pkin(3) * t125 - (t232 * t257 - t234 * t258) * t242 + t244;
t269 = mrSges(4,2) * t232;
t250 = mrSges(4,3) * qJDD(1) + t242 * (-mrSges(4,1) * t234 + t269);
t123 = m(4) * t172 - t250 * t232 + t125;
t260 = -t237 * t129 + t240 * t130;
t124 = m(4) * t173 + t250 * t234 + t260;
t261 = -t232 * t123 + t234 * t124;
t115 = m(3) * t194 - t242 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t261;
t193 = t235 * t208 - t233 * t209;
t255 = qJDD(3) - t193;
t178 = -qJDD(1) * pkin(2) - t242 * qJ(3) + t255;
t225 = t232 ^ 2;
t171 = (-pkin(2) - t270) * qJDD(1) + (-qJ(3) + (-t225 - t226) * pkin(6)) * t242 + t255;
t149 = -t190 * pkin(4) - t198 * pkin(7) + t200 * t197 + t171;
t254 = m(6) * t149 - t158 * mrSges(6,1) + t159 * mrSges(6,2) - t182 * t175 + t183 * t176;
t246 = m(5) * t171 - t190 * mrSges(5,1) + t191 * mrSges(5,2) - t199 * t195 + t200 * t196 + t254;
t245 = -m(4) * t178 + mrSges(4,1) * t263 - t246 + (t225 * t242 + t268) * mrSges(4,3);
t136 = (mrSges(3,1) - t269) * qJDD(1) - t242 * mrSges(3,2) + m(3) * t193 + t245;
t112 = t233 * t115 + t235 * t136;
t117 = t234 * t123 + t232 * t124;
t256 = Ifges(4,5) * t232 + Ifges(4,6) * t234;
t266 = t242 * t256;
t262 = t235 * t115 - t233 * t136;
t160 = Ifges(6,5) * t183 + Ifges(6,6) * t182 + Ifges(6,3) * t227;
t133 = -mrSges(6,1) * t149 + mrSges(6,3) * t145 + Ifges(6,4) * t159 + Ifges(6,2) * t158 + Ifges(6,6) * t224 - t183 * t160 + t227 * t162;
t134 = mrSges(6,2) * t149 - mrSges(6,3) * t144 + Ifges(6,1) * t159 + Ifges(6,4) * t158 + Ifges(6,5) * t224 + t182 * t160 - t227 * t161;
t179 = Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * qJD(4);
t118 = -mrSges(5,1) * t171 + mrSges(5,3) * t152 + Ifges(5,4) * t191 + Ifges(5,2) * t190 + Ifges(5,6) * qJDD(4) - pkin(4) * t254 + pkin(7) * t259 + qJD(4) * t181 + t239 * t133 + t236 * t134 - t200 * t179;
t119 = mrSges(5,2) * t171 - mrSges(5,3) * t151 + Ifges(5,1) * t191 + Ifges(5,4) * t190 + Ifges(5,5) * qJDD(4) - pkin(7) * t132 - qJD(4) * t180 - t236 * t133 + t239 * t134 + t199 * t179;
t106 = -mrSges(4,1) * t178 + mrSges(4,3) * t173 - pkin(3) * t246 + pkin(6) * t260 + t257 * qJDD(1) + t240 * t118 + t237 * t119 - t232 * t266;
t108 = mrSges(4,2) * t178 - mrSges(4,3) * t172 - pkin(6) * t125 + t258 * qJDD(1) - t237 * t118 + t240 * t119 + t234 * t266;
t249 = -mrSges(3,2) * t194 + qJ(3) * t261 + t234 * t106 + t232 * t108 + pkin(2) * (-qJDD(1) * t269 + t245) + mrSges(3,1) * t193 + Ifges(3,3) * qJDD(1);
t247 = mrSges(2,1) * t212 - mrSges(2,2) * t211 + Ifges(2,3) * qJDD(1) + pkin(1) * t112 + t249;
t110 = m(2) * t212 + qJDD(1) * mrSges(2,1) - t242 * mrSges(2,2) + t112;
t109 = m(2) * t211 - t242 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t262;
t104 = (Ifges(3,6) - t256) * qJDD(1) + t242 * Ifges(3,5) - mrSges(3,1) * t231 + mrSges(3,3) * t194 - pkin(2) * t117 + t271;
t103 = mrSges(3,2) * t231 - mrSges(3,3) * t193 + Ifges(3,5) * qJDD(1) - t242 * Ifges(3,6) - qJ(3) * t117 - t232 * t106 + t234 * t108;
t102 = -mrSges(2,2) * g(1) - mrSges(2,3) * t212 + Ifges(2,5) * qJDD(1) - t242 * Ifges(2,6) - qJ(2) * t112 + t235 * t103 - t233 * t104;
t101 = Ifges(2,6) * qJDD(1) + t242 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t211 + t233 * t103 + t235 * t104 - pkin(1) * (m(3) * t231 + t117) + qJ(2) * t262;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t247, t102, t103, t108, t119, t134; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t238 * t102 - t241 * t101 - pkin(5) * (t241 * t109 - t238 * t110), t101, t104, t106, t118, t133; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t241 * t102 - t238 * t101 + pkin(5) * (-t238 * t109 - t241 * t110), t247, t249, t256 * qJDD(1) - t271, -t244, -t248;];
m_new = t1;
