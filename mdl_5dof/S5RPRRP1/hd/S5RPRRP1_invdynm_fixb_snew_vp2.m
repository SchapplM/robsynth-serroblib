% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:19
% EndTime: 2019-12-05 17:59:26
% DurationCPUTime: 2.66s
% Computational Cost: add. (25008->264), mult. (49902->312), div. (0->0), fcn. (28771->6), ass. (0->98)
t234 = sin(qJ(4));
t235 = sin(qJ(3));
t237 = cos(qJ(4));
t238 = cos(qJ(3));
t199 = (-t234 * t238 - t235 * t237) * qJD(1);
t264 = qJD(1) * qJD(3);
t208 = -qJDD(1) * t235 - t238 * t264;
t259 = t235 * t264;
t209 = qJDD(1) * t238 - t259;
t164 = qJD(4) * t199 + t208 * t234 + t209 * t237;
t200 = (-t234 * t235 + t237 * t238) * qJD(1);
t174 = -mrSges(6,1) * t199 + mrSges(6,2) * t200;
t236 = sin(qJ(1));
t239 = cos(qJ(1));
t214 = t236 * g(1) - g(2) * t239;
t240 = qJD(1) ^ 2;
t253 = -t240 * qJ(2) + qJDD(2) - t214;
t274 = -pkin(1) - pkin(6);
t189 = qJDD(1) * t274 + t253;
t177 = g(3) * t235 + t189 * t238;
t144 = (-t209 - t259) * pkin(7) + (-t235 * t238 * t240 + qJDD(3)) * pkin(3) + t177;
t178 = -g(3) * t238 + t189 * t235;
t265 = qJD(1) * t238;
t213 = qJD(3) * pkin(3) - pkin(7) * t265;
t231 = t235 ^ 2;
t145 = -pkin(3) * t231 * t240 + pkin(7) * t208 - qJD(3) * t213 + t178;
t139 = t144 * t237 - t234 * t145;
t224 = qJDD(3) + qJDD(4);
t225 = qJD(3) + qJD(4);
t131 = -0.2e1 * qJD(5) * t200 + (t199 * t225 - t164) * qJ(5) + (t199 * t200 + t224) * pkin(4) + t139;
t183 = -mrSges(6,2) * t225 + mrSges(6,3) * t199;
t262 = m(6) * t131 + mrSges(6,1) * t224 + t183 * t225;
t128 = -t164 * mrSges(6,3) - t200 * t174 + t262;
t140 = t144 * t234 + t145 * t237;
t163 = -qJD(4) * t200 + t208 * t237 - t209 * t234;
t169 = Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * t225;
t170 = Ifges(6,1) * t200 + Ifges(6,4) * t199 + Ifges(6,5) * t225;
t171 = Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * t225;
t185 = pkin(4) * t225 - qJ(5) * t200;
t195 = t199 ^ 2;
t134 = -pkin(4) * t195 + t163 * qJ(5) + 0.2e1 * qJD(5) * t199 - t185 * t225 + t140;
t168 = Ifges(6,4) * t200 + Ifges(6,2) * t199 + Ifges(6,6) * t225;
t251 = -mrSges(6,1) * t131 + mrSges(6,2) * t134 - Ifges(6,5) * t164 - Ifges(6,6) * t163 - Ifges(6,3) * t224 - t168 * t200;
t279 = mrSges(5,1) * t139 - mrSges(5,2) * t140 + Ifges(5,5) * t164 + Ifges(5,6) * t163 + Ifges(5,3) * t224 + pkin(4) * t128 + t200 * t169 - t251 - (t170 + t171) * t199;
t175 = -mrSges(5,1) * t199 + mrSges(5,2) * t200;
t184 = -mrSges(5,2) * t225 + mrSges(5,3) * t199;
t123 = m(5) * t139 + t224 * mrSges(5,1) + t225 * t184 + (-t174 - t175) * t200 + (-mrSges(5,3) - mrSges(6,3)) * t164 + t262;
t186 = mrSges(6,1) * t225 - mrSges(6,3) * t200;
t187 = mrSges(5,1) * t225 - mrSges(5,3) * t200;
t261 = m(6) * t134 + t163 * mrSges(6,3) + t174 * t199;
t126 = m(5) * t140 + t163 * mrSges(5,3) + t199 * t175 + (-t186 - t187) * t225 + (-mrSges(5,2) - mrSges(6,2)) * t224 + t261;
t118 = t123 * t237 + t126 * t234;
t197 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t238 - Ifges(4,2) * t235) * qJD(1);
t198 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t238 - Ifges(4,4) * t235) * qJD(1);
t266 = qJD(1) * t235;
t278 = mrSges(4,1) * t177 - mrSges(4,2) * t178 + Ifges(4,5) * t209 + Ifges(4,6) * t208 + Ifges(4,3) * qJDD(3) + pkin(3) * t118 + t197 * t265 + t198 * t266 + t279;
t207 = (mrSges(4,1) * t235 + mrSges(4,2) * t238) * qJD(1);
t211 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t266;
t115 = m(4) * t177 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t209 + qJD(3) * t211 - t207 * t265 + t118;
t212 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t265;
t257 = -t123 * t234 + t126 * t237;
t116 = m(4) * t178 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t208 - qJD(3) * t212 - t207 * t266 + t257;
t110 = t238 * t115 + t235 * t116;
t194 = -qJDD(1) * pkin(1) + t253;
t277 = mrSges(3,1) * t194 + pkin(2) * t110 + t278;
t215 = -t239 * g(1) - t236 * g(2);
t254 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t215;
t188 = t240 * t274 + t254;
t150 = -t208 * pkin(3) + t213 * t265 + (-pkin(7) * t231 + t274) * t240 + t254;
t137 = -t163 * pkin(4) - t195 * qJ(5) + t200 * t185 + qJDD(5) + t150;
t260 = m(6) * t137 + t164 * mrSges(6,2) + t186 * t200;
t246 = m(5) * t150 + t164 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t163 + t200 * t187 - (t183 + t184) * t199 + t260;
t276 = -m(4) * t188 + t208 * mrSges(4,1) - t209 * mrSges(4,2) - t211 * t266 - t212 * t265 - t246;
t273 = mrSges(2,1) - mrSges(3,2);
t271 = Ifges(2,5) - Ifges(3,4);
t270 = -Ifges(2,6) + Ifges(3,5);
t111 = -t115 * t235 + t116 * t238;
t252 = -mrSges(6,1) * t137 + mrSges(6,3) * t134 + Ifges(6,4) * t164 + Ifges(6,2) * t163 + Ifges(6,6) * t224 + t170 * t225;
t166 = Ifges(6,5) * t200 + Ifges(6,6) * t199 + Ifges(6,3) * t225;
t250 = mrSges(6,2) * t137 - mrSges(6,3) * t131 + Ifges(6,1) * t164 + Ifges(6,4) * t163 + Ifges(6,5) * t224 + t166 * t199;
t249 = -m(3) * t194 + mrSges(3,3) * t240 - t110;
t167 = Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * t225;
t112 = Ifges(5,4) * t164 + Ifges(5,2) * t163 + Ifges(5,6) * t224 + t225 * t171 - mrSges(5,1) * t150 + mrSges(5,3) * t140 - pkin(4) * (-t163 * mrSges(6,1) - t199 * t183 + t260) + qJ(5) * (-t224 * mrSges(6,2) - t225 * t186 + t261) + (-t167 - t166) * t200 + t252;
t113 = mrSges(5,2) * t150 - mrSges(5,3) * t139 + Ifges(5,1) * t164 + Ifges(5,4) * t163 + Ifges(5,5) * t224 - qJ(5) * t128 + t199 * t167 + (-t168 - t169) * t225 + t250;
t196 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t238 - Ifges(4,6) * t235) * qJD(1);
t104 = -mrSges(4,1) * t188 + mrSges(4,3) * t178 + Ifges(4,4) * t209 + Ifges(4,2) * t208 + Ifges(4,6) * qJDD(3) - pkin(3) * t246 + pkin(7) * t257 + qJD(3) * t198 + t237 * t112 + t234 * t113 - t196 * t265;
t106 = mrSges(4,2) * t188 - mrSges(4,3) * t177 + Ifges(4,1) * t209 + Ifges(4,4) * t208 + Ifges(4,5) * qJDD(3) - pkin(7) * t118 - qJD(3) * t197 - t112 * t234 + t113 * t237 - t196 * t266;
t192 = pkin(1) * t240 - t254;
t248 = mrSges(3,2) * t194 - mrSges(3,3) * t192 + Ifges(3,1) * qJDD(1) - pkin(6) * t110 - t104 * t235 + t106 * t238;
t247 = -mrSges(3,1) * t192 - pkin(2) * t276 - pkin(6) * t111 - t238 * t104 - t235 * t106;
t243 = -m(3) * t192 + mrSges(3,2) * t240 + qJDD(1) * mrSges(3,3) - t276;
t244 = -mrSges(2,2) * t215 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t249) + qJ(2) * t243 + mrSges(2,1) * t214 + Ifges(2,3) * qJDD(1) + t248;
t119 = m(2) * t215 - t240 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t243;
t109 = -m(3) * g(3) + t111;
t107 = m(2) * t214 - t240 * mrSges(2,2) + qJDD(1) * t273 + t249;
t103 = -qJ(2) * t109 + t270 * t240 + t271 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t214 + t277;
t102 = mrSges(2,3) * t215 - pkin(1) * t109 + g(3) * t273 - qJDD(1) * t270 + t240 * t271 + t247;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t239 * t103 - t236 * t102 - pkin(5) * (t107 * t239 + t119 * t236), t103, t248, t106, t113, -t168 * t225 + t250; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t236 * t103 + t239 * t102 + pkin(5) * (-t107 * t236 + t119 * t239), t102, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t240 * Ifges(3,5) - t277, t104, t112, -t200 * t166 + t252; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t244, t244, mrSges(3,2) * g(3) + t240 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t247, t278, t279, -t199 * t170 - t251;];
m_new = t1;
