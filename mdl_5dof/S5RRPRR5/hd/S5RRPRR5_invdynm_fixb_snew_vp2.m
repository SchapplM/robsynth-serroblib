% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:33:45
% EndTime: 2019-12-05 18:33:51
% DurationCPUTime: 5.72s
% Computational Cost: add. (131840->248), mult. (182364->312), div. (0->0), fcn. (122991->10), ass. (0->112)
t229 = qJD(1) + qJD(2);
t223 = t229 ^ 2;
t237 = sin(qJ(1));
t241 = cos(qJ(1));
t214 = t241 * g(2) + t237 * g(3);
t208 = qJDD(1) * pkin(1) + t214;
t213 = t237 * g(2) - g(3) * t241;
t242 = qJD(1) ^ 2;
t209 = -pkin(1) * t242 + t213;
t236 = sin(qJ(2));
t240 = cos(qJ(2));
t194 = t236 * t208 + t240 * t209;
t225 = qJDD(1) + qJDD(2);
t191 = -pkin(2) * t223 + qJ(3) * t225 + t194;
t232 = sin(pkin(9));
t233 = cos(pkin(9));
t265 = qJD(3) * t229;
t263 = -t233 * g(1) - 0.2e1 * t232 * t265;
t270 = pkin(3) * t233;
t169 = (-pkin(7) * t225 + t223 * t270 - t191) * t232 + t263;
t174 = -t232 * g(1) + (t191 + 0.2e1 * t265) * t233;
t267 = t225 * t233;
t227 = t233 ^ 2;
t268 = t223 * t227;
t170 = -pkin(3) * t268 + pkin(7) * t267 + t174;
t235 = sin(qJ(4));
t239 = cos(qJ(4));
t151 = t239 * t169 - t235 * t170;
t251 = t232 * t239 + t233 * t235;
t250 = -t232 * t235 + t233 * t239;
t199 = t250 * t229;
t264 = t199 * qJD(4);
t190 = t225 * t251 + t264;
t200 = t251 * t229;
t146 = (-t190 + t264) * pkin(8) + (t199 * t200 + qJDD(4)) * pkin(4) + t151;
t152 = t235 * t169 + t239 * t170;
t189 = -t200 * qJD(4) + t225 * t250;
t197 = qJD(4) * pkin(4) - pkin(8) * t200;
t198 = t199 ^ 2;
t147 = -pkin(4) * t198 + pkin(8) * t189 - qJD(4) * t197 + t152;
t234 = sin(qJ(5));
t238 = cos(qJ(5));
t144 = t146 * t238 - t147 * t234;
t180 = t199 * t238 - t200 * t234;
t159 = qJD(5) * t180 + t189 * t234 + t190 * t238;
t181 = t199 * t234 + t200 * t238;
t165 = -mrSges(6,1) * t180 + mrSges(6,2) * t181;
t228 = qJD(4) + qJD(5);
t175 = -mrSges(6,2) * t228 + mrSges(6,3) * t180;
t224 = qJDD(4) + qJDD(5);
t141 = m(6) * t144 + mrSges(6,1) * t224 - mrSges(6,3) * t159 - t165 * t181 + t175 * t228;
t145 = t146 * t234 + t147 * t238;
t158 = -qJD(5) * t181 + t189 * t238 - t190 * t234;
t176 = mrSges(6,1) * t228 - mrSges(6,3) * t181;
t142 = m(6) * t145 - mrSges(6,2) * t224 + mrSges(6,3) * t158 + t165 * t180 - t176 * t228;
t132 = t141 * t238 + t142 * t234;
t187 = -mrSges(5,1) * t199 + mrSges(5,2) * t200;
t195 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t199;
t129 = m(5) * t151 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t190 + qJD(4) * t195 - t187 * t200 + t132;
t196 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t200;
t259 = -t141 * t234 + t142 * t238;
t130 = m(5) * t152 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t189 - qJD(4) * t196 + t187 * t199 + t259;
t125 = t129 * t239 + t130 * t235;
t173 = -t232 * t191 + t263;
t178 = Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * qJD(4);
t179 = Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * qJD(4);
t161 = Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t228;
t162 = Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t228;
t248 = -mrSges(6,1) * t144 + mrSges(6,2) * t145 - Ifges(6,5) * t159 - Ifges(6,6) * t158 - Ifges(6,3) * t224 - t181 * t161 + t180 * t162;
t244 = -mrSges(5,1) * t151 + mrSges(5,2) * t152 - Ifges(5,5) * t190 - Ifges(5,6) * t189 - Ifges(5,3) * qJDD(4) - pkin(4) * t132 - t200 * t178 + t199 * t179 + t248;
t257 = Ifges(4,4) * t232 + Ifges(4,2) * t233;
t258 = Ifges(4,1) * t232 + Ifges(4,4) * t233;
t271 = -mrSges(4,1) * t173 + mrSges(4,2) * t174 - pkin(3) * t125 - (t232 * t257 - t233 * t258) * t223 + t244;
t269 = mrSges(4,2) * t232;
t256 = Ifges(4,5) * t232 + Ifges(4,6) * t233;
t266 = t223 * t256;
t253 = mrSges(4,3) * t225 + (-mrSges(4,1) * t233 + t269) * t223;
t123 = m(4) * t173 - t232 * t253 + t125;
t260 = -t235 * t129 + t130 * t239;
t124 = m(4) * t174 + t233 * t253 + t260;
t261 = -t123 * t232 + t124 * t233;
t115 = m(3) * t194 - mrSges(3,1) * t223 - mrSges(3,2) * t225 + t261;
t193 = t240 * t208 - t236 * t209;
t255 = qJDD(3) - t193;
t188 = -t225 * pkin(2) - t223 * qJ(3) + t255;
t226 = t232 ^ 2;
t172 = (-pkin(2) - t270) * t225 + (-qJ(3) + (-t226 - t227) * pkin(7)) * t223 + t255;
t149 = -t189 * pkin(4) - t198 * pkin(8) + t200 * t197 + t172;
t254 = m(6) * t149 - t158 * mrSges(6,1) + t159 * mrSges(6,2) - t180 * t175 + t181 * t176;
t246 = m(5) * t172 - t189 * mrSges(5,1) + t190 * mrSges(5,2) - t199 * t195 + t200 * t196 + t254;
t245 = -m(4) * t188 + mrSges(4,1) * t267 - t246 + (t223 * t226 + t268) * mrSges(4,3);
t136 = t245 + (mrSges(3,1) - t269) * t225 - t223 * mrSges(3,2) + m(3) * t193;
t112 = t115 * t236 + t136 * t240;
t117 = t123 * t233 + t124 * t232;
t262 = t115 * t240 - t136 * t236;
t160 = Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t228;
t133 = -mrSges(6,1) * t149 + mrSges(6,3) * t145 + Ifges(6,4) * t159 + Ifges(6,2) * t158 + Ifges(6,6) * t224 - t160 * t181 + t162 * t228;
t134 = mrSges(6,2) * t149 - mrSges(6,3) * t144 + Ifges(6,1) * t159 + Ifges(6,4) * t158 + Ifges(6,5) * t224 + t160 * t180 - t161 * t228;
t177 = Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * qJD(4);
t118 = -mrSges(5,1) * t172 + mrSges(5,3) * t152 + Ifges(5,4) * t190 + Ifges(5,2) * t189 + Ifges(5,6) * qJDD(4) - pkin(4) * t254 + pkin(8) * t259 + qJD(4) * t179 + t238 * t133 + t234 * t134 - t200 * t177;
t119 = mrSges(5,2) * t172 - mrSges(5,3) * t151 + Ifges(5,1) * t190 + Ifges(5,4) * t189 + Ifges(5,5) * qJDD(4) - pkin(8) * t132 - qJD(4) * t178 - t133 * t234 + t134 * t238 + t177 * t199;
t106 = -mrSges(4,1) * t188 + mrSges(4,3) * t174 - pkin(3) * t246 + pkin(7) * t260 + t239 * t118 + t235 * t119 + t225 * t257 - t232 * t266;
t108 = mrSges(4,2) * t188 - mrSges(4,3) * t173 - pkin(7) * t125 - t235 * t118 + t239 * t119 + t225 * t258 + t233 * t266;
t249 = -mrSges(3,2) * t194 + qJ(3) * t261 + t233 * t106 + t232 * t108 + pkin(2) * (-t225 * t269 + t245) + mrSges(3,1) * t193 + Ifges(3,3) * t225;
t247 = mrSges(2,1) * t214 - mrSges(2,2) * t213 + Ifges(2,3) * qJDD(1) + pkin(1) * t112 + t249;
t110 = m(2) * t214 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t242 + t112;
t109 = m(2) * t213 - mrSges(2,1) * t242 - qJDD(1) * mrSges(2,2) + t262;
t104 = mrSges(3,1) * g(1) + (Ifges(3,6) - t256) * t225 + t223 * Ifges(3,5) + mrSges(3,3) * t194 - pkin(2) * t117 + t271;
t103 = -mrSges(3,2) * g(1) - mrSges(3,3) * t193 + Ifges(3,5) * t225 - Ifges(3,6) * t223 - qJ(3) * t117 - t106 * t232 + t108 * t233;
t102 = -mrSges(2,2) * g(1) - mrSges(2,3) * t214 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t242 - pkin(6) * t112 + t103 * t240 - t104 * t236;
t101 = Ifges(2,6) * qJDD(1) + t242 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t213 + t236 * t103 + t240 * t104 - pkin(1) * (-m(3) * g(1) + t117) + pkin(6) * t262;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t247, t102, t103, t108, t119, t134; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t237 * t102 - t241 * t101 - pkin(5) * (t109 * t241 - t110 * t237), t101, t104, t106, t118, t133; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t241 * t102 - t237 * t101 + pkin(5) * (-t109 * t237 - t110 * t241), t247, t249, t225 * t256 - t271, -t244, -t248;];
m_new = t1;
