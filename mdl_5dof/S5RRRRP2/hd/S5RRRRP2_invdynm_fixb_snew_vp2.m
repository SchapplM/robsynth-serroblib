% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:38
% EndTime: 2019-12-05 18:47:45
% DurationCPUTime: 3.88s
% Computational Cost: add. (63584->264), mult. (80564->325), div. (0->0), fcn. (47650->8), ass. (0->102)
t237 = qJD(1) + qJD(2);
t240 = sin(qJ(4));
t241 = sin(qJ(3));
t244 = cos(qJ(4));
t245 = cos(qJ(3));
t205 = (-t240 * t241 + t244 * t245) * t237;
t235 = qJDD(1) + qJDD(2);
t267 = qJD(3) * t237;
t211 = t241 * t235 + t245 * t267;
t212 = t245 * t235 - t241 * t267;
t174 = t205 * qJD(4) + t244 * t211 + t240 * t212;
t206 = (t240 * t245 + t241 * t244) * t237;
t188 = -t205 * mrSges(6,1) + t206 * mrSges(6,2);
t243 = sin(qJ(1));
t247 = cos(qJ(1));
t223 = t247 * g(2) + t243 * g(3);
t216 = qJDD(1) * pkin(1) + t223;
t222 = t243 * g(2) - t247 * g(3);
t248 = qJD(1) ^ 2;
t217 = -t248 * pkin(1) + t222;
t242 = sin(qJ(2));
t246 = cos(qJ(2));
t194 = t242 * t216 + t246 * t217;
t233 = t237 ^ 2;
t191 = -t233 * pkin(2) + t235 * pkin(7) + t194;
t269 = t241 * t191;
t272 = pkin(3) * t233;
t154 = qJDD(3) * pkin(3) - t211 * pkin(8) - t269 + (pkin(8) * t267 + t241 * t272 - g(1)) * t245;
t177 = -t241 * g(1) + t245 * t191;
t271 = t237 * t241;
t220 = qJD(3) * pkin(3) - pkin(8) * t271;
t239 = t245 ^ 2;
t155 = t212 * pkin(8) - qJD(3) * t220 - t239 * t272 + t177;
t149 = t244 * t154 - t240 * t155;
t234 = qJDD(3) + qJDD(4);
t236 = qJD(3) + qJD(4);
t141 = -0.2e1 * qJD(5) * t206 + (t205 * t236 - t174) * qJ(5) + (t205 * t206 + t234) * pkin(4) + t149;
t196 = -t236 * mrSges(6,2) + t205 * mrSges(6,3);
t266 = m(6) * t141 + t234 * mrSges(6,1) + t236 * t196;
t138 = -t174 * mrSges(6,3) - t206 * t188 + t266;
t150 = t240 * t154 + t244 * t155;
t173 = -t206 * qJD(4) - t240 * t211 + t244 * t212;
t181 = Ifges(5,4) * t206 + Ifges(5,2) * t205 + Ifges(5,6) * t236;
t182 = Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t236;
t183 = Ifges(5,1) * t206 + Ifges(5,4) * t205 + Ifges(5,5) * t236;
t198 = t236 * pkin(4) - t206 * qJ(5);
t201 = t205 ^ 2;
t144 = -t201 * pkin(4) + t173 * qJ(5) + 0.2e1 * qJD(5) * t205 - t236 * t198 + t150;
t180 = Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t236;
t256 = -mrSges(6,1) * t141 + mrSges(6,2) * t144 - Ifges(6,5) * t174 - Ifges(6,6) * t173 - Ifges(6,3) * t234 - t206 * t180;
t275 = mrSges(5,1) * t149 - mrSges(5,2) * t150 + Ifges(5,5) * t174 + Ifges(5,6) * t173 + Ifges(5,3) * t234 + pkin(4) * t138 + t206 * t181 - t256 + (-t183 - t182) * t205;
t189 = -t205 * mrSges(5,1) + t206 * mrSges(5,2);
t197 = -t236 * mrSges(5,2) + t205 * mrSges(5,3);
t133 = m(5) * t149 + t234 * mrSges(5,1) + t236 * t197 + (-t188 - t189) * t206 + (-mrSges(5,3) - mrSges(6,3)) * t174 + t266;
t199 = t236 * mrSges(6,1) - t206 * mrSges(6,3);
t200 = t236 * mrSges(5,1) - t206 * mrSges(5,3);
t265 = m(6) * t144 + t173 * mrSges(6,3) + t205 * t188;
t136 = m(5) * t150 + t173 * mrSges(5,3) + t205 * t189 + (-t199 - t200) * t236 + (-mrSges(5,2) - mrSges(6,2)) * t234 + t265;
t128 = t244 * t133 + t240 * t136;
t176 = -t245 * g(1) - t269;
t203 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t241 + Ifges(4,2) * t245) * t237;
t204 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t241 + Ifges(4,4) * t245) * t237;
t274 = mrSges(4,1) * t176 - mrSges(4,2) * t177 + Ifges(4,5) * t211 + Ifges(4,6) * t212 + Ifges(4,3) * qJDD(3) + pkin(3) * t128 + (t241 * t203 - t245 * t204) * t237 + t275;
t270 = t237 * t245;
t210 = (-mrSges(4,1) * t245 + mrSges(4,2) * t241) * t237;
t219 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t270;
t126 = m(4) * t176 + qJDD(3) * mrSges(4,1) - t211 * mrSges(4,3) + qJD(3) * t219 - t210 * t271 + t128;
t218 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t271;
t261 = -t240 * t133 + t244 * t136;
t127 = m(4) * t177 - qJDD(3) * mrSges(4,2) + t212 * mrSges(4,3) - qJD(3) * t218 + t210 * t270 + t261;
t262 = -t241 * t126 + t245 * t127;
t118 = m(3) * t194 - t233 * mrSges(3,1) - t235 * mrSges(3,2) + t262;
t193 = t246 * t216 - t242 * t217;
t258 = -t235 * pkin(2) - t193;
t190 = -t233 * pkin(7) + t258;
t156 = -t212 * pkin(3) + t220 * t271 + (-pkin(8) * t239 - pkin(7)) * t233 + t258;
t147 = -t173 * pkin(4) - t201 * qJ(5) + t206 * t198 + qJDD(5) + t156;
t260 = m(6) * t147 - t173 * mrSges(6,1) + t174 * mrSges(6,2) - t205 * t196 + t206 * t199;
t252 = m(5) * t156 - t173 * mrSges(5,1) + t174 * mrSges(5,2) - t205 * t197 + t206 * t200 + t260;
t250 = -m(4) * t190 + t212 * mrSges(4,1) - t211 * mrSges(4,2) - t218 * t271 + t219 * t270 - t252;
t130 = m(3) * t193 + t235 * mrSges(3,1) - t233 * mrSges(3,2) + t250;
t115 = t242 * t118 + t246 * t130;
t120 = t245 * t126 + t241 * t127;
t263 = t246 * t118 - t242 * t130;
t257 = -mrSges(6,1) * t147 + mrSges(6,3) * t144 + Ifges(6,4) * t174 + Ifges(6,2) * t173 + Ifges(6,6) * t234 + t236 * t182;
t178 = Ifges(6,5) * t206 + Ifges(6,6) * t205 + Ifges(6,3) * t236;
t255 = mrSges(6,2) * t147 - mrSges(6,3) * t141 + Ifges(6,1) * t174 + Ifges(6,4) * t173 + Ifges(6,5) * t234 + t205 * t178;
t179 = Ifges(5,5) * t206 + Ifges(5,6) * t205 + Ifges(5,3) * t236;
t121 = Ifges(5,4) * t174 + Ifges(5,2) * t173 + Ifges(5,6) * t234 + t236 * t183 - mrSges(5,1) * t156 + mrSges(5,3) * t150 - pkin(4) * t260 + qJ(5) * (-t234 * mrSges(6,2) - t236 * t199 + t265) + (-t179 - t178) * t206 + t257;
t122 = mrSges(5,2) * t156 - mrSges(5,3) * t149 + Ifges(5,1) * t174 + Ifges(5,4) * t173 + Ifges(5,5) * t234 - qJ(5) * t138 + t205 * t179 + (-t180 - t181) * t236 + t255;
t202 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t241 + Ifges(4,6) * t245) * t237;
t109 = -mrSges(4,1) * t190 + mrSges(4,3) * t177 + Ifges(4,4) * t211 + Ifges(4,2) * t212 + Ifges(4,6) * qJDD(3) - pkin(3) * t252 + pkin(8) * t261 + qJD(3) * t204 + t244 * t121 + t240 * t122 - t202 * t271;
t111 = mrSges(4,2) * t190 - mrSges(4,3) * t176 + Ifges(4,1) * t211 + Ifges(4,4) * t212 + Ifges(4,5) * qJDD(3) - pkin(8) * t128 - qJD(3) * t203 - t240 * t121 + t244 * t122 + t202 * t270;
t254 = mrSges(3,1) * t193 - mrSges(3,2) * t194 + Ifges(3,3) * t235 + pkin(2) * t250 + pkin(7) * t262 + t245 * t109 + t241 * t111;
t253 = mrSges(2,1) * t223 - mrSges(2,2) * t222 + Ifges(2,3) * qJDD(1) + pkin(1) * t115 + t254;
t113 = m(2) * t223 + qJDD(1) * mrSges(2,1) - t248 * mrSges(2,2) + t115;
t112 = m(2) * t222 - t248 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t263;
t107 = mrSges(3,1) * g(1) + mrSges(3,3) * t194 + t233 * Ifges(3,5) + Ifges(3,6) * t235 - pkin(2) * t120 - t274;
t106 = -mrSges(3,2) * g(1) - mrSges(3,3) * t193 + Ifges(3,5) * t235 - t233 * Ifges(3,6) - pkin(7) * t120 - t241 * t109 + t245 * t111;
t105 = -mrSges(2,2) * g(1) - mrSges(2,3) * t223 + Ifges(2,5) * qJDD(1) - t248 * Ifges(2,6) - pkin(6) * t115 + t246 * t106 - t242 * t107;
t104 = Ifges(2,6) * qJDD(1) + t248 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t222 + t242 * t106 + t246 * t107 - pkin(1) * (-m(3) * g(1) + t120) + pkin(6) * t263;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t253, t105, t106, t111, t122, -t236 * t180 + t255; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t243 * t105 - t247 * t104 - pkin(5) * (t247 * t112 - t243 * t113), t104, t107, t109, t121, -t206 * t178 + t257; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t247 * t105 - t243 * t104 + pkin(5) * (-t243 * t112 - t247 * t113), t253, t254, t274, t275, -t205 * t182 - t256;];
m_new = t1;
