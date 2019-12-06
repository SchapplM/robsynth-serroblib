% Calculate vector of cutting torques with Newton-Euler for
% S5PRPRP1
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:23
% EndTime: 2019-12-05 15:28:28
% DurationCPUTime: 2.67s
% Computational Cost: add. (28215->233), mult. (62889->285), div. (0->0), fcn. (40544->8), ass. (0->101)
t229 = qJD(2) ^ 2;
t221 = sin(pkin(8));
t257 = qJD(2) * t221;
t223 = cos(pkin(8));
t267 = qJD(2) * t223;
t222 = sin(pkin(7));
t224 = cos(pkin(7));
t201 = t222 * g(1) - t224 * g(2);
t202 = -t224 * g(1) - t222 * g(2);
t226 = sin(qJ(2));
t227 = cos(qJ(2));
t184 = t226 * t201 + t227 * t202;
t181 = -t229 * pkin(2) + qJDD(2) * qJ(3) + t184;
t220 = -g(3) + qJDD(1);
t254 = qJD(2) * qJD(3);
t258 = t223 * t220 - 0.2e1 * t221 * t254;
t264 = pkin(3) * t223;
t148 = (-pkin(6) * qJDD(2) + t229 * t264 - t181) * t221 + t258;
t154 = t221 * t220 + (t181 + 0.2e1 * t254) * t223;
t252 = qJDD(2) * t223;
t213 = t223 ^ 2;
t261 = t213 * t229;
t149 = -pkin(3) * t261 + pkin(6) * t252 + t154;
t225 = sin(qJ(4));
t265 = cos(qJ(4));
t145 = t225 * t148 + t265 * t149;
t250 = t223 * t265;
t253 = qJDD(2) * t221;
t238 = t265 * t221 + t223 * t225;
t192 = t238 * qJD(2);
t255 = t192 * qJD(4);
t179 = -qJDD(2) * t250 + t225 * t253 + t255;
t188 = qJD(4) * mrSges(5,1) - t192 * mrSges(5,3);
t191 = -qJD(2) * t250 + t225 * t257;
t166 = t191 * pkin(4) - t192 * qJ(5);
t228 = qJD(4) ^ 2;
t138 = -t228 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t191 * t166 + t145;
t189 = -qJD(4) * mrSges(6,1) + t192 * mrSges(6,2);
t251 = m(6) * t138 + qJDD(4) * mrSges(6,3) + qJD(4) * t189;
t167 = t191 * mrSges(6,1) - t192 * mrSges(6,3);
t259 = -t191 * mrSges(5,1) - t192 * mrSges(5,2) - t167;
t263 = -mrSges(5,3) - mrSges(6,2);
t129 = m(5) * t145 - qJDD(4) * mrSges(5,2) - qJD(4) * t188 + t263 * t179 + t259 * t191 + t251;
t144 = t265 * t148 - t225 * t149;
t256 = t191 * qJD(4);
t180 = t238 * qJDD(2) - t256;
t187 = -qJD(4) * mrSges(5,2) - t191 * mrSges(5,3);
t140 = -qJDD(4) * pkin(4) - t228 * qJ(5) + t192 * t166 + qJDD(5) - t144;
t190 = -t191 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t246 = -m(6) * t140 + qJDD(4) * mrSges(6,1) + qJD(4) * t190;
t130 = m(5) * t144 + qJDD(4) * mrSges(5,1) + qJD(4) * t187 + t263 * t180 + t259 * t192 + t246;
t122 = t225 * t129 + t265 * t130;
t153 = -t221 * t181 + t258;
t159 = Ifges(5,4) * t192 - Ifges(5,2) * t191 + Ifges(5,6) * qJD(4);
t161 = Ifges(5,1) * t192 - Ifges(5,4) * t191 + Ifges(5,5) * qJD(4);
t156 = Ifges(6,5) * t192 + Ifges(6,6) * qJD(4) + Ifges(6,3) * t191;
t160 = Ifges(6,1) * t192 + Ifges(6,4) * qJD(4) + Ifges(6,5) * t191;
t235 = mrSges(6,1) * t140 - mrSges(6,3) * t138 - Ifges(6,4) * t180 - Ifges(6,2) * qJDD(4) - Ifges(6,6) * t179 + t192 * t156 - t191 * t160;
t231 = mrSges(5,2) * t145 - t191 * t161 - qJ(5) * (-t179 * mrSges(6,2) - t191 * t167 + t251) - pkin(4) * (-t180 * mrSges(6,2) - t192 * t167 + t246) - mrSges(5,1) * t144 - t192 * t159 + Ifges(5,6) * t179 - Ifges(5,5) * t180 - Ifges(5,3) * qJDD(4) + t235;
t243 = Ifges(4,4) * t221 + Ifges(4,2) * t223;
t244 = Ifges(4,1) * t221 + Ifges(4,4) * t223;
t266 = -mrSges(4,1) * t153 + mrSges(4,2) * t154 - pkin(3) * t122 - (t243 * t257 - t244 * t267) * qJD(2) + t231;
t262 = mrSges(4,2) * t221;
t239 = mrSges(4,3) * qJDD(2) + t229 * (-mrSges(4,1) * t223 + t262);
t120 = m(4) * t153 - t239 * t221 + t122;
t247 = t265 * t129 - t225 * t130;
t121 = m(4) * t154 + t239 * t223 + t247;
t248 = -t221 * t120 + t223 * t121;
t112 = m(3) * t184 - t229 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t248;
t183 = t227 * t201 - t226 * t202;
t241 = qJDD(3) - t183;
t178 = -qJDD(2) * pkin(2) - t229 * qJ(3) + t241;
t212 = t221 ^ 2;
t152 = (-pkin(2) - t264) * qJDD(2) + (-qJ(3) + (-t212 - t213) * pkin(6)) * t229 + t241;
t142 = -0.2e1 * qJD(5) * t192 + (-t180 + t256) * qJ(5) + (t179 + t255) * pkin(4) + t152;
t131 = m(6) * t142 + t179 * mrSges(6,1) - t180 * mrSges(6,3) - t192 * t189 + t191 * t190;
t233 = m(5) * t152 + t179 * mrSges(5,1) + t180 * mrSges(5,2) + t191 * t187 + t192 * t188 + t131;
t232 = -m(4) * t178 + mrSges(4,1) * t252 - t233 + (t212 * t229 + t261) * mrSges(4,3);
t124 = t232 + m(3) * t183 + (mrSges(3,1) - t262) * qJDD(2) - t229 * mrSges(3,2);
t109 = t226 * t112 + t227 * t124;
t114 = t223 * t120 + t221 * t121;
t158 = Ifges(6,4) * t192 + Ifges(6,2) * qJD(4) + Ifges(6,6) * t191;
t260 = -Ifges(5,5) * t192 + Ifges(5,6) * t191 - Ifges(5,3) * qJD(4) - t158;
t249 = t227 * t112 - t226 * t124;
t245 = -mrSges(6,1) * t142 + mrSges(6,2) * t138;
t242 = Ifges(4,5) * t221 + Ifges(4,6) * t223;
t237 = mrSges(6,2) * t140 - mrSges(6,3) * t142 + Ifges(6,1) * t180 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t179 + qJD(4) * t156;
t115 = -mrSges(5,1) * t152 + mrSges(5,3) * t145 - pkin(4) * t131 + t260 * t192 + (Ifges(5,4) - Ifges(6,5)) * t180 + (-Ifges(5,2) - Ifges(6,3)) * t179 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t160 + t161) * qJD(4) + t245;
t116 = mrSges(5,2) * t152 - mrSges(5,3) * t144 + Ifges(5,1) * t180 - Ifges(5,4) * t179 + Ifges(5,5) * qJDD(4) - qJ(5) * t131 - qJD(4) * t159 + t260 * t191 + t237;
t197 = t242 * qJD(2);
t103 = -mrSges(4,1) * t178 + mrSges(4,3) * t154 - pkin(3) * t233 + pkin(6) * t247 + t243 * qJDD(2) + t265 * t115 + t225 * t116 - t197 * t257;
t105 = mrSges(4,2) * t178 - mrSges(4,3) * t153 - pkin(6) * t122 + t244 * qJDD(2) - t225 * t115 + t265 * t116 + t197 * t267;
t236 = -mrSges(3,2) * t184 + qJ(3) * t248 + t223 * t103 + t221 * t105 + pkin(2) * (-mrSges(4,2) * t253 + t232) + mrSges(3,1) * t183 + Ifges(3,3) * qJDD(2);
t234 = mrSges(2,1) * t201 - mrSges(2,2) * t202 + pkin(1) * t109 + t236;
t107 = m(2) * t202 + t249;
t106 = m(2) * t201 + t109;
t101 = -pkin(2) * t114 - mrSges(3,1) * t220 + mrSges(3,3) * t184 + t229 * Ifges(3,5) + (Ifges(3,6) - t242) * qJDD(2) + t266;
t100 = mrSges(3,2) * t220 - mrSges(3,3) * t183 + Ifges(3,5) * qJDD(2) - t229 * Ifges(3,6) - qJ(3) * t114 - t221 * t103 + t223 * t105;
t99 = mrSges(2,2) * t220 - mrSges(2,3) * t201 - pkin(5) * t109 + t227 * t100 - t226 * t101;
t98 = -mrSges(2,1) * t220 + mrSges(2,3) * t202 + t226 * t100 + t227 * t101 - pkin(1) * (m(3) * t220 + t114) + pkin(5) * t249;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t224 * t99 - t222 * t98 - qJ(1) * (t224 * t106 + t222 * t107), t99, t100, t105, t116, -t191 * t158 + t237; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t222 * t99 + t224 * t98 + qJ(1) * (-t222 * t106 + t224 * t107), t98, t101, t103, t115, -t235; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t234, t234, t236, t242 * qJDD(2) - t266, -t231, Ifges(6,5) * t180 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t179 - qJD(4) * t160 + t192 * t158 - t245;];
m_new = t1;
