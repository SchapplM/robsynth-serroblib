% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:58
% EndTime: 2020-01-03 12:09:07
% DurationCPUTime: 6.12s
% Computational Cost: add. (147781->269), mult. (197708->345), div. (0->0), fcn. (125769->10), ass. (0->111)
t227 = qJDD(1) + qJDD(2);
t236 = sin(qJ(3));
t240 = cos(qJ(3));
t229 = qJD(1) + qJD(2);
t258 = qJD(3) * t229;
t209 = t236 * t227 + t240 * t258;
t238 = sin(qJ(1));
t242 = cos(qJ(1));
t221 = -t242 * g(2) - t238 * g(3);
t214 = qJDD(1) * pkin(1) + t221;
t220 = -t238 * g(2) + t242 * g(3);
t243 = qJD(1) ^ 2;
t215 = -t243 * pkin(1) + t220;
t237 = sin(qJ(2));
t241 = cos(qJ(2));
t193 = t237 * t214 + t241 * t215;
t225 = t229 ^ 2;
t185 = -t225 * pkin(2) + t227 * pkin(7) + t193;
t259 = t236 * t185;
t262 = pkin(3) * t225;
t168 = qJDD(3) * pkin(3) - t209 * qJ(4) - t259 + (qJ(4) * t258 + t236 * t262 - g(1)) * t240;
t175 = -t236 * g(1) + t240 * t185;
t210 = t240 * t227 - t236 * t258;
t261 = t229 * t236;
t216 = qJD(3) * pkin(3) - qJ(4) * t261;
t232 = t240 ^ 2;
t169 = t210 * qJ(4) - qJD(3) * t216 - t232 * t262 + t175;
t233 = sin(pkin(9));
t234 = cos(pkin(9));
t201 = (t233 * t240 + t234 * t236) * t229;
t148 = -0.2e1 * qJD(4) * t201 + t234 * t168 - t233 * t169;
t190 = t234 * t209 + t233 * t210;
t200 = (-t233 * t236 + t234 * t240) * t229;
t145 = (qJD(3) * t200 - t190) * pkin(8) + (t200 * t201 + qJDD(3)) * pkin(4) + t148;
t149 = 0.2e1 * qJD(4) * t200 + t233 * t168 + t234 * t169;
t189 = -t233 * t209 + t234 * t210;
t196 = qJD(3) * pkin(4) - t201 * pkin(8);
t199 = t200 ^ 2;
t146 = -t199 * pkin(4) + t189 * pkin(8) - qJD(3) * t196 + t149;
t235 = sin(qJ(5));
t239 = cos(qJ(5));
t143 = t239 * t145 - t235 * t146;
t179 = t239 * t200 - t235 * t201;
t158 = t179 * qJD(5) + t235 * t189 + t239 * t190;
t180 = t235 * t200 + t239 * t201;
t164 = -t179 * mrSges(6,1) + t180 * mrSges(6,2);
t228 = qJD(3) + qJD(5);
t172 = -t228 * mrSges(6,2) + t179 * mrSges(6,3);
t226 = qJDD(3) + qJDD(5);
t140 = m(6) * t143 + t226 * mrSges(6,1) - t158 * mrSges(6,3) - t180 * t164 + t228 * t172;
t144 = t235 * t145 + t239 * t146;
t157 = -t180 * qJD(5) + t239 * t189 - t235 * t190;
t173 = t228 * mrSges(6,1) - t180 * mrSges(6,3);
t141 = m(6) * t144 - t226 * mrSges(6,2) + t157 * mrSges(6,3) + t179 * t164 - t228 * t173;
t131 = t239 * t140 + t235 * t141;
t183 = -t200 * mrSges(5,1) + t201 * mrSges(5,2);
t194 = -qJD(3) * mrSges(5,2) + t200 * mrSges(5,3);
t128 = m(5) * t148 + qJDD(3) * mrSges(5,1) - t190 * mrSges(5,3) + qJD(3) * t194 - t201 * t183 + t131;
t195 = qJD(3) * mrSges(5,1) - t201 * mrSges(5,3);
t254 = -t235 * t140 + t239 * t141;
t129 = m(5) * t149 - qJDD(3) * mrSges(5,2) + t189 * mrSges(5,3) - qJD(3) * t195 + t200 * t183 + t254;
t124 = t234 * t128 + t233 * t129;
t174 = -t240 * g(1) - t259;
t203 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t236 + Ifges(4,2) * t240) * t229;
t204 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t236 + Ifges(4,4) * t240) * t229;
t177 = Ifges(5,4) * t201 + Ifges(5,2) * t200 + Ifges(5,6) * qJD(3);
t178 = Ifges(5,1) * t201 + Ifges(5,4) * t200 + Ifges(5,5) * qJD(3);
t160 = Ifges(6,4) * t180 + Ifges(6,2) * t179 + Ifges(6,6) * t228;
t161 = Ifges(6,1) * t180 + Ifges(6,4) * t179 + Ifges(6,5) * t228;
t249 = -mrSges(6,1) * t143 + mrSges(6,2) * t144 - Ifges(6,5) * t158 - Ifges(6,6) * t157 - Ifges(6,3) * t226 - t180 * t160 + t179 * t161;
t246 = -mrSges(5,1) * t148 + mrSges(5,2) * t149 - Ifges(5,5) * t190 - Ifges(5,6) * t189 - Ifges(5,3) * qJDD(3) - pkin(4) * t131 - t201 * t177 + t200 * t178 + t249;
t263 = mrSges(4,1) * t174 - mrSges(4,2) * t175 + Ifges(4,5) * t209 + Ifges(4,6) * t210 + Ifges(4,3) * qJDD(3) + pkin(3) * t124 + (t236 * t203 - t240 * t204) * t229 - t246;
t260 = t229 * t240;
t208 = (-mrSges(4,1) * t240 + mrSges(4,2) * t236) * t229;
t218 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t260;
t122 = m(4) * t174 + qJDD(3) * mrSges(4,1) - t209 * mrSges(4,3) + qJD(3) * t218 - t208 * t261 + t124;
t217 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t261;
t255 = -t233 * t128 + t234 * t129;
t123 = m(4) * t175 - qJDD(3) * mrSges(4,2) + t210 * mrSges(4,3) - qJD(3) * t217 + t208 * t260 + t255;
t256 = -t236 * t122 + t240 * t123;
t114 = m(3) * t193 - t225 * mrSges(3,1) - t227 * mrSges(3,2) + t256;
t192 = t241 * t214 - t237 * t215;
t251 = -t227 * pkin(2) - t192;
t184 = -t225 * pkin(7) + t251;
t170 = -t210 * pkin(3) + qJDD(4) + t216 * t261 + (-qJ(4) * t232 - pkin(7)) * t225 + t251;
t151 = -t189 * pkin(4) - t199 * pkin(8) + t201 * t196 + t170;
t253 = m(6) * t151 - t157 * mrSges(6,1) + t158 * mrSges(6,2) - t179 * t172 + t180 * t173;
t247 = m(5) * t170 - t189 * mrSges(5,1) + t190 * mrSges(5,2) - t200 * t194 + t201 * t195 + t253;
t245 = -m(4) * t184 + t210 * mrSges(4,1) - t209 * mrSges(4,2) - t217 * t261 + t218 * t260 - t247;
t135 = m(3) * t192 + t227 * mrSges(3,1) - t225 * mrSges(3,2) + t245;
t111 = t237 * t114 + t241 * t135;
t116 = t240 * t122 + t236 * t123;
t257 = t241 * t114 - t237 * t135;
t159 = Ifges(6,5) * t180 + Ifges(6,6) * t179 + Ifges(6,3) * t228;
t132 = -mrSges(6,1) * t151 + mrSges(6,3) * t144 + Ifges(6,4) * t158 + Ifges(6,2) * t157 + Ifges(6,6) * t226 - t180 * t159 + t228 * t161;
t133 = mrSges(6,2) * t151 - mrSges(6,3) * t143 + Ifges(6,1) * t158 + Ifges(6,4) * t157 + Ifges(6,5) * t226 + t179 * t159 - t228 * t160;
t176 = Ifges(5,5) * t201 + Ifges(5,6) * t200 + Ifges(5,3) * qJD(3);
t117 = -mrSges(5,1) * t170 + mrSges(5,3) * t149 + Ifges(5,4) * t190 + Ifges(5,2) * t189 + Ifges(5,6) * qJDD(3) - pkin(4) * t253 + pkin(8) * t254 + qJD(3) * t178 + t239 * t132 + t235 * t133 - t201 * t176;
t118 = mrSges(5,2) * t170 - mrSges(5,3) * t148 + Ifges(5,1) * t190 + Ifges(5,4) * t189 + Ifges(5,5) * qJDD(3) - pkin(8) * t131 - qJD(3) * t177 - t235 * t132 + t239 * t133 + t200 * t176;
t202 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t236 + Ifges(4,6) * t240) * t229;
t105 = -mrSges(4,1) * t184 + mrSges(4,3) * t175 + Ifges(4,4) * t209 + Ifges(4,2) * t210 + Ifges(4,6) * qJDD(3) - pkin(3) * t247 + qJ(4) * t255 + qJD(3) * t204 + t234 * t117 + t233 * t118 - t202 * t261;
t107 = mrSges(4,2) * t184 - mrSges(4,3) * t174 + Ifges(4,1) * t209 + Ifges(4,4) * t210 + Ifges(4,5) * qJDD(3) - qJ(4) * t124 - qJD(3) * t203 - t233 * t117 + t234 * t118 + t202 * t260;
t250 = mrSges(3,1) * t192 - mrSges(3,2) * t193 + Ifges(3,3) * t227 + pkin(2) * t245 + pkin(7) * t256 + t240 * t105 + t236 * t107;
t248 = mrSges(2,1) * t221 - mrSges(2,2) * t220 + Ifges(2,3) * qJDD(1) + pkin(1) * t111 + t250;
t109 = m(2) * t221 + qJDD(1) * mrSges(2,1) - t243 * mrSges(2,2) + t111;
t108 = m(2) * t220 - t243 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t257;
t103 = mrSges(3,1) * g(1) + mrSges(3,3) * t193 + t225 * Ifges(3,5) + Ifges(3,6) * t227 - pkin(2) * t116 - t263;
t102 = -mrSges(3,2) * g(1) - mrSges(3,3) * t192 + Ifges(3,5) * t227 - t225 * Ifges(3,6) - pkin(7) * t116 - t236 * t105 + t240 * t107;
t101 = -mrSges(2,2) * g(1) - mrSges(2,3) * t221 + Ifges(2,5) * qJDD(1) - t243 * Ifges(2,6) - pkin(6) * t111 + t241 * t102 - t237 * t103;
t100 = Ifges(2,6) * qJDD(1) + t243 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t220 + t237 * t102 + t241 * t103 - pkin(1) * (-m(3) * g(1) + t116) + pkin(6) * t257;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t248, t101, t102, t107, t118, t133; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t238 * t101 + t242 * t100 - pkin(5) * (-t242 * t108 + t238 * t109), t100, t103, t105, t117, t132; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t242 * t101 + t238 * t100 + pkin(5) * (t238 * t108 + t242 * t109), t248, t250, t263, -t246, -t249;];
m_new = t1;
