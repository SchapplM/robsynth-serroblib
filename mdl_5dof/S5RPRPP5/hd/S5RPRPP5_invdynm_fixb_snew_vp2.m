% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:10
% EndTime: 2019-12-31 18:16:12
% DurationCPUTime: 1.34s
% Computational Cost: add. (7842->269), mult. (15175->308), div. (0->0), fcn. (5734->4), ass. (0->94)
t284 = -2 * qJD(1);
t235 = sin(qJ(1));
t237 = cos(qJ(1));
t210 = -t237 * g(1) - t235 * g(2);
t239 = qJD(1) ^ 2;
t153 = (t239 * pkin(1)) - qJDD(1) * qJ(2) + (qJD(2) * t284) - t210;
t145 = -(t239 * pkin(6)) - t153;
t234 = sin(qJ(3));
t236 = cos(qJ(3));
t263 = qJD(1) * qJD(3);
t195 = t234 * qJDD(1) + t236 * t263;
t196 = t236 * qJDD(1) - t234 * t263;
t270 = t196 * qJ(4);
t271 = qJ(4) * t234;
t134 = t195 * pkin(3) - t270 + (-0.2e1 * qJD(4) * t236 + (pkin(3) * t236 + t271) * qJD(3)) * qJD(1) + t145;
t264 = qJD(1) * t236;
t205 = -(qJD(3) * mrSges(6,1)) - mrSges(6,3) * t264;
t207 = -(qJD(3) * mrSges(5,1)) + mrSges(5,2) * t264;
t265 = qJD(1) * t234;
t208 = -mrSges(5,2) * t265 + (qJD(3) * mrSges(5,3));
t204 = -(qJD(3) * pkin(4)) - qJ(5) * t264;
t230 = t234 ^ 2;
t280 = -pkin(3) - pkin(4);
t281 = 0.2e1 * qJD(4);
t126 = t270 + qJDD(5) + (-qJ(5) * t230 + pkin(6)) * t239 + t280 * t195 + (-qJD(3) * t271 + (-pkin(3) * qJD(3) + t204 + t281) * t236) * qJD(1) + t153;
t202 = (qJD(3) * mrSges(6,2)) + mrSges(6,3) * t265;
t259 = m(6) * t126 - t195 * mrSges(6,1) - t202 * t265;
t117 = m(5) * t134 + t195 * mrSges(5,1) - (mrSges(6,2) + mrSges(5,3)) * t196 + t208 * t265 - (t205 + t207) * t264 - t259;
t203 = -(qJD(3) * mrSges(4,2)) - mrSges(4,3) * t265;
t206 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t264;
t283 = -m(4) * t145 - t195 * mrSges(4,1) - t196 * mrSges(4,2) - t203 * t265 - t206 * t264 - t117;
t279 = pkin(4) * t239;
t278 = mrSges(2,1) - mrSges(3,2);
t276 = -mrSges(4,3) - mrSges(5,2);
t275 = Ifges(2,5) - Ifges(3,4);
t274 = (-Ifges(2,6) + Ifges(3,5));
t273 = Ifges(5,6) - Ifges(6,6);
t165 = (Ifges(5,2) * qJD(3)) + (Ifges(5,4) * t236 + Ifges(5,6) * t234) * qJD(1);
t269 = -(Ifges(4,3) * qJD(3)) - (Ifges(4,5) * t236 - Ifges(4,6) * t234) * qJD(1) - t165;
t167 = -(Ifges(6,5) * qJD(3)) + (Ifges(6,1) * t236 + Ifges(6,4) * t234) * qJD(1);
t168 = Ifges(5,4) * qJD(3) + (Ifges(5,1) * t236 + Ifges(5,5) * t234) * qJD(1);
t268 = -t167 - t168;
t266 = Ifges(6,5) * qJDD(3);
t209 = t235 * g(1) - t237 * g(2);
t255 = -t239 * qJ(2) + qJDD(2) - t209;
t146 = (-pkin(1) - pkin(6)) * qJDD(1) + t255;
t142 = -t236 * g(3) + t234 * t146;
t191 = (pkin(3) * t234 - qJ(4) * t236) * qJD(1);
t238 = qJD(3) ^ 2;
t252 = -t238 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t281 + t142;
t131 = -t230 * t279 + t195 * qJ(5) + qJD(3) * t204 + ((2 * qJD(5)) - t191) * t265 + t252;
t193 = (-mrSges(6,1) * t234 + mrSges(6,2) * t236) * qJD(1);
t123 = m(6) * t131 + (qJDD(3) * mrSges(6,2)) + t195 * mrSges(6,3) + qJD(3) * t205 + t193 * t265;
t137 = -t191 * t265 + t252;
t251 = m(5) * t137 + (qJDD(3) * mrSges(5,3)) + qJD(3) * t207 + t123;
t192 = (mrSges(5,1) * t234 - mrSges(5,3) * t236) * qJD(1);
t261 = qJD(1) * (-t192 - (mrSges(4,1) * t234 + mrSges(4,2) * t236) * qJD(1));
t115 = m(4) * t142 - (qJDD(3) * mrSges(4,2)) - qJD(3) * t206 + t276 * t195 + t234 * t261 + t251;
t141 = t234 * g(3) + t236 * t146;
t258 = -t238 * qJ(4) + t191 * t264 + qJDD(4);
t132 = -t196 * qJ(5) + ((qJD(5) * t284) - t146) * t236 + t280 * qJDD(3) + (-qJ(5) * t263 + t236 * t279 - g(3)) * t234 + t258;
t124 = m(6) * t132 - (qJDD(3) * mrSges(6,1)) - t196 * mrSges(6,3) - qJD(3) * t202 - t193 * t264;
t139 = -qJDD(3) * pkin(3) - t141 + t258;
t248 = -m(5) * t139 + (qJDD(3) * mrSges(5,1)) + qJD(3) * t208 - t124;
t116 = m(4) * t141 + (qJDD(3) * mrSges(4,1)) + qJD(3) * t203 + t276 * t196 + t236 * t261 + t248;
t110 = t236 * t115 - t234 * t116;
t109 = t234 * t115 + t236 * t116;
t257 = mrSges(6,1) * t126 - mrSges(6,3) * t131 - Ifges(6,4) * t196 - Ifges(6,2) * t195;
t164 = -(Ifges(6,6) * qJD(3)) + (Ifges(6,4) * t236 + Ifges(6,2) * t234) * qJD(1);
t256 = mrSges(6,1) * t132 - mrSges(6,2) * t131 + Ifges(6,5) * t196 + Ifges(6,6) * t195 - (Ifges(6,3) * qJDD(3)) + t164 * t264 - t167 * t265;
t161 = -(Ifges(6,3) * qJD(3)) + (Ifges(6,5) * t236 + Ifges(6,6) * t234) * qJD(1);
t253 = mrSges(6,2) * t126 - mrSges(6,3) * t132 + Ifges(6,1) * t196 + Ifges(6,4) * t195 + qJD(3) * t164 + t161 * t265;
t157 = -qJDD(1) * pkin(1) + t255;
t250 = -m(3) * t157 + (t239 * mrSges(3,3)) - t109;
t169 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t236 - Ifges(4,4) * t234) * qJD(1);
t246 = mrSges(5,1) * t134 - mrSges(5,2) * t137 + pkin(4) * (-t196 * mrSges(6,2) - t205 * t264 - t259) + qJ(5) * t123 - t257;
t103 = (Ifges(4,4) - Ifges(5,5)) * t196 + (-Ifges(4,2) - Ifges(5,3)) * t195 + (Ifges(4,6) - t273) * qJDD(3) + (t169 - t268) * qJD(3) - t246 + (t161 + t269) * t264 + mrSges(4,3) * t142 - mrSges(4,1) * t145 - pkin(3) * t117;
t166 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t236 - Ifges(4,2) * t234) * qJD(1);
t162 = Ifges(5,6) * qJD(3) + (Ifges(5,5) * t236 + Ifges(5,3) * t234) * qJD(1);
t244 = mrSges(5,2) * t139 - mrSges(5,3) * t134 + Ifges(5,1) * t196 + Ifges(5,4) * qJDD(3) + Ifges(5,5) * t195 - qJ(5) * t124 + qJD(3) * t162 + t253;
t105 = (Ifges(4,5) - Ifges(6,5)) * qJDD(3) + t244 + t269 * t265 - Ifges(4,4) * t195 + Ifges(4,1) * t196 - qJD(3) * t166 - mrSges(4,3) * t141 + mrSges(4,2) * t145 - qJ(4) * t117;
t249 = mrSges(3,2) * t157 - mrSges(3,3) * t153 + Ifges(3,1) * qJDD(1) - pkin(6) * t109 - t234 * t103 + t236 * t105;
t247 = -mrSges(3,1) * t153 - pkin(2) * t283 - pkin(6) * t110 - t236 * t103 - t234 * t105;
t242 = -m(3) * t153 + t239 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t283;
t245 = -mrSges(2,2) * t210 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t250) + qJ(2) * t242 + mrSges(2,1) * t209 + Ifges(2,3) * qJDD(1) + t249;
t243 = mrSges(5,1) * t139 - mrSges(5,3) * t137 - Ifges(5,4) * t196 - (Ifges(5,2) * qJDD(3)) - Ifges(5,6) * t195 + pkin(4) * t124 + t162 * t264 - t168 * t265 + t256;
t241 = mrSges(4,2) * t142 - qJ(4) * (-t195 * mrSges(5,2) - t192 * t265 + t251) - pkin(3) * (-t196 * mrSges(5,2) - t192 * t264 + t248) - mrSges(4,1) * t141 - t169 * t265 - t166 * t264 + Ifges(4,6) * t195 - Ifges(4,5) * t196 - (Ifges(4,3) * qJDD(3)) + t243;
t240 = -mrSges(3,1) * t157 - pkin(2) * t109 + t241;
t111 = m(2) * t210 - t239 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t242;
t108 = -m(3) * g(3) + t110;
t106 = m(2) * t209 - t239 * mrSges(2,2) + t278 * qJDD(1) + t250;
t102 = (t274 * t239) + t275 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - t240 - mrSges(2,3) * t209 - qJ(2) * t108;
t101 = mrSges(2,3) * t210 - pkin(1) * t108 + t278 * g(3) - t274 * qJDD(1) + t275 * t239 + t247;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t237 * t102 - t235 * t101 - pkin(5) * (t237 * t106 + t235 * t111), t102, t249, t105, -t165 * t265 + t244 - t266, t253 - t266; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t235 * t102 + t237 * t101 + pkin(5) * (-t235 * t106 + t237 * t111), t101, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t239 * Ifges(3,5)) + t240, t103, -t243, -(Ifges(6,6) * qJDD(3)) - qJD(3) * t167 - t161 * t264 - t257; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t245, t245, mrSges(3,2) * g(3) + t239 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t247, -t241, Ifges(5,5) * t196 + Ifges(5,3) * t195 + t273 * qJDD(3) + t268 * qJD(3) + (-t161 + t165) * t264 + t246, t256;];
m_new = t1;
