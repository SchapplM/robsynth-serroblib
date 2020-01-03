% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR15_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:42
% EndTime: 2019-12-31 18:36:49
% DurationCPUTime: 3.79s
% Computational Cost: add. (47310->266), mult. (97754->327), div. (0->0), fcn. (58218->8), ass. (0->108)
t228 = sin(qJ(1));
t231 = cos(qJ(1));
t213 = -t231 * g(1) - t228 * g(2);
t259 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t213;
t258 = (-pkin(1) - pkin(6));
t257 = mrSges(2,1) - mrSges(3,2);
t256 = -Ifges(3,4) + Ifges(2,5);
t255 = (Ifges(3,5) - Ifges(2,6));
t233 = qJD(1) ^ 2;
t177 = (t258 * t233) - t259;
t227 = sin(qJ(3));
t230 = cos(qJ(3));
t252 = qJD(1) * qJD(3);
t249 = t230 * t252;
t207 = t227 * qJDD(1) + t249;
t250 = t227 * t252;
t208 = t230 * qJDD(1) - t250;
t158 = (-t208 + t250) * qJ(4) + (t207 + t249) * pkin(3) + t177;
t212 = t228 * g(1) - t231 * g(2);
t245 = -t233 * qJ(2) + qJDD(2) - t212;
t183 = t258 * qJDD(1) + t245;
t173 = -t230 * g(3) + t227 * t183;
t205 = (pkin(3) * t227 - qJ(4) * t230) * qJD(1);
t232 = qJD(3) ^ 2;
t253 = t227 * qJD(1);
t161 = -t232 * pkin(3) + qJDD(3) * qJ(4) - t205 * t253 + t173;
t224 = sin(pkin(8));
t225 = cos(pkin(8));
t254 = qJD(1) * t230;
t200 = t224 * qJD(3) + t225 * t254;
t142 = -0.2e1 * qJD(4) * t200 + t225 * t158 - t224 * t161;
t181 = t224 * qJDD(3) + t225 * t208;
t199 = t225 * qJD(3) - t224 * t254;
t140 = (t199 * t253 - t181) * pkin(7) + (t199 * t200 + t207) * pkin(4) + t142;
t143 = 0.2e1 * qJD(4) * t199 + t224 * t158 + t225 * t161;
t180 = t225 * qJDD(3) - t224 * t208;
t182 = pkin(4) * t253 - t200 * pkin(7);
t198 = t199 ^ 2;
t141 = -t198 * pkin(4) + t180 * pkin(7) - t182 * t253 + t143;
t226 = sin(qJ(5));
t229 = cos(qJ(5));
t138 = t229 * t140 - t226 * t141;
t168 = t229 * t199 - t226 * t200;
t150 = t168 * qJD(5) + t226 * t180 + t229 * t181;
t169 = t226 * t199 + t229 * t200;
t155 = -t168 * mrSges(6,1) + t169 * mrSges(6,2);
t214 = qJD(5) + t253;
t162 = -t214 * mrSges(6,2) + t168 * mrSges(6,3);
t204 = qJDD(5) + t207;
t133 = m(6) * t138 + t204 * mrSges(6,1) - t150 * mrSges(6,3) - t169 * t155 + t214 * t162;
t139 = t226 * t140 + t229 * t141;
t149 = -t169 * qJD(5) + t229 * t180 - t226 * t181;
t163 = t214 * mrSges(6,1) - t169 * mrSges(6,3);
t134 = m(6) * t139 - t204 * mrSges(6,2) + t149 * mrSges(6,3) + t168 * t155 - t214 * t163;
t126 = t229 * t133 + t226 * t134;
t170 = -t199 * mrSges(5,1) + t200 * mrSges(5,2);
t178 = -mrSges(5,2) * t253 + t199 * mrSges(5,3);
t124 = m(5) * t142 + t207 * mrSges(5,1) - t181 * mrSges(5,3) - t200 * t170 + t178 * t253 + t126;
t179 = mrSges(5,1) * t253 - t200 * mrSges(5,3);
t247 = -t226 * t133 + t229 * t134;
t125 = m(5) * t143 - t207 * mrSges(5,2) + t180 * mrSges(5,3) + t199 * t170 - t179 * t253 + t247;
t119 = t225 * t124 + t224 * t125;
t206 = (mrSges(4,1) * t227 + mrSges(4,2) * t230) * qJD(1);
t211 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t254;
t248 = -t224 * t124 + t225 * t125;
t117 = m(4) * t173 - qJDD(3) * mrSges(4,2) - t207 * mrSges(4,3) - qJD(3) * t211 - t206 * t253 + t248;
t172 = t227 * g(3) + t230 * t183;
t210 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t253;
t160 = -qJDD(3) * pkin(3) - t232 * qJ(4) + t205 * t254 + qJDD(4) - t172;
t144 = -t180 * pkin(4) - t198 * pkin(7) + t200 * t182 + t160;
t243 = m(6) * t144 - t149 * mrSges(6,1) + t150 * mrSges(6,2) - t168 * t162 + t169 * t163;
t235 = -m(5) * t160 + t180 * mrSges(5,1) - t181 * mrSges(5,2) + t199 * t178 - t200 * t179 - t243;
t129 = m(4) * t172 + qJDD(3) * mrSges(4,1) - t208 * mrSges(4,3) + qJD(3) * t210 - t206 * t254 + t235;
t110 = t230 * t117 - t227 * t129;
t109 = t227 * t117 + t230 * t129;
t188 = -qJDD(1) * pkin(1) + t245;
t244 = -m(3) * t188 + (t233 * mrSges(3,3)) - t109;
t115 = -m(4) * t177 - t207 * mrSges(4,1) - t208 * mrSges(4,2) - t210 * t253 - t211 * t254 - t119;
t152 = Ifges(6,4) * t169 + Ifges(6,2) * t168 + Ifges(6,6) * t214;
t153 = Ifges(6,1) * t169 + Ifges(6,4) * t168 + Ifges(6,5) * t214;
t242 = -mrSges(6,1) * t138 + mrSges(6,2) * t139 - Ifges(6,5) * t150 - Ifges(6,6) * t149 - Ifges(6,3) * t204 - t169 * t152 + t168 * t153;
t151 = Ifges(6,5) * t169 + Ifges(6,6) * t168 + Ifges(6,3) * t214;
t127 = -mrSges(6,1) * t144 + mrSges(6,3) * t139 + Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t204 - t169 * t151 + t214 * t153;
t128 = mrSges(6,2) * t144 - mrSges(6,3) * t138 + Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t204 + t168 * t151 - t214 * t152;
t164 = Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * t253;
t166 = Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * t253;
t106 = -mrSges(5,1) * t160 + mrSges(5,3) * t143 + Ifges(5,4) * t181 + Ifges(5,2) * t180 + Ifges(5,6) * t207 - pkin(4) * t243 + pkin(7) * t247 + t229 * t127 + t226 * t128 - t200 * t164 + t166 * t253;
t165 = Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * t253;
t112 = mrSges(5,2) * t160 - mrSges(5,3) * t142 + Ifges(5,1) * t181 + Ifges(5,4) * t180 + Ifges(5,5) * t207 - pkin(7) * t126 - t226 * t127 + t229 * t128 + t199 * t164 - t165 * t253;
t192 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t230 - Ifges(4,6) * t227) * qJD(1);
t193 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t230 - Ifges(4,2) * t227) * qJD(1);
t102 = mrSges(4,2) * t177 - mrSges(4,3) * t172 + Ifges(4,1) * t208 - Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3) - qJ(4) * t119 - qJD(3) * t193 - t224 * t106 + t225 * t112 - t192 * t253;
t194 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t230 - Ifges(4,4) * t227) * qJD(1);
t234 = -mrSges(5,1) * t142 + mrSges(5,2) * t143 - Ifges(5,5) * t181 - Ifges(5,6) * t180 - pkin(4) * t126 - t200 * t165 + t199 * t166 + t242;
t103 = t234 + Ifges(4,6) * qJDD(3) + (-Ifges(4,2) - Ifges(5,3)) * t207 + Ifges(4,4) * t208 + qJD(3) * t194 + mrSges(4,3) * t173 - mrSges(4,1) * t177 - pkin(3) * t119 - t192 * t254;
t186 = t233 * pkin(1) + t259;
t241 = mrSges(3,2) * t188 - mrSges(3,3) * t186 + Ifges(3,1) * qJDD(1) - pkin(6) * t109 + t230 * t102 - t227 * t103;
t240 = -mrSges(3,1) * t186 - pkin(2) * t115 - pkin(6) * t110 - t227 * t102 - t230 * t103;
t239 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,5) * t208 - Ifges(4,6) * t207 + Ifges(4,3) * qJDD(3) + pkin(3) * t235 + qJ(4) * t248 + t225 * t106 + t224 * t112 + t193 * t254 + t194 * t253;
t238 = -m(3) * t186 + t233 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t115;
t237 = -mrSges(2,2) * t213 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t244) + qJ(2) * t238 + mrSges(2,1) * t212 + Ifges(2,3) * qJDD(1) + t241;
t236 = mrSges(3,1) * t188 + pkin(2) * t109 + t239;
t113 = m(2) * t213 - t233 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t238;
t108 = -m(3) * g(3) + t110;
t104 = m(2) * t212 - t233 * mrSges(2,2) + t257 * qJDD(1) + t244;
t100 = (t255 * t233) + t256 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t212 + t236 - qJ(2) * t108;
t99 = mrSges(2,3) * t213 - pkin(1) * t108 + t257 * g(3) - t255 * qJDD(1) + t256 * t233 + t240;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t231 * t100 - t228 * t99 - pkin(5) * (t231 * t104 + t228 * t113), t100, t241, t102, t112, t128; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t228 * t100 + t231 * t99 + pkin(5) * (-t228 * t104 + t231 * t113), t99, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t233 * Ifges(3,5)) - t236, t103, t106, t127; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t237, t237, mrSges(3,2) * g(3) + t233 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t240, t239, Ifges(5,3) * t207 - t234, -t242;];
m_new = t1;
