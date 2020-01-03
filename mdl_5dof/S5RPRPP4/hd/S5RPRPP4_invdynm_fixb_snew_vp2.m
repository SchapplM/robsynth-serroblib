% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:20
% EndTime: 2019-12-31 18:14:23
% DurationCPUTime: 2.05s
% Computational Cost: add. (21273->265), mult. (45199->313), div. (0->0), fcn. (25040->6), ass. (0->98)
t229 = sin(qJ(1));
t231 = cos(qJ(1));
t208 = -t231 * g(1) - t229 * g(2);
t246 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t208;
t265 = -2 * qJD(4);
t264 = -pkin(1) - pkin(6);
t263 = mrSges(2,1) - mrSges(3,2);
t262 = -mrSges(5,3) - mrSges(6,2);
t261 = Ifges(2,5) - Ifges(3,4);
t260 = (-Ifges(2,6) + Ifges(3,5));
t207 = t229 * g(1) - t231 * g(2);
t233 = qJD(1) ^ 2;
t245 = -t233 * qJ(2) + qJDD(2) - t207;
t181 = qJDD(1) * t264 + t245;
t228 = sin(qJ(3));
t230 = cos(qJ(3));
t169 = t228 * g(3) + t230 * t181;
t255 = qJD(1) * qJD(3);
t252 = t228 * t255;
t202 = qJDD(1) * t230 - t252;
t140 = (-t202 - t252) * qJ(4) + (-t228 * t230 * t233 + qJDD(3)) * pkin(3) + t169;
t170 = -g(3) * t230 + t228 * t181;
t201 = -qJDD(1) * t228 - t230 * t255;
t256 = qJD(1) * t230;
t205 = qJD(3) * pkin(3) - qJ(4) * t256;
t223 = t228 ^ 2;
t141 = -pkin(3) * t223 * t233 + qJ(4) * t201 - qJD(3) * t205 + t170;
t226 = sin(pkin(7));
t227 = cos(pkin(7));
t188 = (t226 * t230 + t227 * t228) * qJD(1);
t137 = t140 * t226 + t141 * t227 + t188 * t265;
t166 = -t201 * t227 + t202 * t226;
t257 = qJD(1) * t228;
t189 = -t226 * t257 + t227 * t256;
t178 = (qJD(3) * mrSges(5,1)) - mrSges(5,3) * t189;
t155 = pkin(4) * t188 - qJ(5) * t189;
t232 = qJD(3) ^ 2;
t130 = -pkin(4) * t232 + qJDD(3) * qJ(5) + (2 * qJD(5) * qJD(3)) - t155 * t188 + t137;
t179 = -(qJD(3) * mrSges(6,1)) + mrSges(6,2) * t189;
t253 = m(6) * t130 + qJDD(3) * mrSges(6,3) + qJD(3) * t179;
t156 = mrSges(6,1) * t188 - mrSges(6,3) * t189;
t258 = -mrSges(5,1) * t188 - mrSges(5,2) * t189 - t156;
t121 = m(5) * t137 - qJDD(3) * mrSges(5,2) - qJD(3) * t178 + t166 * t262 + t188 * t258 + t253;
t248 = -t227 * t140 + t226 * t141;
t136 = t189 * t265 - t248;
t167 = t201 * t226 + t202 * t227;
t177 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t188;
t132 = -qJDD(3) * pkin(4) - t232 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t155) * t189 + t248;
t180 = -mrSges(6,2) * t188 + qJD(3) * mrSges(6,3);
t250 = -m(6) * t132 + qJDD(3) * mrSges(6,1) + qJD(3) * t180;
t122 = m(5) * t136 + qJDD(3) * mrSges(5,1) + qJD(3) * t177 + t167 * t262 + t189 * t258 + t250;
t114 = t121 * t226 + t122 * t227;
t149 = Ifges(6,4) * t189 + (Ifges(6,2) * qJD(3)) + Ifges(6,6) * t188;
t259 = -Ifges(5,5) * t189 + Ifges(5,6) * t188 - (Ifges(5,3) * qJD(3)) - t149;
t200 = (mrSges(4,1) * t228 + mrSges(4,2) * t230) * qJD(1);
t204 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t257;
t111 = m(4) * t169 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t202 + qJD(3) * t204 - t200 * t256 + t114;
t206 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t256;
t251 = t121 * t227 - t122 * t226;
t112 = m(4) * t170 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t201 - qJD(3) * t206 - t200 * t257 + t251;
t107 = -t111 * t228 + t112 * t230;
t143 = -t201 * pkin(3) + qJDD(4) + t205 * t256 + (-qJ(4) * t223 + t264) * t233 + t246;
t134 = -0.2e1 * qJD(5) * t189 + (qJD(3) * t188 - t167) * qJ(5) + (qJD(3) * t189 + t166) * pkin(4) + t143;
t249 = -mrSges(6,1) * t134 + mrSges(6,2) * t130;
t106 = t230 * t111 + t228 * t112;
t147 = Ifges(6,5) * t189 + Ifges(6,6) * qJD(3) + Ifges(6,3) * t188;
t244 = mrSges(6,2) * t132 - mrSges(6,3) * t134 + Ifges(6,1) * t167 + Ifges(6,4) * qJDD(3) + Ifges(6,5) * t166 + qJD(3) * t147;
t124 = m(6) * t134 + t166 * mrSges(6,1) - t167 * mrSges(6,3) - t179 * t189 + t188 * t180;
t187 = -qJDD(1) * pkin(1) + t245;
t243 = -m(3) * t187 + (t233 * mrSges(3,3)) - t106;
t151 = Ifges(6,1) * t189 + Ifges(6,4) * qJD(3) + Ifges(6,5) * t188;
t242 = mrSges(6,1) * t132 - mrSges(6,3) * t130 - Ifges(6,4) * t167 - Ifges(6,2) * qJDD(3) - Ifges(6,6) * t166 + t189 * t147 - t188 * t151;
t152 = Ifges(5,1) * t189 - Ifges(5,4) * t188 + Ifges(5,5) * qJD(3);
t108 = -mrSges(5,1) * t143 + mrSges(5,3) * t137 - pkin(4) * t124 + t259 * t189 + (Ifges(5,4) - Ifges(6,5)) * t167 + (-Ifges(5,2) - Ifges(6,3)) * t166 + (Ifges(5,6) - Ifges(6,6)) * qJDD(3) + (t151 + t152) * qJD(3) + t249;
t150 = Ifges(5,4) * t189 - Ifges(5,2) * t188 + Ifges(5,6) * qJD(3);
t109 = mrSges(5,2) * t143 - mrSges(5,3) * t136 + Ifges(5,1) * t167 - Ifges(5,4) * t166 + Ifges(5,5) * qJDD(3) - qJ(5) * t124 - qJD(3) * t150 + t188 * t259 + t244;
t176 = t233 * t264 + t246;
t190 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t230 - Ifges(4,6) * t228) * qJD(1);
t192 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t230 - Ifges(4,4) * t228) * qJD(1);
t239 = m(5) * t143 + t166 * mrSges(5,1) + t167 * mrSges(5,2) + t177 * t188 + t189 * t178 + t124;
t100 = -mrSges(4,1) * t176 + mrSges(4,3) * t170 + Ifges(4,4) * t202 + Ifges(4,2) * t201 + Ifges(4,6) * qJDD(3) - pkin(3) * t239 + qJ(4) * t251 + qJD(3) * t192 + t227 * t108 + t226 * t109 - t190 * t256;
t191 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t230 - Ifges(4,2) * t228) * qJD(1);
t102 = mrSges(4,2) * t176 - mrSges(4,3) * t169 + Ifges(4,1) * t202 + Ifges(4,4) * t201 + Ifges(4,5) * qJDD(3) - qJ(4) * t114 - qJD(3) * t191 - t108 * t226 + t109 * t227 - t190 * t257;
t184 = t233 * pkin(1) - t246;
t241 = mrSges(3,2) * t187 - mrSges(3,3) * t184 + Ifges(3,1) * qJDD(1) - pkin(6) * t106 - t100 * t228 + t102 * t230;
t117 = -m(4) * t176 + mrSges(4,1) * t201 - t202 * mrSges(4,2) - t204 * t257 - t206 * t256 - t239;
t240 = -mrSges(3,1) * t184 - pkin(2) * t117 - pkin(6) * t107 - t230 * t100 - t228 * t102;
t236 = -m(3) * t184 + t233 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t117;
t238 = -mrSges(2,2) * t208 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t243) + qJ(2) * t236 + mrSges(2,1) * t207 + Ifges(2,3) * qJDD(1) + t241;
t237 = mrSges(5,2) * t137 - t188 * t152 - qJ(5) * (-t166 * mrSges(6,2) - t188 * t156 + t253) - pkin(4) * (-t167 * mrSges(6,2) - t189 * t156 + t250) - mrSges(5,1) * t136 - t189 * t150 + Ifges(5,6) * t166 - Ifges(5,5) * t167 - Ifges(5,3) * qJDD(3) + t242;
t235 = -mrSges(4,1) * t169 + mrSges(4,2) * t170 - Ifges(4,5) * t202 - Ifges(4,6) * t201 - Ifges(4,3) * qJDD(3) - pkin(3) * t114 - t191 * t256 - t192 * t257 + t237;
t234 = -mrSges(3,1) * t187 - pkin(2) * t106 + t235;
t115 = m(2) * t208 - mrSges(2,1) * t233 - qJDD(1) * mrSges(2,2) + t236;
t105 = -m(3) * g(3) + t107;
t103 = m(2) * t207 - t233 * mrSges(2,2) + qJDD(1) * t263 + t243;
t99 = (t260 * t233) + t261 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t207 - qJ(2) * t105 - t234;
t98 = mrSges(2,3) * t208 - pkin(1) * t105 + g(3) * t263 - qJDD(1) * t260 + t233 * t261 + t240;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t231 * t99 - t229 * t98 - pkin(5) * (t103 * t231 + t115 * t229), t99, t241, t102, t109, -t149 * t188 + t244; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t229 * t99 + t231 * t98 + pkin(5) * (-t103 * t229 + t115 * t231), t98, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t233 * Ifges(3,5)) + t234, t100, t108, -t242; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t238, t238, mrSges(3,2) * g(3) + t233 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t240, -t235, -t237, Ifges(6,5) * t167 + Ifges(6,6) * qJDD(3) + Ifges(6,3) * t166 - qJD(3) * t151 + t189 * t149 - t249;];
m_new = t1;
