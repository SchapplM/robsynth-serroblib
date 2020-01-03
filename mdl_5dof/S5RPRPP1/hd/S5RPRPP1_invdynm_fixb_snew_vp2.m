% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:42
% EndTime: 2019-12-31 18:08:48
% DurationCPUTime: 3.01s
% Computational Cost: add. (34996->264), mult. (73586->324), div. (0->0), fcn. (42162->8), ass. (0->99)
t231 = sin(qJ(1));
t233 = cos(qJ(1));
t212 = t231 * g(1) - t233 * g(2);
t203 = qJDD(1) * pkin(1) + t212;
t213 = -t233 * g(1) - t231 * g(2);
t235 = qJD(1) ^ 2;
t205 = -t235 * pkin(1) + t213;
t228 = sin(pkin(7));
t229 = cos(pkin(7));
t181 = t228 * t203 + t229 * t205;
t163 = -t235 * pkin(2) + qJDD(1) * pkin(6) + t181;
t226 = -g(3) + qJDD(2);
t230 = sin(qJ(3));
t232 = cos(qJ(3));
t152 = -t230 * t163 + t232 * t226;
t254 = qJD(1) * qJD(3);
t252 = t232 * t254;
t206 = t230 * qJDD(1) + t252;
t147 = (-t206 + t252) * qJ(4) + (t230 * t232 * t235 + qJDD(3)) * pkin(3) + t152;
t153 = t232 * t163 + t230 * t226;
t207 = t232 * qJDD(1) - t230 * t254;
t256 = qJD(1) * t230;
t209 = qJD(3) * pkin(3) - qJ(4) * t256;
t225 = t232 ^ 2;
t148 = -t225 * t235 * pkin(3) + t207 * qJ(4) - qJD(3) * t209 + t153;
t227 = sin(pkin(8));
t255 = qJD(1) * t232;
t259 = cos(pkin(8));
t191 = t227 * t256 - t259 * t255;
t261 = -2 * qJD(4);
t144 = t227 * t147 + t259 * t148 + t191 * t261;
t182 = t227 * t206 - t259 * t207;
t192 = (t227 * t232 + t259 * t230) * qJD(1);
t187 = qJD(3) * mrSges(5,1) - t192 * mrSges(5,3);
t167 = t191 * pkin(4) - t192 * qJ(5);
t234 = qJD(3) ^ 2;
t137 = -t234 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t191 * t167 + t144;
t188 = -qJD(3) * mrSges(6,1) + t192 * mrSges(6,2);
t253 = m(6) * t137 + qJDD(3) * mrSges(6,3) + qJD(3) * t188;
t168 = t191 * mrSges(6,1) - t192 * mrSges(6,3);
t257 = -t191 * mrSges(5,1) - t192 * mrSges(5,2) - t168;
t260 = -mrSges(5,3) - mrSges(6,2);
t127 = m(5) * t144 - qJDD(3) * mrSges(5,2) - qJD(3) * t187 + t260 * t182 + t257 * t191 + t253;
t244 = t259 * t147 - t227 * t148;
t143 = t192 * t261 + t244;
t183 = t259 * t206 + t227 * t207;
t186 = -qJD(3) * mrSges(5,2) - t191 * mrSges(5,3);
t139 = -qJDD(3) * pkin(4) - t234 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t167) * t192 - t244;
t189 = -t191 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t248 = -m(6) * t139 + qJDD(3) * mrSges(6,1) + qJD(3) * t189;
t128 = m(5) * t143 + qJDD(3) * mrSges(5,1) + qJD(3) * t186 + t260 * t183 + t257 * t192 + t248;
t121 = t227 * t127 + t259 * t128;
t197 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t230 + Ifges(4,2) * t232) * qJD(1);
t198 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t230 + Ifges(4,4) * t232) * qJD(1);
t159 = Ifges(5,4) * t192 - Ifges(5,2) * t191 + Ifges(5,6) * qJD(3);
t161 = Ifges(5,1) * t192 - Ifges(5,4) * t191 + Ifges(5,5) * qJD(3);
t156 = Ifges(6,5) * t192 + Ifges(6,6) * qJD(3) + Ifges(6,3) * t191;
t160 = Ifges(6,1) * t192 + Ifges(6,4) * qJD(3) + Ifges(6,5) * t191;
t241 = mrSges(6,1) * t139 - mrSges(6,3) * t137 - Ifges(6,4) * t183 - Ifges(6,2) * qJDD(3) - Ifges(6,6) * t182 + t192 * t156 - t191 * t160;
t238 = mrSges(5,2) * t144 - t191 * t161 - qJ(5) * (-t182 * mrSges(6,2) - t191 * t168 + t253) - pkin(4) * (-t183 * mrSges(6,2) - t192 * t168 + t248) - mrSges(5,1) * t143 - t192 * t159 + Ifges(5,6) * t182 - Ifges(5,5) * t183 - Ifges(5,3) * qJDD(3) + t241;
t262 = mrSges(4,1) * t152 - mrSges(4,2) * t153 + Ifges(4,5) * t206 + Ifges(4,6) * t207 + Ifges(4,3) * qJDD(3) + pkin(3) * t121 + (t230 * t197 - t232 * t198) * qJD(1) - t238;
t204 = (-mrSges(4,1) * t232 + mrSges(4,2) * t230) * qJD(1);
t211 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t255;
t117 = m(4) * t152 + qJDD(3) * mrSges(4,1) - t206 * mrSges(4,3) + qJD(3) * t211 - t204 * t256 + t121;
t210 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t256;
t249 = t259 * t127 - t227 * t128;
t118 = m(4) * t153 - qJDD(3) * mrSges(4,2) + t207 * mrSges(4,3) - qJD(3) * t210 + t204 * t255 + t249;
t250 = -t230 * t117 + t232 * t118;
t111 = m(3) * t181 - t235 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t250;
t180 = t229 * t203 - t228 * t205;
t245 = -qJDD(1) * pkin(2) - t180;
t162 = -t235 * pkin(6) + t245;
t149 = -t207 * pkin(3) + qJDD(4) + t209 * t256 + (-qJ(4) * t225 - pkin(6)) * t235 + t245;
t141 = -0.2e1 * qJD(5) * t192 + (qJD(3) * t191 - t183) * qJ(5) + (qJD(3) * t192 + t182) * pkin(4) + t149;
t134 = m(6) * t141 + t182 * mrSges(6,1) - t183 * mrSges(6,3) - t192 * t188 + t191 * t189;
t239 = m(5) * t149 + t182 * mrSges(5,1) + t183 * mrSges(5,2) + t191 * t186 + t192 * t187 + t134;
t237 = -m(4) * t162 + t207 * mrSges(4,1) - t206 * mrSges(4,2) - t210 * t256 + t211 * t255 - t239;
t123 = m(3) * t180 + qJDD(1) * mrSges(3,1) - t235 * mrSges(3,2) + t237;
t108 = t228 * t111 + t229 * t123;
t113 = t232 * t117 + t230 * t118;
t158 = Ifges(6,4) * t192 + Ifges(6,2) * qJD(3) + Ifges(6,6) * t191;
t258 = -Ifges(5,5) * t192 + Ifges(5,6) * t191 - Ifges(5,3) * qJD(3) - t158;
t251 = t229 * t111 - t228 * t123;
t247 = -mrSges(6,1) * t141 + mrSges(6,2) * t137;
t243 = mrSges(6,2) * t139 - mrSges(6,3) * t141 + Ifges(6,1) * t183 + Ifges(6,4) * qJDD(3) + Ifges(6,5) * t182 + qJD(3) * t156;
t119 = -mrSges(5,1) * t149 + mrSges(5,3) * t144 - pkin(4) * t134 + t258 * t192 + (Ifges(5,4) - Ifges(6,5)) * t183 + (-Ifges(5,2) - Ifges(6,3)) * t182 + (Ifges(5,6) - Ifges(6,6)) * qJDD(3) + (t160 + t161) * qJD(3) + t247;
t120 = mrSges(5,2) * t149 - mrSges(5,3) * t143 + Ifges(5,1) * t183 - Ifges(5,4) * t182 + Ifges(5,5) * qJDD(3) - qJ(5) * t134 - qJD(3) * t159 + t258 * t191 + t243;
t196 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t230 + Ifges(4,6) * t232) * qJD(1);
t102 = -mrSges(4,1) * t162 + mrSges(4,3) * t153 + Ifges(4,4) * t206 + Ifges(4,2) * t207 + Ifges(4,6) * qJDD(3) - pkin(3) * t239 + qJ(4) * t249 + qJD(3) * t198 + t259 * t119 + t227 * t120 - t196 * t256;
t104 = mrSges(4,2) * t162 - mrSges(4,3) * t152 + Ifges(4,1) * t206 + Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3) - qJ(4) * t121 - qJD(3) * t197 - t227 * t119 + t259 * t120 + t196 * t255;
t242 = mrSges(3,1) * t180 - mrSges(3,2) * t181 + Ifges(3,3) * qJDD(1) + pkin(2) * t237 + pkin(6) * t250 + t232 * t102 + t230 * t104;
t240 = mrSges(2,1) * t212 - mrSges(2,2) * t213 + Ifges(2,3) * qJDD(1) + pkin(1) * t108 + t242;
t106 = m(2) * t213 - t235 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t251;
t105 = m(2) * t212 + qJDD(1) * mrSges(2,1) - t235 * mrSges(2,2) + t108;
t100 = -mrSges(3,1) * t226 + mrSges(3,3) * t181 + t235 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t113 - t262;
t99 = mrSges(3,2) * t226 - mrSges(3,3) * t180 + Ifges(3,5) * qJDD(1) - t235 * Ifges(3,6) - pkin(6) * t113 - t230 * t102 + t232 * t104;
t98 = -mrSges(2,2) * g(3) - mrSges(2,3) * t212 + Ifges(2,5) * qJDD(1) - t235 * Ifges(2,6) - qJ(2) * t108 - t228 * t100 + t229 * t99;
t97 = Ifges(2,6) * qJDD(1) + t235 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t213 + t228 * t99 + t229 * t100 - pkin(1) * (m(3) * t226 + t113) + qJ(2) * t251;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t233 * t98 - t231 * t97 - pkin(5) * (t233 * t105 + t231 * t106), t98, t99, t104, t120, -t191 * t158 + t243; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t231 * t98 + t233 * t97 + pkin(5) * (-t231 * t105 + t233 * t106), t97, t100, t102, t119, -t241; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t240, t240, t242, t262, -t238, Ifges(6,5) * t183 + Ifges(6,6) * qJDD(3) + Ifges(6,3) * t182 - qJD(3) * t160 + t192 * t158 - t247;];
m_new = t1;
