% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:28
% EndTime: 2019-12-31 20:49:34
% DurationCPUTime: 3.38s
% Computational Cost: add. (55364->266), mult. (73586->326), div. (0->0), fcn. (42162->8), ass. (0->102)
t218 = qJDD(1) + qJDD(2);
t228 = sin(qJ(3));
t231 = cos(qJ(3));
t219 = qJD(1) + qJD(2);
t253 = qJD(3) * t219;
t201 = t228 * t218 + t231 * t253;
t230 = sin(qJ(1));
t233 = cos(qJ(1));
t213 = t230 * g(1) - t233 * g(2);
t206 = qJDD(1) * pkin(1) + t213;
t214 = -t233 * g(1) - t230 * g(2);
t235 = qJD(1) ^ 2;
t207 = -t235 * pkin(1) + t214;
t229 = sin(qJ(2));
t232 = cos(qJ(2));
t183 = t229 * t206 + t232 * t207;
t217 = t219 ^ 2;
t169 = -t217 * pkin(2) + t218 * pkin(7) + t183;
t256 = t228 * t169;
t261 = pkin(3) * t217;
t147 = qJDD(3) * pkin(3) - t201 * qJ(4) - t256 + (qJ(4) * t253 + t228 * t261 - g(3)) * t231;
t153 = -t228 * g(3) + t231 * t169;
t202 = t231 * t218 - t228 * t253;
t258 = t219 * t228;
t208 = qJD(3) * pkin(3) - qJ(4) * t258;
t226 = t231 ^ 2;
t148 = t202 * qJ(4) - qJD(3) * t208 - t226 * t261 + t153;
t227 = sin(pkin(8));
t257 = t219 * t231;
t259 = cos(pkin(8));
t191 = t227 * t258 - t259 * t257;
t262 = -2 * qJD(4);
t144 = t227 * t147 + t259 * t148 + t191 * t262;
t179 = t227 * t201 - t259 * t202;
t192 = (t227 * t231 + t259 * t228) * t219;
t187 = qJD(3) * mrSges(5,1) - t192 * mrSges(5,3);
t165 = t191 * pkin(4) - t192 * qJ(5);
t234 = qJD(3) ^ 2;
t137 = -t234 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t191 * t165 + t144;
t188 = -qJD(3) * mrSges(6,1) + t192 * mrSges(6,2);
t252 = m(6) * t137 + qJDD(3) * mrSges(6,3) + qJD(3) * t188;
t166 = t191 * mrSges(6,1) - t192 * mrSges(6,3);
t254 = -t191 * mrSges(5,1) - t192 * mrSges(5,2) - t166;
t260 = -mrSges(5,3) - mrSges(6,2);
t128 = m(5) * t144 - qJDD(3) * mrSges(5,2) - qJD(3) * t187 + t260 * t179 + t254 * t191 + t252;
t244 = t259 * t147 - t227 * t148;
t143 = t192 * t262 + t244;
t180 = t259 * t201 + t227 * t202;
t186 = -qJD(3) * mrSges(5,2) - t191 * mrSges(5,3);
t139 = -qJDD(3) * pkin(4) - t234 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t165) * t192 - t244;
t189 = -t191 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t248 = -m(6) * t139 + qJDD(3) * mrSges(6,1) + qJD(3) * t189;
t129 = m(5) * t143 + qJDD(3) * mrSges(5,1) + qJD(3) * t186 + t260 * t180 + t254 * t192 + t248;
t121 = t227 * t128 + t259 * t129;
t152 = -t231 * g(3) - t256;
t195 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t228 + Ifges(4,2) * t231) * t219;
t196 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t228 + Ifges(4,4) * t231) * t219;
t158 = Ifges(5,4) * t192 - Ifges(5,2) * t191 + Ifges(5,6) * qJD(3);
t160 = Ifges(5,1) * t192 - Ifges(5,4) * t191 + Ifges(5,5) * qJD(3);
t155 = Ifges(6,5) * t192 + Ifges(6,6) * qJD(3) + Ifges(6,3) * t191;
t159 = Ifges(6,1) * t192 + Ifges(6,4) * qJD(3) + Ifges(6,5) * t191;
t241 = mrSges(6,1) * t139 - mrSges(6,3) * t137 - Ifges(6,4) * t180 - Ifges(6,2) * qJDD(3) - Ifges(6,6) * t179 + t192 * t155 - t191 * t159;
t238 = mrSges(5,2) * t144 - t191 * t160 - qJ(5) * (-t179 * mrSges(6,2) - t191 * t166 + t252) - pkin(4) * (-t180 * mrSges(6,2) - t192 * t166 + t248) - mrSges(5,1) * t143 - t192 * t158 + Ifges(5,6) * t179 - Ifges(5,5) * t180 - Ifges(5,3) * qJDD(3) + t241;
t263 = mrSges(4,1) * t152 - mrSges(4,2) * t153 + Ifges(4,5) * t201 + Ifges(4,6) * t202 + Ifges(4,3) * qJDD(3) + pkin(3) * t121 + (t228 * t195 - t231 * t196) * t219 - t238;
t200 = (-mrSges(4,1) * t231 + mrSges(4,2) * t228) * t219;
t210 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t257;
t117 = m(4) * t152 + qJDD(3) * mrSges(4,1) - t201 * mrSges(4,3) + qJD(3) * t210 - t200 * t258 + t121;
t209 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t258;
t249 = t259 * t128 - t227 * t129;
t118 = m(4) * t153 - qJDD(3) * mrSges(4,2) + t202 * mrSges(4,3) - qJD(3) * t209 + t200 * t257 + t249;
t250 = -t228 * t117 + t231 * t118;
t111 = m(3) * t183 - t217 * mrSges(3,1) - t218 * mrSges(3,2) + t250;
t182 = t232 * t206 - t229 * t207;
t245 = -t218 * pkin(2) - t182;
t168 = -t217 * pkin(7) + t245;
t149 = -t202 * pkin(3) + qJDD(4) + t208 * t258 + (-qJ(4) * t226 - pkin(7)) * t217 + t245;
t141 = -0.2e1 * qJD(5) * t192 + (qJD(3) * t191 - t180) * qJ(5) + (qJD(3) * t192 + t179) * pkin(4) + t149;
t134 = m(6) * t141 + t179 * mrSges(6,1) - t180 * mrSges(6,3) - t192 * t188 + t191 * t189;
t239 = m(5) * t149 + t179 * mrSges(5,1) + t180 * mrSges(5,2) + t191 * t186 + t192 * t187 + t134;
t237 = -m(4) * t168 + t202 * mrSges(4,1) - t201 * mrSges(4,2) - t209 * t258 + t210 * t257 - t239;
t123 = m(3) * t182 + t218 * mrSges(3,1) - t217 * mrSges(3,2) + t237;
t108 = t229 * t111 + t232 * t123;
t113 = t231 * t117 + t228 * t118;
t157 = Ifges(6,4) * t192 + Ifges(6,2) * qJD(3) + Ifges(6,6) * t191;
t255 = -Ifges(5,5) * t192 + Ifges(5,6) * t191 - Ifges(5,3) * qJD(3) - t157;
t251 = t232 * t111 - t229 * t123;
t247 = -mrSges(6,1) * t141 + mrSges(6,2) * t137;
t243 = mrSges(6,2) * t139 - mrSges(6,3) * t141 + Ifges(6,1) * t180 + Ifges(6,4) * qJDD(3) + Ifges(6,5) * t179 + qJD(3) * t155;
t119 = -mrSges(5,1) * t149 + mrSges(5,3) * t144 - pkin(4) * t134 + t255 * t192 + (Ifges(5,4) - Ifges(6,5)) * t180 + (-Ifges(5,2) - Ifges(6,3)) * t179 + (Ifges(5,6) - Ifges(6,6)) * qJDD(3) + (t159 + t160) * qJD(3) + t247;
t120 = mrSges(5,2) * t149 - mrSges(5,3) * t143 + Ifges(5,1) * t180 - Ifges(5,4) * t179 + Ifges(5,5) * qJDD(3) - qJ(5) * t134 - qJD(3) * t158 + t255 * t191 + t243;
t194 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t228 + Ifges(4,6) * t231) * t219;
t102 = -mrSges(4,1) * t168 + mrSges(4,3) * t153 + Ifges(4,4) * t201 + Ifges(4,2) * t202 + Ifges(4,6) * qJDD(3) - pkin(3) * t239 + qJ(4) * t249 + qJD(3) * t196 + t259 * t119 + t227 * t120 - t194 * t258;
t104 = mrSges(4,2) * t168 - mrSges(4,3) * t152 + Ifges(4,1) * t201 + Ifges(4,4) * t202 + Ifges(4,5) * qJDD(3) - qJ(4) * t121 - qJD(3) * t195 - t227 * t119 + t259 * t120 + t194 * t257;
t242 = mrSges(3,1) * t182 - mrSges(3,2) * t183 + Ifges(3,3) * t218 + pkin(2) * t237 + pkin(7) * t250 + t231 * t102 + t228 * t104;
t240 = mrSges(2,1) * t213 - mrSges(2,2) * t214 + Ifges(2,3) * qJDD(1) + pkin(1) * t108 + t242;
t106 = m(2) * t214 - t235 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t251;
t105 = m(2) * t213 + qJDD(1) * mrSges(2,1) - t235 * mrSges(2,2) + t108;
t100 = mrSges(3,1) * g(3) + mrSges(3,3) * t183 + t217 * Ifges(3,5) + Ifges(3,6) * t218 - pkin(2) * t113 - t263;
t99 = -mrSges(3,2) * g(3) - mrSges(3,3) * t182 + Ifges(3,5) * t218 - t217 * Ifges(3,6) - pkin(7) * t113 - t228 * t102 + t231 * t104;
t98 = -mrSges(2,2) * g(3) - mrSges(2,3) * t213 + Ifges(2,5) * qJDD(1) - t235 * Ifges(2,6) - pkin(6) * t108 - t229 * t100 + t232 * t99;
t97 = Ifges(2,6) * qJDD(1) + t235 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t214 + t229 * t99 + t232 * t100 - pkin(1) * (-m(3) * g(3) + t113) + pkin(6) * t251;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t233 * t98 - t230 * t97 - pkin(5) * (t233 * t105 + t230 * t106), t98, t99, t104, t120, -t191 * t157 + t243; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t230 * t98 + t233 * t97 + pkin(5) * (-t230 * t105 + t233 * t106), t97, t100, t102, t119, -t241; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t240, t240, t242, t263, -t238, Ifges(6,5) * t180 + Ifges(6,6) * qJDD(3) + Ifges(6,3) * t179 - qJD(3) * t159 + t192 * t157 - t247;];
m_new = t1;
