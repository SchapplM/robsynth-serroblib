% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:50:45
% EndTime: 2019-12-31 21:50:50
% DurationCPUTime: 3.50s
% Computational Cost: add. (60028->266), mult. (75710->324), div. (0->0), fcn. (44282->8), ass. (0->102)
t224 = qJDD(1) + qJDD(2);
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t226 = qJD(1) + qJD(2);
t253 = qJD(3) * t226;
t200 = t230 * t224 + t233 * t253;
t232 = sin(qJ(1));
t235 = cos(qJ(1));
t212 = t232 * g(1) - t235 * g(2);
t205 = qJDD(1) * pkin(1) + t212;
t213 = -t235 * g(1) - t232 * g(2);
t236 = qJD(1) ^ 2;
t206 = -t236 * pkin(1) + t213;
t231 = sin(qJ(2));
t234 = cos(qJ(2));
t183 = t231 * t205 + t234 * t206;
t222 = t226 ^ 2;
t180 = -t222 * pkin(2) + t224 * pkin(7) + t183;
t256 = t230 * t180;
t260 = pkin(3) * t222;
t147 = qJDD(3) * pkin(3) - t200 * pkin(8) - t256 + (pkin(8) * t253 + t230 * t260 - g(3)) * t233;
t165 = -t230 * g(3) + t233 * t180;
t201 = t233 * t224 - t230 * t253;
t258 = t226 * t230;
t209 = qJD(3) * pkin(3) - pkin(8) * t258;
t228 = t233 ^ 2;
t148 = t201 * pkin(8) - qJD(3) * t209 - t228 * t260 + t165;
t229 = sin(qJ(4));
t261 = cos(qJ(4));
t144 = t229 * t147 + t261 * t148;
t194 = (t229 * t233 + t261 * t230) * t226;
t161 = t194 * qJD(4) + t229 * t200 - t261 * t201;
t225 = qJD(3) + qJD(4);
t187 = t225 * mrSges(5,1) - t194 * mrSges(5,3);
t257 = t226 * t233;
t193 = t229 * t258 - t261 * t257;
t223 = qJDD(3) + qJDD(4);
t176 = t193 * pkin(4) - t194 * qJ(5);
t221 = t225 ^ 2;
t139 = -t221 * pkin(4) + t223 * qJ(5) + 0.2e1 * qJD(5) * t225 - t193 * t176 + t144;
t188 = -t225 * mrSges(6,1) + t194 * mrSges(6,2);
t252 = m(6) * t139 + t223 * mrSges(6,3) + t225 * t188;
t177 = t193 * mrSges(6,1) - t194 * mrSges(6,3);
t254 = -t193 * mrSges(5,1) - t194 * mrSges(5,2) - t177;
t259 = -mrSges(5,3) - mrSges(6,2);
t127 = m(5) * t144 - t223 * mrSges(5,2) + t259 * t161 - t225 * t187 + t254 * t193 + t252;
t143 = t261 * t147 - t229 * t148;
t162 = -t193 * qJD(4) + t261 * t200 + t229 * t201;
t186 = -t225 * mrSges(5,2) - t193 * mrSges(5,3);
t141 = -t223 * pkin(4) - t221 * qJ(5) + t194 * t176 + qJDD(5) - t143;
t189 = -t193 * mrSges(6,2) + t225 * mrSges(6,3);
t248 = -m(6) * t141 + t223 * mrSges(6,1) + t225 * t189;
t129 = m(5) * t143 + t223 * mrSges(5,1) + t259 * t162 + t225 * t186 + t254 * t194 + t248;
t121 = t229 * t127 + t261 * t129;
t164 = -t233 * g(3) - t256;
t191 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t230 + Ifges(4,2) * t233) * t226;
t192 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t230 + Ifges(4,4) * t233) * t226;
t169 = Ifges(5,4) * t194 - Ifges(5,2) * t193 + Ifges(5,6) * t225;
t171 = Ifges(5,1) * t194 - Ifges(5,4) * t193 + Ifges(5,5) * t225;
t166 = Ifges(6,5) * t194 + Ifges(6,6) * t225 + Ifges(6,3) * t193;
t170 = Ifges(6,1) * t194 + Ifges(6,4) * t225 + Ifges(6,5) * t193;
t242 = mrSges(6,1) * t141 - mrSges(6,3) * t139 - Ifges(6,4) * t162 - Ifges(6,2) * t223 - Ifges(6,6) * t161 + t194 * t166 - t193 * t170;
t239 = mrSges(5,2) * t144 - t193 * t171 - qJ(5) * (-t161 * mrSges(6,2) - t193 * t177 + t252) - pkin(4) * (-t162 * mrSges(6,2) - t194 * t177 + t248) - mrSges(5,1) * t143 - t194 * t169 + Ifges(5,6) * t161 - Ifges(5,5) * t162 - Ifges(5,3) * t223 + t242;
t262 = mrSges(4,1) * t164 - mrSges(4,2) * t165 + Ifges(4,5) * t200 + Ifges(4,6) * t201 + Ifges(4,3) * qJDD(3) + pkin(3) * t121 + (t230 * t191 - t233 * t192) * t226 - t239;
t199 = (-mrSges(4,1) * t233 + mrSges(4,2) * t230) * t226;
t208 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t257;
t119 = m(4) * t164 + qJDD(3) * mrSges(4,1) - t200 * mrSges(4,3) + qJD(3) * t208 - t199 * t258 + t121;
t207 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t258;
t249 = t261 * t127 - t229 * t129;
t120 = m(4) * t165 - qJDD(3) * mrSges(4,2) + t201 * mrSges(4,3) - qJD(3) * t207 + t199 * t257 + t249;
t250 = -t230 * t119 + t233 * t120;
t111 = m(3) * t183 - t222 * mrSges(3,1) - t224 * mrSges(3,2) + t250;
t182 = t234 * t205 - t231 * t206;
t245 = -t224 * pkin(2) - t182;
t179 = -t222 * pkin(7) + t245;
t149 = -t201 * pkin(3) + t209 * t258 + (-pkin(8) * t228 - pkin(7)) * t222 + t245;
t136 = -0.2e1 * qJD(5) * t194 + (t193 * t225 - t162) * qJ(5) + (t194 * t225 + t161) * pkin(4) + t149;
t130 = m(6) * t136 + t161 * mrSges(6,1) - t162 * mrSges(6,3) - t194 * t188 + t193 * t189;
t240 = m(5) * t149 + t161 * mrSges(5,1) + t162 * mrSges(5,2) + t193 * t186 + t194 * t187 + t130;
t238 = -m(4) * t179 + t201 * mrSges(4,1) - t200 * mrSges(4,2) - t207 * t258 + t208 * t257 - t240;
t123 = m(3) * t182 + t224 * mrSges(3,1) - t222 * mrSges(3,2) + t238;
t108 = t231 * t111 + t234 * t123;
t113 = t233 * t119 + t230 * t120;
t168 = Ifges(6,4) * t194 + Ifges(6,2) * t225 + Ifges(6,6) * t193;
t255 = -Ifges(5,5) * t194 + Ifges(5,6) * t193 - Ifges(5,3) * t225 - t168;
t251 = t234 * t111 - t231 * t123;
t247 = -mrSges(6,1) * t136 + mrSges(6,2) * t139;
t244 = mrSges(6,2) * t141 - mrSges(6,3) * t136 + Ifges(6,1) * t162 + Ifges(6,4) * t223 + Ifges(6,5) * t161 + t225 * t166;
t114 = -mrSges(5,1) * t149 + mrSges(5,3) * t144 - pkin(4) * t130 + (t170 + t171) * t225 + (Ifges(5,6) - Ifges(6,6)) * t223 + t255 * t194 + (Ifges(5,4) - Ifges(6,5)) * t162 + (-Ifges(5,2) - Ifges(6,3)) * t161 + t247;
t115 = mrSges(5,2) * t149 - mrSges(5,3) * t143 + Ifges(5,1) * t162 - Ifges(5,4) * t161 + Ifges(5,5) * t223 - qJ(5) * t130 - t225 * t169 + t255 * t193 + t244;
t190 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t230 + Ifges(4,6) * t233) * t226;
t102 = -mrSges(4,1) * t179 + mrSges(4,3) * t165 + Ifges(4,4) * t200 + Ifges(4,2) * t201 + Ifges(4,6) * qJDD(3) - pkin(3) * t240 + pkin(8) * t249 + qJD(3) * t192 + t261 * t114 + t229 * t115 - t190 * t258;
t104 = mrSges(4,2) * t179 - mrSges(4,3) * t164 + Ifges(4,1) * t200 + Ifges(4,4) * t201 + Ifges(4,5) * qJDD(3) - pkin(8) * t121 - qJD(3) * t191 - t229 * t114 + t261 * t115 + t190 * t257;
t243 = mrSges(3,1) * t182 - mrSges(3,2) * t183 + Ifges(3,3) * t224 + pkin(2) * t238 + pkin(7) * t250 + t233 * t102 + t230 * t104;
t241 = mrSges(2,1) * t212 - mrSges(2,2) * t213 + Ifges(2,3) * qJDD(1) + pkin(1) * t108 + t243;
t106 = m(2) * t213 - t236 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t251;
t105 = m(2) * t212 + qJDD(1) * mrSges(2,1) - t236 * mrSges(2,2) + t108;
t100 = mrSges(3,1) * g(3) + mrSges(3,3) * t183 + t222 * Ifges(3,5) + Ifges(3,6) * t224 - pkin(2) * t113 - t262;
t99 = -mrSges(3,2) * g(3) - mrSges(3,3) * t182 + Ifges(3,5) * t224 - t222 * Ifges(3,6) - pkin(7) * t113 - t230 * t102 + t233 * t104;
t98 = -mrSges(2,2) * g(3) - mrSges(2,3) * t212 + Ifges(2,5) * qJDD(1) - t236 * Ifges(2,6) - pkin(6) * t108 - t231 * t100 + t234 * t99;
t97 = Ifges(2,6) * qJDD(1) + t236 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t213 + t231 * t99 + t234 * t100 - pkin(1) * (-m(3) * g(3) + t113) + pkin(6) * t251;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t235 * t98 - t232 * t97 - pkin(5) * (t235 * t105 + t232 * t106), t98, t99, t104, t115, -t193 * t168 + t244; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t232 * t98 + t235 * t97 + pkin(5) * (-t232 * t105 + t235 * t106), t97, t100, t102, t114, -t242; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t241, t241, t243, t262, -t239, Ifges(6,5) * t162 + Ifges(6,6) * t223 + Ifges(6,3) * t161 + t194 * t168 - t225 * t170 - t247;];
m_new = t1;
