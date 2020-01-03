% Calculate vector of cutting torques with Newton-Euler for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:49
% EndTime: 2019-12-31 17:36:51
% DurationCPUTime: 2.00s
% Computational Cost: add. (18644->213), mult. (40470->273), div. (0->0), fcn. (24553->8), ass. (0->99)
t213 = sin(pkin(7));
t215 = cos(pkin(7));
t192 = t213 * g(1) - t215 * g(2);
t193 = -t215 * g(1) - t213 * g(2);
t217 = sin(qJ(2));
t219 = cos(qJ(2));
t169 = t217 * t192 + t219 * t193;
t220 = qJD(2) ^ 2;
t260 = -t220 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t169;
t212 = sin(pkin(8));
t207 = t212 ^ 2;
t214 = cos(pkin(8));
t208 = t214 ^ 2;
t248 = t208 * t220;
t259 = t207 * t220 + t248;
t211 = -g(3) + qJDD(1);
t149 = t214 * t211 - t260 * t212;
t250 = qJ(4) * t212;
t233 = -pkin(3) * t214 - t250;
t182 = t233 * qJD(2);
t244 = qJD(2) * t212;
t143 = t182 * t244 + qJDD(4) - t149;
t137 = (-pkin(4) * t214 * t220 - pkin(6) * qJDD(2)) * t212 + t143;
t150 = t212 * t211 + t260 * t214;
t243 = qJD(2) * t214;
t145 = t182 * t243 + t150;
t240 = qJDD(2) * t214;
t138 = -pkin(4) * t248 - pkin(6) * t240 + t145;
t216 = sin(qJ(5));
t218 = cos(qJ(5));
t135 = t218 * t137 - t216 * t138;
t231 = -t212 * t216 - t214 * t218;
t176 = t231 * qJD(2);
t232 = t212 * t218 - t214 * t216;
t177 = t232 * qJD(2);
t158 = -t176 * mrSges(6,1) + t177 * mrSges(6,2);
t165 = t176 * qJD(5) + t232 * qJDD(2);
t170 = -qJD(5) * mrSges(6,2) + t176 * mrSges(6,3);
t130 = m(6) * t135 + qJDD(5) * mrSges(6,1) - t165 * mrSges(6,3) + qJD(5) * t170 - t177 * t158;
t136 = t216 * t137 + t218 * t138;
t164 = -t177 * qJD(5) + t231 * qJDD(2);
t171 = qJD(5) * mrSges(6,1) - t177 * mrSges(6,3);
t131 = m(6) * t136 - qJDD(5) * mrSges(6,2) + t164 * mrSges(6,3) - qJD(5) * t171 + t176 * t158;
t121 = -t216 * t130 + t218 * t131;
t183 = (-mrSges(5,1) * t214 - mrSges(5,3) * t212) * qJD(2);
t119 = m(5) * t145 + mrSges(5,2) * t240 + t183 * t243 + t121;
t234 = Ifges(5,5) * t212 - Ifges(5,3) * t214;
t185 = t234 * qJD(2);
t189 = (Ifges(5,1) * t212 - Ifges(5,5) * t214) * qJD(2);
t120 = t218 * t130 + t216 * t131;
t152 = Ifges(6,4) * t177 + Ifges(6,2) * t176 + Ifges(6,6) * qJD(5);
t153 = Ifges(6,1) * t177 + Ifges(6,4) * t176 + Ifges(6,5) * qJD(5);
t230 = mrSges(6,1) * t135 - mrSges(6,2) * t136 + Ifges(6,5) * t165 + Ifges(6,6) * t164 + Ifges(6,3) * qJDD(5) + t177 * t152 - t176 * t153;
t241 = qJDD(2) * t212;
t224 = -mrSges(5,1) * t143 + mrSges(5,3) * t145 + Ifges(5,4) * t241 - pkin(4) * t120 - t230;
t228 = -m(5) * t143 - t120;
t252 = Ifges(4,1) * t212;
t258 = ((t185 - (Ifges(4,4) * t212 + Ifges(4,2) * t214) * qJD(2)) * t212 + (t189 + (Ifges(4,4) * t214 + t252) * qJD(2)) * t214) * qJD(2) - mrSges(4,1) * t149 + mrSges(4,2) * t150 - pkin(3) * ((-qJDD(2) * mrSges(5,2) - qJD(2) * t183) * t212 + t228) - qJ(4) * t119 - t224;
t251 = Ifges(4,5) * t212;
t257 = (Ifges(4,6) - Ifges(5,6)) * t214 + t251;
t255 = Ifges(4,4) - Ifges(5,5);
t253 = mrSges(4,2) * t212;
t247 = t220 * qJ(3);
t184 = (-mrSges(4,1) * t214 + t253) * qJD(2);
t116 = m(4) * t149 + ((-mrSges(5,2) - mrSges(4,3)) * qJDD(2) + (-t183 - t184) * qJD(2)) * t212 + t228;
t117 = m(4) * t150 + (qJDD(2) * mrSges(4,3) + qJD(2) * t184) * t214 + t119;
t236 = -t212 * t116 + t214 * t117;
t110 = m(3) * t169 - t220 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t236;
t168 = t219 * t192 - t217 * t193;
t238 = -qJDD(3) + t168;
t229 = -0.2e1 * qJD(4) * t244 - t238;
t141 = (qJ(3) + (-t207 - t208) * pkin(6)) * t220 + (t250 + pkin(2) + (pkin(3) + pkin(4)) * t214) * qJDD(2) - t229;
t132 = -m(6) * t141 + t164 * mrSges(6,1) - t165 * mrSges(6,2) + t176 * t170 - t177 * t171;
t148 = -t247 + (-pkin(2) + t233) * qJDD(2) + t229;
t128 = m(5) * t148 - mrSges(5,1) * t240 - t259 * mrSges(5,2) - mrSges(5,3) * t241 + t132;
t163 = -qJDD(2) * pkin(2) - t238 - t247;
t222 = -m(4) * t163 + mrSges(4,1) * t240 + t259 * mrSges(4,3) - t128;
t126 = -t220 * mrSges(3,2) + m(3) * t168 + t222 + (mrSges(3,1) - t253) * qJDD(2);
t107 = t217 * t110 + t219 * t126;
t112 = t214 * t116 + t212 * t117;
t237 = t219 * t110 - t217 * t126;
t186 = (Ifges(4,6) * t214 + t251) * qJD(2);
t187 = (Ifges(5,4) * t212 - Ifges(5,6) * t214) * qJD(2);
t151 = Ifges(6,5) * t177 + Ifges(6,6) * t176 + Ifges(6,3) * qJD(5);
t123 = -mrSges(6,1) * t141 + mrSges(6,3) * t136 + Ifges(6,4) * t165 + Ifges(6,2) * t164 + Ifges(6,6) * qJDD(5) + qJD(5) * t153 - t177 * t151;
t124 = mrSges(6,2) * t141 - mrSges(6,3) * t135 + Ifges(6,1) * t165 + Ifges(6,4) * t164 + Ifges(6,5) * qJDD(5) - qJD(5) * t152 + t176 * t151;
t223 = mrSges(5,1) * t148 - mrSges(5,2) * t145 + pkin(4) * t132 + pkin(6) * t121 + t218 * t123 + t216 * t124;
t101 = -mrSges(4,1) * t163 + mrSges(4,3) * t150 - pkin(3) * t128 + (-t186 - t187) * t244 + ((Ifges(4,2) + Ifges(5,3)) * t214 + t255 * t212) * qJDD(2) - t223;
t226 = mrSges(5,2) * t143 - mrSges(5,3) * t148 + Ifges(5,1) * t241 - pkin(6) * t120 - t216 * t123 + t218 * t124 + t187 * t243;
t103 = t186 * t243 + mrSges(4,2) * t163 - mrSges(4,3) * t149 - qJ(4) * t128 + (t255 * t214 + t252) * qJDD(2) + t226;
t227 = -mrSges(3,2) * t169 + qJ(3) * t236 + t214 * t101 + t212 * t103 + pkin(2) * (-mrSges(4,2) * t241 + t222) + mrSges(3,1) * t168 + Ifges(3,3) * qJDD(2);
t225 = mrSges(2,1) * t192 - mrSges(2,2) * t193 + pkin(1) * t107 + t227;
t105 = m(2) * t193 + t237;
t104 = m(2) * t192 + t107;
t99 = t220 * Ifges(3,5) - mrSges(3,1) * t211 + mrSges(3,3) * t169 - pkin(2) * t112 + (Ifges(3,6) - t257) * qJDD(2) + t258;
t98 = mrSges(3,2) * t211 - mrSges(3,3) * t168 + Ifges(3,5) * qJDD(2) - t220 * Ifges(3,6) - qJ(3) * t112 - t212 * t101 + t214 * t103;
t97 = mrSges(2,2) * t211 - mrSges(2,3) * t192 - pkin(5) * t107 - t217 * t99 + t219 * t98;
t96 = -mrSges(2,1) * t211 + mrSges(2,3) * t193 + t217 * t98 + t219 * t99 - pkin(1) * (m(3) * t211 + t112) + pkin(5) * t237;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t215 * t97 - t213 * t96 - qJ(1) * (t215 * t104 + t213 * t105), t97, t98, t103, -Ifges(5,5) * t240 + t226, t124; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t213 * t97 + t215 * t96 + qJ(1) * (-t213 * t104 + t215 * t105), t96, t99, t101, t224 + (-t212 * t185 - t214 * t189) * qJD(2) - Ifges(5,6) * t240, t123; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t225, t225, t227, t257 * qJDD(2) - t258, t234 * qJDD(2) + t187 * t244 + t223, t230;];
m_new = t1;
