% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% m [6x1]
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
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:38
% EndTime: 2022-01-20 10:19:43
% DurationCPUTime: 2.28s
% Computational Cost: add. (32487->221), mult. (41797->272), div. (0->0), fcn. (20254->8), ass. (0->91)
t193 = qJD(1) + qJD(2);
t203 = sin(qJ(4));
t206 = cos(qJ(4));
t167 = (-mrSges(6,1) * t206 + mrSges(6,2) * t203) * t193;
t192 = qJDD(1) + qJDD(2);
t228 = qJD(4) * t193;
t223 = t206 * t228;
t169 = t203 * t192 + t223;
t205 = sin(qJ(1));
t208 = cos(qJ(1));
t183 = t205 * g(1) - t208 * g(2);
t175 = qJDD(1) * pkin(1) + t183;
t184 = -t208 * g(1) - t205 * g(2);
t209 = qJD(1) ^ 2;
t176 = -t209 * pkin(1) + t184;
t204 = sin(qJ(2));
t207 = cos(qJ(2));
t144 = t207 * t175 - t204 * t176;
t141 = t192 * pkin(2) + t144;
t145 = t204 * t175 + t207 * t176;
t191 = t193 ^ 2;
t142 = -t191 * pkin(2) + t145;
t201 = sin(pkin(8));
t202 = cos(pkin(8));
t137 = t201 * t141 + t202 * t142;
t134 = -t191 * pkin(3) + t192 * pkin(7) + t137;
t200 = -g(3) + qJDD(3);
t186 = t206 * t200;
t227 = qJD(5) * t193;
t235 = pkin(4) * t191;
t126 = qJDD(4) * pkin(4) + t186 + (-t169 + t223) * qJ(5) + (t206 * t235 - t134 - 0.2e1 * t227) * t203;
t232 = t193 * t206;
t180 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t232;
t226 = m(6) * t126 + qJDD(4) * mrSges(6,1) + qJD(4) * t180;
t233 = t193 * t203;
t121 = -t169 * mrSges(6,3) - t167 * t233 + t226;
t130 = -t203 * t134 + t186;
t131 = t206 * t134 + t203 * t200;
t153 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t203 + Ifges(5,2) * t206) * t193;
t154 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t203 + Ifges(6,4) * t206) * t193;
t155 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t203 + Ifges(5,4) * t206) * t193;
t170 = t206 * t192 - t203 * t228;
t177 = qJD(4) * pkin(4) - qJ(5) * t233;
t199 = t206 ^ 2;
t127 = t170 * qJ(5) - qJD(4) * t177 - t199 * t235 + 0.2e1 * t206 * t227 + t131;
t152 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t203 + Ifges(6,2) * t206) * t193;
t216 = -mrSges(6,1) * t126 + mrSges(6,2) * t127 - Ifges(6,5) * t169 - Ifges(6,6) * t170 - Ifges(6,3) * qJDD(4) - t152 * t233;
t237 = mrSges(5,1) * t130 - mrSges(5,2) * t131 + Ifges(5,5) * t169 + Ifges(5,6) * t170 + Ifges(5,3) * qJDD(4) + pkin(4) * t121 - (-t203 * t153 + (t154 + t155) * t206) * t193 - t216;
t234 = -mrSges(5,2) - mrSges(6,2);
t168 = (-mrSges(5,1) * t206 + mrSges(5,2) * t203) * t193;
t181 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t232;
t119 = m(5) * t130 + qJDD(4) * mrSges(5,1) + qJD(4) * t181 + (-t167 - t168) * t233 + (-mrSges(5,3) - mrSges(6,3)) * t169 + t226;
t225 = m(6) * t127 + t170 * mrSges(6,3) + t167 * t232;
t178 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t233;
t229 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t233 - t178;
t120 = m(5) * t131 + t170 * mrSges(5,3) + t229 * qJD(4) + t234 * qJDD(4) + t168 * t232 + t225;
t220 = -t203 * t119 + t206 * t120;
t110 = m(4) * t137 - t191 * mrSges(4,1) - t192 * mrSges(4,2) + t220;
t136 = t202 * t141 - t201 * t142;
t218 = -t192 * pkin(3) - t136;
t133 = -t191 * pkin(7) + t218;
t129 = t177 * t233 - t170 * pkin(4) + qJDD(5) + (-qJ(5) * t199 - pkin(7)) * t191 + t218;
t219 = -m(6) * t129 + t170 * mrSges(6,1) + t180 * t232;
t212 = -m(5) * t133 + t170 * mrSges(5,1) + t234 * t169 + t181 * t232 + t229 * t233 + t219;
t114 = m(4) * t136 + t192 * mrSges(4,1) - t191 * mrSges(4,2) + t212;
t103 = t201 * t110 + t202 * t114;
t100 = m(3) * t144 + t192 * mrSges(3,1) - t191 * mrSges(3,2) + t103;
t221 = t202 * t110 - t201 * t114;
t101 = m(3) * t145 - t191 * mrSges(3,1) - t192 * mrSges(3,2) + t221;
t95 = t207 * t100 + t204 * t101;
t112 = t206 * t119 + t203 * t120;
t224 = m(4) * t200 + t112;
t222 = -t204 * t100 + t207 * t101;
t217 = -mrSges(6,1) * t129 + mrSges(6,3) * t127 + Ifges(6,4) * t169 + Ifges(6,2) * t170 + Ifges(6,6) * qJDD(4) + qJD(4) * t154;
t150 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t203 + Ifges(6,6) * t206) * t193;
t215 = mrSges(6,2) * t129 - mrSges(6,3) * t126 + Ifges(6,1) * t169 + Ifges(6,4) * t170 + Ifges(6,5) * qJDD(4) + t150 * t232;
t151 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t203 + Ifges(5,6) * t206) * t193;
t105 = Ifges(5,4) * t169 + Ifges(5,2) * t170 + Ifges(5,6) * qJDD(4) + qJD(4) * t155 - mrSges(5,1) * t133 + mrSges(5,3) * t131 - pkin(4) * (t169 * mrSges(6,2) - t219) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t178 + t225) + (-pkin(4) * t178 - t150 - t151) * t233 + t217;
t107 = t151 * t232 + mrSges(5,2) * t133 - mrSges(5,3) * t130 + Ifges(5,1) * t169 + Ifges(5,4) * t170 + Ifges(5,5) * qJDD(4) - qJ(5) * t121 + (-t152 - t153) * qJD(4) + t215;
t214 = mrSges(4,1) * t136 - mrSges(4,2) * t137 + Ifges(4,3) * t192 + pkin(3) * t212 + pkin(7) * t220 + t206 * t105 + t203 * t107;
t213 = mrSges(3,1) * t144 - mrSges(3,2) * t145 + Ifges(3,3) * t192 + pkin(2) * t103 + t214;
t211 = mrSges(2,1) * t183 - mrSges(2,2) * t184 + Ifges(2,3) * qJDD(1) + pkin(1) * t95 + t213;
t96 = -mrSges(4,1) * t200 + mrSges(4,3) * t137 + t191 * Ifges(4,5) + Ifges(4,6) * t192 - pkin(3) * t112 - t237;
t93 = m(2) * t184 - t209 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t222;
t92 = m(2) * t183 + qJDD(1) * mrSges(2,1) - t209 * mrSges(2,2) + t95;
t91 = mrSges(4,2) * t200 - mrSges(4,3) * t136 + Ifges(4,5) * t192 - t191 * Ifges(4,6) - pkin(7) * t112 - t203 * t105 + t206 * t107;
t90 = -mrSges(3,2) * g(3) - mrSges(3,3) * t144 + Ifges(3,5) * t192 - t191 * Ifges(3,6) - qJ(3) * t103 - t201 * t96 + t202 * t91;
t89 = mrSges(3,1) * g(3) + mrSges(3,3) * t145 + t191 * Ifges(3,5) + Ifges(3,6) * t192 - pkin(2) * t224 + qJ(3) * t221 + t201 * t91 + t202 * t96;
t88 = -mrSges(2,2) * g(3) - mrSges(2,3) * t183 + Ifges(2,5) * qJDD(1) - t209 * Ifges(2,6) - pkin(6) * t95 - t204 * t89 + t207 * t90;
t87 = Ifges(2,6) * qJDD(1) + t209 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t184 + t204 * t90 + t207 * t89 - pkin(1) * (-m(3) * g(3) + t224) + pkin(6) * t222;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t208 * t88 - t205 * t87 - pkin(5) * (t205 * t93 + t208 * t92), t88, t90, t91, t107, -qJD(4) * t152 + t215; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t205 * t88 + t208 * t87 + pkin(5) * (-t205 * t92 + t208 * t93), t87, t89, t96, t105, -t150 * t233 + t217; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t211, t211, t213, t214, t237, -t154 * t232 - t216;];
m_new = t1;
