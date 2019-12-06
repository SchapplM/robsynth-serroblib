% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:45:41
% EndTime: 2019-12-05 16:45:45
% DurationCPUTime: 1.79s
% Computational Cost: add. (24305->212), mult. (30185->264), div. (0->0), fcn. (15841->8), ass. (0->84)
t191 = sin(pkin(8));
t219 = cos(pkin(8));
t177 = -t219 * g(1) - t191 * g(2);
t190 = -g(3) + qJDD(1);
t194 = sin(qJ(2));
t197 = cos(qJ(2));
t145 = -t194 * t177 + t197 * t190;
t141 = qJDD(2) * pkin(2) + t145;
t146 = t197 * t177 + t194 * t190;
t199 = qJD(2) ^ 2;
t142 = -t199 * pkin(2) + t146;
t193 = sin(qJ(3));
t196 = cos(qJ(3));
t136 = t193 * t141 + t196 * t142;
t185 = qJD(2) + qJD(3);
t183 = t185 ^ 2;
t184 = qJDD(2) + qJDD(3);
t133 = -t183 * pkin(3) + t184 * pkin(7) + t136;
t192 = sin(qJ(4));
t176 = t191 * g(1) - t219 * g(2);
t195 = cos(qJ(4));
t216 = t195 * t176;
t129 = -t192 * t133 - t216;
t130 = t195 * t133 - t192 * t176;
t147 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t192 - Ifges(6,3) * t195) * t185;
t150 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t192 + Ifges(5,2) * t195) * t185;
t163 = (-mrSges(6,1) * t195 - mrSges(6,3) * t192) * t185;
t213 = qJD(4) * t185;
t165 = t192 * t184 + t195 * t213;
t166 = t195 * t184 - t192 * t213;
t162 = (-pkin(4) * t195 - qJ(5) * t192) * t185;
t198 = qJD(4) ^ 2;
t217 = t185 * t195;
t126 = -t198 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t162 * t217 + t130;
t128 = -qJDD(4) * pkin(4) - t198 * qJ(5) + t216 + qJDD(5) + (t162 * t185 + t133) * t192;
t206 = -mrSges(6,1) * t128 + mrSges(6,3) * t126 + Ifges(6,4) * t165 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t166;
t174 = mrSges(6,2) * t217 + qJD(4) * mrSges(6,3);
t208 = -m(6) * t128 + qJDD(4) * mrSges(6,1) + qJD(4) * t174;
t218 = t185 * t192;
t172 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t218;
t209 = m(6) * t126 + qJDD(4) * mrSges(6,3) + qJD(4) * t172 + t163 * t217;
t151 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t192 - Ifges(6,5) * t195) * t185;
t214 = t151 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t192 + Ifges(5,4) * t195) * t185;
t223 = -(t214 * t195 + (t147 - t150) * t192) * t185 + mrSges(5,1) * t129 - mrSges(5,2) * t130 + Ifges(5,5) * t165 + Ifges(5,6) * t166 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t165 * mrSges(6,2) - t163 * t218 + t208) + qJ(5) * (t166 * mrSges(6,2) + t209) + t206;
t221 = m(3) + m(4);
t220 = mrSges(5,3) + mrSges(6,2);
t164 = (-mrSges(5,1) * t195 + mrSges(5,2) * t192) * t185;
t171 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t218;
t118 = m(5) * t130 - qJDD(4) * mrSges(5,2) - qJD(4) * t171 + t164 * t217 + t220 * t166 + t209;
t173 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t217;
t119 = m(5) * t129 + qJDD(4) * mrSges(5,1) + qJD(4) * t173 + (-t163 - t164) * t218 - t220 * t165 + t208;
t210 = t195 * t118 - t192 * t119;
t105 = m(4) * t136 - t183 * mrSges(4,1) - t184 * mrSges(4,2) + t210;
t135 = t196 * t141 - t193 * t142;
t132 = -t184 * pkin(3) - t183 * pkin(7) - t135;
t123 = -t166 * pkin(4) - t165 * qJ(5) + (-0.2e1 * qJD(5) * t192 + (pkin(4) * t192 - qJ(5) * t195) * qJD(4)) * t185 + t132;
t120 = m(6) * t123 - t166 * mrSges(6,1) - t165 * mrSges(6,3) - t172 * t218 - t174 * t217;
t202 = -m(5) * t132 + t166 * mrSges(5,1) - t165 * mrSges(5,2) - t171 * t218 + t173 * t217 - t120;
t112 = m(4) * t135 + t184 * mrSges(4,1) - t183 * mrSges(4,2) + t202;
t98 = t193 * t105 + t196 * t112;
t109 = t192 * t118 + t195 * t119;
t96 = m(3) * t145 + qJDD(2) * mrSges(3,1) - t199 * mrSges(3,2) + t98;
t211 = t196 * t105 - t193 * t112;
t97 = m(3) * t146 - t199 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t211;
t212 = -t194 * t96 + t197 * t97;
t207 = -mrSges(6,1) * t123 + mrSges(6,2) * t126;
t148 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t192 + Ifges(5,6) * t195) * t185;
t149 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t192 - Ifges(6,6) * t195) * t185;
t101 = -mrSges(5,1) * t132 + mrSges(5,3) * t130 - pkin(4) * t120 + (-t148 - t149) * t218 + (Ifges(5,2) + Ifges(6,3)) * t166 + (Ifges(5,4) - Ifges(6,5)) * t165 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t214 * qJD(4) + t207;
t203 = mrSges(6,2) * t128 - mrSges(6,3) * t123 + Ifges(6,1) * t165 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t166 + qJD(4) * t147 + t149 * t217;
t102 = mrSges(5,2) * t132 - mrSges(5,3) * t129 + Ifges(5,1) * t165 + Ifges(5,4) * t166 + Ifges(5,5) * qJDD(4) - qJ(5) * t120 - qJD(4) * t150 + t148 * t217 + t203;
t93 = -mrSges(4,2) * t176 - mrSges(4,3) * t135 + Ifges(4,5) * t184 - t183 * Ifges(4,6) - pkin(7) * t109 - t192 * t101 + t195 * t102;
t94 = mrSges(4,1) * t176 + mrSges(4,3) * t136 + t183 * Ifges(4,5) + Ifges(4,6) * t184 - pkin(3) * t109 - t223;
t87 = Ifges(3,6) * qJDD(2) + t199 * Ifges(3,5) + mrSges(3,1) * t176 + mrSges(3,3) * t146 + t193 * t93 + t196 * t94 - pkin(2) * (-m(4) * t176 + t109) + pkin(6) * t211;
t89 = -mrSges(3,2) * t176 - mrSges(3,3) * t145 + Ifges(3,5) * qJDD(2) - t199 * Ifges(3,6) - pkin(6) * t98 - t193 * t94 + t196 * t93;
t205 = -mrSges(2,2) * t177 + pkin(5) * t212 + t194 * t89 + t197 * t87 + pkin(1) * (t221 * t176 - t109) + mrSges(2,1) * t176;
t204 = mrSges(4,1) * t135 - mrSges(4,2) * t136 + Ifges(4,3) * t184 + pkin(3) * t202 + pkin(7) * t210 + t195 * t101 + t192 * t102;
t201 = mrSges(3,1) * t145 - mrSges(3,2) * t146 + Ifges(3,3) * qJDD(2) + pkin(2) * t98 + t204;
t106 = (m(2) + t221) * t176 - t109;
t92 = t194 * t97 + t197 * t96;
t90 = m(2) * t177 + t212;
t85 = -mrSges(2,1) * t190 + mrSges(2,3) * t177 - pkin(1) * t92 - t201;
t84 = mrSges(2,2) * t190 - mrSges(2,3) * t176 - pkin(5) * t92 - t194 * t87 + t197 * t89;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t219 * t84 - t191 * t85 - qJ(1) * (t219 * t106 + t191 * t90), t84, t89, t93, t102, t203; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t191 * t84 + t219 * t85 + qJ(1) * (-t191 * t106 + t219 * t90), t85, t87, t94, t101, (-t192 * t147 - t195 * t151) * t185 + t206; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t205, t205, t201, t204, t223, Ifges(6,5) * t165 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t166 - qJD(4) * t151 + t149 * t218 - t207;];
m_new = t1;
