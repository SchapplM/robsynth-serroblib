% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRP1
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:42
% EndTime: 2019-12-05 15:06:45
% DurationCPUTime: 1.54s
% Computational Cost: add. (16408->199), mult. (28837->249), div. (0->0), fcn. (16415->8), ass. (0->85)
t187 = sin(qJ(4));
t189 = cos(qJ(4));
t162 = (-mrSges(6,1) * t189 + mrSges(6,2) * t187) * qJD(3);
t210 = qJD(3) * qJD(4);
t206 = t189 * t210;
t164 = t187 * qJDD(3) + t206;
t184 = sin(pkin(7));
t186 = cos(pkin(7));
t170 = -t186 * g(1) - t184 * g(2);
t182 = -g(3) + qJDD(1);
t183 = sin(pkin(8));
t185 = cos(pkin(8));
t137 = -t183 * t170 + t185 * t182;
t138 = t185 * t170 + t183 * t182;
t188 = sin(qJ(3));
t190 = cos(qJ(3));
t133 = t188 * t137 + t190 * t138;
t191 = qJD(3) ^ 2;
t130 = -t191 * pkin(3) + qJDD(3) * pkin(6) + t133;
t169 = t184 * g(1) - t186 * g(2);
t168 = qJDD(2) - t169;
t150 = t189 * t168;
t209 = qJD(3) * qJD(5);
t217 = pkin(4) * t191;
t122 = qJDD(4) * pkin(4) + t150 + (-t164 + t206) * qJ(5) + (t189 * t217 - t130 - 0.2e1 * t209) * t187;
t211 = qJD(3) * t189;
t174 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t211;
t208 = m(6) * t122 + qJDD(4) * mrSges(6,1) + qJD(4) * t174;
t212 = qJD(3) * t187;
t117 = -t164 * mrSges(6,3) - t162 * t212 + t208;
t126 = -t187 * t130 + t150;
t127 = t189 * t130 + t187 * t168;
t146 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t187 + Ifges(5,2) * t189) * qJD(3);
t147 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t187 + Ifges(6,4) * t189) * qJD(3);
t148 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t187 + Ifges(5,4) * t189) * qJD(3);
t165 = t189 * qJDD(3) - t187 * t210;
t171 = qJD(4) * pkin(4) - qJ(5) * t212;
t181 = t189 ^ 2;
t123 = t165 * qJ(5) - qJD(4) * t171 - t181 * t217 + 0.2e1 * t189 * t209 + t127;
t145 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t187 + Ifges(6,2) * t189) * qJD(3);
t197 = -mrSges(6,1) * t122 + mrSges(6,2) * t123 - Ifges(6,5) * t164 - Ifges(6,6) * t165 - Ifges(6,3) * qJDD(4) - t145 * t212;
t219 = mrSges(5,1) * t126 - mrSges(5,2) * t127 + Ifges(5,5) * t164 + Ifges(5,6) * t165 + Ifges(5,3) * qJDD(4) + pkin(4) * t117 - (-t187 * t146 + (t147 + t148) * t189) * qJD(3) - t197;
t216 = -mrSges(5,2) - mrSges(6,2);
t163 = (-mrSges(5,1) * t189 + mrSges(5,2) * t187) * qJD(3);
t175 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t211;
t113 = m(5) * t126 + qJDD(4) * mrSges(5,1) + qJD(4) * t175 + (-mrSges(5,3) - mrSges(6,3)) * t164 + (-t162 - t163) * t212 + t208;
t207 = m(6) * t123 + t165 * mrSges(6,3) + t162 * t211;
t172 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t212;
t213 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t212 - t172;
t114 = m(5) * t127 + t165 * mrSges(5,3) + t213 * qJD(4) + t216 * qJDD(4) + t163 * t211 + t207;
t203 = -t187 * t113 + t189 * t114;
t102 = m(4) * t133 - t191 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t203;
t132 = t190 * t137 - t188 * t138;
t200 = -qJDD(3) * pkin(3) - t132;
t129 = -t191 * pkin(6) + t200;
t125 = t171 * t212 - t165 * pkin(4) + qJDD(5) + (-qJ(5) * t181 - pkin(6)) * t191 + t200;
t202 = -m(6) * t125 + t165 * mrSges(6,1) + t174 * t211;
t193 = -m(5) * t129 + t165 * mrSges(5,1) + t216 * t164 + t175 * t211 + t213 * t212 + t202;
t109 = m(4) * t132 + qJDD(3) * mrSges(4,1) - t191 * mrSges(4,2) + t193;
t95 = t188 * t102 + t190 * t109;
t106 = t189 * t113 + t187 * t114;
t93 = m(3) * t137 + t95;
t204 = t190 * t102 - t188 * t109;
t94 = m(3) * t138 + t204;
t205 = -t183 * t93 + t185 * t94;
t201 = (-m(3) - m(4)) * t168 - t106;
t143 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t187 + Ifges(6,6) * t189) * qJD(3);
t144 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t187 + Ifges(5,6) * t189) * qJD(3);
t198 = -mrSges(6,1) * t125 + mrSges(6,3) * t123 + Ifges(6,4) * t164 + Ifges(6,2) * t165 + Ifges(6,6) * qJDD(4) + qJD(4) * t147;
t97 = Ifges(5,4) * t164 + Ifges(5,2) * t165 + Ifges(5,6) * qJDD(4) + qJD(4) * t148 - mrSges(5,1) * t129 + mrSges(5,3) * t127 - pkin(4) * (t164 * mrSges(6,2) - t202) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t172 + t207) + (-pkin(4) * t172 - t143 - t144) * t212 + t198;
t196 = mrSges(6,2) * t125 - mrSges(6,3) * t122 + Ifges(6,1) * t164 + Ifges(6,4) * t165 + Ifges(6,5) * qJDD(4) + t143 * t211;
t99 = t144 * t211 + mrSges(5,2) * t129 - mrSges(5,3) * t126 + Ifges(5,1) * t164 + Ifges(5,4) * t165 + Ifges(5,5) * qJDD(4) - qJ(5) * t117 + (-t145 - t146) * qJD(4) + t196;
t87 = mrSges(4,2) * t168 - mrSges(4,3) * t132 + Ifges(4,5) * qJDD(3) - t191 * Ifges(4,6) - pkin(6) * t106 - t187 * t97 + t189 * t99;
t91 = -mrSges(4,1) * t168 + mrSges(4,3) * t133 + t191 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t106 - t219;
t84 = -mrSges(3,1) * t168 + mrSges(3,3) * t138 + t188 * t87 + t190 * t91 - pkin(2) * (m(4) * t168 + t106) + pkin(5) * t204;
t86 = mrSges(3,2) * t168 - mrSges(3,3) * t137 - pkin(5) * t95 - t188 * t91 + t190 * t87;
t199 = mrSges(2,1) * t169 - mrSges(2,2) * t170 + pkin(1) * t201 + qJ(2) * t205 + t183 * t86 + t185 * t84;
t195 = mrSges(4,1) * t132 - mrSges(4,2) * t133 + Ifges(4,3) * qJDD(3) + pkin(3) * t193 + pkin(6) * t203 + t187 * t99 + t189 * t97;
t194 = mrSges(3,1) * t137 - mrSges(3,2) * t138 + pkin(2) * t95 + t195;
t103 = m(2) * t169 + t201;
t90 = t183 * t94 + t185 * t93;
t88 = m(2) * t170 + t205;
t82 = -mrSges(2,1) * t182 + mrSges(2,3) * t170 - pkin(1) * t90 - t194;
t81 = mrSges(2,2) * t182 - mrSges(2,3) * t169 - qJ(2) * t90 - t183 * t84 + t185 * t86;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t186 * t81 - t184 * t82 - qJ(1) * (t186 * t103 + t184 * t88), t81, t86, t87, t99, -qJD(4) * t145 + t196; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t184 * t81 + t186 * t82 + qJ(1) * (-t184 * t103 + t186 * t88), t82, t84, t91, t97, -t143 * t212 + t198; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t199, t199, t194, t195, t219, -t147 * t211 - t197;];
m_new = t1;
