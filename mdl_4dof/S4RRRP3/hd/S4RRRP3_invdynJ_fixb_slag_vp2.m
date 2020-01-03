% Calculate vector of inverse dynamics joint torques for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:05
% DurationCPUTime: 1.86s
% Computational Cost: add. (1119->222), mult. (1768->290), div. (0->0), fcn. (769->8), ass. (0->109)
t218 = Ifges(5,4) + Ifges(4,5);
t217 = Ifges(5,6) - Ifges(4,6);
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t103 = sin(qJ(2));
t182 = pkin(1) * qJD(1);
t154 = t103 * t182;
t98 = qJD(1) + qJD(2);
t66 = pkin(6) * t98 + t154;
t176 = t102 * t66;
t39 = -qJD(3) * pkin(3) + qJD(4) + t176;
t169 = t105 * t66;
t46 = qJD(3) * qJ(4) + t169;
t120 = -t102 * t46 + t105 * t39;
t222 = mrSges(4,1) + mrSges(5,1);
t221 = -mrSges(5,2) - mrSges(4,3);
t129 = t105 * mrSges(5,1) + t102 * mrSges(5,3);
t69 = -t105 * mrSges(4,1) + t102 * mrSges(4,2);
t220 = t69 - t129;
t219 = mrSges(3,2) + t221;
t183 = Ifges(5,5) * t105;
t127 = Ifges(5,1) * t102 - t183;
t173 = t102 * t98;
t166 = t105 * t98;
t77 = Ifges(4,4) * t166;
t216 = Ifges(4,1) * t173 + t218 * qJD(3) + t127 * t98 + t77;
t215 = t217 * t102 + t218 * t105;
t161 = qJD(3) * t102;
t106 = cos(qJ(2));
t181 = pkin(1) * qJD(2);
t148 = qJD(1) * t181;
t164 = pkin(1) * qJDD(1);
t59 = t103 * t164 + t106 * t148;
t97 = qJDD(1) + qJDD(2);
t45 = pkin(6) * t97 + t59;
t30 = t105 * t45;
t13 = -t161 * t66 + t30;
t160 = qJD(3) * t105;
t14 = -t102 * t45 - t160 * t66;
t121 = -t102 * t14 + t105 * t13;
t10 = -qJDD(3) * pkin(3) + qJDD(4) - t14;
t5 = qJDD(3) * qJ(4) + t30 + (qJD(4) - t176) * qJD(3);
t214 = t10 * t102 + t105 * t5;
t58 = -t103 * t148 + t106 * t164;
t44 = -pkin(2) * t97 - t58;
t52 = -t105 * t97 + t161 * t98;
t53 = t102 * t97 + t160 * t98;
t213 = m(4) * t44 + mrSges(4,1) * t52 + mrSges(4,2) * t53;
t128 = t102 * mrSges(5,1) - t105 * mrSges(5,3);
t130 = mrSges(4,1) * t102 + mrSges(4,2) * t105;
t153 = t106 * t182;
t163 = qJ(4) * t102;
t123 = pkin(3) * t105 + t163;
t68 = -pkin(2) - t123;
t22 = t68 * t98 - t153;
t67 = -pkin(2) * t98 - t153;
t212 = t120 * mrSges(5,2) + t22 * t128 + t67 * t130;
t65 = mrSges(5,2) * t166 + qJD(3) * mrSges(5,3);
t208 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t166 + t65;
t209 = t222 * qJD(3) + t221 * t173;
t211 = t102 * t209 - t105 * t208;
t33 = -qJDD(3) * mrSges(5,1) + t53 * mrSges(5,2);
t34 = -mrSges(5,2) * t52 + qJDD(3) * mrSges(5,3);
t210 = m(4) * t121 + m(5) * (qJD(3) * t120 + t214) + (-qJDD(3) * mrSges(4,1) + mrSges(4,3) * t53 + t33) * t102 + (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t52 + t34) * t105 - t161 * t208 - t160 * t209;
t101 = qJ(1) + qJ(2);
t94 = sin(t101);
t95 = cos(t101);
t207 = g(1) * t95 + g(2) * t94;
t206 = t219 * t95 + (-m(5) * t68 + mrSges(3,1) - t220) * t94;
t167 = t105 * t95;
t205 = -t222 * t167 + t219 * t94 + (-mrSges(3,1) + (mrSges(4,2) - mrSges(5,3)) * t102) * t95;
t191 = pkin(1) * t106;
t193 = pkin(1) * t103;
t200 = m(4) * pkin(1) * (t103 * t67 + (t102 ^ 2 + t105 ^ 2) * t66 * t106) + (-mrSges(3,1) * t193 - mrSges(3,2) * t191) * t98;
t197 = pkin(2) * t94;
t194 = t102 / 0.2e1;
t104 = sin(qJ(1));
t192 = pkin(1) * t104;
t107 = cos(qJ(1));
t96 = t107 * pkin(1);
t187 = t95 * pkin(2) + t94 * pkin(6);
t186 = Ifges(4,4) * t102;
t185 = Ifges(4,4) * t105;
t184 = Ifges(5,5) * t102;
t162 = qJD(3) * t98;
t159 = qJD(4) * t102;
t152 = t103 * t181;
t143 = -t162 / 0.2e1;
t87 = t95 * pkin(6);
t141 = t87 - t192;
t136 = pkin(3) * t167 + t95 * t163 + t187;
t126 = Ifges(4,2) * t105 + t186;
t122 = pkin(3) * t102 - qJ(4) * t105;
t115 = t102 * (Ifges(4,1) * t105 - t186);
t114 = t105 * (Ifges(5,3) * t102 + t183);
t113 = (t102 * t39 + t105 * t46) * t106;
t48 = pkin(3) * t161 - qJ(4) * t160 - t159;
t2 = pkin(3) * t52 - qJ(4) * t53 - t159 * t98 + t44;
t76 = Ifges(5,5) * t173;
t40 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t166 + t76;
t41 = Ifges(4,6) * qJD(3) + t126 * t98;
t108 = -t2 * t129 + t105 * (Ifges(4,4) * t53 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t105 * (Ifges(5,5) * t53 + Ifges(5,6) * qJDD(3)) / 0.2e1 + Ifges(3,3) * t97 - t59 * mrSges(3,2) + t44 * t69 + t58 * mrSges(3,1) - t41 * t161 / 0.2e1 + t115 * t162 / 0.2e1 + t114 * t143 + (Ifges(4,1) * t102 + t127 + t185) * t53 / 0.2e1 + ((Ifges(4,1) + Ifges(5,1)) * t53 + t218 * qJDD(3)) * t194 + (t218 * t102 - t217 * t105) * qJDD(3) / 0.2e1 + (t98 * (Ifges(5,1) * t105 + t184) + t40) * t161 / 0.2e1 + t121 * mrSges(4,3) + t214 * mrSges(5,2) + (-t126 / 0.2e1 + t184 / 0.2e1 + (Ifges(5,5) - Ifges(4,4)) * t194 + (-Ifges(4,2) / 0.2e1 - Ifges(5,3)) * t105) * t52 + (t98 * (-Ifges(4,2) * t102 + t185) + t216) * t160 / 0.2e1 + (t212 + t215 * qJD(3) / 0.2e1) * qJD(3);
t57 = t68 - t191;
t51 = t122 * t98;
t50 = t69 * t98;
t49 = t129 * t98;
t29 = t48 + t152;
t15 = mrSges(5,1) * t52 - mrSges(5,3) * t53;
t1 = [m(3) * (t103 * t59 + t106 * t58) * pkin(1) + t108 + t50 * t152 + t57 * t15 - t29 * t49 + m(5) * (t113 * t181 + t2 * t57 + t22 * t29) + Ifges(2,3) * qJDD(1) + (mrSges(3,1) * t191 - mrSges(3,2) * t193) * t97 + t200 * qJD(2) + (-mrSges(2,1) * t107 + t104 * mrSges(2,2) - m(5) * (t96 + t136) - m(4) * (t96 + t187) - m(3) * t96 + t205) * g(2) + (t104 * mrSges(2,1) + mrSges(2,2) * t107 + m(3) * t192 - m(4) * (t141 - t197) - m(5) * t141 + t206) * g(1) + t213 * (-pkin(2) - t191) - t211 * t106 * t181 + t210 * (pkin(6) + t193); t108 + t68 * t15 - t48 * t49 + (-t50 + t49) * t154 + (-(t103 * t22 + t113) * t182 + t2 * t68 + t22 * t48) * m(5) - t200 * qJD(1) + (-m(4) * t187 - m(5) * t136 + t205) * g(2) + (-m(5) * t87 - m(4) * (t87 - t197) + t206) * g(1) + t211 * t153 - t213 * pkin(2) + t210 * pkin(6); qJD(4) * t65 + t51 * t49 - pkin(3) * t33 + qJ(4) * t34 - t10 * mrSges(5,1) - t13 * mrSges(4,2) + t14 * mrSges(4,1) + t5 * mrSges(5,3) + t41 * t173 / 0.2e1 + t218 * t53 + t217 * t52 + t208 * t176 + t209 * t169 - (Ifges(5,1) * t166 + t40 + t76) * t173 / 0.2e1 + t215 * t143 + (Ifges(5,2) + Ifges(4,3)) * qJDD(3) + t220 * g(3) + (-t10 * pkin(3) - g(3) * t123 + t5 * qJ(4) + t46 * qJD(4) - t120 * t66 - t22 * t51) * m(5) - (-Ifges(4,2) * t173 + t216 + t77) * t166 / 0.2e1 + ((t114 / 0.2e1 - t115 / 0.2e1) * t98 - t212) * t98 + (m(5) * t122 + t128 + t130) * t207; -t49 * t173 - qJD(3) * t65 + (g(3) * t105 - t46 * qJD(3) - t102 * t207 + t22 * t173 + t10) * m(5) + t33;];
tau = t1;
