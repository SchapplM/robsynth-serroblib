% Calculate vector of inverse dynamics joint torques for
% S5RPPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:54:58
% EndTime: 2019-12-31 17:55:05
% DurationCPUTime: 4.22s
% Computational Cost: add. (1521->275), mult. (3044->333), div. (0->0), fcn. (1797->8), ass. (0->120)
t192 = mrSges(5,1) + mrSges(6,1);
t193 = mrSges(5,2) - mrSges(6,3);
t91 = sin(pkin(7));
t92 = cos(pkin(7));
t134 = t91 ^ 2 + t92 ^ 2;
t118 = t134 * mrSges(4,3);
t187 = Ifges(5,1) + Ifges(6,1);
t185 = Ifges(5,5) + Ifges(6,4);
t94 = -pkin(1) - qJ(3);
t191 = -qJD(1) * qJD(3) + qJDD(1) * t94;
t190 = Ifges(5,4) - Ifges(6,5);
t87 = pkin(7) + qJ(4);
t79 = sin(t87);
t80 = cos(t87);
t189 = t192 * t79 + t193 * t80;
t188 = -m(5) - m(6);
t184 = Ifges(6,6) - Ifges(5,6);
t97 = cos(qJ(4));
t141 = t97 * t92;
t95 = sin(qJ(4));
t180 = -t91 * t95 + t141;
t53 = t91 * t97 + t92 * t95;
t46 = t53 * qJD(1);
t29 = -qJD(4) * t46 + qJDD(1) * t180;
t15 = -qJDD(4) * mrSges(6,1) + t29 * mrSges(6,2);
t183 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t29 + t15;
t132 = qJD(1) * t91;
t122 = t95 * t132;
t130 = qJD(4) * t97;
t123 = t92 * t130;
t30 = qJD(1) * t123 - qJD(4) * t122 + qJDD(1) * t53;
t17 = -mrSges(6,2) * t30 + qJDD(4) * mrSges(6,3);
t182 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t30 + t17;
t152 = Ifges(6,5) * t46;
t45 = Ifges(5,4) * t46;
t47 = qJD(1) * t141 - t122;
t181 = t185 * qJD(4) + t187 * t47 + t152 - t45;
t146 = t46 * mrSges(5,3);
t147 = t46 * mrSges(6,2);
t41 = qJD(4) * mrSges(6,3) - t147;
t138 = -qJD(4) * mrSges(5,2) - t146 + t41;
t144 = t47 * mrSges(5,3);
t145 = t47 * mrSges(6,2);
t137 = qJD(4) * t192 - t144 - t145;
t89 = qJD(1) * qJD(2);
t62 = -qJDD(1) * qJ(2) - t89;
t96 = sin(qJ(1));
t98 = cos(qJ(1));
t179 = -g(1) * t96 + g(2) * t98;
t61 = qJD(1) * t94 + qJD(2);
t112 = -pkin(6) * qJD(1) + t61;
t42 = t112 * t91;
t43 = t112 * t92;
t13 = t42 * t97 + t43 * t95;
t55 = qJDD(2) + t191;
t111 = -pkin(6) * qJDD(1) + t55;
t36 = t111 * t91;
t37 = t111 * t92;
t5 = -qJD(4) * t13 - t36 * t95 + t37 * t97;
t148 = t42 * t95;
t12 = t43 * t97 - t148;
t7 = -qJD(4) * pkin(4) + qJD(5) - t12;
t175 = -m(6) * t7 + t137;
t10 = qJD(4) * qJ(5) + t13;
t174 = -m(6) * t10 - t138;
t124 = t43 * t130 + t97 * t36 + t95 * t37;
t131 = qJD(4) * t95;
t4 = -t131 * t42 + t124;
t48 = -t130 * t91 - t131 * t92;
t49 = -t131 * t91 + t123;
t173 = -t12 * t48 - t13 * t49 - t180 * t5 - t4 * t53;
t1 = qJDD(4) * qJ(5) + (qJD(5) - t148) * qJD(4) + t124;
t3 = -qJDD(4) * pkin(4) + qJDD(5) - t5;
t172 = -t1 * t53 - t10 * t49 + t180 * t3 + t48 * t7;
t171 = mrSges(3,2) - mrSges(2,1) - mrSges(6,2) - mrSges(4,3) - mrSges(5,3);
t104 = pkin(4) * t79 - qJ(5) * t80;
t110 = mrSges(4,1) * t91 + mrSges(4,2) * t92;
t170 = -m(6) * t104 + mrSges(2,2) - mrSges(3,3) - t110 - t189;
t77 = qJD(1) * qJ(2) + qJD(3);
t58 = pkin(3) * t132 + t77;
t169 = m(4) * t77 + m(5) * t58 + mrSges(5,1) * t46 + mrSges(5,2) * t47 + t110 * qJD(1);
t167 = -t46 / 0.2e1;
t166 = t46 / 0.2e1;
t164 = t47 / 0.2e1;
t162 = t180 / 0.2e1;
t81 = t91 * pkin(3);
t154 = -pkin(6) + t94;
t153 = Ifges(5,4) * t47;
t127 = qJDD(1) * t92;
t128 = qJDD(1) * t91;
t136 = mrSges(4,1) * t128 + mrSges(4,2) * t127;
t135 = t98 * pkin(1) + t96 * qJ(2);
t70 = qJ(2) + t81;
t129 = qJDD(1) * pkin(1);
t126 = -m(4) + t188;
t59 = qJDD(3) - t62;
t83 = t98 * qJ(2);
t119 = -pkin(1) * t96 + t83;
t117 = t134 * t61;
t116 = t134 * t55;
t114 = t30 * mrSges(5,1) + t29 * mrSges(5,2);
t113 = t30 * mrSges(6,1) - t29 * mrSges(6,3);
t50 = pkin(3) * t128 + t59;
t56 = t154 * t91;
t57 = t154 * t92;
t32 = t56 * t97 + t57 * t95;
t31 = t56 * t95 - t97 * t57;
t99 = qJD(1) ^ 2;
t93 = -pkin(6) - qJ(3);
t78 = qJDD(2) - t129;
t44 = Ifges(6,5) * t47;
t24 = mrSges(6,1) * t46 - mrSges(6,3) * t47;
t23 = pkin(4) * t47 + qJ(5) * t46;
t22 = pkin(4) * t53 - qJ(5) * t180 + t70;
t19 = -t46 * Ifges(5,2) + Ifges(5,6) * qJD(4) + t153;
t18 = Ifges(6,6) * qJD(4) + t46 * Ifges(6,3) + t44;
t11 = pkin(4) * t46 - qJ(5) * t47 + t58;
t6 = pkin(4) * t49 - qJ(5) * t48 - qJD(5) * t180 + qJD(2);
t2 = pkin(4) * t30 - qJ(5) * t29 - qJD(5) * t47 + t50;
t8 = [(Ifges(3,1) + Ifges(2,3)) * qJDD(1) + qJ(2) * t136 + m(3) * (-pkin(1) * t78 + (-t62 + t89) * qJ(2)) - (Ifges(4,4) * t92 - Ifges(4,2) * t91) * t128 + t22 * t113 + t70 * t114 + t59 * t110 + m(4) * (qJ(2) * t59 - qJD(3) * t117 + t116 * t94) + (t187 * t48 - t190 * t49) * t164 + (t180 * t187 - t190 * t53) * t29 / 0.2e1 + ((-m(4) - m(3)) * t135 + t188 * (t96 * t81 - t93 * t98 + t135) + (-m(4) * qJ(3) + t171) * t98 + t170 * t96) * g(2) + (-m(3) * t119 - m(4) * t83 + t188 * (t98 * t81 + t96 * t93 + t119) + t170 * t98 + (-m(4) * t94 - t171) * t96) * g(1) + (-m(5) * t5 + m(6) * t3 + t183) * t31 + (t184 * t49 + t185 * t48) * qJD(4) / 0.2e1 + (t185 * qJDD(4) + t187 * t29) * t162 + t181 * t48 / 0.2e1 + (m(5) * t4 + m(6) * t1 + t182) * t32 + (-m(5) * t12 - t175) * (t180 * qJD(3) + qJD(4) * t32) + t172 * mrSges(6,2) + t173 * mrSges(5,3) + (m(5) * t13 - t174) * (-qJD(3) * t53 - qJD(4) * t31) + t169 * qJD(2) + m(6) * (t11 * t6 + t2 * t22) + (t58 * mrSges(5,1) + t11 * mrSges(6,1) + t18 / 0.2e1 - t19 / 0.2e1 + Ifges(6,3) * t166 - Ifges(5,2) * t167) * t49 - t53 * (Ifges(5,4) * t29 + Ifges(5,6) * qJDD(4)) / 0.2e1 - 0.2e1 * t62 * mrSges(3,3) + t6 * t24 + (-t129 + t78) * mrSges(3,2) + ((Ifges(5,2) + Ifges(6,3)) * t53 + t190 * (-t180 / 0.2e1 - t162)) * t30 + (m(5) * t70 + mrSges(5,1) * t53 + mrSges(5,2) * t180) * t50 + (t180 * t185 + t184 * t53) * qJDD(4) / 0.2e1 + t2 * (mrSges(6,1) * t53 - mrSges(6,3) * t180) + (t58 * mrSges(5,2) - t11 * mrSges(6,3) + Ifges(5,4) * t167 + Ifges(6,5) * t166) * t48 + (Ifges(4,1) * t92 - Ifges(4,4) * t91) * t127 + (-t55 - t191) * t118 + t53 * (Ifges(6,5) * t29 + Ifges(6,6) * qJDD(4)) / 0.2e1; (-m(3) * qJ(2) - mrSges(3,3)) * t99 + t182 * t53 - t183 * t180 + t138 * t49 + t137 * t48 + (mrSges(3,2) - t118) * qJDD(1) + m(3) * t78 - m(6) * t172 - m(5) * t173 + m(4) * t116 + (-m(6) * t11 - t169 - t24) * qJD(1) + t179 * (m(3) - t126); -t99 * t118 + t137 * t47 + t138 * t46 + t113 + t114 + t136 + (g(1) * t98 + g(2) * t96) * t126 + (t10 * t46 - t47 * t7 + t2) * m(6) + (t12 * t47 + t13 * t46 + t50) * m(5) + (qJD(1) * t117 + t59) * m(4); t189 * g(3) + (Ifges(5,3) + Ifges(6,2)) * qJDD(4) + t184 * t30 + t185 * t29 - (t184 * t47 - t185 * t46) * qJD(4) / 0.2e1 - (-t187 * t46 - t153 + t18 + t44) * t47 / 0.2e1 + (-Ifges(5,2) * t47 + t181 - t45) * t166 + ((pkin(4) * m(6) + t192) * t80 + (qJ(5) * m(6) - t193) * t79) * t179 + (-t146 + t174) * t12 + (t144 + t175) * t13 - t58 * (mrSges(5,1) * t47 - mrSges(5,2) * t46) - t11 * (mrSges(6,1) * t47 + mrSges(6,3) * t46) + qJD(5) * t41 - t23 * t24 - pkin(4) * t15 + qJ(5) * t17 + t1 * mrSges(6,3) - t3 * mrSges(6,1) - t4 * mrSges(5,2) + t5 * mrSges(5,1) + (-pkin(4) * t3 + g(3) * t104 + qJ(5) * t1 + qJD(5) * t10 - t11 * t23) * m(6) + t10 * t145 + t7 * t147 + t19 * t164 + (Ifges(6,3) * t47 - t152) * t167; -qJD(4) * t41 + t47 * t24 + (-g(3) * t79 - t10 * qJD(4) + t11 * t47 - t179 * t80 + t3) * m(6) + t15;];
tau = t8;
