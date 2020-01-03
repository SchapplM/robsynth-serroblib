% Calculate vector of inverse dynamics joint torques for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:12:57
% DurationCPUTime: 1.61s
% Computational Cost: add. (1130->249), mult. (1752->323), div. (0->0), fcn. (778->8), ass. (0->115)
t179 = Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1;
t107 = sin(qJ(3));
t110 = cos(qJ(3));
t161 = mrSges(4,2) + mrSges(5,2);
t162 = mrSges(4,1) + mrSges(5,1);
t178 = t107 * t161 - t110 * t162 - mrSges(3,1);
t177 = Ifges(4,1) + Ifges(5,1);
t175 = Ifges(5,5) + Ifges(4,5);
t174 = Ifges(4,2) + Ifges(5,2);
t173 = Ifges(5,6) + Ifges(4,6);
t172 = (Ifges(4,4) + Ifges(5,4)) * t110;
t171 = mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t170 = m(4) * pkin(2);
t111 = cos(qJ(2));
t165 = pkin(1) * t111;
t109 = sin(qJ(1));
t164 = g(1) * t109;
t138 = qJD(3) * t110;
t102 = qJD(1) + qJD(2);
t108 = sin(qJ(2));
t155 = pkin(1) * qJD(1);
t68 = pkin(6) * t102 + t108 * t155;
t132 = t68 * t138;
t101 = qJDD(1) + qJDD(2);
t154 = pkin(1) * qJD(2);
t133 = qJD(1) * t154;
t144 = pkin(1) * qJDD(1);
t60 = t108 * t144 + t111 * t133;
t46 = pkin(6) * t101 + t60;
t11 = -t107 * t46 - t132;
t163 = t11 * mrSges(4,3);
t106 = -qJ(4) - pkin(6);
t105 = qJ(1) + qJ(2);
t97 = sin(t105);
t98 = cos(t105);
t160 = t98 * pkin(2) + t97 * pkin(6);
t159 = Ifges(4,4) * t107;
t157 = Ifges(5,4) * t107;
t153 = qJD(3) * pkin(3);
t139 = qJD(3) * t107;
t10 = t110 * t46 - t139 * t68;
t152 = t10 * t110;
t151 = t107 * t11;
t91 = pkin(1) * t108 + pkin(6);
t149 = t110 * t91;
t145 = -qJ(4) - t91;
t143 = t102 * t107;
t142 = t102 * t110;
t141 = t107 * t111;
t140 = t110 * t111;
t96 = t110 * qJD(4);
t136 = t111 * t155;
t135 = t111 * t154;
t134 = pkin(3) * t139;
t92 = pkin(3) * t110 + pkin(2);
t131 = t102 * t139;
t52 = t101 * t110 - t131;
t53 = t101 * t107 + t102 * t138;
t13 = -t52 * mrSges(5,1) + t53 * mrSges(5,2);
t128 = -t106 * t97 + t98 * t92;
t127 = qJD(3) * t106;
t126 = qJ(4) * t102 + t68;
t125 = qJD(3) * t145;
t59 = -t108 * t133 + t111 * t144;
t73 = -t110 * mrSges(4,1) + t107 * mrSges(4,2);
t122 = -t110 * mrSges(5,1) + t107 * mrSges(5,2);
t121 = Ifges(4,2) * t110 + t159;
t120 = Ifges(5,2) * t110 + t157;
t29 = t126 * t107;
t23 = -t29 + t153;
t30 = t126 * t110;
t119 = -t107 * t30 - t110 * t23;
t65 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t143;
t67 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t142;
t118 = t107 * t67 + t110 * t65;
t45 = -pkin(2) * t101 - t59;
t69 = -pkin(2) * t102 - t136;
t116 = m(4) * (t108 * t69 + (t107 ^ 2 + t110 ^ 2) * t111 * t68);
t115 = t171 * t97 + t178 * t98;
t114 = (-m(4) * pkin(6) + m(5) * t106 + t171) * t98 + (m(5) * t92 + t170 - t178) * t97;
t12 = -pkin(3) * t52 + qJDD(4) + t45;
t3 = qJ(4) * t52 + t102 * t96 + t10;
t36 = -t102 * t92 + qJD(4) - t136;
t41 = Ifges(5,6) * qJD(3) + t102 * t120;
t42 = Ifges(4,6) * qJD(3) + t102 * t121;
t79 = Ifges(5,4) * t142;
t43 = Ifges(5,1) * t143 + Ifges(5,5) * qJD(3) + t79;
t80 = Ifges(4,4) * t142;
t44 = Ifges(4,1) * t143 + Ifges(4,5) * qJD(3) + t80;
t113 = t3 * t110 * mrSges(5,3) + t59 * mrSges(3,1) - t60 * mrSges(3,2) + mrSges(4,3) * t152 + Ifges(3,3) * t101 + t12 * t122 + t45 * t73 + (-t173 * t107 + t175 * t110) * qJD(3) ^ 2 / 0.2e1 - (t42 + t41) * t139 / 0.2e1 + (t177 * t110 - t157 - t159) * t131 / 0.2e1 + (t36 * (mrSges(5,1) * t107 + mrSges(5,2) * t110) + t69 * (mrSges(4,1) * t107 + mrSges(4,2) * t110)) * qJD(3) + (t121 / 0.2e1 + t120 / 0.2e1 + t107 * t179 + t174 * t110 / 0.2e1) * t52 + (t177 * t107 + t172 / 0.2e1 + t110 * t179) * t53 + (t44 + t43 + (-t174 * t107 + t172) * t102) * t138 / 0.2e1 + (t175 * t107 + t173 * t110) * qJDD(3);
t112 = cos(qJ(1));
t100 = t112 * pkin(1);
t99 = t110 * qJ(4);
t93 = -pkin(2) - t165;
t74 = pkin(6) * t110 + t99;
t72 = t106 * t107;
t71 = -t92 - t165;
t66 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t142;
t64 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t143;
t63 = t108 * t154 + t134;
t58 = t99 + t149;
t57 = t145 * t107;
t51 = t73 * t102;
t50 = t122 * t102;
t49 = -qJD(4) * t107 + t110 * t127;
t48 = t107 * t127 + t96;
t35 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t53;
t34 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t53;
t33 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t52;
t32 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t52;
t17 = (-qJD(4) - t135) * t107 + t110 * t125;
t16 = t107 * t125 + t110 * t135 + t96;
t14 = -mrSges(4,1) * t52 + mrSges(4,2) * t53;
t2 = -t132 + qJDD(3) * pkin(3) - qJ(4) * t53 + (-qJD(4) * t102 - t46) * t107;
t1 = [(-mrSges(2,1) * t112 + t109 * mrSges(2,2) - m(5) * (t100 + t128) - m(4) * (t100 + t160) + t115) * g(2) + (mrSges(5,3) * t119 - t118 * t91) * qJD(3) + (-t2 * mrSges(5,3) - t91 * t35 - t163) * t107 + t113 + m(4) * (t10 * t149 - t91 * t151 + t45 * t93) + m(5) * (t12 * t71 + t16 * t30 + t17 * t23 + t2 * t57 + t3 * t58 + t36 * t63) + t93 * t14 + t57 * t34 + t58 * t32 + t63 * t50 + t17 * t64 + t16 * t66 + t71 * t13 + (t109 * mrSges(2,1) + mrSges(2,2) * t112 + t114) * g(1) + ((t111 * mrSges(3,1) - t108 * mrSges(3,2)) * t101 + (m(4) + m(5)) * t164 + (-g(2) * t112 + t108 * t60 + t111 * t59 + t164) * m(3) + (t108 * t51 - t65 * t141 + t67 * t140 + t116 + (-t108 * mrSges(3,1) - t111 * mrSges(3,2)) * t102) * qJD(2)) * pkin(1) + t33 * t149 + Ifges(2,3) * qJDD(1); m(5) * (-t12 * t92 + t36 * t134 + t2 * t72 + t23 * t49 + t3 * t74 + t30 * t48) + t113 + (t119 * qJD(3) - t107 * t2) * mrSges(5,3) + t114 * g(1) + (-m(4) * t160 - m(5) * t128 + t115) * g(2) + (t50 * t153 - t163) * t107 - t92 * t13 + t72 * t34 + t74 * t32 + t49 * t64 + t48 * t66 - pkin(2) * t14 + ((t102 * mrSges(3,1) - t50 - t51) * t108 + (t102 * mrSges(3,2) + (-t66 - t67) * t110 + (t64 + t65) * t107) * t111 - t116 - m(5) * (t108 * t36 + t140 * t30 - t141 * t23)) * t155 + (-t107 * t35 + t110 * t33 + m(4) * (-t151 + t152) - t118 * qJD(3)) * pkin(6) - t45 * t170; t11 * mrSges(4,1) + t2 * mrSges(5,1) - t10 * mrSges(4,2) - t3 * mrSges(5,2) + t29 * t66 + t118 * t68 + t175 * t53 + t173 * t52 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (t122 + t73) * g(3) + (t34 + (-g(3) * t110 + t2) * m(5)) * pkin(3) + (t64 - m(5) * (-t23 - t29)) * t30 + ((-t69 * mrSges(4,2) - t36 * mrSges(5,2) + t23 * mrSges(5,3) - t43 / 0.2e1 - t44 / 0.2e1 - t79 / 0.2e1 - t80 / 0.2e1 + (-Ifges(5,5) / 0.2e1 - Ifges(4,5) / 0.2e1) * qJD(3)) * t110 + (-t69 * mrSges(4,1) - t36 * mrSges(5,1) + t30 * mrSges(5,3) + t41 / 0.2e1 + t42 / 0.2e1 + t179 * t143 + (Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1) * qJD(3) + (-m(5) * t36 - t50) * pkin(3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(4,1) / 0.2e1) * t142) * t107) * t102 + (g(1) * t98 + g(2) * t97) * (t161 * t110 + (m(5) * pkin(3) + t162) * t107); (t107 * t64 - t110 * t66) * t102 + (t12 - g(1) * t97 + g(2) * t98 - (-t107 * t23 + t110 * t30) * t102) * m(5) + t13;];
tau = t1;
