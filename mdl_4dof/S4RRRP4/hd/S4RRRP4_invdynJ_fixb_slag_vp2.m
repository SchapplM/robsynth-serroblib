% Calculate vector of inverse dynamics joint torques for
% S4RRRP4
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
% Datum: 2019-12-31 17:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:15:10
% EndTime: 2019-12-31 17:15:18
% DurationCPUTime: 4.47s
% Computational Cost: add. (1715->289), mult. (4075->376), div. (0->0), fcn. (2529->8), ass. (0->131)
t206 = Ifges(4,4) + Ifges(5,4);
t126 = sin(qJ(2));
t129 = cos(qJ(2));
t106 = -mrSges(3,1) * t129 + mrSges(3,2) * t126;
t124 = qJ(2) + qJ(3);
t119 = sin(t124);
t120 = cos(t124);
t171 = mrSges(4,2) + mrSges(5,2);
t172 = mrSges(4,1) + mrSges(5,1);
t142 = t119 * t171 - t172 * t120;
t214 = -t106 - t142;
t207 = Ifges(4,1) + Ifges(5,1);
t205 = Ifges(4,5) + Ifges(5,5);
t204 = Ifges(4,2) + Ifges(5,2);
t203 = Ifges(4,6) + Ifges(5,6);
t125 = sin(qJ(3));
t128 = cos(qJ(3));
t91 = -t125 * t126 + t128 * t129;
t79 = t91 * qJD(1);
t213 = t206 * t79;
t92 = t125 * t129 + t126 * t128;
t80 = t92 * qJD(1);
t212 = t206 * t80;
t123 = qJD(2) + qJD(3);
t211 = t203 * t123 + t204 * t79 + t212;
t210 = t205 * t123 + t207 * t80 + t213;
t157 = qJD(1) * qJD(2);
t101 = qJDD(1) * t129 - t126 * t157;
t127 = sin(qJ(1));
t130 = cos(qJ(1));
t201 = g(1) * t130 + g(2) * t127;
t114 = pkin(2) * t129 + pkin(1);
t107 = t114 * qJD(1);
t131 = -pkin(6) - pkin(5);
t109 = t131 * t129;
t98 = qJD(1) * t109;
t81 = t125 * t98;
t108 = t131 * t126;
t97 = qJD(1) * t108;
t87 = qJD(2) * pkin(2) + t97;
t51 = t128 * t87 + t81;
t73 = t80 * qJ(4);
t26 = t51 - t73;
t17 = pkin(3) * t123 + t26;
t61 = -pkin(3) * t79 + qJD(4) - t107;
t208 = mrSges(4,2) * t107 - mrSges(5,2) * t61 + t51 * mrSges(4,3) + t17 * mrSges(5,3);
t184 = t126 / 0.2e1;
t164 = qJDD(1) * pkin(1);
t89 = t101 * pkin(5);
t102 = qJDD(1) * t126 + t129 * t157;
t90 = t102 * pkin(5);
t202 = t126 * t90 + t129 * t89;
t168 = qJ(4) * t79;
t84 = t128 * t98;
t52 = t125 * t87 - t84;
t27 = t52 + t168;
t146 = t52 * mrSges(4,3) + t27 * mrSges(5,3);
t112 = pkin(3) * t120;
t200 = m(4) * t114 + mrSges(2,1) + m(5) * (t112 + t114) + m(3) * pkin(1) + t214;
t199 = mrSges(2,2) + m(5) * (-qJ(4) + t131) - mrSges(5,3) + m(4) * t131 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t191 = t80 / 0.2e1;
t188 = pkin(3) * t80;
t182 = pkin(2) * t126;
t179 = g(3) * t129;
t56 = t128 * t97 + t81;
t64 = t125 * t108 - t128 * t109;
t170 = Ifges(3,4) * t126;
t169 = Ifges(3,4) * t129;
t166 = t129 * Ifges(3,2);
t163 = qJD(1) * t126;
t162 = qJD(1) * t129;
t161 = qJD(2) * t126;
t160 = qJD(2) * t129;
t159 = qJD(3) * t125;
t158 = qJD(3) * t128;
t155 = pkin(2) * t161;
t153 = qJD(2) * t131;
t138 = t91 * qJD(3);
t36 = qJD(1) * t138 + t101 * t125 + t102 * t128;
t139 = t92 * qJD(3);
t37 = -qJD(1) * t139 + t101 * t128 - t102 * t125;
t150 = -t37 * mrSges(5,1) + t36 * mrSges(5,2);
t149 = t171 * t120;
t55 = -t125 * t97 + t84;
t63 = t128 * t108 + t109 * t125;
t145 = mrSges(3,1) * t126 + mrSges(3,2) * t129;
t144 = t166 + t170;
t143 = Ifges(3,5) * t129 - Ifges(3,6) * t126;
t74 = -pkin(2) * t101 - t164;
t141 = pkin(1) * t145;
t60 = qJDD(2) * pkin(2) - pkin(6) * t102 - t90;
t62 = pkin(6) * t101 + t89;
t7 = t125 * t60 + t128 * t62 + t87 * t158 + t159 * t98;
t100 = t129 * t153;
t99 = t126 * t153;
t11 = t125 * t100 + t108 * t158 + t109 * t159 + t128 * t99;
t140 = t126 * (Ifges(3,1) * t129 - t170);
t8 = -qJD(3) * t52 - t125 * t62 + t128 * t60;
t12 = -qJD(3) * t64 + t128 * t100 - t125 * t99;
t121 = qJDD(2) + qJDD(3);
t2 = pkin(3) * t121 - qJ(4) * t36 - qJD(4) * t80 + t8;
t3 = qJ(4) * t37 + qJD(4) * t79 + t7;
t132 = -t7 * mrSges(4,2) - t3 * mrSges(5,2) + (-t61 * t80 + t2) * mrSges(5,1) + (t107 * t80 + t8) * mrSges(4,1) + t208 * t79 + t203 * t37 + t205 * t36 - (t207 * t79 - t212) * t80 / 0.2e1 + t211 * t191 - (-t203 * t80 + t205 * t79) * t123 / 0.2e1 + (Ifges(5,3) + Ifges(4,3)) * t121 - (-t204 * t80 + t210 + t213) * t79 / 0.2e1;
t116 = Ifges(3,4) * t162;
t105 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t162;
t104 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t163;
t78 = Ifges(3,1) * t163 + Ifges(3,5) * qJD(2) + t116;
t77 = Ifges(3,6) * qJD(2) + qJD(1) * t144;
t70 = -pkin(3) * t91 - t114;
t69 = mrSges(4,1) * t123 - t80 * mrSges(4,3);
t68 = mrSges(5,1) * t123 - t80 * mrSges(5,3);
t67 = -mrSges(4,2) * t123 + mrSges(4,3) * t79;
t66 = -mrSges(5,2) * t123 + mrSges(5,3) * t79;
t65 = pkin(2) * t163 + t188;
t58 = -qJD(2) * t92 - t139;
t57 = qJD(2) * t91 + t138;
t50 = -mrSges(4,1) * t79 + mrSges(4,2) * t80;
t49 = -mrSges(5,1) * t79 + mrSges(5,2) * t80;
t46 = -pkin(3) * t58 + t155;
t45 = qJ(4) * t91 + t64;
t44 = -qJ(4) * t92 + t63;
t35 = -t73 + t56;
t34 = t55 - t168;
t21 = -mrSges(4,2) * t121 + mrSges(4,3) * t37;
t20 = -mrSges(5,2) * t121 + mrSges(5,3) * t37;
t19 = mrSges(4,1) * t121 - mrSges(4,3) * t36;
t18 = mrSges(5,1) * t121 - mrSges(5,3) * t36;
t10 = -pkin(3) * t37 + qJDD(4) + t74;
t5 = -qJ(4) * t57 - qJD(4) * t92 + t12;
t4 = qJ(4) * t58 + qJD(4) * t91 + t11;
t1 = [t211 * t58 / 0.2e1 + t210 * t57 / 0.2e1 - t208 * t57 + (-mrSges(4,1) * t74 - mrSges(5,1) * t10 + mrSges(4,3) * t7 + mrSges(5,3) * t3 + t121 * t203 + t204 * t37 + t206 * t36) * t91 + (m(3) * t164 + mrSges(3,1) * t101 - mrSges(3,2) * t102) * pkin(1) + m(5) * (t10 * t70 + t17 * t5 + t2 * t44 + t27 * t4 + t3 * t45 + t46 * t61) + (mrSges(4,1) * t107 - mrSges(5,1) * t61 + t146) * t58 + (t129 * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t101) - t126 * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t102) - t104 * t160 - t105 * t161 + m(3) * t202) * pkin(5) + t202 * mrSges(3,3) - t114 * (-mrSges(4,1) * t37 + mrSges(4,2) * t36) + t129 * (Ifges(3,4) * t102 + Ifges(3,2) * t101) / 0.2e1 + t63 * t19 + t64 * t21 + t4 * t66 + t11 * t67 + t5 * t68 + t12 * t69 + t46 * t49 + t44 * t18 + t45 * t20 + (t127 * t200 + t130 * t199) * g(1) + (t127 * t199 - t130 * t200) * g(2) + (mrSges(4,2) * t74 + mrSges(5,2) * t10 - mrSges(4,3) * t8 - mrSges(5,3) * t2 + t121 * t205 + t206 * t37 + t207 * t36) * t92 + (0.2e1 * Ifges(3,5) * t184 + Ifges(3,6) * t129) * qJDD(2) + (t203 * t58 + t205 * t57) * t123 / 0.2e1 + (t204 * t58 + t206 * t57) * t79 / 0.2e1 + (t206 * t58 + t207 * t57) * t191 + qJD(2) ^ 2 * t143 / 0.2e1 + t101 * t144 / 0.2e1 + (Ifges(3,1) * t102 + Ifges(3,4) * t101) * t184 + (t129 * (-Ifges(3,2) * t126 + t169) + t140) * t157 / 0.2e1 + m(4) * (-t107 * t155 + t11 * t52 - t114 * t74 + t12 * t51 + t63 * t8 + t64 * t7) + t70 * t150 - t141 * t157 + t78 * t160 / 0.2e1 - t77 * t161 / 0.2e1 + t50 * t155 - t106 * t164 + t102 * (t126 * Ifges(3,1) + t169) / 0.2e1 + Ifges(2,3) * qJDD(1); t146 * t80 + t132 - m(4) * (t51 * t55 + t52 * t56) + (-m(5) * t112 - t214) * g(3) + Ifges(3,6) * t101 + Ifges(3,5) * t102 - t89 * mrSges(3,2) - t90 * mrSges(3,1) - t65 * t49 - t35 * t66 - t56 * t67 - t34 * t68 - t55 * t69 + (-qJD(2) * t143 / 0.2e1 + t77 * t184 + (t166 * t184 - t140 / 0.2e1 + t141) * qJD(1) + (t104 * t129 + t105 * t126) * pkin(5) + (m(4) * t107 - t50) * t182 - (t116 + t78) * t129 / 0.2e1) * qJD(1) + (-m(5) * t179 + (t126 * t201 + t52 * t158 - t51 * t159 - t179) * m(4) + (t19 + (m(5) * t27 + t66 + t67) * qJD(3) + t8 * m(4)) * t128 + (m(5) * t3 + t20 + t21 + (-m(5) * t17 - t68 - t69) * qJD(3) + t7 * m(4)) * t125) * pkin(2) - m(5) * (t17 * t34 + t27 * t35 + t61 * t65) + Ifges(3,3) * qJDD(2) + t201 * (-m(5) * (-pkin(3) * t119 - t182) + t119 * t172 + t145 + t149) + (m(5) * t2 + t18) * (pkin(2) * t128 + pkin(3)); t132 + t142 * g(3) + (-pkin(3) * t49 + t146) * t80 - t26 * t66 - t51 * t67 + t27 * t68 + t52 * t69 + pkin(3) * t18 + (-g(3) * t112 - t61 * t188 - (-t17 + t26) * t27 + t2 * pkin(3)) * m(5) + t201 * (t149 + (m(5) * pkin(3) + t172) * t119); -t79 * t66 + t80 * t68 + (-g(1) * t127 + g(2) * t130 + t17 * t80 - t27 * t79 + t10) * m(5) + t150;];
tau = t1;
