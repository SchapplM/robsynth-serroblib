% Calculate vector of inverse dynamics joint torques for
% S5PRRRP1
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:48
% EndTime: 2019-12-05 16:39:55
% DurationCPUTime: 2.03s
% Computational Cost: add. (1435->283), mult. (2127->355), div. (0->0), fcn. (980->8), ass. (0->123)
t189 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t114 = sin(qJ(4));
t116 = cos(qJ(4));
t167 = mrSges(5,2) + mrSges(6,2);
t168 = mrSges(5,1) + mrSges(6,1);
t188 = t114 * t167 - t168 * t116 - mrSges(4,1);
t187 = Ifges(5,1) + Ifges(6,1);
t185 = Ifges(6,5) + Ifges(5,5);
t184 = Ifges(5,2) + Ifges(6,2);
t183 = Ifges(6,6) + Ifges(5,6);
t110 = qJDD(2) + qJDD(3);
t115 = sin(qJ(3));
t117 = cos(qJ(3));
t158 = pkin(2) * qJD(3);
t138 = qJD(2) * t158;
t152 = pkin(2) * qJDD(2);
t62 = t115 * t152 + t117 * t138;
t48 = pkin(7) * t110 + t62;
t182 = qJD(1) * qJD(4) + t48;
t181 = (Ifges(5,4) + Ifges(6,4)) * t116;
t180 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t112 = qJD(2) + qJD(3);
t151 = t112 * t114;
t67 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t151;
t68 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t151;
t150 = t112 * t116;
t69 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t150;
t70 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t150;
t179 = -(t67 + t68) * t114 + (t69 + t70) * t116;
t178 = m(5) * pkin(3);
t175 = m(5) + m(6);
t172 = pkin(2) * t117;
t111 = pkin(8) + qJ(2);
t102 = sin(t111);
t171 = g(1) * t102;
t146 = qJD(4) * t114;
t159 = pkin(2) * qJD(2);
t71 = pkin(7) * t112 + t115 * t159;
t5 = t114 * qJDD(1) + t182 * t116 - t71 * t146;
t170 = t116 * t5;
t105 = t116 * qJDD(1);
t147 = qJD(1) * t114;
t42 = t116 * t71 + t147;
t6 = -qJD(4) * t42 - t114 * t48 + t105;
t169 = t6 * t114;
t113 = -qJ(5) - pkin(7);
t106 = qJ(3) + t111;
t94 = sin(t106);
t95 = cos(t106);
t164 = t95 * pkin(3) + t94 * pkin(7);
t163 = Ifges(5,4) * t114;
t161 = Ifges(6,4) * t114;
t97 = pkin(2) * t115 + pkin(7);
t155 = t116 * t97;
t153 = -qJ(5) - t97;
t149 = t114 * t117;
t148 = t116 * t117;
t145 = qJD(4) * t116;
t107 = t116 * qJD(5);
t144 = m(2) + m(3) + m(4);
t141 = t117 * t159;
t140 = t117 * t158;
t139 = pkin(4) * t146;
t98 = t116 * pkin(4) + pkin(3);
t137 = t112 * t146;
t54 = t110 * t116 - t137;
t55 = t110 * t114 + t112 * t145;
t13 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t134 = -t113 * t94 + t95 * t98;
t133 = qJD(4) * t113;
t132 = qJ(5) * t112 + t71;
t131 = qJD(4) * t153;
t61 = -t115 * t138 + t117 * t152;
t79 = -t116 * mrSges(5,1) + t114 * mrSges(5,2);
t128 = -t116 * mrSges(6,1) + t114 * mrSges(6,2);
t127 = Ifges(5,2) * t116 + t163;
t126 = Ifges(6,2) * t116 + t161;
t108 = t116 * qJD(1);
t20 = -t132 * t114 + t108;
t16 = qJD(4) * pkin(4) + t20;
t21 = t132 * t116 + t147;
t125 = -t21 * t114 - t16 * t116;
t124 = -t114 * t16 + t116 * t21;
t41 = -t114 * t71 + t108;
t123 = -t42 * t114 - t41 * t116;
t47 = -t110 * pkin(3) - t61;
t72 = -t112 * pkin(3) - t141;
t121 = m(5) * (t115 * t72 + t42 * t148 - t41 * t149);
t120 = t180 * t94 + t188 * t95;
t119 = (-m(5) * pkin(7) + m(6) * t113 + t180) * t95 + (m(6) * t98 + t178 - t188) * t94;
t12 = -t54 * pkin(4) + qJDD(5) + t47;
t3 = qJ(5) * t54 + t112 * t107 + t5;
t36 = -t98 * t112 + qJD(5) - t141;
t43 = Ifges(6,6) * qJD(4) + t126 * t112;
t44 = Ifges(5,6) * qJD(4) + t127 * t112;
t81 = Ifges(6,4) * t150;
t45 = Ifges(6,1) * t151 + Ifges(6,5) * qJD(4) + t81;
t82 = Ifges(5,4) * t150;
t46 = Ifges(5,1) * t151 + Ifges(5,5) * qJD(4) + t82;
t118 = t3 * t116 * mrSges(6,3) + t61 * mrSges(4,1) - t62 * mrSges(4,2) + mrSges(5,3) * t170 + Ifges(4,3) * t110 + t12 * t128 + t47 * t79 + (-t183 * t114 + t185 * t116) * qJD(4) ^ 2 / 0.2e1 - (t44 + t43) * t146 / 0.2e1 + (t187 * t116 - t161 - t163) * t137 / 0.2e1 + (t36 * (mrSges(6,1) * t114 + mrSges(6,2) * t116) + t72 * (mrSges(5,1) * t114 + mrSges(5,2) * t116)) * qJD(4) + (t127 / 0.2e1 + t126 / 0.2e1 + t114 * t189 + t184 * t116 / 0.2e1) * t54 + (t187 * t114 + t181 / 0.2e1 + t116 * t189) * t55 + (t46 + t45 + (-t184 * t114 + t181) * t112) * t145 / 0.2e1 + (t185 * t114 + t183 * t116) * qJDD(4);
t109 = t116 * qJ(5);
t103 = cos(t111);
t99 = -pkin(3) - t172;
t93 = pkin(2) * t103;
t80 = pkin(7) * t116 + t109;
t78 = t113 * t114;
t77 = -t98 - t172;
t66 = t115 * t158 + t139;
t60 = t109 + t155;
t59 = t153 * t114;
t53 = t79 * t112;
t52 = t128 * t112;
t51 = -t114 * qJD(5) + t116 * t133;
t50 = t114 * t133 + t107;
t35 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t55;
t34 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t55;
t33 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t54;
t32 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t54;
t18 = (-qJD(5) - t140) * t114 + t116 * t131;
t17 = t114 * t131 + t116 * t140 + t107;
t14 = -mrSges(5,1) * t54 + mrSges(5,2) * t55;
t2 = -t71 * t145 + qJDD(4) * pkin(4) - qJ(5) * t55 + t105 + (-qJD(5) * t112 - t182) * t114;
t1 = [(t34 + t35) * t116 + (t32 + t33) * t114 + t179 * qJD(4) + m(5) * (t114 * t5 + t116 * t6 + (-t114 * t41 + t116 * t42) * qJD(4)) + m(6) * (t124 * qJD(4) + t114 * t3 + t116 * t2) + t144 * qJDD(1) + (-t144 - t175) * g(3); m(6) * (t12 * t77 + t16 * t18 + t17 * t21 + t2 * t59 + t3 * t60 + t36 * t66) + ((t117 * mrSges(4,1) - t115 * mrSges(4,2)) * t110 + t175 * t171 + (-g(2) * t103 + t115 * t62 + t117 * t61 + t171) * m(4) + (t115 * t53 + t121 - t68 * t149 + t70 * t148 + (-t115 * mrSges(4,1) - t117 * mrSges(4,2)) * t112) * qJD(3)) * pkin(2) + (t125 * mrSges(6,3) + t123 * mrSges(5,3) + (m(5) * t123 - t114 * t70 - t116 * t68) * t97) * qJD(4) + (mrSges(3,1) * t102 + mrSges(3,2) * t103 + t119) * g(1) + (-m(6) * (t134 + t93) - m(5) * (t93 + t164) - mrSges(3,1) * t103 + mrSges(3,2) * t102 + t120) * g(2) + (-t6 * mrSges(5,3) - t2 * mrSges(6,3) - t97 * t35) * t114 + t99 * t14 + t66 * t52 + t18 * t67 + t17 * t69 + t77 * t13 + t59 * t34 + t60 * t32 + m(5) * (t5 * t155 - t97 * t169 + t47 * t99) + Ifges(3,3) * qJDD(2) + t118 + t33 * t155; (-t68 * t145 - t70 * t146 - t114 * t35 + t116 * t33 + m(5) * (-t41 * t145 - t42 * t146 - t169 + t170)) * pkin(7) + ((t112 * mrSges(4,1) - t52 - t53) * t115 + (t112 * mrSges(4,2) - t179) * t117 - m(6) * (t115 * t36 + t21 * t148 - t16 * t149) - t121) * t159 + (t123 * qJD(4) - t169) * mrSges(5,3) + (t125 * qJD(4) - t2 * t114) * mrSges(6,3) + t119 * g(1) + (-m(5) * t164 - m(6) * t134 + t120) * g(2) + t52 * t139 - t98 * t13 + t51 * t67 + t50 * t69 + t78 * t34 + t80 * t32 - pkin(3) * t14 + m(6) * (-t12 * t98 + t36 * t139 + t16 * t51 + t2 * t78 + t21 * t50 + t3 * t80) + t118 - t47 * t178; t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t20 * t69 - t41 * t70 + t42 * t68 + t185 * t55 + t183 * t54 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t128 + t79) * g(3) + (t34 + (-g(3) * t116 + t2) * m(6)) * pkin(4) + (t67 - m(6) * (-t16 + t20)) * t21 + ((t41 * mrSges(5,3) + t16 * mrSges(6,3) - t72 * mrSges(5,2) - t36 * mrSges(6,2) - t45 / 0.2e1 - t46 / 0.2e1 - t81 / 0.2e1 - t82 / 0.2e1 + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4)) * t116 + (t42 * mrSges(5,3) + t21 * mrSges(6,3) - t72 * mrSges(5,1) - t36 * mrSges(6,1) + t43 / 0.2e1 + t44 / 0.2e1 + t189 * t151 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t36 - t52) * pkin(4) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t150) * t114) * t112 + (g(1) * t95 + g(2) * t94) * (t167 * t116 + (m(6) * pkin(4) + t168) * t114); (t114 * t67 - t116 * t69) * t112 + (-g(1) * t94 + g(2) * t95 - t124 * t112 + t12) * m(6) + t13;];
tau = t1;
