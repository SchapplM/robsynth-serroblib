% Calculate vector of inverse dynamics joint torques for
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:47
% EndTime: 2019-12-05 15:04:53
% DurationCPUTime: 2.05s
% Computational Cost: add. (1007->229), mult. (2220->338), div. (0->0), fcn. (1638->12), ass. (0->116)
t72 = cos(pkin(7));
t74 = sin(qJ(3));
t115 = t72 * t74;
t69 = sin(pkin(7));
t76 = cos(qJ(3));
t120 = t69 * t76;
t71 = cos(pkin(8));
t148 = -t71 * t115 + t120;
t68 = sin(pkin(8));
t109 = qJD(1) * t68;
t45 = t76 * qJD(2) - t74 * t109;
t147 = qJD(3) * t45;
t73 = sin(qJ(5));
t136 = t73 / 0.2e1;
t75 = cos(qJ(5));
t89 = -mrSges(6,1) * t75 + mrSges(6,2) * t73;
t80 = m(6) * pkin(4) - t89;
t146 = -t80 - mrSges(5,1);
t133 = g(3) * t68;
t66 = qJ(3) + pkin(9);
t63 = cos(t66);
t145 = t63 * t133;
t47 = t89 * qJD(3);
t144 = -qJD(3) * mrSges(5,1) + t47;
t102 = t73 * qJD(3);
t51 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t102;
t101 = t75 * qJD(3);
t52 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t101;
t86 = -t73 * t51 + t75 * t52;
t143 = -qJD(3) * mrSges(5,2) + t86;
t114 = t72 * t76;
t121 = t69 * t74;
t142 = -t71 * t121 - t114;
t46 = qJD(2) * t74 + t109 * t76;
t141 = qJD(3) * t46;
t99 = qJD(3) * qJD(5);
t48 = qJDD(3) * t75 - t73 * t99;
t34 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t48;
t49 = qJDD(3) * t73 + t75 * t99;
t35 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t49;
t140 = t75 * t34 - t73 * t35;
t67 = sin(pkin(9));
t70 = cos(pkin(9));
t43 = t67 * t74 - t70 * t76;
t103 = qJDD(3) * pkin(3);
t100 = qJDD(1) * t68;
t23 = t76 * qJDD(2) - t100 * t74 - t141;
t19 = t23 + t103;
t22 = qJDD(2) * t74 + t76 * t100 + t147;
t6 = t67 * t19 + t70 * t22;
t4 = qJDD(3) * pkin(6) + t6;
t54 = -qJDD(1) * t71 + qJDD(4);
t36 = t70 * t46;
t38 = qJD(3) * pkin(3) + t45;
t14 = t67 * t38 + t36;
t12 = qJD(3) * pkin(6) + t14;
t55 = -qJD(1) * t71 + qJD(4);
t9 = -t12 * t73 + t55 * t75;
t1 = qJD(5) * t9 + t4 * t75 + t54 * t73;
t10 = t12 * t75 + t55 * t73;
t2 = -qJD(5) * t10 - t4 * t73 + t54 * t75;
t139 = t1 * t75 - t2 * t73;
t127 = t46 * t67;
t13 = t38 * t70 - t127;
t11 = -qJD(3) * pkin(4) - t13;
t138 = t11 * (mrSges(6,1) * t73 + mrSges(6,2) * t75) + qJD(5) * (Ifges(6,5) * t75 - Ifges(6,6) * t73) / 0.2e1;
t137 = m(5) * pkin(3);
t134 = pkin(3) * t74;
t130 = Ifges(6,4) * t73;
t129 = Ifges(6,4) * t75;
t128 = Ifges(6,2) * t75;
t126 = t63 * t72;
t124 = t68 * t73;
t123 = t68 * t75;
t122 = t69 * t71;
t118 = t71 * t73;
t117 = t71 * t75;
t62 = sin(t66);
t116 = t72 * t62;
t106 = qJD(3) * t68;
t105 = qJD(5) * t73;
t104 = qJD(5) * t75;
t95 = m(3) + m(4) + m(5) + m(6);
t93 = t10 * t75 - t73 * t9;
t92 = t10 * t73 + t9 * t75;
t91 = -t76 * mrSges(4,1) + t74 * mrSges(4,2);
t90 = -mrSges(4,1) * t74 - mrSges(4,2) * t76;
t88 = t128 + t130;
t5 = t19 * t70 - t22 * t67;
t33 = t43 * t68;
t21 = -t33 * t75 - t118;
t20 = t33 * t73 - t117;
t44 = t67 * t76 + t70 * t74;
t85 = t148 * pkin(3);
t82 = t73 * (Ifges(6,1) * t75 - t130);
t81 = t90 * t68;
t79 = g(3) * t71 + (-g(1) * t72 - g(2) * t69) * t68;
t78 = -qJD(5) * t92 + t139;
t77 = qJD(3) ^ 2;
t61 = t71 ^ 2 * qJDD(1);
t60 = Ifges(6,4) * t101;
t59 = -pkin(3) * t70 - pkin(4);
t42 = Ifges(6,1) * t102 + Ifges(6,5) * qJD(5) + t60;
t41 = Ifges(6,6) * qJD(5) + qJD(3) * t88;
t40 = t43 * qJD(3);
t39 = t44 * qJD(3);
t32 = t44 * t68;
t30 = t126 * t71 + t62 * t69;
t28 = t122 * t63 - t116;
t26 = t44 * t106;
t25 = t43 * t106;
t24 = -mrSges(6,1) * t48 + mrSges(6,2) * t49;
t8 = -qJD(5) * t21 + t26 * t73;
t7 = qJD(5) * t20 - t26 * t75;
t3 = -qJDD(3) * pkin(4) - t5;
t15 = [m(2) * qJDD(1) + t20 * t35 + t21 * t34 + t32 * t24 - t25 * t47 + t8 * t51 + t7 * t52 + t91 * t77 * t68 + (mrSges(5,1) * t25 + mrSges(5,2) * t26) * qJD(3) + (-mrSges(5,1) * t32 + mrSges(5,2) * t33 + t81) * qJDD(3) + (-m(2) - t95) * g(3) + m(3) * (qJDD(1) * t68 ^ 2 + t61) + m(4) * (t61 + (t22 * t76 - t23 * t74 + (-t45 * t76 - t46 * t74) * qJD(3)) * t68) + m(5) * (t13 * t25 - t14 * t26 - t32 * t5 - t33 * t6 - t54 * t71) + m(6) * (t1 * t21 + t10 * t7 - t11 * t25 + t2 * t20 + t3 * t32 + t8 * t9); m(3) * qJDD(2) + t43 * t24 + t90 * t77 + t144 * t39 - t143 * t40 + ((-t75 * t51 - t73 * t52) * qJD(5) + t140) * t44 + m(6) * (t11 * t39 + t3 * t43 - t40 * t93 + t44 * t78) + m(4) * (t22 * t74 + t23 * t76 + (-t45 * t74 + t46 * t76) * qJD(3)) + m(5) * (-t13 * t39 - t14 * t40 - t43 * t5 + t44 * t6) + (-mrSges(5,1) * t43 - mrSges(5,2) * t44 - t91) * qJDD(3) + (-g(1) * t69 + g(2) * t72) * t95; (0.2e1 * Ifges(6,5) * t136 + Ifges(6,6) * t75) * qJDD(5) + (-g(1) * t30 - g(2) * t28 - t10 * t105 - t104 * t9 + t139 - t145) * mrSges(6,3) + (-t67 * t103 + t145 - t6) * mrSges(5,2) + (-m(6) * (t142 * pkin(3) + pkin(6) * t28) - (-t120 * t71 + t115) * mrSges(4,2) + t28 * mrSges(5,2) + (-mrSges(4,1) - t137) * t142 + t146 * (-t122 * t62 - t126)) * g(2) + (m(6) * t78 - t104 * t51 - t105 * t52 + t140) * (pkin(3) * t67 + pkin(6)) + (t23 + t141) * mrSges(4,1) + (-m(5) * t14 - m(6) * t93 - t143) * (t45 * t70 - t127) + (m(5) * t13 - m(6) * t11 - t144) * (t45 * t67 + t36) + t138 * qJD(5) + (Ifges(6,1) * t49 + Ifges(6,4) * t48) * t136 + (t75 * (-Ifges(6,2) * t73 + t129) + t82) * t99 / 0.2e1 + (m(6) * t59 + t89) * t3 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (-m(6) * (pkin(6) * t63 - t134) + t80 * t62 + m(5) * t134) * t133 + t75 * (Ifges(6,4) * t49 + Ifges(6,2) * t48) / 0.2e1 + t59 * t24 + (-t22 + t147) * mrSges(4,2) + (-m(6) * (pkin(6) * t30 + t85) - m(5) * t85 + t30 * mrSges(5,2) - t148 * mrSges(4,1) - (-t114 * t71 - t121) * mrSges(4,2) + t146 * (-t116 * t71 + t63 * t69)) * g(1) + (t70 * t103 + t62 * t133 + t5) * mrSges(5,1) + (t5 * t70 + t6 * t67) * t137 - g(3) * t81 + t48 * t88 / 0.2e1 + t42 * t104 / 0.2e1 - t41 * t105 / 0.2e1 + t49 * (Ifges(6,1) * t73 + t129) / 0.2e1; t73 * t34 + t75 * t35 + t86 * qJD(5) + (qJD(5) * t93 + t1 * t73 + t2 * t75 + t79) * m(6) + (t54 + t79) * m(5); Ifges(6,5) * t49 + Ifges(6,6) * t48 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t9 * t52 + t10 * t51 - g(1) * ((t123 * t72 - t30 * t73) * mrSges(6,1) + (-t124 * t72 - t30 * t75) * mrSges(6,2)) - g(2) * ((t123 * t69 - t28 * t73) * mrSges(6,1) + (-t124 * t69 - t28 * t75) * mrSges(6,2)) - g(3) * ((-t124 * t63 - t117) * mrSges(6,1) + (-t123 * t63 + t118) * mrSges(6,2)) + (t41 * t136 + (-t82 / 0.2e1 + t128 * t136) * qJD(3) + t92 * mrSges(6,3) - (t42 + t60) * t75 / 0.2e1 - t138) * qJD(3);];
tau = t15;
