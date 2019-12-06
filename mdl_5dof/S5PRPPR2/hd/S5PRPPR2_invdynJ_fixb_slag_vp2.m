% Calculate vector of inverse dynamics joint torques for
% S5PRPPR2
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:58
% EndTime: 2019-12-05 15:24:08
% DurationCPUTime: 3.29s
% Computational Cost: add. (1346->242), mult. (2881->324), div. (0->0), fcn. (2026->14), ass. (0->120)
t145 = m(5) + m(6);
t161 = m(4) + t145;
t162 = pkin(2) * t161;
t92 = sin(pkin(9));
t95 = cos(pkin(9));
t151 = mrSges(5,3) * (t92 ^ 2 + t95 ^ 2);
t102 = cos(qJ(2));
t128 = qJD(1) * t102;
t100 = sin(qJ(2));
t129 = qJD(1) * t100;
t93 = sin(pkin(8));
t72 = t93 * t129;
t96 = cos(pkin(8));
t51 = t128 * t96 - t72;
t160 = qJD(4) - t51;
t113 = -mrSges(5,1) * t95 + mrSges(5,2) * t92;
t76 = pkin(4) * t95 + pkin(3);
t90 = pkin(9) + qJ(5);
t82 = sin(t90);
t84 = cos(t90);
t159 = m(5) * pkin(3) + m(6) * t76 + t84 * mrSges(6,1) - t82 * mrSges(6,2) + mrSges(4,1) - t113;
t158 = -m(5) * qJ(4) + m(6) * (-pkin(6) - qJ(4)) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t94 = sin(pkin(7));
t97 = cos(pkin(7));
t155 = g(1) * t97 + g(2) * t94;
t101 = cos(qJ(5));
t73 = pkin(2) * t93 + qJ(4);
t141 = pkin(6) + t73;
t56 = t141 * t92;
t57 = t141 * t95;
t99 = sin(qJ(5));
t29 = t101 * t57 - t56 * t99;
t65 = t101 * t92 + t95 * t99;
t154 = -qJD(5) * t29 - t160 * t65;
t28 = -t101 * t56 - t57 * t99;
t63 = t101 * t95 - t92 * t99;
t153 = qJD(5) * t28 + t160 * t63;
t50 = t63 * qJD(2);
t53 = t65 * qJD(2);
t152 = mrSges(6,1) * t50 - mrSges(6,2) * t53 - t113 * qJD(2);
t91 = qJ(2) + pkin(8);
t83 = sin(t91);
t150 = t155 * t83;
t54 = t63 * qJD(5);
t124 = qJD(2) * qJD(4);
t130 = qJDD(2) * pkin(2);
t125 = qJD(1) * qJD(2);
t120 = t100 * t125;
t68 = t102 * qJDD(1) - t120;
t60 = t68 + t130;
t119 = t102 * t125;
t69 = qJDD(1) * t100 + t119;
t35 = t93 * t60 + t96 * t69;
t18 = qJDD(2) * qJ(4) + t124 + t35;
t79 = t95 * qJDD(3);
t12 = -t18 * t92 + t79;
t13 = t92 * qJDD(3) + t95 * t18;
t112 = -t12 * t92 + t13 * t95;
t146 = t53 / 0.2e1;
t144 = pkin(2) * t96;
t126 = qJDD(2) * t95;
t127 = qJDD(2) * t92;
t61 = -mrSges(5,1) * t126 + mrSges(5,2) * t127;
t32 = qJD(2) * t54 + qJDD(2) * t65;
t55 = t65 * qJD(5);
t33 = -qJD(2) * t55 + qJDD(2) * t63;
t7 = -t33 * mrSges(6,1) + t32 * mrSges(6,2);
t140 = t61 + t7;
t139 = Ifges(6,4) * t53;
t85 = cos(t91);
t135 = t85 * t94;
t134 = t85 * t97;
t71 = qJD(2) * pkin(2) + t128;
t45 = t129 * t96 + t93 * t71;
t41 = qJD(2) * qJ(4) + t45;
t37 = t92 * qJD(3) + t95 * t41;
t132 = pkin(6) * qJD(2);
t131 = qJD(2) * t51;
t34 = t60 * t96 - t93 * t69;
t44 = t71 * t96 - t72;
t117 = -g(1) * t94 + g(2) * t97;
t116 = qJD(4) - t44;
t115 = qJDD(4) - t34;
t114 = -mrSges(4,2) + t151;
t81 = t95 * qJD(3);
t111 = -(-t41 * t92 + t81) * t92 + t37 * t95;
t26 = t81 + (-t41 - t132) * t92;
t27 = t132 * t95 + t37;
t6 = t101 * t27 + t26 * t99;
t5 = t101 * t26 - t27 * t99;
t110 = t102 * mrSges(3,1) - t100 * mrSges(3,2);
t109 = mrSges(3,1) * t100 + mrSges(3,2) * t102;
t64 = t100 * t96 + t102 * t93;
t62 = t100 * t93 - t96 * t102;
t103 = qJD(2) ^ 2;
t77 = -pkin(3) - t144;
t70 = -t76 - t144;
t52 = t62 * qJD(2);
t49 = t64 * qJD(2);
t48 = t64 * qJD(1);
t47 = Ifges(6,4) * t50;
t43 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t53;
t42 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t50;
t40 = -qJD(2) * pkin(3) + t116;
t38 = -qJD(2) * t76 + t116;
t25 = -qJDD(2) * pkin(3) + t115;
t24 = Ifges(6,1) * t53 + Ifges(6,5) * qJD(5) + t47;
t23 = Ifges(6,2) * t50 + Ifges(6,6) * qJD(5) + t139;
t22 = t63 * t64;
t21 = t65 * t64;
t20 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t33;
t19 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t32;
t14 = -qJDD(2) * t76 + t115;
t9 = pkin(6) * t126 + t13;
t8 = t79 + (-pkin(6) * qJDD(2) - t18) * t92;
t4 = t52 * t65 - t64 * t54;
t3 = -t52 * t63 - t55 * t64;
t2 = -qJD(5) * t6 + t101 * t8 - t9 * t99;
t1 = qJD(5) * t5 + t101 * t9 + t8 * t99;
t10 = [m(2) * qJDD(1) - t21 * t19 + t22 * t20 + t3 * t42 + t4 * t43 + t140 * t62 - t152 * t49 - t109 * t103 + (-t49 * mrSges(4,1) - t114 * t52) * qJD(2) + (-m(2) - m(3) - t161) * g(3) + m(3) * (t100 * t69 + t102 * t68) + m(5) * (-t111 * t52 + t112 * t64 + t25 * t62 + t40 * t49) + m(6) * (t1 * t22 + t14 * t62 - t2 * t21 + t3 * t6 + t38 * t49 + t4 * t5) + m(4) * (-t34 * t62 + t35 * t64 - t44 * t49 - t45 * t52) + (-t62 * mrSges(4,1) + t114 * t64 + t110) * qJDD(2); (t68 + t120) * mrSges(3,1) + (qJDD(2) * t73 + t124 - t131) * t151 + (mrSges(6,2) * t14 - mrSges(6,3) * t2 + Ifges(6,1) * t32 + Ifges(6,4) * t33 + Ifges(6,5) * qJDD(5)) * t65 + (-mrSges(6,1) * t14 + mrSges(6,3) * t1 + Ifges(6,4) * t32 + Ifges(6,2) * t33 + Ifges(6,6) * qJDD(5)) * t63 + (-t69 + t119) * mrSges(3,2) + ((t34 * t96 + t35 * t93) * pkin(2) + t44 * t48 - t45 * t51) * m(4) + t152 * t48 + t153 * t42 + (t1 * t29 + t14 * t70 + t153 * t6 + t154 * t5 + t2 * t28 - t38 * t48) * m(6) + t154 * t43 + (Ifges(5,4) * t92 + Ifges(5,2) * t95) * t126 + (Ifges(5,1) * t92 + Ifges(5,4) * t95) * t127 + t112 * mrSges(5,3) + (Ifges(6,1) * t54 - Ifges(6,4) * t55) * t146 + t50 * (Ifges(6,4) * t54 - Ifges(6,2) * t55) / 0.2e1 + t38 * (mrSges(6,1) * t55 + mrSges(6,2) * t54) + qJD(5) * (Ifges(6,5) * t54 - Ifges(6,6) * t55) / 0.2e1 + (-t5 * t54 - t6 * t55) * mrSges(6,3) + t25 * t113 + (qJD(2) * t48 + t130 * t96 + t34) * mrSges(4,1) + t77 * t61 + t70 * t7 + t54 * t24 / 0.2e1 - t55 * t23 / 0.2e1 + t28 * t19 + t29 * t20 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (-t130 * t93 + t131 - t35) * mrSges(4,2) + (-t102 * t162 + t158 * t83 - t159 * t85 - t110) * g(3) + t155 * (t100 * t162 + t158 * t85 + t159 * t83 + t109) + (t160 * t111 + t112 * t73 + t25 * t77 - t40 * t48) * m(5); t63 * t19 + t65 * t20 + t54 * t42 - t55 * t43 + (t1 * t65 + t2 * t63 - t5 * t55 + t54 * t6 + t117) * m(6) + (t12 * t95 + t13 * t92 + t117) * m(5) + (qJDD(3) + t117) * m(4); -t103 * t151 + t145 * t85 * g(3) - t50 * t42 + t53 * t43 + t140 + (t5 * t53 - t50 * t6 + t14 - t150) * m(6) + (-qJD(2) * t111 - t150 + t25) * m(5); Ifges(6,5) * t32 + Ifges(6,6) * t33 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t38 * (mrSges(6,1) * t53 + mrSges(6,2) * t50) - t53 * (Ifges(6,1) * t50 - t139) / 0.2e1 + t23 * t146 - qJD(5) * (Ifges(6,5) * t50 - Ifges(6,6) * t53) / 0.2e1 - t5 * t42 + t6 * t43 - g(1) * ((-t134 * t82 + t84 * t94) * mrSges(6,1) + (-t134 * t84 - t82 * t94) * mrSges(6,2)) - g(2) * ((-t135 * t82 - t84 * t97) * mrSges(6,1) + (-t135 * t84 + t82 * t97) * mrSges(6,2)) - g(3) * (-mrSges(6,1) * t82 - mrSges(6,2) * t84) * t83 + (t5 * t50 + t53 * t6) * mrSges(6,3) - (-Ifges(6,2) * t53 + t24 + t47) * t50 / 0.2e1;];
tau = t10;
