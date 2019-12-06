% Calculate vector of inverse dynamics joint torques for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:00
% EndTime: 2019-12-05 15:03:04
% DurationCPUTime: 1.59s
% Computational Cost: add. (723->169), mult. (1488->225), div. (0->0), fcn. (944->10), ass. (0->84)
t96 = -mrSges(4,1) + mrSges(5,2);
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t70 = t51 * mrSges(6,1) + t53 * mrSges(6,2);
t131 = -t70 + mrSges(4,2);
t48 = sin(pkin(7));
t50 = cos(pkin(7));
t122 = g(1) * t50 + g(2) * t48;
t47 = sin(pkin(8));
t49 = cos(pkin(8));
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t22 = t47 * t52 - t54 * t49;
t16 = t22 * qJD(1);
t65 = t16 + qJD(4);
t130 = -t51 / 0.2e1;
t117 = t53 / 0.2e1;
t82 = qJD(1) * qJD(3);
t129 = qJDD(1) * t52 + t54 * t82;
t128 = qJDD(1) * t54 - t52 * t82;
t80 = qJD(3) * qJD(5);
t25 = qJDD(3) * t53 - t51 * t80;
t14 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t25;
t26 = -qJDD(3) * t51 - t53 * t80;
t15 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t26;
t92 = mrSges(6,3) * qJD(3);
t29 = -qJD(5) * mrSges(6,2) - t51 * t92;
t30 = qJD(5) * mrSges(6,1) - t53 * t92;
t57 = (t53 * t29 - t51 * t30) * qJD(5) + t53 * t14 + t51 * t15;
t44 = pkin(8) + qJ(3);
t43 = cos(t44);
t114 = g(3) * t43;
t126 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4);
t125 = t51 * t29 + t53 * t30;
t55 = -pkin(3) - pkin(6);
t7 = t128 * t49 - t129 * t47;
t60 = qJDD(4) - t7;
t3 = qJDD(3) * t55 + t60;
t11 = qJD(3) * t55 + t65;
t8 = -qJD(2) * t51 + t11 * t53;
t1 = qJD(5) * t8 + qJDD(2) * t53 + t3 * t51;
t9 = qJD(2) * t53 + t11 * t51;
t2 = -qJD(5) * t9 - qJDD(2) * t51 + t3 * t53;
t123 = t1 * t51 + t2 * t53;
t110 = Ifges(6,4) * t53;
t111 = Ifges(6,4) * t51;
t121 = (-Ifges(6,2) * t53 - t111) * t130 + (-Ifges(6,1) * t51 - t110) * t117;
t23 = t47 * t54 + t49 * t52;
t17 = t23 * qJD(1);
t13 = qJD(3) * qJ(4) + t17;
t71 = mrSges(6,1) * t53 - mrSges(6,2) * t51;
t119 = t13 * t71 + qJD(5) * (-Ifges(6,5) * t51 - Ifges(6,6) * t53) / 0.2e1;
t6 = t128 * t47 + t129 * t49;
t4 = -t6 + t126;
t118 = t65 * t13 + (-t122 * t43 - t4) * qJ(4);
t107 = t48 * t51;
t106 = t48 * t53;
t105 = t50 * t51;
t104 = t50 * t53;
t95 = mrSges(4,2) - mrSges(5,3);
t91 = qJD(5) * t51;
t90 = qJD(5) * t53;
t88 = t16 * qJD(3);
t87 = t17 * qJD(3);
t84 = qJDD(3) * mrSges(5,2);
t83 = m(3) + m(4) + m(5);
t78 = m(6) + t83;
t75 = t51 * t9 + t53 * t8;
t74 = -t8 * t51 + t9 * t53;
t18 = t22 * qJD(3);
t73 = -t13 * t18 - t23 * t4;
t69 = Ifges(6,1) * t53 - t111;
t68 = -Ifges(6,2) * t51 + t110;
t59 = qJD(5) * t74 + t123;
t42 = sin(t44);
t58 = -t13 * qJD(3) - t122 * t42 + t114;
t24 = t70 * qJD(3);
t21 = Ifges(6,5) * qJD(5) + qJD(3) * t69;
t20 = Ifges(6,6) * qJD(5) + qJD(3) * t68;
t19 = t23 * qJD(3);
t12 = -qJD(3) * pkin(3) + t65;
t10 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t5 = -qJDD(3) * pkin(3) + t60;
t27 = [-t18 * t24 + t125 * t19 + (-qJDD(3) * t95 + t10) * t23 + (t18 * t95 + t19 * t96) * qJD(3) + (-m(2) - t78) * g(3) + (qJDD(3) * t96 + t57) * t22 + m(6) * (t19 * t75 + t22 * t59 + t73) + m(4) * (t16 * t19 - t17 * t18 - t22 * t7 + t23 * t6) + m(5) * (t12 * t19 + t22 * t5 + t73) + (m(2) + m(3) * (t47 ^ 2 + t49 ^ 2)) * qJDD(1); -t29 * t91 - t30 * t90 + t53 * t15 - t51 * t14 + m(6) * (-qJD(5) * t75 + t1 * t53 - t2 * t51) + t83 * qJDD(2) + (-g(1) * t48 + g(2) * t50) * t78; (-t88 - t6) * mrSges(4,2) + (t87 + t7) * mrSges(4,1) + t26 * t68 / 0.2e1 + t25 * t69 / 0.2e1 - t4 * t70 + (-t17 * t75 + t118) * m(6) + ((-m(6) - m(5)) * (t43 * pkin(3) + t42 * qJ(4)) + (-m(6) * pkin(6) + t96) * t43 + t131 * t42) * g(3) + ((-mrSges(5,3) + t131) * t43 + (m(5) * pkin(3) - m(6) * t55 + mrSges(6,3) - t96) * t42) * t122 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + (m(6) * t59 + t57) * t55 + (-t87 + t5) * mrSges(5,2) + (Ifges(6,1) * t25 + Ifges(6,4) * t26) * t117 - t20 * t90 / 0.2e1 - t21 * t91 / 0.2e1 - pkin(3) * t84 + qJ(4) * t10 + (t8 * t91 - t9 * t90 - t114 - t123) * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t117 - Ifges(6,6) * t51) * qJDD(5) + (-g(3) * t42 - t126 - t4 + t88) * mrSges(5,3) + t65 * t24 - t125 * t17 + t119 * qJD(5) + t121 * t80 + (-pkin(3) * t5 - t12 * t17 + t118) * m(5) + (Ifges(6,4) * t25 + Ifges(6,2) * t26) * t130; t84 + (-qJD(3) * mrSges(5,3) - t24) * qJD(3) + (t58 + t59) * m(6) + (t5 + t58) * m(5) + t57; Ifges(6,5) * t25 + Ifges(6,6) * t26 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t8 * t29 + t9 * t30 - g(1) * ((t104 * t42 - t107) * mrSges(6,1) + (-t105 * t42 - t106) * mrSges(6,2)) - g(2) * ((t106 * t42 + t105) * mrSges(6,1) + (-t107 * t42 + t104) * mrSges(6,2)) + t71 * t114 + (t51 * t21 / 0.2e1 + t20 * t117 - t121 * qJD(3) + t74 * mrSges(6,3) - t119) * qJD(3);];
tau = t27;
