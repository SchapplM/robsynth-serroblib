% Calculate vector of inverse dynamics joint torques for
% S5PRPPR3
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:18
% EndTime: 2019-12-05 15:26:25
% DurationCPUTime: 2.16s
% Computational Cost: add. (850->202), mult. (1611->258), div. (0->0), fcn. (959->10), ass. (0->103)
t132 = -m(6) - m(5);
t138 = m(4) - t132;
t102 = -mrSges(4,1) + mrSges(5,2);
t55 = sin(qJ(5));
t57 = cos(qJ(5));
t75 = mrSges(6,1) * t55 + mrSges(6,2) * t57;
t137 = -t75 + mrSges(4,2);
t52 = sin(pkin(7));
t54 = cos(pkin(7));
t136 = g(1) * t54 + g(2) * t52;
t51 = sin(pkin(8));
t56 = sin(qJ(2));
t98 = qJD(1) * t56;
t38 = t51 * t98;
t53 = cos(pkin(8));
t58 = cos(qJ(2));
t97 = qJD(1) * t58;
t20 = t53 * t97 - t38;
t135 = -t20 + qJD(4);
t134 = -t55 / 0.2e1;
t123 = t57 / 0.2e1;
t50 = qJ(2) + pkin(8);
t47 = cos(t50);
t118 = g(3) * t47;
t99 = mrSges(6,3) * qJD(2);
t35 = -qJD(5) * mrSges(6,2) - t55 * t99;
t36 = qJD(5) * mrSges(6,1) - t57 * t99;
t71 = t55 * t35 + t57 * t36;
t88 = qJD(2) * qJD(5);
t29 = qJDD(2) * t57 - t55 * t88;
t16 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t29;
t30 = -qJDD(2) * t55 - t57 * t88;
t17 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t30;
t129 = t57 * t16 + t55 * t17;
t128 = -t51 * t56 + t53 * t58;
t122 = -pkin(3) - pkin(6);
t90 = qJD(1) * qJD(2);
t86 = t56 * t90;
t31 = t58 * qJDD(1) - t86;
t94 = qJDD(2) * pkin(2);
t25 = t31 + t94;
t85 = t58 * t90;
t32 = qJDD(1) * t56 + t85;
t8 = t25 * t53 - t51 * t32;
t79 = qJDD(4) - t8;
t3 = qJDD(2) * t122 + t79;
t37 = qJD(2) * pkin(2) + t97;
t14 = t37 * t53 - t38;
t81 = qJD(4) - t14;
t10 = qJD(2) * t122 + t81;
t6 = -qJD(3) * t55 + t10 * t57;
t1 = qJD(5) * t6 + qJDD(3) * t57 + t3 * t55;
t7 = qJD(3) * t57 + t10 * t55;
t2 = -qJD(5) * t7 - qJDD(3) * t55 + t3 * t57;
t127 = t1 * t55 + t2 * t57;
t105 = t57 * Ifges(6,4);
t115 = Ifges(6,4) * t55;
t126 = (-Ifges(6,2) * t57 - t115) * t134 + (-Ifges(6,1) * t55 - t105) * t123;
t15 = t37 * t51 + t53 * t98;
t13 = qJD(2) * qJ(4) + t15;
t76 = mrSges(6,1) * t57 - mrSges(6,2) * t55;
t125 = t13 * t76 + qJD(5) * (-Ifges(6,5) * t55 - Ifges(6,6) * t57) / 0.2e1;
t89 = qJD(2) * qJD(4);
t9 = t25 * t51 + t32 * t53;
t4 = qJDD(2) * qJ(4) + t89 + t9;
t42 = pkin(2) * t51 + qJ(4);
t124 = -t136 * qJ(4) * t47 + t135 * t13 + t4 * t42;
t49 = t58 * pkin(2);
t112 = t52 * t55;
t111 = t52 * t57;
t109 = t54 * t55;
t108 = t54 * t57;
t101 = -mrSges(4,2) + mrSges(5,3);
t96 = qJD(5) * t55;
t95 = qJD(5) * t57;
t27 = t51 * t58 + t53 * t56;
t18 = t27 * qJD(1);
t93 = t18 * qJD(2);
t92 = t20 * qJD(2);
t91 = qJDD(2) * mrSges(5,2);
t45 = -pkin(2) * t53 - pkin(3);
t84 = -g(1) * t52 + g(2) * t54;
t83 = t55 * t7 + t57 * t6;
t82 = -t6 * t55 + t7 * t57;
t21 = t128 * qJD(2);
t80 = t13 * t21 + t27 * t4;
t78 = t58 * mrSges(3,1) - t56 * mrSges(3,2);
t77 = mrSges(3,1) * t56 + mrSges(3,2) * t58;
t74 = Ifges(6,1) * t57 - t115;
t73 = -t55 * Ifges(6,2) + t105;
t63 = qJD(5) * t82 + t127;
t46 = sin(t50);
t62 = -t13 * qJD(2) - t136 * t46 + t118;
t61 = (t57 * t35 - t55 * t36) * qJD(5) + t129;
t59 = qJD(2) ^ 2;
t28 = t75 * qJD(2);
t23 = Ifges(6,5) * qJD(5) + qJD(2) * t74;
t22 = Ifges(6,6) * qJD(5) + qJD(2) * t73;
t19 = t27 * qJD(2);
t12 = -qJD(2) * pkin(3) + t81;
t11 = -mrSges(6,1) * t30 + mrSges(6,2) * t29;
t5 = -qJDD(2) * pkin(3) + t79;
t24 = [m(2) * qJDD(1) + t27 * t11 + t21 * t28 - t77 * t59 + t71 * t19 - t61 * t128 + (t101 * t21 + t102 * t19) * qJD(2) + (-m(2) - m(3) - t138) * g(3) + m(6) * (-t128 * t63 + t19 * t83 + t80) + m(3) * (t31 * t58 + t32 * t56) + m(4) * (t128 * t8 - t14 * t19 + t15 * t21 + t27 * t9) + m(5) * (t12 * t19 - t128 * t5 + t80) + (t101 * t27 - t102 * t128 + t78) * qJDD(2); (t5 - t93) * mrSges(5,2) + (-t51 * t94 - t9 + t92) * mrSges(4,2) + (t53 * t94 + t8 + t93) * mrSges(4,1) + (t14 * t18 - t15 * t20 + (t51 * t9 + t53 * t8) * pkin(2)) * m(4) + (0.2e1 * Ifges(6,5) * t123 - Ifges(6,6) * t55) * qJDD(5) + (t6 * t96 - t7 * t95 - t118 - t127) * mrSges(6,3) + (t63 * m(6) + t35 * t95 - t36 * t96 + t129) * (-pkin(6) + t45) - t71 * t18 + t125 * qJD(5) + t126 * t88 + (-t12 * t18 + t45 * t5 + t124) * m(5) + (-t83 * t18 + t124) * m(6) + (Ifges(3,3) + Ifges(5,1) + Ifges(4,3)) * qJDD(2) + (Ifges(6,1) * t29 + Ifges(6,4) * t30) * t123 + (-g(3) * t46 + qJDD(2) * t42 + t4 + t89 - t92) * mrSges(5,3) + (t85 - t32) * mrSges(3,2) + (t86 + t31) * mrSges(3,1) + t42 * t11 + (-m(4) * t49 - t78 + t132 * (t47 * pkin(3) + t46 * qJ(4) + t49) + (-m(6) * pkin(6) + t102) * t47 + t137 * t46) * g(3) + t136 * (t77 + (-mrSges(5,3) + t137) * t47 + (m(5) * pkin(3) - m(6) * t122 + mrSges(6,3) - t102) * t46 + t138 * pkin(2) * t56) + t30 * t73 / 0.2e1 + t29 * t74 / 0.2e1 + t4 * t75 - t22 * t95 / 0.2e1 - t23 * t96 / 0.2e1 + (Ifges(6,4) * t29 + Ifges(6,2) * t30) * t134 + t45 * t91 + t135 * t28; -t55 * t16 + t57 * t17 - t71 * qJD(5) + (-qJD(5) * t83 + t1 * t57 - t2 * t55 + t84) * m(6) + (m(4) + m(5)) * (qJDD(3) + t84); t91 - t59 * mrSges(5,3) - qJD(2) * t28 + (t62 + t63) * m(6) + (t5 + t62) * m(5) + t61; Ifges(6,5) * t29 + Ifges(6,6) * t30 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t6 * t35 + t7 * t36 - g(1) * ((t108 * t46 - t112) * mrSges(6,1) + (-t109 * t46 - t111) * mrSges(6,2)) - g(2) * ((t111 * t46 + t109) * mrSges(6,1) + (-t112 * t46 + t108) * mrSges(6,2)) + t76 * t118 + (t55 * t23 / 0.2e1 + t22 * t123 - t126 * qJD(2) + t82 * mrSges(6,3) - t125) * qJD(2);];
tau = t24;
