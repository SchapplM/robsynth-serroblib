% Calculate vector of inverse dynamics joint torques for
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:50
% EndTime: 2019-12-05 14:57:56
% DurationCPUTime: 1.21s
% Computational Cost: add. (853->186), mult. (1937->289), div. (0->0), fcn. (1527->12), ass. (0->93)
t61 = sin(qJ(5));
t109 = t61 / 0.2e1;
t63 = cos(qJ(5));
t73 = -t63 * mrSges(6,1) + t61 * mrSges(6,2);
t115 = m(6) * pkin(4) + mrSges(5,1) - t73;
t114 = -m(6) * pkin(6) + mrSges(5,2) - mrSges(6,3);
t87 = qJD(4) * t61;
t45 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t87;
t85 = qJD(4) * t63;
t46 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t85;
t69 = -t61 * t45 + t63 * t46;
t113 = -qJD(4) * mrSges(5,2) + t69;
t79 = qJD(4) * qJD(5);
t42 = qJDD(4) * t63 - t61 * t79;
t27 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t42;
t43 = qJDD(4) * t61 + t63 * t79;
t28 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t43;
t112 = t63 * t27 - t61 * t28;
t55 = sin(pkin(9));
t58 = cos(pkin(9));
t62 = sin(qJ(4));
t64 = cos(qJ(4));
t40 = t55 * t64 + t62 * t58;
t32 = t40 * qJD(4);
t56 = sin(pkin(8));
t81 = qJDD(1) * t56;
t35 = qJDD(2) * t58 - t55 * t81;
t36 = qJDD(2) * t55 + t58 * t81;
t89 = qJD(1) * t56;
t37 = qJD(2) * t58 - t55 * t89;
t38 = qJD(2) * t55 + t58 * t89;
t84 = qJD(4) * t64;
t86 = qJD(4) * t62;
t5 = t62 * t35 + t64 * t36 + t37 * t84 - t38 * t86;
t3 = qJDD(4) * pkin(6) + t5;
t59 = cos(pkin(8));
t48 = -qJDD(1) * t59 + qJDD(3);
t14 = t62 * t37 + t38 * t64;
t12 = qJD(4) * pkin(6) + t14;
t49 = -qJD(1) * t59 + qJD(3);
t9 = -t12 * t61 + t49 * t63;
t1 = t9 * qJD(5) + t3 * t63 + t48 * t61;
t10 = t12 * t63 + t49 * t61;
t2 = -t10 * qJD(5) - t3 * t61 + t48 * t63;
t111 = t1 * t63 - t2 * t61;
t6 = -t14 * qJD(4) + t35 * t64 - t62 * t36;
t13 = t37 * t64 - t62 * t38;
t11 = -qJD(4) * pkin(4) - t13;
t110 = t11 * (mrSges(6,1) * t61 + mrSges(6,2) * t63) + qJD(5) * (Ifges(6,5) * t63 - Ifges(6,6) * t61) / 0.2e1;
t104 = Ifges(6,4) * t61;
t103 = Ifges(6,4) * t63;
t102 = Ifges(6,2) * t63;
t101 = t48 * t59;
t54 = pkin(9) + qJ(4);
t53 = cos(t54);
t60 = cos(pkin(7));
t100 = t53 * t60;
t99 = t56 * t61;
t98 = t56 * t63;
t57 = sin(pkin(7));
t97 = t57 * t59;
t96 = t59 * t61;
t95 = t59 * t63;
t52 = sin(t54);
t94 = t60 * t52;
t83 = qJD(5) * t61;
t82 = qJD(5) * t63;
t80 = m(4) + m(5) + m(6);
t78 = m(3) + t80;
t75 = t10 * t63 - t61 * t9;
t74 = t10 * t61 + t9 * t63;
t72 = t102 + t104;
t39 = t55 * t62 - t64 * t58;
t25 = t39 * t56;
t16 = -t25 * t63 - t96;
t15 = t25 * t61 - t95;
t67 = t61 * (Ifges(6,1) * t63 - t104);
t65 = -t74 * qJD(5) + t111;
t51 = Ifges(6,4) * t85;
t41 = t73 * qJD(4);
t34 = Ifges(6,1) * t87 + Ifges(6,5) * qJD(5) + t51;
t33 = Ifges(6,6) * qJD(5) + t72 * qJD(4);
t31 = t39 * qJD(4);
t24 = t40 * t56;
t23 = t59 * t100 + t52 * t57;
t21 = t53 * t97 - t94;
t19 = (-t55 * t86 + t58 * t84) * t56;
t18 = t56 * t32;
t17 = -mrSges(6,1) * t42 + mrSges(6,2) * t43;
t8 = -t16 * qJD(5) + t18 * t61;
t7 = t15 * qJD(5) - t18 * t63;
t4 = -qJDD(4) * pkin(4) - t6;
t20 = [t15 * t28 + t16 * t27 + t24 * t17 + t19 * t41 + t8 * t45 + t7 * t46 + (qJD(4) * t18 + qJDD(4) * t25) * mrSges(5,2) + (-qJD(4) * t19 - qJDD(4) * t24) * mrSges(5,1) + (-m(2) - t78) * g(3) + m(4) * (-t101 + (-t35 * t55 + t36 * t58) * t56) + m(5) * (-t13 * t19 - t14 * t18 - t24 * t6 - t25 * t5 - t101) + m(6) * (t1 * t16 + t10 * t7 + t11 * t19 + t15 * t2 + t24 * t4 + t8 * t9) + (m(2) + m(3) * (t56 ^ 2 + t59 ^ 2)) * qJDD(1); m(3) * qJDD(2) + t39 * t17 + t32 * t41 + (-qJD(4) * t32 - qJDD(4) * t39) * mrSges(5,1) - t113 * t31 + (-qJDD(4) * mrSges(5,2) + (-t63 * t45 - t61 * t46) * qJD(5) + t112) * t40 + m(4) * (t35 * t58 + t36 * t55) + m(5) * (-t13 * t32 - t14 * t31 - t39 * t6 + t40 * t5) + m(6) * (t11 * t32 - t75 * t31 + t39 * t4 + t65 * t40) + (-g(1) * t57 + g(2) * t60) * t78; t61 * t27 + t63 * t28 + t69 * qJD(5) + m(6) * (t75 * qJD(5) + t1 * t61 + t2 * t63) + 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1) * t48 + (g(3) * t59 + (-g(1) * t60 - g(2) * t57) * t56) * t80; (Ifges(6,1) * t43 + Ifges(6,4) * t42) * t109 + t43 * (Ifges(6,1) * t61 + t103) / 0.2e1 + t34 * t82 / 0.2e1 - t33 * t83 / 0.2e1 + t42 * t72 / 0.2e1 + t4 * t73 + t63 * (Ifges(6,4) * t43 + Ifges(6,2) * t42) / 0.2e1 - t5 * mrSges(5,2) + t6 * mrSges(5,1) + Ifges(5,3) * qJDD(4) + (t63 * (-Ifges(6,2) * t61 + t103) + t67) * t79 / 0.2e1 + (t114 * t53 + t115 * t52) * g(3) * t56 + t110 * qJD(5) + (t114 * t21 - t115 * (-t52 * t97 - t100)) * g(2) + (t114 * t23 - t115 * (t53 * t57 - t59 * t94)) * g(1) + (-m(6) * t4 - t17) * pkin(4) + (0.2e1 * Ifges(6,5) * t109 + Ifges(6,6) * t63) * qJDD(5) + (-m(6) * t11 + qJD(4) * mrSges(5,1) - t41) * t14 + (-m(6) * t75 - t113) * t13 + (m(6) * t65 - t45 * t82 - t46 * t83 + t112) * pkin(6) + (-t10 * t83 - t9 * t82 + t111) * mrSges(6,3); Ifges(6,5) * t43 + Ifges(6,6) * t42 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t9 * t46 + t10 * t45 - g(1) * ((-t23 * t61 + t60 * t98) * mrSges(6,1) + (-t23 * t63 - t60 * t99) * mrSges(6,2)) - g(2) * ((-t21 * t61 + t57 * t98) * mrSges(6,1) + (-t21 * t63 - t57 * t99) * mrSges(6,2)) - g(3) * ((-t53 * t99 - t95) * mrSges(6,1) + (-t53 * t98 + t96) * mrSges(6,2)) + (t33 * t109 + (-t67 / 0.2e1 + t102 * t109) * qJD(4) + t74 * mrSges(6,3) - (t34 + t51) * t63 / 0.2e1 - t110) * qJD(4);];
tau = t20;
