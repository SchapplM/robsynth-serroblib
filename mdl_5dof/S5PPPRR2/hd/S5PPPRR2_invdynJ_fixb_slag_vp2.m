% Calculate vector of inverse dynamics joint torques for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:26
% EndTime: 2019-12-05 14:59:30
% DurationCPUTime: 1.34s
% Computational Cost: add. (875->205), mult. (1998->318), div. (0->0), fcn. (1580->10), ass. (0->101)
t56 = sin(pkin(9));
t59 = cos(pkin(9));
t57 = sin(pkin(8));
t96 = qJD(1) * t57;
t42 = qJD(2) * t56 + t59 * t96;
t60 = cos(pkin(8));
t52 = -qJD(1) * t60 + qJD(3);
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t22 = t42 * t65 + t52 * t63;
t12 = qJD(4) * pkin(6) + t22;
t41 = -qJD(2) * t59 + t56 * t96;
t62 = sin(qJ(5));
t64 = cos(qJ(5));
t10 = t12 * t64 + t41 * t62;
t87 = qJDD(1) * t57;
t40 = qJDD(2) * t56 + t59 * t87;
t51 = -qJDD(1) * t60 + qJDD(3);
t92 = qJD(4) * t65;
t94 = qJD(4) * t63;
t5 = t65 * t40 - t42 * t94 + t63 * t51 + t52 * t92;
t3 = qJDD(4) * pkin(6) + t5;
t39 = -t59 * qJDD(2) + t56 * t87;
t9 = -t12 * t62 + t41 * t64;
t1 = t9 * qJD(5) + t3 * t64 + t39 * t62;
t2 = -t10 * qJD(5) - t3 * t62 + t39 * t64;
t77 = t1 * t64 - t2 * t62;
t90 = qJD(5) * t64;
t91 = qJD(5) * t62;
t121 = -t10 * t91 - t9 * t90 + t77;
t114 = t62 / 0.2e1;
t21 = -t42 * t63 + t52 * t65;
t120 = -t21 * qJD(4) + t5;
t88 = t22 * qJD(4);
t6 = -t40 * t63 + t51 * t65 - t88;
t119 = t6 + t88;
t95 = qJD(4) * t62;
t48 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t95;
t93 = qJD(4) * t64;
t49 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t93;
t117 = -t62 * t48 + t64 * t49;
t85 = qJD(4) * qJD(5);
t45 = qJDD(4) * t64 - t62 * t85;
t26 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t45;
t46 = qJDD(4) * t62 + t64 * t85;
t27 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t46;
t116 = t64 * t26 - t62 * t27;
t11 = -qJD(4) * pkin(4) - t21;
t115 = t11 * (mrSges(6,1) * t62 + mrSges(6,2) * t64) + qJD(5) * (Ifges(6,5) * t64 - Ifges(6,6) * t62) / 0.2e1;
t110 = Ifges(6,4) * t62;
t109 = Ifges(6,4) * t64;
t108 = Ifges(6,2) * t64;
t107 = t39 * t59;
t106 = t56 * t57;
t105 = t56 * t65;
t104 = t57 * t63;
t103 = t57 * t65;
t58 = sin(pkin(7));
t102 = t58 * t60;
t61 = cos(pkin(7));
t101 = t60 * t61;
t86 = m(4) + m(5) + m(6);
t84 = t59 * t103;
t82 = t56 * t94;
t80 = m(3) + t86;
t79 = m(6) * pkin(6) + mrSges(6,3);
t76 = t10 * t62 + t9 * t64;
t75 = -mrSges(6,1) * t64 + mrSges(6,2) * t62;
t74 = t108 + t110;
t34 = -t60 * t63 + t84;
t17 = t64 * t106 - t34 * t62;
t18 = t62 * t106 + t34 * t64;
t36 = t64 * t105 - t59 * t62;
t35 = -t62 * t105 - t59 * t64;
t33 = t59 * t104 + t60 * t65;
t23 = -mrSges(6,1) * t45 + mrSges(6,2) * t46;
t66 = qJD(4) ^ 2;
t71 = qJDD(4) * mrSges(5,1) - mrSges(5,2) * t66 - t23;
t69 = t62 * (Ifges(6,1) * t64 - t110);
t68 = m(6) * pkin(4) - t75;
t43 = t75 * qJD(4);
t67 = -mrSges(5,1) * t66 - qJDD(4) * mrSges(5,2) + qJD(4) * t43;
t54 = Ifges(6,4) * t93;
t38 = Ifges(6,1) * t95 + Ifges(6,5) * qJD(5) + t54;
t37 = Ifges(6,6) * qJD(5) + t74 * qJD(4);
t32 = t59 * t101 + t56 * t58;
t31 = t56 * t101 - t58 * t59;
t30 = t59 * t102 - t56 * t61;
t29 = t56 * t102 + t59 * t61;
t25 = qJD(4) * t84 - t60 * t94;
t24 = t33 * qJD(4);
t20 = -t36 * qJD(5) + t62 * t82;
t19 = t35 * qJD(5) - t64 * t82;
t16 = t61 * t104 + t32 * t65;
t15 = t61 * t103 - t32 * t63;
t14 = t58 * t104 + t30 * t65;
t13 = t58 * t103 - t30 * t63;
t8 = t17 * qJD(5) - t24 * t64;
t7 = -t18 * qJD(5) + t24 * t62;
t4 = -qJDD(4) * pkin(4) - t6;
t28 = [t17 * t27 + t18 * t26 + t33 * t23 + t25 * t43 + t7 * t48 + t8 * t49 + (qJD(4) * t24 - qJDD(4) * t34) * mrSges(5,2) + (-qJD(4) * t25 - qJDD(4) * t33) * mrSges(5,1) + (-m(2) - t80) * g(3) + m(4) * (-t51 * t60 + (t39 * t56 + t40 * t59) * t57) + m(5) * (t39 * t106 - t21 * t25 - t22 * t24 - t33 * t6 + t34 * t5) + m(6) * (t1 * t18 + t10 * t8 + t11 * t25 + t17 * t2 + t33 * t4 + t7 * t9) + (m(2) + m(3) * (t57 ^ 2 + t60 ^ 2)) * qJDD(1); m(3) * qJDD(2) + t19 * t49 + t20 * t48 + t36 * t26 + t35 * t27 + (-t71 * t63 + t67 * t65) * t56 + m(4) * (t40 * t56 - t107) + m(5) * (-t107 + (t5 * t65 - t6 * t63 + (-t21 * t65 - t22 * t63) * qJD(4)) * t56) + m(6) * (t1 * t36 + t10 * t19 + t2 * t35 + t20 * t9 + (t11 * t92 + t4 * t63) * t56) + (-g(1) * t58 + g(2) * t61) * t80; m(4) * t51 + (t117 * qJD(4) + m(5) * t119 + m(6) * (t10 * t93 - t9 * t95 - t4) + t71) * t65 + ((-t64 * t48 - t62 * t49) * qJD(5) + m(5) * t120 + m(6) * (qJD(4) * t11 + t121) + t67 + t116) * t63 + (g(3) * t60 + (-g(1) * t61 - g(2) * t58) * t57) * t86; (Ifges(6,1) * t46 + Ifges(6,4) * t45) * t114 + t46 * (Ifges(6,1) * t62 + t109) / 0.2e1 + t38 * t90 / 0.2e1 - t37 * t91 / 0.2e1 - g(1) * (t68 * t15 + t79 * t16) - g(2) * (t68 * t13 + t79 * t14) - g(3) * (-t68 * t33 + t79 * t34) + t64 * (Ifges(6,4) * t46 + Ifges(6,2) * t45) / 0.2e1 + t45 * t74 / 0.2e1 + t4 * t75 - t22 * t43 - pkin(4) * t23 + Ifges(5,3) * qJDD(4) + (t64 * (-Ifges(6,2) * t62 + t109) + t69) * t85 / 0.2e1 - t117 * t21 + t115 * qJD(5) + (g(1) * t16 + g(2) * t14 + g(3) * t34 - t120) * mrSges(5,2) + (-g(1) * t15 - g(2) * t13 + g(3) * t33 + t119) * mrSges(5,1) + (-t11 * t22 - (t10 * t64 - t62 * t9) * t21 - pkin(4) * t4) * m(6) + (0.2e1 * Ifges(6,5) * t114 + Ifges(6,6) * t64) * qJDD(5) + (m(6) * (-t76 * qJD(5) + t77) - t48 * t90 - t49 * t91 + t116) * pkin(6) + t121 * mrSges(6,3); Ifges(6,5) * t46 + Ifges(6,6) * t45 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t9 * t49 + t10 * t48 - g(1) * ((-t16 * t62 + t31 * t64) * mrSges(6,1) + (-t16 * t64 - t31 * t62) * mrSges(6,2)) - g(2) * ((-t14 * t62 + t29 * t64) * mrSges(6,1) + (-t14 * t64 - t29 * t62) * mrSges(6,2)) - g(3) * (mrSges(6,1) * t17 - mrSges(6,2) * t18) + (t37 * t114 + (-t69 / 0.2e1 + t108 * t114) * qJD(4) + t76 * mrSges(6,3) - (t38 + t54) * t64 / 0.2e1 - t115) * qJD(4);];
tau = t28;
