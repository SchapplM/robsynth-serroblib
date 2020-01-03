% Calculate vector of inverse dynamics joint torques for
% S4PRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:22:58
% EndTime: 2019-12-31 16:23:01
% DurationCPUTime: 1.34s
% Computational Cost: add. (540->161), mult. (1175->232), div. (0->0), fcn. (733->10), ass. (0->82)
t48 = sin(qJ(4));
t50 = cos(qJ(4));
t63 = -t50 * mrSges(5,1) + t48 * mrSges(5,2);
t110 = -mrSges(4,1) + t63;
t109 = m(4) + m(5);
t108 = m(5) * pkin(3) - t110;
t100 = t48 / 0.2e1;
t43 = qJ(2) + pkin(7);
t40 = sin(t43);
t107 = g(3) * t40;
t106 = t110 * qJD(2);
t78 = qJD(2) * t48;
t33 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t78;
t77 = qJD(2) * t50;
t34 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t77;
t59 = -t48 * t33 + t50 * t34;
t105 = -qJD(2) * mrSges(4,2) + t59;
t72 = qJD(2) * qJD(4);
t29 = qJDD(2) * t50 - t48 * t72;
t15 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t29;
t30 = qJDD(2) * t48 + t50 * t72;
t16 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t30;
t104 = t50 * t15 - t48 * t16;
t51 = cos(qJ(2));
t49 = sin(qJ(2));
t73 = qJD(1) * qJD(2);
t71 = t49 * t73;
t31 = t51 * qJDD(1) - t71;
t74 = qJDD(2) * pkin(2);
t24 = t31 + t74;
t70 = t51 * t73;
t32 = qJDD(1) * t49 + t70;
t44 = sin(pkin(7));
t46 = cos(pkin(7));
t6 = t44 * t24 + t46 * t32;
t4 = qJDD(2) * pkin(5) + t6;
t35 = qJD(2) * pkin(2) + qJD(1) * t51;
t81 = qJD(1) * t49;
t13 = t44 * t35 + t46 * t81;
t11 = qJD(2) * pkin(5) + t13;
t7 = qJD(3) * t50 - t11 * t48;
t1 = qJD(4) * t7 + qJDD(3) * t48 + t4 * t50;
t8 = qJD(3) * t48 + t11 * t50;
t2 = -qJD(4) * t8 + qJDD(3) * t50 - t4 * t48;
t103 = t1 * t50 - t2 * t48;
t12 = t35 * t46 - t44 * t81;
t10 = -qJD(2) * pkin(3) - t12;
t62 = mrSges(5,1) * t48 + mrSges(5,2) * t50;
t101 = t10 * t62 + qJD(4) * (Ifges(5,5) * t50 - Ifges(5,6) * t48) / 0.2e1;
t97 = pkin(2) * t51;
t92 = Ifges(5,4) * t48;
t91 = Ifges(5,4) * t50;
t90 = Ifges(5,2) * t50;
t45 = sin(pkin(6));
t89 = t45 * t48;
t88 = t45 * t50;
t47 = cos(pkin(6));
t87 = t47 * t48;
t86 = t47 * t50;
t76 = qJD(4) * t48;
t75 = qJD(4) * t50;
t68 = -g(1) * t45 + g(2) * t47;
t67 = t48 * t8 + t50 * t7;
t66 = -t48 * t7 + t50 * t8;
t65 = t51 * mrSges(3,1) - t49 * mrSges(3,2);
t64 = mrSges(3,1) * t49 + mrSges(3,2) * t51;
t61 = t90 + t92;
t5 = t24 * t46 - t32 * t44;
t26 = t44 * t51 + t46 * t49;
t25 = t44 * t49 - t46 * t51;
t57 = t48 * (Ifges(5,1) * t50 - t92);
t54 = -qJD(4) * t67 + t103;
t41 = cos(t43);
t39 = Ifges(5,4) * t77;
t38 = -pkin(2) * t46 - pkin(3);
t22 = Ifges(5,1) * t78 + Ifges(5,5) * qJD(4) + t39;
t21 = Ifges(5,6) * qJD(4) + qJD(2) * t61;
t20 = t25 * qJD(2);
t18 = t26 * qJD(2);
t9 = -mrSges(5,1) * t29 + mrSges(5,2) * t30;
t3 = -qJDD(2) * pkin(3) - t5;
t14 = [m(2) * qJDD(1) + t25 * t9 - t64 * qJD(2) ^ 2 + t106 * t18 - t105 * t20 + ((-t50 * t33 - t48 * t34) * qJD(4) + t104) * t26 + m(3) * (t31 * t51 + t32 * t49) + m(4) * (-t12 * t18 - t13 * t20 - t25 * t5 + t26 * t6) + m(5) * (t10 * t18 - t20 * t66 + t25 * t3 + t26 * t54) + (-mrSges(4,1) * t25 - mrSges(4,2) * t26 + t65) * qJDD(2) + (-m(2) - m(3) - t109) * g(3); (Ifges(5,1) * t30 + Ifges(5,4) * t29) * t100 + t50 * (Ifges(5,4) * t30 + Ifges(5,2) * t29) / 0.2e1 + m(4) * (t44 * t6 + t46 * t5) * pkin(2) + t30 * (Ifges(5,1) * t48 + t91) / 0.2e1 + t22 * t75 / 0.2e1 - t21 * t76 / 0.2e1 + t29 * t61 / 0.2e1 + t38 * t9 + (t50 * (-Ifges(5,2) * t48 + t91) + t57) * t72 / 0.2e1 + (m(5) * t38 + t63) * t3 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t101 * qJD(4) + (-t44 * t74 - t6) * mrSges(4,2) + (t70 - t32) * mrSges(3,2) + (t46 * t74 + t5) * mrSges(4,1) + (t71 + t31) * mrSges(3,1) + (0.2e1 * Ifges(5,5) * t100 + Ifges(5,6) * t50) * qJDD(4) + (-m(5) * (pkin(5) * t40 + t97) - m(4) * t97 + t40 * mrSges(4,2) - t65 - t108 * t41) * g(3) + (m(5) * t54 - t33 * t75 - t34 * t76 + t104) * (pkin(2) * t44 + pkin(5)) + (-t7 * t75 - t76 * t8 + t103 - t107) * mrSges(5,3) + (g(1) * t47 + g(2) * t45) * (t64 + t109 * pkin(2) * t49 + (-m(5) * pkin(5) + mrSges(4,2) - mrSges(5,3)) * t41 + t108 * t40) + ((m(4) * t12 - m(5) * t10 - t106) * t26 - (-m(4) * t13 - m(5) * t66 - t105) * t25) * qJD(1); t48 * t15 + t50 * t16 + t59 * qJD(4) + (qJD(4) * t66 + t1 * t48 + t2 * t50 + t68) * m(5) + (qJDD(3) + t68) * m(4); Ifges(5,5) * t30 + Ifges(5,6) * t29 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t7 * t34 + t8 * t33 - g(1) * ((-t41 * t87 + t88) * mrSges(5,1) + (-t41 * t86 - t89) * mrSges(5,2)) - g(2) * ((-t41 * t89 - t86) * mrSges(5,1) + (-t41 * t88 + t87) * mrSges(5,2)) + t62 * t107 + (t21 * t100 + (-t57 / 0.2e1 + t90 * t100) * qJD(2) + t67 * mrSges(5,3) - (t22 + t39) * t50 / 0.2e1 - t101) * qJD(2);];
tau = t14;
