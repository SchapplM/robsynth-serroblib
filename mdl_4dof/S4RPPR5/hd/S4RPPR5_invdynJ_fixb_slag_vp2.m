% Calculate vector of inverse dynamics joint torques for
% S4RPPR5
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:45
% DurationCPUTime: 1.23s
% Computational Cost: add. (567->165), mult. (1050->223), div. (0->0), fcn. (500->6), ass. (0->74)
t53 = -pkin(1) - pkin(2);
t34 = qJD(1) * t53 + qJD(2);
t47 = sin(pkin(6));
t48 = cos(pkin(6));
t80 = qJ(2) * qJD(1);
t13 = t47 * t34 + t48 * t80;
t11 = -qJD(1) * pkin(5) + t13;
t49 = sin(qJ(4));
t51 = cos(qJ(4));
t3 = qJD(3) * t51 - t11 * t49;
t4 = qJD(3) * t49 + t11 * t51;
t33 = qJDD(1) * t53 + qJDD(2);
t76 = qJD(1) * qJD(2);
t35 = qJDD(1) * qJ(2) + t76;
t9 = t47 * t33 + t48 * t35;
t7 = -qJDD(1) * pkin(5) + t9;
t1 = qJD(4) * t3 + qJDD(3) * t49 + t51 * t7;
t2 = -qJD(4) * t4 + qJDD(3) * t51 - t49 * t7;
t67 = t1 * t51 - t2 * t49;
t81 = qJD(4) * t51;
t82 = qJD(4) * t49;
t111 = -t3 * t81 - t4 * t82 + t67;
t75 = qJD(1) * qJD(4);
t24 = -qJDD(1) * t51 + t49 * t75;
t25 = -qJDD(1) * t49 - t51 * t75;
t8 = t33 * t48 - t35 * t47;
t6 = qJDD(1) * pkin(3) - t8;
t110 = m(5) * t6 - mrSges(5,1) * t24 + mrSges(5,2) * t25;
t32 = mrSges(5,1) * t51 - mrSges(5,2) * t49;
t92 = t34 * t48;
t109 = m(4) * (-(-t47 * t80 + t92) * t47 + t13 * t48) + t47 * t32 * qJD(1);
t100 = -t49 / 0.2e1;
t108 = -mrSges(3,1) - mrSges(2,1);
t107 = -mrSges(3,3) + mrSges(2,2);
t106 = m(5) * pkin(3) + mrSges(4,1) + t32;
t84 = mrSges(5,3) * qJD(1);
t30 = qJD(4) * mrSges(5,1) + t49 * t84;
t31 = -qJD(4) * mrSges(5,2) - t51 * t84;
t104 = -t49 * t30 + t51 * t31;
t14 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t24;
t15 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t25;
t103 = t51 * t14 - t49 * t15;
t65 = -t3 * t49 + t4 * t51;
t85 = qJ(2) * t47;
t10 = -t92 + (pkin(3) + t85) * qJD(1);
t64 = mrSges(5,1) * t49 + mrSges(5,2) * t51;
t102 = -t10 * t64 + qJD(4) * (-Ifges(5,5) * t51 + Ifges(5,6) * t49) / 0.2e1;
t50 = sin(qJ(1));
t52 = cos(qJ(1));
t18 = -t50 * t47 - t52 * t48;
t19 = t47 * t52 - t50 * t48;
t101 = -g(1) * t18 - g(2) * t19;
t95 = Ifges(5,4) * t49;
t94 = Ifges(5,4) * t51;
t93 = t10 * t47;
t27 = t48 * qJ(2) + t47 * t53;
t86 = t52 * pkin(1) + t50 * qJ(2);
t79 = qJDD(1) * mrSges(3,1);
t78 = qJDD(1) * mrSges(4,1);
t77 = qJDD(1) * mrSges(4,2);
t73 = t52 * pkin(2) + t86;
t44 = t52 * qJ(2);
t70 = -pkin(1) * t50 + t44;
t66 = -t3 * t51 - t4 * t49;
t63 = -Ifges(5,1) * t49 - t94;
t62 = -t51 * Ifges(5,2) - t95;
t26 = t48 * t53 - t85;
t57 = t49 * (-Ifges(5,1) * t51 + t95);
t56 = t51 * (Ifges(5,2) * t49 - t94);
t54 = qJD(1) ^ 2;
t41 = -qJDD(1) * pkin(1) + qJDD(2);
t17 = Ifges(5,5) * qJD(4) + qJD(1) * t63;
t16 = Ifges(5,6) * qJD(4) + qJD(1) * t62;
t5 = [(Ifges(4,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) + (Ifges(5,1) * t25 + Ifges(5,4) * t24) * t100 + (t48 * t76 + t9) * mrSges(4,2) + (t47 * t76 - t8) * mrSges(4,1) + m(4) * (t26 * t8 + t27 * t9) + t27 * t77 + pkin(1) * t79 - t51 * (Ifges(5,4) * t25 + Ifges(5,2) * t24) / 0.2e1 + m(3) * (-pkin(1) * t41 + (t35 + t76) * qJ(2)) - t26 * t78 - t17 * t81 / 0.2e1 + t16 * t82 / 0.2e1 + t24 * t62 / 0.2e1 + t25 * t63 / 0.2e1 + t110 * (pkin(3) - t26) + (t101 - t111) * mrSges(5,3) + (m(5) * (t48 * t65 + t93) + t104 * t48 + t109) * qJD(2) + t6 * t32 + 0.2e1 * t35 * mrSges(3,3) - t41 * mrSges(3,1) - (t57 + t56) * t75 / 0.2e1 + t102 * qJD(4) + (-t30 * t81 - t31 * t82 + m(5) * (qJD(4) * t66 + t67) + t103) * (-pkin(5) + t27) + (-m(5) * (pkin(5) * t18 + t70) - m(4) * t44 + t18 * mrSges(4,2) - m(3) * t70 + t107 * t52 + (-m(4) * t53 + m(5) * pkin(2) - t108) * t50 - t106 * t19) * g(1) + (-m(5) * (pkin(5) * t19 + t73) - m(3) * t86 - m(4) * t73 + t19 * mrSges(4,2) + t108 * t52 + t107 * t50 + t106 * t18) * g(2) + (0.2e1 * Ifges(5,5) * t100 - Ifges(5,6) * t51) * qJDD(4); -t79 - t54 * mrSges(3,3) + (m(4) * t8 - t54 * mrSges(4,2) - t110 - t78) * t48 + (-t54 * mrSges(4,1) + t77 + (-t51 * t30 - t49 * t31) * qJD(4) + m(4) * t9 + m(5) * t111 + t103) * t47 + (-t54 * qJ(2) + t41) * m(3) + (-m(5) * t93 + (-m(5) * t65 - t104) * t48 - t109) * qJD(1) + (m(5) + m(3) + m(4)) * (-g(1) * t50 + g(2) * t52); t49 * t14 + t51 * t15 + t104 * qJD(4) + (qJD(4) * t65 + t1 * t49 + t2 * t51 + g(3)) * m(5) + (qJDD(3) + g(3)) * m(4); t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t25 + Ifges(5,6) * t24 + Ifges(5,3) * qJDD(4) + g(3) * t32 - t3 * t31 + t4 * t30 + (t51 * t17 / 0.2e1 + t16 * t100 + (t57 / 0.2e1 + t56 / 0.2e1) * qJD(1) + t66 * mrSges(5,3) - t102) * qJD(1) + t101 * t64;];
tau = t5;
