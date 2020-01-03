% Calculate vector of inverse dynamics joint torques for
% S4RPPR7
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:36
% EndTime: 2019-12-31 16:41:38
% DurationCPUTime: 1.57s
% Computational Cost: add. (695->170), mult. (1405->220), div. (0->0), fcn. (789->8), ass. (0->80)
t59 = sin(pkin(6));
t60 = cos(pkin(6));
t88 = t59 ^ 2 + t60 ^ 2;
t81 = t88 * mrSges(4,3);
t62 = -pkin(1) - qJ(3);
t113 = -qJD(1) * qJD(3) + qJDD(1) * t62;
t63 = sin(qJ(4));
t65 = cos(qJ(4));
t91 = t60 * t65;
t111 = -t59 * t63 + t91;
t26 = t111 * qJD(4);
t66 = cos(qJ(1));
t100 = g(2) * t66;
t29 = -t65 * t59 - t63 * t60;
t24 = t29 * qJD(1);
t87 = qJD(1) * t59;
t25 = qJD(1) * t91 - t63 * t87;
t72 = t59 * mrSges(4,1) + t60 * mrSges(4,2);
t112 = mrSges(5,1) * t24 - mrSges(5,2) * t25 - t72 * qJD(1);
t57 = qJD(1) * qJD(2);
t38 = -qJDD(1) * qJ(2) - t57;
t64 = sin(qJ(1));
t73 = -g(1) * t64 + t100;
t109 = -m(4) - m(5) - m(3);
t108 = mrSges(2,1) + mrSges(4,3) - mrSges(3,2);
t102 = pkin(3) * t59;
t55 = pkin(6) + qJ(4);
t48 = sin(t55);
t49 = cos(t55);
t70 = -t48 * mrSges(5,1) - t49 * mrSges(5,2);
t106 = -m(5) * t102 + mrSges(2,2) - mrSges(3,3) + t70 - t72;
t32 = qJDD(2) + t113;
t75 = -pkin(5) * qJDD(1) + t32;
t17 = t75 * t59;
t18 = t75 * t60;
t37 = qJD(1) * t62 + qJD(2);
t76 = -pkin(5) * qJD(1) + t37;
t21 = t76 * t59;
t22 = t76 * t60;
t5 = -t21 * t63 + t22 * t65;
t1 = qJD(4) * t5 + t17 * t65 + t18 * t63;
t6 = t21 * t65 + t22 * t63;
t2 = -qJD(4) * t6 - t17 * t63 + t18 * t65;
t27 = t29 * qJD(4);
t105 = t1 * t29 - t111 * t2 - t26 * t6 - t27 * t5;
t103 = t25 / 0.2e1;
t95 = -pkin(5) + t62;
t94 = Ifges(5,4) * t25;
t84 = qJDD(1) * t60;
t85 = qJDD(1) * t59;
t90 = mrSges(4,1) * t85 + mrSges(4,2) * t84;
t86 = qJDD(1) * pkin(1);
t46 = qJD(1) * qJ(2) + qJD(3);
t36 = qJDD(3) - t38;
t80 = t88 * t37;
t79 = t88 * t32;
t13 = qJD(1) * t27 + qJDD(1) * t111;
t14 = -qJD(1) * t26 + qJDD(1) * t29;
t77 = -t14 * mrSges(5,1) + t13 * mrSges(5,2);
t74 = -g(1) * t66 - g(2) * t64;
t33 = t95 * t59;
t34 = t95 * t60;
t16 = t33 * t65 + t34 * t63;
t15 = -t33 * t63 + t34 * t65;
t67 = qJD(1) ^ 2;
t61 = -pkin(5) - qJ(3);
t47 = qJDD(2) - t86;
t41 = qJ(2) + t102;
t35 = pkin(3) * t87 + t46;
t28 = pkin(3) * t85 + t36;
t23 = Ifges(5,4) * t24;
t20 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t25;
t19 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t24;
t10 = Ifges(5,1) * t25 + Ifges(5,5) * qJD(4) + t23;
t9 = Ifges(5,2) * t24 + Ifges(5,6) * qJD(4) + t94;
t8 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t14;
t7 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t13;
t4 = -qJD(3) * t111 - qJD(4) * t16;
t3 = qJD(3) * t29 + qJD(4) * t15;
t11 = [(-t86 + t47) * mrSges(3,2) + (-t32 - t113) * t81 + (Ifges(5,1) * t27 - Ifges(5,4) * t26) * t103 + (Ifges(4,1) * t60 - Ifges(4,4) * t59) * t84 + (-mrSges(5,1) * t28 + Ifges(5,4) * t13 + Ifges(5,2) * t14 + Ifges(5,6) * qJDD(4)) * t29 + (mrSges(5,2) * t28 + Ifges(5,1) * t13 + Ifges(5,4) * t14 + Ifges(5,5) * qJDD(4)) * t111 + t36 * t72 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + m(4) * (qJ(2) * t36 - qJD(3) * t80 + t62 * t79) + qJ(2) * t90 - (Ifges(4,4) * t60 - Ifges(4,2) * t59) * t85 + m(3) * (-pkin(1) * t47 + (-t38 + t57) * qJ(2)) + m(5) * (t1 * t16 + t15 * t2 + t28 * t41 + t3 * t6 + t4 * t5) + ((-m(5) * (-pkin(1) + t61) + mrSges(5,3) - m(4) * t62 + m(3) * pkin(1) + t108) * t64 + (t109 * qJ(2) + t106) * t66) * g(1) + (t105 - t100) * mrSges(5,3) + (t109 * (t66 * pkin(1) + t64 * qJ(2)) + (-m(4) * qJ(3) + m(5) * t61 - t108) * t66 + t106 * t64) * g(2) + (m(4) * t46 + m(5) * t35 - t112) * qJD(2) + t35 * (t26 * mrSges(5,1) + t27 * mrSges(5,2)) - 0.2e1 * t38 * mrSges(3,3) - t26 * t9 / 0.2e1 + t24 * (Ifges(5,4) * t27 - Ifges(5,2) * t26) / 0.2e1 + qJD(4) * (Ifges(5,5) * t27 - Ifges(5,6) * t26) / 0.2e1 + t27 * t10 / 0.2e1 + t15 * t7 + t16 * t8 + t3 * t19 + t4 * t20 + t41 * t77; -t67 * mrSges(3,3) + t26 * t19 + t27 * t20 - t29 * t8 + t111 * t7 + t112 * qJD(1) + (mrSges(3,2) - t81) * qJDD(1) + (-t35 * qJD(1) - t105 + t73) * m(5) + (-t46 * qJD(1) + t73 + t79) * m(4) + (-t67 * qJ(2) + t47 + t73) * m(3); -t24 * t19 + t25 * t20 - t67 * t81 + t77 + t90 + (-t6 * t24 + t5 * t25 + t28 + t74) * m(5) + (qJD(1) * t80 + t36 + t74) * m(4); Ifges(5,5) * t13 + Ifges(5,6) * t14 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t35 * (mrSges(5,1) * t25 + mrSges(5,2) * t24) - t25 * (Ifges(5,1) * t24 - t94) / 0.2e1 + t9 * t103 - qJD(4) * (Ifges(5,5) * t24 - Ifges(5,6) * t25) / 0.2e1 - t5 * t19 + t6 * t20 - g(3) * t70 + (t24 * t5 + t25 * t6) * mrSges(5,3) + t73 * (mrSges(5,1) * t49 - mrSges(5,2) * t48) - (-Ifges(5,2) * t25 + t10 + t23) * t24 / 0.2e1;];
tau = t11;
