% Calculate vector of inverse dynamics joint torques for
% S4RPRP3
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:32
% EndTime: 2019-12-31 16:42:36
% DurationCPUTime: 2.11s
% Computational Cost: add. (579->205), mult. (1232->264), div. (0->0), fcn. (581->8), ass. (0->87)
t114 = Ifges(5,4) + Ifges(4,4);
t52 = sin(pkin(6));
t39 = pkin(1) * t52 + pkin(5);
t30 = t39 * qJDD(1);
t118 = qJD(2) * qJD(3) + t30;
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t38 = -mrSges(4,1) * t57 + mrSges(4,2) * t55;
t65 = -t57 * mrSges(5,1) + t55 * mrSges(5,2);
t117 = t38 + t65;
t106 = m(3) + m(5) + m(4);
t116 = pkin(1) * t106 + mrSges(2,1);
t115 = Ifges(5,1) + Ifges(4,1);
t113 = Ifges(5,5) + Ifges(4,5);
t112 = Ifges(5,2) + Ifges(4,2);
t111 = Ifges(5,6) + Ifges(4,6);
t32 = t39 * qJD(1);
t82 = qJD(2) * t55;
t12 = t32 * t57 + t82;
t69 = qJ(4) * qJD(1) + t32;
t8 = t57 * t69 + t82;
t110 = -t12 * mrSges(4,3) - t8 * mrSges(5,3);
t49 = t57 * qJD(2);
t11 = -t32 * t55 + t49;
t7 = -t55 * t69 + t49;
t5 = qJD(3) * pkin(3) + t7;
t109 = -t11 * mrSges(4,3) - t5 * mrSges(5,3);
t108 = t114 * t57;
t81 = qJD(3) * t55;
t3 = t55 * qJDD(2) + t118 * t57 - t32 * t81;
t48 = t57 * qJDD(2);
t4 = -t12 * qJD(3) - t30 * t55 + t48;
t107 = t3 * t57 - t4 * t55;
t42 = pkin(3) * t57 + pkin(2);
t105 = -m(4) * pkin(2) - m(5) * t42 - mrSges(3,1) + t117;
t104 = -m(4) * pkin(5) + m(5) * (-qJ(4) - pkin(5)) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t53 = cos(pkin(6));
t99 = pkin(1) * t53;
t93 = Ifges(4,4) * t55;
t91 = Ifges(5,4) * t55;
t86 = qJ(4) + t39;
t84 = qJD(1) * t55;
t83 = qJD(1) * t57;
t80 = qJD(3) * t57;
t79 = qJD(1) * qJD(3);
t78 = qJD(1) * qJD(4);
t75 = pkin(3) * t81;
t40 = -pkin(2) - t99;
t27 = qJDD(1) * t57 - t55 * t79;
t28 = qJDD(1) * t55 + t57 * t79;
t72 = -t27 * mrSges(5,1) + t28 * mrSges(5,2);
t70 = qJD(3) * t86;
t66 = -t5 * t55 + t57 * t8;
t64 = Ifges(4,2) * t57 + t93;
t63 = Ifges(5,2) * t57 + t91;
t31 = t40 * qJDD(1);
t29 = -t42 - t99;
t58 = cos(qJ(1));
t56 = sin(qJ(1));
t51 = qJ(1) + pkin(6);
t46 = cos(t51);
t45 = sin(t51);
t44 = Ifges(4,4) * t83;
t43 = Ifges(5,4) * t83;
t37 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t83;
t36 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t83;
t35 = t40 * qJD(1);
t34 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t84;
t33 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t84;
t26 = t65 * qJD(1);
t23 = t86 * t57;
t22 = t86 * t55;
t21 = Ifges(4,1) * t84 + Ifges(4,5) * qJD(3) + t44;
t20 = Ifges(5,1) * t84 + Ifges(5,5) * qJD(3) + t43;
t19 = Ifges(4,6) * qJD(3) + qJD(1) * t64;
t18 = Ifges(5,6) * qJD(3) + qJD(1) * t63;
t17 = qJD(1) * t29 + qJD(4);
t16 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t28;
t15 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t28;
t14 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t27;
t13 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t27;
t10 = -qJD(4) * t55 - t57 * t70;
t9 = qJD(4) * t57 - t55 * t70;
t6 = -pkin(3) * t27 + qJDD(4) + t31;
t2 = qJ(4) * t27 + t57 * t78 + t3;
t1 = -t32 * t80 + qJDD(3) * pkin(3) - qJ(4) * t28 + t48 + (-t78 - t118) * t55;
t24 = [t115 * t55 * t28 + (t114 * t55 + t63 + t64) * t27 / 0.2e1 + ((-t112 * t55 + t108) * t57 + (t115 * t57 - t91 - t93) * t55) * t79 / 0.2e1 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t53 * mrSges(3,1) - 0.2e1 * t52 * mrSges(3,2) + m(3) * (t52 ^ 2 + t53 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t112 * t27 + t114 * t28) * t57 / 0.2e1 + (t58 * mrSges(2,2) + t104 * t46 - t105 * t45 + t116 * t56) * g(1) + (t56 * mrSges(2,2) + t104 * t45 + t105 * t46 - t116 * t58) * g(2) + t108 * t28 / 0.2e1 + t6 * t65 + (t17 * (mrSges(5,1) * t55 + mrSges(5,2) * t57) + t35 * (mrSges(4,1) * t55 + mrSges(4,2) * t57)) * qJD(3) + t26 * t75 + (-t111 * t55 + t113 * t57) * qJD(3) ^ 2 / 0.2e1 + (t111 * t57 + t113 * t55) * qJDD(3) + (-t1 * t55 + t2 * t57) * mrSges(5,3) + t10 * t33 + t9 * t36 + t40 * (-mrSges(4,1) * t27 + mrSges(4,2) * t28) - t22 * t15 + t23 * t13 + (t21 + t20) * t80 / 0.2e1 - (t19 + t18) * t81 / 0.2e1 + m(5) * (-t1 * t22 + t10 * t5 + t17 * t75 + t2 * t23 + t29 * t6 + t8 * t9) + (-t34 * t80 - t37 * t81 + m(4) * ((-t11 * t57 - t12 * t55) * qJD(3) + t107) - t55 * t16 + t57 * t14) * t39 + t107 * mrSges(4,3) + t109 * t80 + t110 * t81 + t29 * t72 + (m(4) * t40 + t38) * t31; m(3) * qJDD(2) + (t15 + t16) * t57 + (t13 + t14) * t55 + ((t36 + t37) * t57 + (-t33 - t34) * t55) * qJD(3) + m(4) * (t3 * t55 + t4 * t57 + (-t11 * t55 + t12 * t57) * qJD(3)) + m(5) * (qJD(3) * t66 + t1 * t57 + t2 * t55) - t106 * g(3); t4 * mrSges(4,1) + t1 * mrSges(5,1) - t3 * mrSges(4,2) - t2 * mrSges(5,2) - t11 * t37 + t12 * t34 - t7 * t36 + t113 * t28 + t111 * t27 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t117 * g(3) + (t15 + (-g(3) * t57 + t1) * m(5)) * pkin(3) + (t33 - m(5) * (-t5 + t7)) * t8 + ((-t35 * mrSges(4,2) - t17 * mrSges(5,2) - t20 / 0.2e1 - t21 / 0.2e1 - t43 / 0.2e1 - t44 / 0.2e1 + (-Ifges(5,5) / 0.2e1 - Ifges(4,5) / 0.2e1) * qJD(3) - t109) * t57 + (-t35 * mrSges(4,1) - t17 * mrSges(5,1) + t18 / 0.2e1 + t19 / 0.2e1 + (Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t84 + (Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1) * qJD(3) + (-m(5) * t17 - t26) * pkin(3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(4,1) / 0.2e1) * t83 - t110) * t55) * qJD(1) + (g(1) * t46 + g(2) * t45) * ((mrSges(4,2) + mrSges(5,2)) * t57 + (m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1)) * t55); (t55 * t33 - t57 * t36) * qJD(1) + (-g(1) * t45 + g(2) * t46 - qJD(1) * t66 + t6) * m(5) + t72;];
tau = t24;
