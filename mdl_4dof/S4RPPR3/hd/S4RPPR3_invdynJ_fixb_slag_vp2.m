% Calculate vector of inverse dynamics joint torques for
% S4RPPR3
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
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:52
% DurationCPUTime: 1.41s
% Computational Cost: add. (692->162), mult. (1507->225), div. (0->0), fcn. (949->12), ass. (0->77)
t61 = sin(pkin(6));
t43 = pkin(1) * t61 + qJ(3);
t36 = qJD(1) * qJD(3) + qJDD(1) * t43;
t60 = sin(pkin(7));
t62 = cos(pkin(7));
t102 = t60 ^ 2 + t62 ^ 2;
t97 = m(3) + m(4) + m(5);
t101 = pkin(1) * t97 + mrSges(2,1);
t48 = t62 * qJDD(2);
t20 = -t36 * t60 + t48;
t21 = t60 * qJDD(2) + t62 * t36;
t99 = -t20 * t60 + t21 * t62;
t45 = pkin(3) * t62 + pkin(2);
t58 = pkin(7) + qJ(4);
t51 = sin(t58);
t53 = cos(t58);
t75 = t53 * mrSges(5,1) - t51 * mrSges(5,2);
t76 = -t62 * mrSges(4,1) + t60 * mrSges(4,2);
t96 = m(4) * pkin(2) + m(5) * t45 + mrSges(3,1) + t75 - t76;
t95 = -m(4) * qJ(3) + m(5) * (-pkin(5) - qJ(3)) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t39 = t60 * t67 + t62 * t65;
t31 = t39 * qJD(1);
t93 = t31 / 0.2e1;
t63 = cos(pkin(6));
t92 = pkin(1) * t63;
t90 = pkin(5) + t43;
t89 = Ifges(5,4) * t31;
t42 = t43 * qJD(1);
t25 = t60 * qJD(2) + t62 * t42;
t86 = pkin(5) * qJD(1);
t84 = qJDD(1) * t60;
t83 = qJDD(1) * t62;
t46 = -pkin(2) - t92;
t38 = -t60 * t65 + t62 * t67;
t32 = t38 * qJD(4);
t14 = qJD(1) * t32 + qJDD(1) * t39;
t33 = t39 * qJD(4);
t15 = -qJD(1) * t33 + qJDD(1) * t38;
t79 = -t15 * mrSges(5,1) + t14 * mrSges(5,2);
t59 = qJ(1) + pkin(6);
t52 = sin(t59);
t54 = cos(t59);
t78 = -g(1) * t52 + g(2) * t54;
t77 = -mrSges(4,1) * t83 + mrSges(4,2) * t84;
t50 = t62 * qJD(2);
t18 = t50 + (-t42 - t86) * t60;
t19 = t62 * t86 + t25;
t3 = t18 * t67 - t19 * t65;
t4 = t18 * t65 + t19 * t67;
t73 = -(-t42 * t60 + t50) * t60 + t25 * t62;
t34 = t90 * t60;
t35 = t90 * t62;
t11 = -t34 * t67 - t35 * t65;
t12 = -t34 * t65 + t35 * t67;
t41 = -t45 - t92;
t68 = cos(qJ(1));
t66 = sin(qJ(1));
t40 = qJDD(1) * t46 + qJDD(3);
t30 = t38 * qJD(1);
t29 = qJD(1) * t41 + qJD(3);
t27 = qJDD(1) * t41 + qJDD(3);
t26 = Ifges(5,4) * t30;
t23 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t31;
t22 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t30;
t17 = pkin(5) * t83 + t21;
t16 = t48 + (-pkin(5) * qJDD(1) - t36) * t60;
t10 = Ifges(5,1) * t31 + Ifges(5,5) * qJD(4) + t26;
t9 = Ifges(5,2) * t30 + Ifges(5,6) * qJD(4) + t89;
t8 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t15;
t7 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t14;
t6 = -qJD(3) * t39 - qJD(4) * t12;
t5 = qJD(3) * t38 + qJD(4) * t11;
t2 = -qJD(4) * t4 + t16 * t67 - t17 * t65;
t1 = qJD(4) * t3 + t16 * t65 + t17 * t67;
t13 = [m(5) * (t1 * t12 + t11 * t2 + t27 * t41 + t3 * t6 + t4 * t5) + t41 * t79 + t40 * t76 + t46 * t77 + (mrSges(5,2) * t27 - mrSges(5,3) * t2 + Ifges(5,1) * t14 + Ifges(5,4) * t15 + Ifges(5,5) * qJDD(4)) * t39 + (-mrSges(5,1) * t27 + mrSges(5,3) * t1 + Ifges(5,4) * t14 + Ifges(5,2) * t15 + Ifges(5,6) * qJDD(4)) * t38 + (t68 * mrSges(2,2) + t101 * t66 + t52 * t96 + t54 * t95) * g(1) + (t66 * mrSges(2,2) - t101 * t68 + t52 * t95 - t54 * t96) * g(2) + t32 * t10 / 0.2e1 - t33 * t9 / 0.2e1 + t5 * t22 + t6 * t23 + t11 * t7 + t12 * t8 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t63 * mrSges(3,1) - 0.2e1 * t61 * mrSges(3,2) + m(3) * (t61 ^ 2 + t63 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (Ifges(4,4) * t60 + Ifges(4,2) * t62) * t83 + (Ifges(4,1) * t60 + Ifges(4,4) * t62) * t84 + (-t3 * t32 - t4 * t33) * mrSges(5,3) + t29 * (mrSges(5,1) * t33 + mrSges(5,2) * t32) + t30 * (Ifges(5,4) * t32 - Ifges(5,2) * t33) / 0.2e1 + qJD(4) * (Ifges(5,5) * t32 - Ifges(5,6) * t33) / 0.2e1 + (Ifges(5,1) * t32 - Ifges(5,4) * t33) * t93 + m(4) * (t73 * qJD(3) + t40 * t46 + t43 * t99) + (t102 * t36 + t99) * mrSges(4,3); m(3) * qJDD(2) + t32 * t22 - t33 * t23 + t38 * t7 + t39 * t8 + m(4) * (t20 * t62 + t21 * t60) + m(5) * (t1 * t39 + t2 * t38 - t3 * t33 + t32 * t4) - t97 * g(3); -t30 * t22 + t31 * t23 - t102 * qJD(1) ^ 2 * mrSges(4,3) + t77 + t79 + (t3 * t31 - t4 * t30 + t27 + t78) * m(5) + (-qJD(1) * t73 + t40 + t78) * m(4); Ifges(5,5) * t14 + Ifges(5,6) * t15 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t29 * (mrSges(5,1) * t31 + mrSges(5,2) * t30) - t31 * (Ifges(5,1) * t30 - t89) / 0.2e1 + t9 * t93 - qJD(4) * (Ifges(5,5) * t30 - Ifges(5,6) * t31) / 0.2e1 - t3 * t22 + t4 * t23 - g(3) * t75 + (t3 * t30 + t31 * t4) * mrSges(5,3) - (-Ifges(5,2) * t31 + t10 + t26) * t30 / 0.2e1 + (g(1) * t54 + g(2) * t52) * (mrSges(5,1) * t51 + mrSges(5,2) * t53);];
tau = t13;
