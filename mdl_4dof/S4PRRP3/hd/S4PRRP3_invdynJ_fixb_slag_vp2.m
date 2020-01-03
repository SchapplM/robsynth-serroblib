% Calculate vector of inverse dynamics joint torques for
% S4PRRP3
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:44
% EndTime: 2019-12-31 16:26:46
% DurationCPUTime: 1.66s
% Computational Cost: add. (471->194), mult. (1028->249), div. (0->0), fcn. (472->4), ass. (0->79)
t104 = Ifges(4,4) + Ifges(5,4);
t45 = sin(qJ(3));
t46 = cos(qJ(3));
t31 = -t46 * mrSges(4,1) + t45 * mrSges(4,2);
t52 = -t46 * mrSges(5,1) + t45 * mrSges(5,2);
t106 = t31 + t52;
t105 = Ifges(5,1) + Ifges(4,1);
t103 = Ifges(5,5) + Ifges(4,5);
t102 = Ifges(5,2) + Ifges(4,2);
t101 = Ifges(5,6) + Ifges(4,6);
t44 = -qJ(4) - pkin(5);
t32 = t44 * t46;
t72 = t45 * qJD(1);
t12 = -qJD(2) * t32 + t72;
t77 = qJD(2) * t46;
t24 = pkin(5) * t77 + t72;
t100 = -t24 * mrSges(4,3) - t12 * mrSges(5,3);
t99 = t104 * t46;
t42 = t46 * qJD(1);
t78 = qJD(2) * t45;
t23 = -pkin(5) * t78 + t42;
t30 = t44 * t45;
t11 = qJD(2) * t30 + t42;
t5 = qJD(3) * pkin(3) + t11;
t98 = -t23 * mrSges(4,3) - t5 * mrSges(5,3);
t69 = qJD(2) * qJD(3);
t62 = t45 * t69;
t70 = qJD(1) * qJD(3);
t71 = qJDD(2) * t46;
t67 = pkin(5) * t71 + t45 * qJDD(1) + t46 * t70;
t3 = -pkin(5) * t62 + t67;
t61 = t46 * t69;
t22 = t45 * qJDD(2) + t61;
t41 = t46 * qJDD(1);
t4 = -pkin(5) * t22 - t45 * t70 + t41;
t97 = t3 * t46 - t4 * t45;
t73 = qJDD(2) * pkin(2);
t34 = pkin(3) * t46 + pkin(2);
t96 = m(4) * pkin(2) + m(5) * t34 + mrSges(3,1) - t106;
t95 = -m(4) * pkin(5) + m(5) * t44 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t90 = m(2) + m(3);
t86 = Ifges(4,4) * t45;
t84 = Ifges(5,4) * t45;
t76 = qJD(3) * t45;
t75 = qJD(3) * t46;
t74 = qJD(4) * t46;
t66 = pkin(3) * t76;
t65 = pkin(5) * t76;
t21 = -t62 + t71;
t58 = -t21 * mrSges(5,1) + t22 * mrSges(5,2);
t57 = qJD(3) * t44;
t53 = t12 * t46 - t45 * t5;
t51 = Ifges(4,2) * t46 + t86;
t50 = Ifges(5,2) * t46 + t84;
t43 = pkin(6) + qJ(2);
t39 = cos(t43);
t38 = sin(t43);
t37 = Ifges(4,4) * t77;
t36 = Ifges(5,4) * t77;
t29 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t77;
t28 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t77;
t27 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t78;
t26 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t78;
t25 = -qJD(2) * t34 + qJD(4);
t20 = t52 * qJD(2);
t18 = Ifges(4,1) * t78 + Ifges(4,5) * qJD(3) + t37;
t17 = Ifges(5,1) * t78 + Ifges(5,5) * qJD(3) + t36;
t16 = Ifges(4,6) * qJD(3) + qJD(2) * t51;
t15 = Ifges(5,6) * qJD(3) + qJD(2) * t50;
t14 = -t45 * qJD(4) + t46 * t57;
t13 = t45 * t57 + t74;
t10 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t22;
t9 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t22;
t8 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t21;
t7 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t21;
t6 = -pkin(3) * t21 + qJDD(4) - t73;
t2 = t21 * qJ(4) + (-t65 + t74) * qJD(2) + t67;
t1 = -pkin(5) * t61 + qJDD(3) * pkin(3) - t22 * qJ(4) + t41 + (-pkin(5) * qJDD(2) - qJD(2) * qJD(4) - t70) * t45;
t19 = [(t10 + t9) * t46 + (t7 + t8) * t45 + t90 * qJDD(1) + ((t28 + t29) * t46 + (-t26 - t27) * t45) * qJD(3) + m(4) * (t3 * t45 + t4 * t46 + (-t23 * t45 + t24 * t46) * qJD(3)) + m(5) * (qJD(3) * t53 + t1 * t46 + t2 * t45) + (-m(4) - m(5) - t90) * g(3); t20 * t66 + (t102 * t21 + t104 * t22) * t46 / 0.2e1 + t99 * t22 / 0.2e1 + (t104 * t45 + t50 + t51) * t21 / 0.2e1 + (-t1 * t45 + t2 * t46) * mrSges(5,3) + t25 * (mrSges(5,1) * t45 + mrSges(5,2) * t46) * qJD(3) + t105 * t45 * t22 + (-t101 * t45 + t103 * t46) * qJD(3) ^ 2 / 0.2e1 + (t101 * t46 + t103 * t45) * qJDD(3) + (-(mrSges(4,1) * t45 + mrSges(4,2) * t46) * t69 + t21 * mrSges(4,1) - t22 * mrSges(4,2) + m(4) * t73) * pkin(2) + (m(4) * ((-t23 * t46 - t24 * t45) * qJD(3) + t97) - t45 * t10 + t46 * t8 - t27 * t75) * pkin(5) + t97 * mrSges(4,3) + t98 * t75 + t100 * t76 + (t18 + t17) * t75 / 0.2e1 - (t16 + t15) * t76 / 0.2e1 + m(5) * (t1 * t30 + t12 * t13 + t14 * t5 - t2 * t32 + t25 * t66 - t34 * t6) + (t38 * t96 + t39 * t95) * g(1) + (t38 * t95 - t39 * t96) * g(2) - t31 * t73 - t29 * t65 - t34 * t58 + t6 * t52 + ((-t102 * t45 + t99) * t46 + (t105 * t46 - t84 - t86) * t45) * t69 / 0.2e1 + Ifges(3,3) * qJDD(2) + t14 * t26 + t13 * t28 + t30 * t9 - t32 * t7; t4 * mrSges(4,1) + t1 * mrSges(5,1) - t3 * mrSges(4,2) - t2 * mrSges(5,2) - t11 * t28 - t23 * t29 + t24 * t27 + t103 * t22 + t101 * t21 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t106 * g(3) + (t9 + (-g(3) * t46 + t1) * m(5)) * pkin(3) + (-m(5) * (t11 - t5) + t26) * t12 + ((-t25 * mrSges(5,2) - t17 / 0.2e1 - t18 / 0.2e1 - t36 / 0.2e1 - t37 / 0.2e1 + pkin(2) * mrSges(4,2) * qJD(2) + (-Ifges(5,5) / 0.2e1 - Ifges(4,5) / 0.2e1) * qJD(3) - t98) * t46 + (-t25 * mrSges(5,1) + t15 / 0.2e1 + t16 / 0.2e1 + (Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1) * qJD(3) + (pkin(2) * mrSges(4,1) + (Ifges(5,4) / 0.2e1 + Ifges(4,4) / 0.2e1) * t45) * qJD(2) + (-m(5) * t25 - t20) * pkin(3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(4,1) / 0.2e1) * t77 - t100) * t45) * qJD(2) + (g(1) * t39 + g(2) * t38) * ((mrSges(4,2) + mrSges(5,2)) * t46 + (m(5) * pkin(3) + mrSges(4,1) + mrSges(5,1)) * t45); (t45 * t26 - t46 * t28) * qJD(2) + (-g(1) * t38 + g(2) * t39 - qJD(2) * t53 + t6) * m(5) + t58;];
tau = t19;
