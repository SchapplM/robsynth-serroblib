% Calculate vector of inverse dynamics joint torques for
% S4RPPR4
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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

function tau = S4RPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:50
% DurationCPUTime: 0.87s
% Computational Cost: add. (412->132), mult. (721->176), div. (0->0), fcn. (307->8), ass. (0->61)
t32 = sin(qJ(4));
t77 = -t32 / 0.2e1;
t34 = cos(qJ(4));
t68 = t34 / 0.2e1;
t76 = -m(4) - m(5);
t29 = qJ(1) + pkin(6);
t27 = cos(t29);
t64 = g(2) * t27;
t75 = mrSges(3,1) - mrSges(4,2);
t51 = qJD(1) * qJD(4);
t13 = qJDD(1) * t34 - t32 * t51;
t5 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t13;
t14 = -qJDD(1) * t32 - t34 * t51;
t6 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t14;
t74 = t32 * t6 + t34 * t5;
t31 = cos(pkin(6));
t25 = -pkin(1) * t31 - pkin(2);
t20 = -pkin(5) + t25;
t10 = t20 * qJD(1) + qJD(3);
t3 = -qJD(2) * t32 + t10 * t34;
t9 = t20 * qJDD(1) + qJDD(3);
t1 = qJD(4) * t3 + qJDD(2) * t34 + t32 * t9;
t4 = qJD(2) * t34 + t10 * t32;
t2 = -qJD(4) * t4 - qJDD(2) * t32 + t34 * t9;
t73 = t1 * t32 + t2 * t34;
t26 = sin(t29);
t72 = -g(1) * t26 + t64;
t46 = mrSges(5,1) * t32 + mrSges(5,2) * t34;
t71 = -mrSges(4,3) - t46 + mrSges(3,2);
t58 = Ifges(5,4) * t34;
t59 = Ifges(5,4) * t32;
t70 = (-Ifges(5,2) * t34 - t59) * t77 + (-Ifges(5,1) * t32 - t58) * t68;
t30 = sin(pkin(6));
t23 = pkin(1) * t30 + qJ(3);
t16 = t23 * qJD(1);
t47 = mrSges(5,1) * t34 - mrSges(5,2) * t32;
t69 = t16 * t47 + qJD(4) * (-Ifges(5,5) * t32 - Ifges(5,6) * t34) / 0.2e1;
t67 = m(3) + m(4);
t33 = sin(qJ(1));
t66 = pkin(1) * t33;
t35 = cos(qJ(1));
t28 = t35 * pkin(1);
t57 = mrSges(5,3) * qJD(1);
t55 = qJD(4) * t32;
t54 = qJD(4) * t34;
t53 = qJDD(1) * mrSges(4,2);
t52 = qJD(1) * qJD(3);
t48 = -t3 * t32 + t34 * t4;
t45 = Ifges(5,1) * t34 - t59;
t44 = -Ifges(5,2) * t32 + t58;
t11 = t23 * qJDD(1) + t52;
t42 = qJD(3) * t16 + t11 * t23;
t38 = -t16 * qJD(1) + t72;
t37 = t48 * qJD(4) + t73;
t18 = qJD(4) * mrSges(5,1) - t34 * t57;
t17 = -qJD(4) * mrSges(5,2) - t32 * t57;
t15 = t25 * qJDD(1) + qJDD(3);
t12 = t46 * qJD(1);
t8 = Ifges(5,5) * qJD(4) + t45 * qJD(1);
t7 = Ifges(5,6) * qJD(4) + t44 * qJD(1);
t19 = [-t8 * t55 / 0.2e1 - t7 * t54 / 0.2e1 + t14 * t44 / 0.2e1 + t13 * t45 / 0.2e1 + t25 * t53 + (Ifges(5,1) * t13 + Ifges(5,4) * t14) * t68 + m(5) * t42 + m(4) * (t15 * t25 + t42) + qJD(3) * t12 + t15 * mrSges(4,2) + t11 * t46 + t23 * (-t14 * mrSges(5,1) + t13 * mrSges(5,2)) + (Ifges(5,4) * t13 + Ifges(5,2) * t14) * t77 + t70 * t51 + t69 * qJD(4) + (t52 + t11) * mrSges(4,3) + (0.2e1 * Ifges(5,5) * t68 - Ifges(5,6) * t32) * qJDD(4) + (t37 * m(5) + t17 * t54 - t18 * t55 + t74) * t20 + (-m(3) * t28 - t35 * mrSges(2,1) + t33 * mrSges(2,2) + t76 * (t27 * pkin(2) + t26 * qJ(3) + t28) + (-m(5) * pkin(5) - t75) * t27 + t71 * t26) * g(2) + (t3 * t55 - t4 * t54 - t64 - t73) * mrSges(5,3) + (m(3) * t66 + t33 * mrSges(2,1) + t35 * mrSges(2,2) + t76 * (t27 * qJ(3) - t66) + t71 * t27 + (m(4) * pkin(2) - m(5) * (-pkin(2) - pkin(5)) + mrSges(5,3) + t75) * t26) * g(1) + (t23 * mrSges(4,3) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t31 * mrSges(3,1) - 0.2e1 * t30 * mrSges(3,2) + m(3) * (t30 ^ 2 + t31 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1); m(5) * (t1 * t34 - t2 * t32 + (-t3 * t34 - t32 * t4) * qJD(4)) - t17 * t55 + t34 * t6 - t18 * t54 - t32 * t5 + t67 * qJDD(2) + (-m(5) - t67) * g(3); t53 + (t34 * t17 - t32 * t18) * qJD(4) + (-mrSges(4,3) * qJD(1) - t12) * qJD(1) + (t37 + t38) * m(5) + (t15 + t38) * m(4) + t74; t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t13 + Ifges(5,6) * t14 + Ifges(5,3) * qJDD(4) + g(3) * t46 - t3 * t17 + t4 * t18 + (t32 * t8 / 0.2e1 + t7 * t68 - t70 * qJD(1) + t48 * mrSges(5,3) - t69) * qJD(1) + t72 * t47;];
tau = t19;
