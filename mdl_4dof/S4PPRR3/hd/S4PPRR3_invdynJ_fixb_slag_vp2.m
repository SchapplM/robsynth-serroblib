% Calculate vector of inverse dynamics joint torques for
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:23
% DurationCPUTime: 0.67s
% Computational Cost: add. (310->125), mult. (716->182), div. (0->0), fcn. (398->6), ass. (0->61)
t33 = sin(qJ(3));
t61 = qJD(2) * t33;
t22 = qJD(3) * pkin(5) + t61;
t32 = sin(qJ(4));
t34 = cos(qJ(4));
t41 = t32 * qJD(1) - t22 * t34;
t35 = cos(qJ(3));
t54 = qJD(2) * qJD(3);
t27 = t35 * t54;
t18 = t33 * qJDD(2) + t27;
t11 = qJDD(3) * pkin(5) + t18;
t6 = -qJD(1) * t34 - t22 * t32;
t1 = qJD(4) * t6 - t32 * qJDD(1) + t34 * t11;
t2 = qJD(4) * t41 - qJDD(1) * t34 - t11 * t32;
t47 = t1 * t34 - t2 * t32;
t56 = qJD(4) * t34;
t57 = qJD(4) * t32;
t82 = t41 * t57 - t6 * t56 + t47;
t31 = cos(pkin(6));
t62 = sin(pkin(6));
t12 = -t31 * t35 - t33 * t62;
t13 = t31 * t33 - t35 * t62;
t81 = -g(1) * t12 - g(2) * t13;
t73 = t32 / 0.2e1;
t80 = m(3) + m(4);
t59 = qJD(3) * t32;
t19 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t59;
t58 = qJD(3) * t34;
t20 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t58;
t76 = -t32 * t19 + t34 * t20;
t53 = qJD(3) * qJD(4);
t15 = qJDD(3) * t34 - t32 * t53;
t4 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t15;
t16 = qJDD(3) * t32 + t34 * t53;
t5 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t16;
t75 = -t32 * t5 + t34 * t4;
t60 = qJD(2) * t35;
t23 = -qJD(3) * pkin(3) - t60;
t44 = mrSges(5,1) * t32 + mrSges(5,2) * t34;
t74 = t23 * t44 + qJD(4) * (Ifges(5,5) * t34 - Ifges(5,6) * t32) / 0.2e1;
t67 = Ifges(5,4) * t32;
t66 = Ifges(5,4) * t34;
t65 = Ifges(5,2) * t34;
t55 = m(2) + t80;
t50 = -m(5) * pkin(5) - mrSges(5,3);
t49 = t33 * t54;
t46 = -t32 * t41 + t34 * t6;
t45 = t32 * t6 + t34 * t41;
t21 = -mrSges(5,1) * t34 + mrSges(5,2) * t32;
t43 = t65 + t67;
t38 = t32 * (Ifges(5,1) * t34 - t67);
t17 = qJDD(2) * t35 - t49;
t37 = m(5) * pkin(3) - t21;
t36 = qJD(3) ^ 2;
t28 = Ifges(5,4) * t58;
t14 = t21 * qJD(3);
t10 = -qJDD(3) * pkin(3) - t17;
t9 = Ifges(5,1) * t59 + Ifges(5,5) * qJD(4) + t28;
t8 = Ifges(5,6) * qJD(4) + qJD(3) * t43;
t3 = -mrSges(5,1) * t15 + mrSges(5,2) * t16;
t7 = [-t20 * t56 + t19 * t57 - t32 * t4 - t34 * t5 + m(5) * (qJD(4) * t45 - t1 * t32 - t2 * t34) + t55 * qJDD(1) + (-m(5) - t55) * g(3); m(3) * qJDD(2) + (qJDD(3) * mrSges(4,1) - t36 * mrSges(4,2) - t3 + t76 * qJD(3) + m(4) * t17 + m(5) * (-t41 * t58 - t59 * t6 - t10)) * t35 + (-t36 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + qJD(3) * t14 + (-t34 * t19 - t32 * t20) * qJD(4) + m(4) * t18 + m(5) * (qJD(3) * t23 + t82) + t75) * t33 + (m(5) + t80) * (-g(1) * t62 + g(2) * t31); t16 * (Ifges(5,1) * t32 + t66) / 0.2e1 - t14 * t61 + t9 * t56 / 0.2e1 - t8 * t57 / 0.2e1 - g(2) * (t12 * t37 + t13 * t50) - g(1) * (t12 * t50 - t13 * t37) + t15 * t43 / 0.2e1 + (Ifges(5,1) * t16 + Ifges(5,4) * t15) * t73 + t34 * (Ifges(5,4) * t16 + Ifges(5,2) * t15) / 0.2e1 + t10 * t21 - pkin(3) * t3 + Ifges(4,3) * qJDD(3) - t76 * t60 + (t34 * (-Ifges(5,2) * t32 + t66) + t38) * t53 / 0.2e1 + t74 * qJD(4) + (-t18 + t27 + t81) * mrSges(4,2) + (g(1) * t13 - g(2) * t12 + t17 + t49) * mrSges(4,1) + (-pkin(3) * t10 - (t23 * t33 - t35 * t45) * qJD(2)) * m(5) + (0.2e1 * Ifges(5,5) * t73 + Ifges(5,6) * t34) * qJDD(4) + (m(5) * (-qJD(4) * t46 + t47) - t19 * t56 - t20 * t57 + t75) * pkin(5) + t82 * mrSges(5,3); t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t16 + Ifges(5,6) * t15 + Ifges(5,3) * qJDD(4) - g(3) * t21 - t41 * t19 - t6 * t20 + (t8 * t73 + (-t38 / 0.2e1 + t65 * t73) * qJD(3) + t46 * mrSges(5,3) - (t9 + t28) * t34 / 0.2e1 - t74) * qJD(3) + t81 * t44;];
tau = t7;
