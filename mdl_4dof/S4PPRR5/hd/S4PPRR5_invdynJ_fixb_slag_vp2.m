% Calculate vector of inverse dynamics joint torques for
% S4PPRR5
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
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:41
% EndTime: 2019-12-31 16:19:43
% DurationCPUTime: 0.85s
% Computational Cost: add. (348->118), mult. (788->174), div. (0->0), fcn. (429->6), ass. (0->60)
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t17 = -qJD(1) * t30 + qJD(2) * t32;
t88 = t17 * qJD(3);
t27 = sin(pkin(6));
t28 = cos(pkin(6));
t83 = -g(1) * t27 + g(2) * t28;
t29 = sin(qJ(4));
t78 = t29 / 0.2e1;
t87 = -m(4) - m(5);
t86 = g(3) * t32;
t6 = t32 * qJDD(1) + t30 * qJDD(2) + t88;
t85 = t6 - t88;
t18 = qJD(1) * t32 + qJD(2) * t30;
t59 = qJD(3) * t18;
t7 = -qJDD(1) * t30 + qJDD(2) * t32 - t59;
t51 = t7 + t59;
t58 = qJD(3) * t29;
t19 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t58;
t31 = cos(qJ(4));
t57 = qJD(3) * t31;
t20 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t57;
t41 = t29 * t19 - t31 * t20;
t53 = qJD(3) * qJD(4);
t15 = qJDD(3) * t31 - t29 * t53;
t16 = qJDD(3) * t29 + t31 * t53;
t84 = t31 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t15) - t29 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t16);
t13 = qJD(3) * pkin(5) + t18;
t4 = qJDD(3) * pkin(5) + t6;
t56 = qJD(4) * t29;
t1 = -t13 * t56 + t31 * t4;
t55 = qJD(4) * t31;
t2 = -t13 * t55 - t29 * t4;
t48 = t1 * t31 - t2 * t29;
t82 = (t29 ^ 2 + t31 ^ 2) * t13;
t12 = -qJD(3) * pkin(3) - t17;
t45 = mrSges(5,1) * t29 + mrSges(5,2) * t31;
t81 = t12 * t45 + qJD(4) * (Ifges(5,5) * t31 - Ifges(5,6) * t29) / 0.2e1;
t3 = -mrSges(5,1) * t15 + mrSges(5,2) * t16;
t33 = qJD(3) ^ 2;
t5 = -qJDD(3) * pkin(3) - t7;
t79 = m(4) * t51 + m(5) * (qJD(3) * t82 - t5) - t41 * qJD(3) + qJDD(3) * mrSges(4,1) - t33 * mrSges(4,2) - t3;
t76 = m(2) + m(3);
t69 = Ifges(5,4) * t29;
t68 = Ifges(5,2) * t31;
t66 = t29 * t30;
t65 = t30 * t31;
t64 = t31 * Ifges(5,4);
t52 = m(5) * pkin(5) + mrSges(5,3);
t46 = -t31 * mrSges(5,1) + t29 * mrSges(5,2);
t44 = t68 + t69;
t42 = t31 * t19 + t29 * t20;
t37 = t29 * (Ifges(5,1) * t31 - t69);
t36 = m(5) * pkin(3) - t46;
t14 = t46 * qJD(3);
t34 = -t33 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + qJD(3) * t14 - t42 * qJD(4) + m(4) * t85 + m(5) * (qJD(3) * t12 + t48) + t84;
t22 = Ifges(5,4) * t57;
t11 = Ifges(5,1) * t58 + Ifges(5,5) * qJD(4) + t22;
t10 = Ifges(5,6) * qJD(4) + t44 * qJD(3);
t8 = [t76 * qJDD(1) + (-t76 + t87) * g(3) - t79 * t30 + t34 * t32; m(3) * qJDD(2) + t79 * t32 + t34 * t30 + (m(3) - t87) * t83; t31 * (Ifges(5,4) * t16 + Ifges(5,2) * t15) / 0.2e1 + t16 * (Ifges(5,1) * t29 + t64) / 0.2e1 + t11 * t55 / 0.2e1 - t10 * t56 / 0.2e1 - g(3) * (-t36 * t30 + t52 * t32) + t15 * t44 / 0.2e1 + t5 * t46 + (Ifges(5,1) * t16 + Ifges(5,4) * t15) * t78 - pkin(3) * t3 - t18 * t14 + Ifges(4,3) * qJDD(3) + (t31 * (-Ifges(5,2) * t29 + t64) + t37) * t53 / 0.2e1 + t41 * t17 + t81 * qJD(4) + t48 * mrSges(5,3) + (-t85 + t86) * mrSges(4,2) + (g(3) * t30 + t51) * mrSges(4,1) + (-pkin(3) * t5 - t12 * t18 - t17 * t82) * m(5) + (0.2e1 * Ifges(5,5) * t78 + t31 * Ifges(5,6)) * qJDD(4) + (m(5) * t48 - t19 * t55 - t20 * t56 + t84) * pkin(5) + ((mrSges(4,1) + t36) * t32 + (-mrSges(4,2) + t52) * t30) * t83; Ifges(5,5) * t16 + Ifges(5,6) * t15 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - g(1) * ((-t27 * t66 + t28 * t31) * mrSges(5,1) + (-t27 * t65 - t28 * t29) * mrSges(5,2)) - g(2) * ((t27 * t31 + t28 * t66) * mrSges(5,1) + (-t27 * t29 + t28 * t65) * mrSges(5,2)) + t45 * t86 + t42 * t13 + (t10 * t78 + (-t37 / 0.2e1 + t68 * t78) * qJD(3) - (t11 + t22) * t31 / 0.2e1 - t81) * qJD(3);];
tau = t8;
