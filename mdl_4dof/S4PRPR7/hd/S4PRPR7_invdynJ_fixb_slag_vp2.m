% Calculate vector of inverse dynamics joint torques for
% S4PRPR7
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:25:37
% EndTime: 2019-12-31 16:25:39
% DurationCPUTime: 1.28s
% Computational Cost: add. (372->144), mult. (770->187), div. (0->0), fcn. (364->6), ass. (0->70)
t71 = mrSges(3,1) - mrSges(4,2);
t101 = -mrSges(5,3) - t71;
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t49 = t30 * mrSges(5,1) + t32 * mrSges(5,2);
t100 = -t49 + mrSges(3,2);
t28 = sin(pkin(6));
t29 = cos(pkin(6));
t94 = g(1) * t29 + g(2) * t28;
t99 = -t30 / 0.2e1;
t88 = t32 / 0.2e1;
t98 = -m(5) - m(4);
t67 = mrSges(5,3) * qJD(2);
t18 = -qJD(4) * mrSges(5,2) - t30 * t67;
t19 = qJD(4) * mrSges(5,1) - t32 * t67;
t96 = -t30 * t18 - t32 * t19;
t45 = t18 * t32 - t19 * t30;
t56 = qJD(2) * qJD(4);
t13 = qJDD(2) * t32 - t30 * t56;
t14 = -qJDD(2) * t30 - t32 * t56;
t95 = t32 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t13) + t30 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t14);
t34 = -pkin(2) - pkin(5);
t33 = cos(qJ(2));
t65 = qJD(1) * t33;
t53 = qJD(3) - t65;
t11 = qJD(2) * t34 + t53;
t31 = sin(qJ(2));
t57 = qJD(1) * qJD(2);
t23 = t31 * t57;
t15 = qJDD(1) * t33 - t23;
t43 = qJDD(3) - t15;
t4 = qJDD(2) * t34 + t43;
t63 = qJD(4) * t30;
t1 = -t11 * t63 + t32 * t4;
t62 = qJD(4) * t32;
t2 = t11 * t62 + t30 * t4;
t52 = t1 * t32 + t2 * t30;
t79 = Ifges(5,4) * t32;
t80 = Ifges(5,4) * t30;
t93 = (-Ifges(5,2) * t32 - t80) * t99 + (-Ifges(5,1) * t30 - t79) * t88;
t92 = (t30 ^ 2 + t32 ^ 2) * t11;
t12 = t49 * qJD(2);
t91 = -qJD(2) * t12 + t45 * qJD(4) + t95;
t66 = qJD(1) * t31;
t20 = qJD(2) * qJ(3) + t66;
t50 = mrSges(5,1) * t32 - mrSges(5,2) * t30;
t90 = t20 * t50 + qJD(4) * (-Ifges(5,5) * t30 - Ifges(5,6) * t32) / 0.2e1;
t55 = qJDD(2) * qJ(3);
t59 = qJDD(1) * t31;
t5 = t55 + t59 + (qJD(3) + t65) * qJD(2);
t89 = qJD(3) * t20 + (-t94 * t33 + t5) * qJ(3);
t85 = g(3) * t33;
t76 = t20 * t33;
t74 = t30 * t31;
t73 = t31 * t32;
t70 = -mrSges(3,2) + mrSges(4,3);
t60 = t20 * qJD(2);
t58 = qJDD(2) * mrSges(4,2);
t54 = t33 * t57;
t48 = Ifges(5,1) * t32 - t80;
t47 = -t30 * Ifges(5,2) + t79;
t37 = -t94 * t31 - t60 + t85;
t35 = qJD(2) ^ 2;
t17 = -qJD(2) * pkin(2) + t53;
t16 = t54 + t59;
t10 = Ifges(5,5) * qJD(4) + qJD(2) * t48;
t9 = Ifges(5,6) * qJD(4) + qJD(2) * t47;
t8 = -qJDD(2) * pkin(2) + t43;
t3 = -mrSges(5,1) * t14 + mrSges(5,2) * t13;
t6 = [m(2) * qJDD(1) + (-m(2) - m(3) + t98) * g(3) + (t3 - t71 * t35 + t70 * qJDD(2) - t96 * qJD(2) + m(3) * t16 + m(4) * (qJD(2) * t17 + t5) + m(5) * (qJD(2) * t92 + t5)) * t31 + (t70 * t35 + t71 * qJDD(2) + m(3) * t15 + m(4) * (-t8 + t60) + m(5) * (-t52 + t60) - t91) * t33; (t34 * t45 + t90) * qJD(4) - pkin(2) * t58 - t9 * t62 / 0.2e1 - t10 * t63 / 0.2e1 + (-t23 + t8) * mrSges(4,2) + (t54 - t16) * mrSges(3,2) + t14 * t47 / 0.2e1 + t13 * t48 / 0.2e1 + t5 * t49 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + ((-mrSges(4,3) + t100) * t33 + (m(4) * pkin(2) - m(5) * t34 - t101) * t31) * t94 + (t98 * (t33 * pkin(2) + t31 * qJ(3)) + t100 * t31 + (-m(5) * pkin(5) + t101) * t33) * g(3) + (t23 + t15) * mrSges(3,1) + (Ifges(5,1) * t13 + Ifges(5,4) * t14) * t88 + (-g(3) * t31 + qJD(2) * qJD(3) + t5 - t54 + t55) * mrSges(4,3) + (0.2e1 * Ifges(5,5) * t88 - Ifges(5,6) * t30) * qJDD(4) + t96 * t66 + t53 * t12 - t52 * mrSges(5,3) + t95 * t34 + (t52 * t34 - (t31 * t92 + t76) * qJD(1) + t89) * m(5) + t93 * t56 + (-pkin(2) * t8 - (t17 * t31 + t76) * qJD(1) + t89) * m(4) + qJ(3) * t3 + (Ifges(5,4) * t13 + Ifges(5,2) * t14) * t99; t58 - t35 * mrSges(4,3) + (t37 + t52) * m(5) + (t8 + t37) * m(4) + t91; Ifges(5,5) * t13 + Ifges(5,6) * t14 + Ifges(5,3) * qJDD(4) - t2 * mrSges(5,2) + t1 * mrSges(5,1) - g(1) * ((-t28 * t30 + t29 * t73) * mrSges(5,1) + (-t28 * t32 - t29 * t74) * mrSges(5,2)) - g(2) * ((t28 * t73 + t29 * t30) * mrSges(5,1) + (-t28 * t74 + t29 * t32) * mrSges(5,2)) + t50 * t85 - t45 * t11 + (t30 * t10 / 0.2e1 + t9 * t88 - t93 * qJD(2) - t90) * qJD(2);];
tau = t6;
