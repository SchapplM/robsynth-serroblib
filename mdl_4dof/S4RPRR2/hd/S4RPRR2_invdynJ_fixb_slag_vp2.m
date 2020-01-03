% Calculate vector of inverse dynamics joint torques for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:09
% DurationCPUTime: 0.78s
% Computational Cost: add. (868->157), mult. (1606->218), div. (0->0), fcn. (818->12), ass. (0->82)
t74 = sin(qJ(4));
t119 = t74 / 0.2e1;
t125 = mrSges(4,2) - mrSges(5,3);
t77 = cos(qJ(4));
t49 = -t77 * mrSges(5,1) + t74 * mrSges(5,2);
t73 = cos(pkin(7));
t63 = pkin(1) * t73 + pkin(2);
t48 = t63 * qJD(1);
t75 = sin(qJ(3));
t108 = t48 * t75;
t78 = cos(qJ(3));
t72 = sin(pkin(7));
t116 = pkin(1) * t72;
t97 = qJD(1) * t116;
t23 = t78 * t97 + t108;
t70 = qJD(1) + qJD(3);
t19 = pkin(6) * t70 + t23;
t14 = qJD(2) * t77 - t19 * t74;
t69 = qJDD(1) + qJDD(3);
t101 = qJD(3) * t78;
t47 = t63 * qJDD(1);
t59 = t78 * t116;
t95 = t75 * t97;
t8 = -qJD(3) * t95 + qJDD(1) * t59 + t48 * t101 + t75 * t47;
t5 = pkin(6) * t69 + t8;
t2 = qJD(4) * t14 + qJDD(2) * t74 + t5 * t77;
t15 = qJD(2) * t74 + t19 * t77;
t3 = -qJD(4) * t15 + qJDD(2) * t77 - t5 * t74;
t124 = t2 * t77 - t3 * t74;
t123 = -mrSges(4,1) + t49;
t22 = t48 * t78 - t95;
t18 = -pkin(3) * t70 - t22;
t94 = mrSges(5,1) * t74 + mrSges(5,2) * t77;
t122 = t18 * t94 + qJD(4) * (Ifges(5,5) * t77 - Ifges(5,6) * t74) / 0.2e1 - (t14 * t77 + t15 * t74) * mrSges(5,3);
t121 = m(3) * pkin(1);
t35 = t75 * t63 + t59;
t90 = -t14 * t74 + t15 * t77;
t9 = -(qJD(1) * t101 + qJDD(1) * t75) * t116 - qJD(3) * t108 + t47 * t78;
t120 = m(5) * pkin(3);
t117 = m(4) + m(5);
t113 = Ifges(5,4) * t74;
t112 = Ifges(5,4) * t77;
t111 = Ifges(5,2) * t77;
t107 = t70 * t74;
t106 = t70 * t77;
t71 = qJ(1) + pkin(7);
t67 = qJ(3) + t71;
t61 = sin(t67);
t62 = cos(t67);
t103 = t62 * pkin(3) + t61 * pkin(6);
t66 = cos(t71);
t79 = cos(qJ(1));
t102 = t79 * pkin(1) + pkin(2) * t66;
t100 = qJD(4) * t74;
t99 = qJD(4) * t77;
t98 = m(3) + t117;
t93 = t111 + t113;
t43 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t107;
t44 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t106;
t89 = -t74 * t43 + t77 * t44;
t34 = -t75 * t116 + t63 * t78;
t86 = t74 * (Ifges(5,1) * t77 - t113);
t84 = t123 * t62 + t125 * t61;
t83 = (t120 - t123) * t61 + (-m(5) * pkin(6) + t125) * t62;
t30 = Ifges(5,6) * qJD(4) + t93 * t70;
t52 = Ifges(5,4) * t106;
t31 = Ifges(5,1) * t107 + Ifges(5,5) * qJD(4) + t52;
t37 = -t70 * t100 + t69 * t77;
t38 = t69 * t74 + t70 * t99;
t6 = -pkin(3) * t69 - t9;
t82 = (Ifges(5,1) * t38 + Ifges(5,4) * t37) * t119 + t77 * (Ifges(5,4) * t38 + Ifges(5,2) * t37) / 0.2e1 + t37 * t93 / 0.2e1 + t38 * (Ifges(5,1) * t74 + t112) / 0.2e1 - t30 * t100 / 0.2e1 + t6 * t49 + Ifges(4,3) * t69 + t9 * mrSges(4,1) + (t31 + t70 * (-Ifges(5,2) * t74 + t112)) * t99 / 0.2e1 + t124 * mrSges(5,3) + (0.2e1 * Ifges(5,5) * t119 + Ifges(5,6) * t77) * qJDD(4) + (t86 * t70 / 0.2e1 + t122) * qJD(4);
t24 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t37;
t25 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t38;
t81 = -t44 * t100 - t43 * t99 - t74 * t25 + t77 * t24 + m(5) * (-t15 * t100 - t14 * t99 + t124);
t76 = sin(qJ(1));
t65 = sin(t71);
t36 = t49 * t70;
t32 = -pkin(3) - t34;
t29 = t35 * qJD(3);
t28 = t34 * qJD(3);
t13 = -mrSges(5,1) * t37 + mrSges(5,2) * t38;
t1 = [(t79 * mrSges(2,2) + t66 * mrSges(3,2) + (t117 * pkin(2) + mrSges(3,1)) * t65 + (t98 * pkin(1) + mrSges(2,1)) * t76 + t83) * g(1) + m(4) * (-t22 * t29 + t23 * t28 + t34 * t9 + t35 * t8) + (t76 * mrSges(2,2) - t66 * mrSges(3,1) + t65 * mrSges(3,2) - m(5) * (t102 + t103) - m(4) * t102 + (-mrSges(2,1) - t121) * t79 + t84) * g(2) + t89 * t28 + t81 * (pkin(6) + t35) + t82 + (-t28 * t70 - t35 * t69 - t8) * mrSges(4,2) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t73 * mrSges(3,1) - 0.2e1 * t72 * mrSges(3,2) + (t72 ^ 2 + t73 ^ 2) * t121) * pkin(1)) * qJDD(1) + (-t29 * t70 + t34 * t69) * mrSges(4,1) + m(5) * (t18 * t29 + t90 * t28 + t6 * t32) + t29 * t36 + t32 * t13; m(5) * (t90 * qJD(4) + t2 * t74 + t3 * t77) + t44 * t99 + t74 * t24 - t43 * t100 + t77 * t25 + (m(3) + m(4)) * qJDD(2) - t98 * g(3); t81 * pkin(6) + t82 + t83 * g(1) + (t70 * mrSges(4,1) - t36) * t23 + (t70 * mrSges(4,2) - t89) * t22 - t6 * t120 + (-m(5) * t103 + t84) * g(2) - m(5) * (t18 * t23 + t90 * t22) - pkin(3) * t13 - t8 * mrSges(4,2); t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t38 + Ifges(5,6) * t37 + Ifges(5,3) * qJDD(4) + g(3) * t49 - t14 * t44 + t15 * t43 + (t30 * t119 + (-t86 / 0.2e1 + t111 * t119) * t70 - (t31 + t52) * t77 / 0.2e1 - t122) * t70 + (g(1) * t62 + g(2) * t61) * t94;];
tau = t1;
