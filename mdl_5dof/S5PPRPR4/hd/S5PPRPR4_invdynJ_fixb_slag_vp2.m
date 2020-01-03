% Calculate vector of inverse dynamics joint torques for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:17
% EndTime: 2019-12-31 17:32:20
% DurationCPUTime: 1.84s
% Computational Cost: add. (894->199), mult. (2017->272), div. (0->0), fcn. (1376->10), ass. (0->92)
t133 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t78 = cos(qJ(3));
t93 = -qJD(2) * t78 + qJD(4);
t72 = sin(pkin(8));
t73 = cos(pkin(8));
t132 = t72 ^ 2 + t73 ^ 2;
t131 = m(3) + m(4);
t130 = -m(5) - m(6);
t75 = sin(qJ(5));
t77 = cos(qJ(5));
t46 = t72 * t75 - t77 * t73;
t39 = t46 * qJD(3);
t48 = t72 * t77 + t73 * t75;
t16 = -qJD(5) * t39 + qJDD(3) * t48;
t124 = t48 * qJD(5);
t17 = -qJD(3) * t124 - qJDD(3) * t46;
t3 = -t17 * mrSges(6,1) + mrSges(6,2) * t16;
t102 = qJDD(3) * t73;
t103 = qJDD(3) * t72;
t44 = -mrSges(5,1) * t102 + mrSges(5,2) * t103;
t129 = t3 + t44;
t112 = pkin(6) + qJ(4);
t53 = t112 * t72;
t54 = t112 * t73;
t22 = -t53 * t77 - t54 * t75;
t128 = qJD(5) * t22 - t46 * t93;
t23 = -t53 * t75 + t54 * t77;
t84 = t48 * t78;
t127 = qJD(2) * t84 - qJD(4) * t48 - qJD(5) * t23;
t76 = sin(qJ(3));
t34 = t46 * t76;
t40 = t48 * qJD(3);
t90 = -mrSges(5,1) * t73 + mrSges(5,2) * t72;
t126 = mrSges(6,1) * t39 + mrSges(6,2) * t40 + t90 * qJD(3);
t125 = mrSges(5,3) * t132;
t104 = qJDD(1) * t73;
t100 = qJD(2) * qJD(3);
t65 = t78 * t100;
t52 = qJDD(2) * t76 + t65;
t32 = t52 + t133;
t20 = -t32 * t72 - t104;
t21 = -qJDD(1) * t72 + t32 * t73;
t87 = -t20 * t72 + t21 * t73;
t63 = pkin(4) * t73 + pkin(3);
t71 = pkin(8) + qJ(5);
t66 = sin(t71);
t67 = cos(t71);
t89 = -mrSges(6,1) * t67 + mrSges(6,2) * t66;
t123 = m(5) * pkin(3) + m(6) * t63 + mrSges(4,1) - t89 - t90;
t122 = m(5) * qJ(4) + m(6) * t112 - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t118 = t40 / 0.2e1;
t117 = Ifges(6,4) * t40;
t110 = pkin(6) * qJD(3);
t109 = cos(pkin(7));
t108 = sin(pkin(7));
t107 = qJD(1) * t73;
t106 = qJD(2) * t76;
t101 = m(2) + t131;
t64 = t76 * t100;
t56 = qJD(3) * qJ(4) + t106;
t36 = -qJD(1) * t72 + t56 * t73;
t51 = qJDD(2) * t78 - t64;
t45 = -t108 * t76 - t109 * t78;
t47 = -t108 * t78 + t109 * t76;
t91 = -g(1) * t47 + g(2) * t45;
t24 = -t107 + (-t56 - t110) * t72;
t25 = t110 * t73 + t36;
t4 = t24 * t77 - t25 * t75;
t5 = t24 * t75 + t25 * t77;
t86 = -(-t56 * t72 - t107) * t72 + t36 * t73;
t85 = qJDD(4) - t51;
t80 = (-qJD(3) * pkin(3) + t93) * t76 + t78 * t86;
t79 = qJD(3) ^ 2;
t43 = -qJD(3) * t63 + t93;
t41 = t46 * qJD(5);
t38 = -qJDD(3) * pkin(3) + t85;
t37 = Ifges(6,4) * t39;
t33 = t48 * t76;
t29 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t40;
t28 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t39;
t26 = -qJDD(3) * t63 + t85;
t19 = pkin(6) * t102 + t21;
t18 = -t104 + (-pkin(6) * qJDD(3) - t32) * t72;
t13 = Ifges(6,1) * t40 + Ifges(6,5) * qJD(5) - t37;
t12 = -Ifges(6,2) * t39 + Ifges(6,6) * qJD(5) + t117;
t11 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t17;
t10 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t16;
t9 = -qJD(3) * t84 + qJD(5) * t34;
t8 = -t124 * t76 - t39 * t78;
t2 = -qJD(5) * t5 + t18 * t77 - t19 * t75;
t1 = qJD(5) * t4 + t18 * t75 + t19 * t77;
t6 = [t46 * t10 - t48 * t11 + t41 * t28 + t124 * t29 + m(5) * (-t20 * t73 - t21 * t72) + m(6) * (-t1 * t48 + t124 * t4 + t2 * t46 + t41 * t5) + t101 * qJDD(1) + (-t101 + t130) * g(3); m(3) * qJDD(2) - t33 * t10 - t34 * t11 + t8 * t28 + t9 * t29 + (mrSges(4,1) * qJDD(3) - mrSges(4,2) * t79 - t129) * t78 + (-t79 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + qJD(3) * t126) * t76 + m(6) * (qJD(3) * t43 * t76 - t1 * t34 - t2 * t33 - t26 * t78 + t4 * t9 + t5 * t8) + m(4) * (t51 * t78 + t52 * t76) + m(5) * (qJD(3) * t80 - t38 * t78 + t76 * t87) + (qJDD(3) * t76 + t78 * t79) * t125 + (-g(1) * t108 + g(2) * t109) * (-t130 + t131); -t126 * t106 + (t87 + t132 * (-t65 + t133)) * mrSges(5,3) + (t122 * t45 + t123 * t47) * g(1) + (t122 * t47 - t123 * t45) * g(2) + (-t124 * t5 + t4 * t41) * mrSges(6,3) - t39 * (-Ifges(6,4) * t41 - Ifges(6,2) * t124) / 0.2e1 + qJD(5) * (-Ifges(6,5) * t41 - Ifges(6,6) * t124) / 0.2e1 - t124 * t12 / 0.2e1 + t43 * (mrSges(6,1) * t124 - mrSges(6,2) * t41) + (-Ifges(6,1) * t41 - Ifges(6,4) * t124) * t118 + t38 * t90 + (mrSges(6,2) * t26 - mrSges(6,3) * t2 + Ifges(6,1) * t16 + Ifges(6,4) * t17 + Ifges(6,5) * qJDD(5)) * t48 + (mrSges(6,1) * t26 - mrSges(6,3) * t1 - Ifges(6,4) * t16 - Ifges(6,2) * t17 - Ifges(6,6) * qJDD(5)) * t46 - t63 * t3 - pkin(3) * t44 - t41 * t13 / 0.2e1 + t22 * t10 + t23 * t11 + (-pkin(3) * t38 + qJ(4) * t87 - qJD(2) * t80 + qJD(4) * t86) * m(5) + (t64 + t51) * mrSges(4,1) + (t65 - t52) * mrSges(4,2) + (Ifges(5,4) * t72 + Ifges(5,2) * t73) * t102 + (Ifges(5,1) * t72 + Ifges(5,4) * t73) * t103 + t127 * t29 + t128 * t28 + (t1 * t23 - t106 * t43 + t127 * t4 + t128 * t5 + t2 * t22 - t26 * t63) * m(6) + Ifges(4,3) * qJDD(3); -t79 * t125 + t39 * t28 + t40 * t29 + (t39 * t5 + t4 * t40 + t26 + t91) * m(6) + (-qJD(3) * t86 + t38 + t91) * m(5) + t129; Ifges(6,5) * t16 + Ifges(6,6) * t17 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t43 * (mrSges(6,1) * t40 - mrSges(6,2) * t39) - t40 * (-Ifges(6,1) * t39 - t117) / 0.2e1 + t12 * t118 - qJD(5) * (-Ifges(6,5) * t39 - Ifges(6,6) * t40) / 0.2e1 - t4 * t28 + t5 * t29 - g(3) * t89 + (-t39 * t4 + t40 * t5) * mrSges(6,3) + (-g(1) * t45 - g(2) * t47) * (mrSges(6,1) * t66 + mrSges(6,2) * t67) + (-Ifges(6,2) * t40 + t13 - t37) * t39 / 0.2e1;];
tau = t6;
