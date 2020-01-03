% Calculate vector of inverse dynamics joint torques for
% S4RPPR6
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:42
% DurationCPUTime: 2.90s
% Computational Cost: add. (612->190), mult. (1510->259), div. (0->0), fcn. (909->6), ass. (0->90)
t62 = sin(pkin(6));
t60 = t62 ^ 2;
t63 = cos(pkin(6));
t61 = t63 ^ 2;
t118 = (mrSges(4,2) + mrSges(3,3)) * (t60 + t61);
t103 = m(4) + m(5);
t115 = m(3) + t103;
t114 = Ifges(3,4) - Ifges(4,5);
t82 = qJDD(1) * qJ(2);
t83 = qJD(1) * qJD(2);
t46 = t82 + t83;
t113 = t46 + t83;
t99 = t62 * mrSges(4,3);
t73 = t63 * mrSges(4,1) + t99;
t74 = mrSges(3,1) * t63 - mrSges(3,2) * t62;
t111 = -t74 - t73 - mrSges(2,1);
t110 = m(5) * pkin(5) + mrSges(2,2) + mrSges(5,3);
t91 = t62 * qJ(3);
t77 = pkin(1) + t91;
t33 = (pkin(2) + pkin(3)) * t63 + t77;
t18 = qJD(1) * t33 - qJD(2);
t66 = cos(qJ(4));
t88 = qJD(1) * t66;
t64 = sin(qJ(4));
t89 = qJD(1) * t64;
t31 = -t62 * t89 - t63 * t88;
t32 = t62 * t88 - t63 * t89;
t109 = m(5) * t18 - mrSges(5,1) * t31 + mrSges(5,2) * t32 + t73 * qJD(1);
t65 = sin(qJ(1));
t67 = cos(qJ(1));
t108 = -g(1) * t67 - g(2) * t65;
t107 = t60 * t82 + t108 + (t46 + t82) * t61;
t42 = -pkin(2) * t63 - t77;
t106 = t42 * qJDD(1);
t104 = t32 / 0.2e1;
t101 = Ifges(5,4) * t32;
t98 = t67 * t63;
t97 = -pkin(5) + qJ(2);
t92 = t61 * qJ(2);
t96 = t113 * t92;
t95 = t67 * pkin(1) + t65 * qJ(2);
t93 = t60 * qJ(2);
t90 = qJD(1) * t62;
t87 = qJD(3) * t62;
t45 = qJ(2) * t90 + qJD(3);
t86 = qJDD(1) * t62;
t35 = t62 * t46 + qJDD(3);
t44 = t97 * t63;
t76 = pkin(2) * t98 + t67 * t91 + t95;
t70 = t62 * t64 + t63 * t66;
t29 = t70 * qJD(4);
t37 = t62 * t66 - t63 * t64;
t10 = -qJD(1) * t29 + qJDD(1) * t37;
t30 = t37 * qJD(4);
t11 = -qJD(1) * t30 - qJDD(1) * t70;
t72 = -t11 * mrSges(5,1) + t10 * mrSges(5,2);
t71 = -mrSges(5,1) * t70 - mrSges(5,2) * t37;
t34 = -pkin(5) * t90 + t45;
t39 = qJD(1) * t44;
t12 = t34 * t66 - t39 * t64;
t13 = t34 * t64 + t39 * t66;
t43 = t97 * t62;
t15 = t43 * t66 - t44 * t64;
t16 = t43 * t64 + t44 * t66;
t69 = -qJD(1) * t87 + qJDD(2);
t68 = qJD(1) ^ 2;
t56 = -qJDD(1) * pkin(1) + qJDD(2);
t55 = t68 * t92;
t53 = mrSges(3,2) * t86;
t28 = qJD(1) * t42 + qJD(2);
t27 = Ifges(5,4) * t31;
t26 = (-pkin(5) * qJDD(1) + t46) * t63;
t25 = t70 * t67;
t24 = t37 * t67;
t23 = t70 * t65;
t22 = t37 * t65;
t21 = -pkin(5) * t86 + t35;
t20 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t32;
t19 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t31;
t17 = t69 + t106;
t14 = qJDD(1) * t33 - t69;
t8 = Ifges(5,1) * t32 + Ifges(5,5) * qJD(4) + t27;
t7 = Ifges(5,2) * t31 + Ifges(5,6) * qJD(4) + t101;
t6 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t11;
t5 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t10;
t4 = qJD(2) * t37 - qJD(4) * t16;
t3 = qJD(2) * t70 + qJD(4) * t15;
t2 = -qJD(4) * t13 + t21 * t66 - t26 * t64;
t1 = qJD(4) * t12 + t21 * t64 + t26 * t66;
t9 = [(t12 * t29 - t13 * t30) * mrSges(5,3) + t18 * (mrSges(5,1) * t30 - mrSges(5,2) * t29) + qJD(4) * (-Ifges(5,5) * t29 - Ifges(5,6) * t30) / 0.2e1 + t31 * (-Ifges(5,4) * t29 - Ifges(5,2) * t30) / 0.2e1 + (-Ifges(5,1) * t29 - Ifges(5,4) * t30) * t104 - t56 * t74 - t14 * t71 + t33 * t72 + (t60 * t46 + t107) * mrSges(3,3) + (t35 * t62 + t107) * mrSges(4,2) + t109 * t87 + (-m(3) * t95 - m(4) * t76 - m(5) * (pkin(3) * t98 + t76) - t25 * mrSges(5,1) - t24 * mrSges(5,2) + t110 * t65 + t111 * t67) * g(2) + m(3) * (-pkin(1) * t56 + t113 * t93 + t96) + (t114 * t63 + (Ifges(3,1) + Ifges(4,1)) * t62) * t86 + (t23 * mrSges(5,1) + t22 * mrSges(5,2) + (m(3) * pkin(1) - m(4) * t42 + m(5) * t33 - t111) * t65 + (-t115 * qJ(2) + t110) * t67) * g(1) + ((pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t63 + t114 * t62) * t63 + Ifges(2,3)) * qJDD(1) + (-mrSges(5,3) * t2 + Ifges(5,1) * t10 + Ifges(5,4) * t11 + Ifges(5,5) * qJDD(4)) * t37 + (-mrSges(5,3) * t1 - Ifges(5,4) * t10 - Ifges(5,2) * t11 - Ifges(5,6) * qJDD(4)) * t70 + t83 * t118 + m(5) * (t1 * t16 + t12 * t4 + t13 * t3 + t14 * t33 + t15 * t2) + m(4) * (t17 * t42 + (qJ(2) * t35 + qJD(2) * t45 - qJD(3) * t28) * t62 + t96) - pkin(1) * t53 - t29 * t8 / 0.2e1 - t30 * t7 / 0.2e1 + t15 * t5 + t16 * t6 + t3 * t19 + t4 * t20 + (-t106 - t17) * t73; t31 * t19 - t32 * t20 + t53 + (-t99 + (-mrSges(3,1) - mrSges(4,1)) * t63) * qJDD(1) - t72 - t68 * t118 + (-g(1) * t65 + g(2) * t67) * t115 + (-t12 * t32 + t13 * t31 - t14) * m(5) + (-t45 * t90 + t17 - t55) * m(4) + (-t68 * t93 - t55 + t56) * m(3); t66 * t5 + t64 * t6 + (t19 * t66 - t20 * t64) * qJD(4) + t103 * t63 * g(3) + m(4) * t35 + m(5) * (t1 * t64 + t2 * t66 + (-t12 * t64 + t13 * t66) * qJD(4)) + (qJDD(1) * mrSges(4,2) + (m(4) * t28 - t109) * qJD(1) + t103 * t108) * t62; Ifges(5,5) * t10 + Ifges(5,6) * t11 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t18 * (mrSges(5,1) * t32 + mrSges(5,2) * t31) - t32 * (Ifges(5,1) * t31 - t101) / 0.2e1 + t7 * t104 - qJD(4) * (Ifges(5,5) * t31 - Ifges(5,6) * t32) / 0.2e1 - t12 * t19 + t13 * t20 - g(1) * (mrSges(5,1) * t24 - mrSges(5,2) * t25) - g(2) * (mrSges(5,1) * t22 - mrSges(5,2) * t23) - g(3) * t71 + (t12 * t31 + t13 * t32) * mrSges(5,3) - (-Ifges(5,2) * t32 + t27 + t8) * t31 / 0.2e1;];
tau = t9;
