% Calculate vector of inverse dynamics joint torques for
% S4PRRP4
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
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:45
% EndTime: 2019-12-31 16:27:49
% DurationCPUTime: 1.88s
% Computational Cost: add. (477->190), mult. (1021->242), div. (0->0), fcn. (457->4), ass. (0->80)
t112 = Ifges(5,4) + Ifges(4,5);
t111 = Ifges(5,6) - Ifges(4,6);
t42 = pkin(6) + qJ(2);
t38 = sin(t42);
t39 = cos(t42);
t105 = g(1) * t39 + g(2) * t38;
t118 = -m(4) - m(5);
t43 = sin(qJ(3));
t77 = qJD(2) * qJD(3);
t67 = t43 * t77;
t44 = cos(qJ(3));
t79 = t44 * qJDD(2);
t21 = t67 - t79;
t9 = -mrSges(5,2) * t21 + qJDD(3) * mrSges(5,3);
t117 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t21 + t9;
t22 = t43 * qJDD(2) + t44 * t77;
t8 = -qJDD(3) * mrSges(5,1) + t22 * mrSges(5,2);
t116 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t22 - t8;
t85 = qJD(2) * t44;
t24 = pkin(5) * t85 + t43 * qJD(1);
t17 = qJD(3) * qJ(4) + t24;
t78 = qJD(1) * qJD(3);
t5 = -pkin(5) * t22 + qJDD(1) * t44 - t43 * t78;
t3 = -qJDD(3) * pkin(3) + qJDD(4) - t5;
t115 = -t17 * qJD(3) + t3;
t61 = t44 * mrSges(4,1) - t43 * mrSges(4,2);
t114 = -mrSges(3,1) - t61;
t113 = mrSges(3,2) - mrSges(5,2);
t37 = Ifges(4,4) * t85;
t90 = Ifges(5,5) * t44;
t57 = Ifges(5,1) * t43 - t90;
t86 = qJD(2) * t43;
t110 = Ifges(4,1) * t86 + qJD(2) * t57 + t112 * qJD(3) + t37;
t70 = mrSges(4,3) * t86;
t71 = mrSges(5,2) * t86;
t109 = t70 + t71 + (-mrSges(4,1) - mrSges(5,1)) * qJD(3);
t69 = mrSges(5,2) * t85;
t29 = qJD(3) * mrSges(5,3) + t69;
t68 = mrSges(4,3) * t85;
t108 = -qJD(3) * mrSges(4,2) + t29 + t68;
t107 = t111 * t43 + t112 * t44;
t75 = pkin(5) * t79 + t43 * qJDD(1) + t44 * t78;
t4 = -pkin(5) * t67 + t75;
t106 = t4 * t44 - t5 * t43;
t102 = t43 / 0.2e1;
t101 = m(2) + m(3);
t72 = pkin(5) * t86;
t2 = qJDD(3) * qJ(4) + (qJD(4) - t72) * qJD(3) + t75;
t96 = t2 * t44;
t93 = Ifges(4,4) * t43;
t92 = Ifges(4,4) * t44;
t91 = Ifges(5,5) * t43;
t84 = qJD(3) * t43;
t83 = qJD(3) * t44;
t82 = qJDD(2) * pkin(2);
t80 = t43 * qJD(4);
t63 = -t77 / 0.2e1;
t60 = mrSges(4,1) * t43 + mrSges(4,2) * t44;
t59 = t44 * mrSges(5,1) + t43 * mrSges(5,3);
t58 = t43 * mrSges(5,1) - t44 * mrSges(5,3);
t56 = Ifges(4,2) * t44 + t93;
t53 = pkin(3) * t44 + t43 * qJ(4);
t52 = pkin(3) * t43 - qJ(4) * t44;
t25 = -pkin(2) - t53;
t51 = pkin(2) * t60;
t23 = qJD(1) * t44 - t72;
t16 = t25 * qJD(2);
t50 = t16 * t58;
t49 = t43 * (Ifges(4,1) * t44 - t93);
t48 = t44 * (Ifges(5,3) * t43 + t90);
t46 = m(5) * t53 + t59;
t36 = Ifges(5,5) * t86;
t20 = t52 * qJD(2);
t19 = t59 * qJD(2);
t13 = Ifges(4,6) * qJD(3) + qJD(2) * t56;
t12 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t85 + t36;
t11 = -qJD(3) * pkin(3) + qJD(4) - t23;
t10 = qJD(3) * t52 - t80;
t1 = pkin(3) * t21 - qJ(4) * t22 - qJD(2) * t80 - t82;
t6 = [t116 * t44 + t117 * t43 + t101 * qJDD(1) + (t108 * t44 + t109 * t43) * qJD(3) + m(4) * (t4 * t43 + t44 * t5 + (-t23 * t43 + t24 * t44) * qJD(3)) + m(5) * (t2 * t43 - t3 * t44 + (t11 * t43 + t17 * t44) * qJD(3)) + (-t101 + t118) * g(3); (t91 / 0.2e1 - t56 / 0.2e1 - pkin(2) * mrSges(4,1) + t25 * mrSges(5,1) + (-Ifges(5,3) - Ifges(4,2) / 0.2e1) * t44 + (Ifges(5,5) - Ifges(4,4)) * t102) * t21 - t44 * (Ifges(5,5) * t22 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t44 * (Ifges(4,4) * t22 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t25 * t22 * mrSges(5,3) - pkin(2) * t22 * mrSges(4,2) + m(5) * (t1 * t25 + t16 * t10) + t61 * t82 + (-t13 / 0.2e1 + t12 / 0.2e1) * t84 - t51 * t77 - t1 * t59 + (t43 * (Ifges(5,1) * t44 + t91) + t44 * (-Ifges(4,2) * t43 + t92) + t49) * t77 / 0.2e1 + (Ifges(4,1) * t43 + t57 + t92) * t22 / 0.2e1 + (t113 * t39 + (m(4) * pkin(2) - m(5) * t25 - t114 + t59) * t38) * g(1) + (-t23 * t83 - t24 * t84 - t105 + t106) * mrSges(4,3) + (t113 * t38 + (t118 * pkin(2) + t114 - t46) * t39) * g(2) + t110 * t83 / 0.2e1 + ((Ifges(5,1) + Ifges(4,1)) * t22 + t112 * qJDD(3)) * t102 + (-t111 * t44 + t112 * t43) * qJDD(3) / 0.2e1 + t48 * t63 - t10 * t19 + (m(4) * pkin(2) ^ 2 + Ifges(3,3)) * qJDD(2) + (t50 + t107 * qJD(3) / 0.2e1) * qJD(3) + (-t116 * t43 + t117 * t44 - t108 * t84 + t109 * t83 + m(5) * (t96 + t3 * t43 + (t11 * t44 - t17 * t43) * qJD(3)) + m(4) * ((-t23 * t44 - t24 * t43) * qJD(3) + t106) + t105 * t118) * pkin(5) + (t11 * t83 + t115 * t43 + t96) * mrSges(5,2); t13 * t86 / 0.2e1 - t11 * t69 - (Ifges(5,1) * t85 + t12 + t36) * t86 / 0.2e1 + (-t46 - t61) * g(3) + (Ifges(5,2) + Ifges(4,3)) * qJDD(3) + t111 * t21 + t112 * t22 + t107 * t63 + (-m(5) * t17 - t108 + t68) * t23 + (-m(5) * t11 - t109 + t70) * t24 - (-Ifges(4,2) * t86 + t110 + t37) * t85 / 0.2e1 + t17 * t71 + (-t3 * pkin(3) + t2 * qJ(4) + t17 * qJD(4) - t16 * t20) * m(5) + t2 * mrSges(5,3) - t3 * mrSges(5,1) - t4 * mrSges(4,2) + t5 * mrSges(4,1) - pkin(3) * t8 + qJ(4) * t9 + t20 * t19 + qJD(4) * t29 + (m(5) * t52 + t58 + t60) * t105 + (-t50 + (-t49 / 0.2e1 + t48 / 0.2e1 + t51) * qJD(2)) * qJD(2); -t19 * t86 - qJD(3) * t29 + (g(3) * t44 - t105 * t43 + t16 * t86 + t115) * m(5) + t8;];
tau = t6;
