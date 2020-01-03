% Calculate vector of cutting forces with Newton-Euler
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:58
% EndTime: 2019-12-31 21:04:00
% DurationCPUTime: 0.67s
% Computational Cost: add. (4951->161), mult. (9671->188), div. (0->0), fcn. (5649->6), ass. (0->72)
t77 = sin(qJ(2));
t102 = qJD(1) * t77;
t110 = cos(qJ(3));
t76 = sin(qJ(3));
t60 = t76 * qJD(2) + t110 * t102;
t79 = cos(qJ(2));
t101 = t79 * qJD(1);
t69 = -qJD(3) + t101;
t47 = t69 * mrSges(5,1) + t60 * mrSges(5,2);
t100 = qJD(1) * qJD(2);
t71 = t77 * t100;
t64 = t79 * qJDD(1) - t71;
t58 = -qJDD(3) + t64;
t82 = qJD(1) ^ 2;
t78 = sin(qJ(1));
t80 = cos(qJ(1));
t97 = t78 * g(1) - t80 * g(2);
t52 = -qJDD(1) * pkin(1) - t82 * pkin(6) - t97;
t95 = t79 * t100;
t63 = t77 * qJDD(1) + t95;
t19 = (-t63 - t95) * pkin(7) + (-t64 + t71) * pkin(2) + t52;
t62 = (-pkin(2) * t79 - pkin(7) * t77) * qJD(1);
t81 = qJD(2) ^ 2;
t92 = -t80 * g(1) - t78 * g(2);
t53 = -t82 * pkin(1) + qJDD(1) * pkin(6) + t92;
t96 = -t77 * g(3) + t79 * t53;
t23 = -t81 * pkin(2) + qJDD(2) * pkin(7) + t62 * t101 + t96;
t106 = t110 * t23 + t76 * t19;
t112 = -2 * qJD(4);
t59 = -t110 * qJD(2) + t76 * t102;
t38 = t59 * pkin(3) - t60 * qJ(4);
t68 = t69 ^ 2;
t88 = -t68 * pkin(3) - t58 * qJ(4) + t69 * t112 - t59 * t38 + t106;
t116 = m(5) * t88 - t58 * mrSges(5,3) - t69 * t47;
t109 = t59 * t69;
t35 = -t59 * qJD(3) + t76 * qJDD(2) + t110 * t63;
t103 = -t79 * g(3) - t77 * t53;
t85 = qJDD(2) * pkin(2) + t81 * pkin(7) - t62 * t102 + t103;
t115 = -(t35 + t109) * qJ(4) - t85;
t65 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t102;
t66 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t101;
t39 = t59 * mrSges(5,1) - t60 * mrSges(5,3);
t105 = -t59 * mrSges(4,1) - t60 * mrSges(4,2) - t39;
t107 = -mrSges(4,3) - mrSges(5,2);
t40 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t43 = t69 * mrSges(4,2) - t59 * mrSges(4,3);
t48 = -t59 * mrSges(5,2) - t69 * mrSges(5,3);
t91 = t110 * t19 - t76 * t23;
t17 = t58 * pkin(3) - t68 * qJ(4) + t60 * t38 + qJDD(4) - t91;
t42 = -t69 * mrSges(6,2) + t59 * mrSges(6,3);
t99 = t69 * t42 + t58 * mrSges(6,1) + m(6) * (-0.2e1 * qJD(5) * t60 + (-t35 + t109) * qJ(5) + (t59 * t60 + t58) * pkin(4) + t17);
t89 = m(5) * t17 + t99;
t7 = m(4) * t91 + (-t43 - t48) * t69 + (-mrSges(4,1) - mrSges(5,1)) * t58 + (t40 + t105) * t60 + (mrSges(6,3) + t107) * t35 - t89;
t34 = t60 * qJD(3) - t110 * qJDD(2) + t76 * t63;
t45 = t69 * mrSges(6,1) - t60 * mrSges(6,3);
t46 = -t69 * mrSges(4,1) - t60 * mrSges(4,3);
t44 = t69 * pkin(4) - t60 * qJ(5);
t57 = t59 ^ 2;
t98 = m(6) * (-t57 * pkin(4) + t34 * qJ(5) + 0.2e1 * qJD(5) * t59 - t69 * t44 + t88) + t59 * t40 + t34 * mrSges(6,3);
t8 = m(4) * t106 + (t46 - t45) * t69 + t105 * t59 + (mrSges(4,2) - mrSges(6,2)) * t58 + t107 * t34 + t98 + t116;
t114 = m(3) * t52 - t64 * mrSges(3,1) + t63 * mrSges(3,2) + (t65 * t77 - t66 * t79) * qJD(1) + t110 * t7 + t76 * t8;
t93 = m(6) * (-t57 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t34 + (pkin(3) * t69 + (2 * qJD(4)) + t44) * t60 - t115) + t35 * mrSges(6,2) - t34 * mrSges(6,1) + t60 * t45 - t59 * t42;
t87 = m(5) * (t60 * t112 + (-t60 * t69 + t34) * pkin(3) + t115) + t34 * mrSges(5,1) + t59 * t48 - t93;
t113 = -m(4) * t85 + t34 * mrSges(4,1) + (t46 - t47) * t60 + (mrSges(4,2) - mrSges(5,3)) * t35 + t59 * t43 + t87;
t61 = (-mrSges(3,1) * t79 + mrSges(3,2) * t77) * qJD(1);
t4 = m(3) * t96 - qJDD(2) * mrSges(3,2) + t64 * mrSges(3,3) - qJD(2) * t65 + t61 * t101 + t110 * t8 - t76 * t7;
t6 = m(3) * t103 + qJDD(2) * mrSges(3,1) - t63 * mrSges(3,3) + qJD(2) * t66 - t61 * t102 - t113;
t111 = t77 * t4 + t79 * t6;
t86 = -t58 * mrSges(6,2) - t69 * t45 + t98;
t2 = m(2) * t97 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t114;
t1 = m(2) * t92 - t82 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t79 * t4 - t77 * t6;
t3 = [-m(1) * g(1) + t80 * t1 - t78 * t2, t1, t4, t8, -t34 * mrSges(5,2) - t59 * t39 + t116 + t86, t86; -m(1) * g(2) + t78 * t1 + t80 * t2, t2, t6, t7, -t35 * mrSges(5,3) - t60 * t47 + t87, -t35 * mrSges(6,3) - t60 * t40 + t99; (-m(1) - m(2)) * g(3) + t111, -m(2) * g(3) + t111, t114, t113, t58 * mrSges(5,1) + t69 * t48 + (t39 - t40) * t60 + (mrSges(5,2) - mrSges(6,3)) * t35 + t89, t93;];
f_new = t3;
