% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:25:24
% EndTime: 2019-05-05 21:25:27
% DurationCPUTime: 0.97s
% Computational Cost: add. (10059->178), mult. (18910->203), div. (0->0), fcn. (10959->8), ass. (0->84)
t123 = cos(qJ(4));
t83 = sin(qJ(3));
t110 = qJD(1) * t83;
t82 = sin(qJ(4));
t62 = t82 * qJD(3) + t123 * t110;
t108 = qJD(1) * qJD(3);
t85 = cos(qJ(3));
t101 = t85 * t108;
t67 = t83 * qJDD(1) + t101;
t38 = t62 * qJD(4) - t123 * qJDD(3) + t82 * t67;
t61 = -t123 * qJD(3) + t82 * t110;
t39 = -t61 * qJD(4) + t82 * qJDD(3) + t123 * t67;
t109 = t85 * qJD(1);
t72 = -qJD(4) + t109;
t51 = t62 * mrSges(6,1) - t72 * mrSges(6,2);
t122 = t61 * t72;
t125 = -2 * qJD(5);
t84 = sin(qJ(1));
t86 = cos(qJ(1));
t103 = t84 * g(1) - t86 * g(2);
t63 = qJDD(1) * pkin(1) + t103;
t88 = qJD(1) ^ 2;
t98 = -t86 * g(1) - t84 * g(2);
t65 = -t88 * pkin(1) + t98;
t80 = sin(pkin(9));
t81 = cos(pkin(9));
t111 = t80 * t63 + t81 * t65;
t32 = -t88 * pkin(2) + qJDD(1) * pkin(7) + t111;
t79 = -g(3) + qJDD(2);
t100 = -t83 * t32 + t85 * t79;
t66 = (-pkin(3) * t85 - pkin(8) * t83) * qJD(1);
t87 = qJD(3) ^ 2;
t25 = -qJDD(3) * pkin(3) - t87 * pkin(8) + t66 * t110 - t100;
t89 = (-t39 - t122) * qJ(5) + t25 + (-t72 * pkin(4) + t125) * t62;
t128 = m(6) * (t38 * pkin(4) + t89) - t39 * mrSges(6,3) - t62 * t51;
t44 = -t61 * mrSges(6,2) - t62 * mrSges(6,3);
t115 = -t61 * mrSges(5,1) - t62 * mrSges(5,2) - t44;
t102 = t83 * t108;
t99 = t81 * t63 - t80 * t65;
t31 = -qJDD(1) * pkin(2) - t88 * pkin(7) - t99;
t68 = t85 * qJDD(1) - t102;
t22 = (-t67 - t101) * pkin(8) + (-t68 + t102) * pkin(3) + t31;
t116 = t85 * t32 + t83 * t79;
t26 = -t87 * pkin(3) + qJDD(3) * pkin(8) + t66 * t109 + t116;
t117 = t123 * t26 + t82 * t22;
t118 = -mrSges(5,3) - mrSges(6,1);
t41 = -t62 * mrSges(7,2) + t61 * mrSges(7,3);
t46 = -t72 * mrSges(5,1) - t62 * mrSges(5,3);
t60 = qJDD(4) - t68;
t47 = t62 * pkin(5) + t72 * qJ(6);
t48 = t62 * mrSges(7,1) + t72 * mrSges(7,3);
t59 = t61 ^ 2;
t42 = t61 * pkin(4) - t62 * qJ(5);
t71 = t72 ^ 2;
t91 = -t71 * pkin(4) + t60 * qJ(5) - t61 * t42 + t117;
t104 = t72 * t48 - m(7) * (-t38 * pkin(5) - t59 * qJ(6) + qJDD(6) + (t125 - t47) * t72 + t91) - t60 * mrSges(7,2);
t94 = m(6) * (0.2e1 * qJD(5) * t72 - t91) + t104;
t10 = m(5) * t117 + (t46 - t51) * t72 + (-mrSges(5,2) + mrSges(6,3)) * t60 + (-t41 + t115) * t61 + (-mrSges(7,1) + t118) * t38 - t94;
t69 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t110;
t70 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t109;
t49 = t61 * mrSges(6,1) + t72 * mrSges(6,3);
t50 = -t61 * mrSges(7,1) - t72 * mrSges(7,2);
t112 = -t49 + t50;
t119 = mrSges(6,2) - mrSges(7,3);
t45 = t72 * mrSges(5,2) - t61 * mrSges(5,3);
t124 = 2 * qJD(6);
t97 = t123 * t22 - t82 * t26;
t18 = -t60 * pkin(4) - t71 * qJ(5) + t62 * t42 + qJDD(5) - t97;
t106 = m(7) * (t72 * t124 + (t61 * t62 - t60) * qJ(6) + (t39 - t122) * pkin(5) + t18) + t62 * t41 + t39 * mrSges(7,1);
t95 = m(6) * t18 + t106;
t9 = m(5) * t97 + t115 * t62 + t118 * t39 + (-t45 - t112) * t72 + (mrSges(5,1) - t119) * t60 - t95;
t127 = m(4) * t31 - t68 * mrSges(4,1) + t67 * mrSges(4,2) + (t69 * t83 - t70 * t85) * qJD(1) + t82 * t10 + t123 * t9;
t105 = m(7) * (-t59 * pkin(5) + t61 * t124 - t62 * t47 + (pkin(4) + qJ(6)) * t38 + t89) + t38 * mrSges(7,3) + t61 * t50;
t126 = m(5) * t25 + (t46 - t48) * t62 + (t45 - t49) * t61 + (mrSges(5,2) - mrSges(7,2)) * t39 + (mrSges(5,1) - mrSges(6,2)) * t38 + t105 + t128;
t64 = (-mrSges(4,1) * t85 + mrSges(4,2) * t83) * qJD(1);
t6 = m(4) * t116 - qJDD(3) * mrSges(4,2) + t68 * mrSges(4,3) - qJD(3) * t69 + t123 * t10 + t64 * t109 - t82 * t9;
t8 = m(4) * t100 + qJDD(3) * mrSges(4,1) - t67 * mrSges(4,3) + qJD(3) * t70 - t64 * t110 - t126;
t107 = m(3) * t79 + t83 * t6 + t85 * t8;
t93 = -t39 * mrSges(7,2) - t62 * t48 + t105;
t4 = m(3) * t99 + qJDD(1) * mrSges(3,1) - t88 * mrSges(3,2) - t127;
t3 = m(3) * t111 - t88 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t85 * t6 - t83 * t8;
t2 = m(2) * t98 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t81 * t3 - t80 * t4;
t1 = m(2) * t103 + qJDD(1) * mrSges(2,1) - t88 * mrSges(2,2) + t80 * t3 + t81 * t4;
t5 = [-m(1) * g(1) - t84 * t1 + t86 * t2, t2, t3, t6, t10, -t38 * mrSges(6,2) - t61 * t49 + t128 + t93, t93; -m(1) * g(2) + t86 * t1 + t84 * t2, t1, t4, t8, t9, -t60 * mrSges(6,3) + t72 * t51 + (t41 + t44) * t61 + (mrSges(6,1) + mrSges(7,1)) * t38 + t94, -t60 * mrSges(7,3) + t72 * t50 + t106; (-m(1) - m(2)) * g(3) + t107, -m(2) * g(3) + t107, t107, t127, t126, t39 * mrSges(6,1) + t112 * t72 + t119 * t60 + t62 * t44 + t95, -t38 * mrSges(7,1) - t61 * t41 - t104;];
f_new  = t5;
