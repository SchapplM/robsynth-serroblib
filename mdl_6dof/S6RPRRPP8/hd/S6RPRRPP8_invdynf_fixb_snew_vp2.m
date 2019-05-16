% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-05-05 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:50:47
% EndTime: 2019-05-05 21:50:49
% DurationCPUTime: 0.78s
% Computational Cost: add. (6690->179), mult. (12604->197), div. (0->0), fcn. (7012->6), ass. (0->82)
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t98 = -t81 * g(1) - t79 * g(2);
t129 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t98;
t123 = cos(qJ(4));
t80 = cos(qJ(3));
t110 = qJD(1) * t80;
t77 = sin(qJ(4));
t63 = t77 * qJD(3) + t123 * t110;
t108 = qJD(1) * qJD(3);
t78 = sin(qJ(3));
t100 = t78 * t108;
t67 = t80 * qJDD(1) - t100;
t33 = t63 * qJD(4) - t123 * qJDD(3) + t77 * t67;
t62 = -t123 * qJD(3) + t77 * t110;
t34 = -t62 * qJD(4) + t77 * qJDD(3) + t123 * t67;
t109 = t78 * qJD(1);
t71 = qJD(4) + t109;
t46 = t63 * mrSges(6,1) + t71 * mrSges(6,2);
t122 = t62 * t71;
t126 = -2 * qJD(5);
t65 = (pkin(3) * t78 - pkin(8) * t80) * qJD(1);
t82 = qJD(3) ^ 2;
t124 = -pkin(1) - pkin(7);
t102 = t79 * g(1) - t81 * g(2);
t83 = qJD(1) ^ 2;
t91 = -t83 * qJ(2) + qJDD(2) - t102;
t51 = t124 * qJDD(1) + t91;
t96 = t78 * g(3) + t80 * t51;
t24 = -qJDD(3) * pkin(3) - t82 * pkin(8) + t65 * t110 - t96;
t84 = (-t34 + t122) * qJ(5) + t24 + (t71 * pkin(4) + t126) * t63;
t128 = m(6) * (t33 * pkin(4) + t84) - t34 * mrSges(6,3) - t63 * t46;
t42 = t63 * pkin(5) - t71 * qJ(6);
t45 = -t62 * mrSges(7,1) + t71 * mrSges(7,2);
t60 = t62 ^ 2;
t105 = m(7) * (-t60 * pkin(5) + 0.2e1 * qJD(6) * t62 - t63 * t42 + (pkin(4) + qJ(6)) * t33 + t84) + t33 * mrSges(7,3) + t62 * t45;
t40 = -t71 * mrSges(5,2) - t62 * mrSges(5,3);
t41 = t71 * mrSges(5,1) - t63 * mrSges(5,3);
t43 = t63 * mrSges(7,1) - t71 * mrSges(7,3);
t44 = t62 * mrSges(6,1) - t71 * mrSges(6,3);
t127 = m(5) * t24 + (t41 - t43) * t63 + (t40 - t44) * t62 + (mrSges(5,2) - mrSges(7,2)) * t34 + (mrSges(5,1) - mrSges(6,2)) * t33 + t105 + t128;
t125 = -m(2) - m(3);
t121 = (mrSges(2,1) - mrSges(3,2));
t119 = -mrSges(2,2) + mrSges(3,3);
t117 = mrSges(6,2) - mrSges(7,3);
t116 = -mrSges(5,3) - mrSges(6,1);
t99 = t80 * t108;
t66 = -t78 * qJDD(1) - t99;
t85 = t124 * t83 - t129;
t21 = (-t67 + t100) * pkin(8) + (-t66 + t99) * pkin(3) + t85;
t101 = -t80 * g(3) + t78 * t51;
t25 = -t82 * pkin(3) + qJDD(3) * pkin(8) - t65 * t109 + t101;
t115 = t123 * t25 + t77 * t21;
t39 = -t62 * mrSges(6,2) - t63 * mrSges(6,3);
t114 = -t62 * mrSges(5,1) - t63 * mrSges(5,2) - t39;
t111 = t44 - t45;
t37 = t62 * pkin(4) - t63 * qJ(5);
t61 = qJDD(4) - t66;
t70 = t71 ^ 2;
t97 = t123 * t21 - t77 * t25;
t17 = -t61 * pkin(4) - t70 * qJ(5) + t63 * t37 + qJDD(5) - t97;
t36 = -t63 * mrSges(7,2) + t62 * mrSges(7,3);
t106 = m(7) * (-0.2e1 * qJD(6) * t71 + (t62 * t63 - t61) * qJ(6) + (t34 + t122) * pkin(5) + t17) + t63 * t36 + t34 * mrSges(7,1);
t87 = -t70 * pkin(4) + t61 * qJ(5) - t62 * t37 + t115;
t104 = m(7) * (-t33 * pkin(5) - t60 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t42) * t71 + t87) + t71 * t43 + t61 * mrSges(7,2);
t64 = (mrSges(4,1) * t78 + mrSges(4,2) * t80) * qJD(1);
t69 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t110;
t95 = m(6) * t17 + t106;
t7 = m(5) * t97 + t114 * t63 + t116 * t34 + (t40 - t111) * t71 + (mrSges(5,1) - t117) * t61 - t95;
t94 = m(6) * (t71 * t126 - t87) - t104;
t9 = m(5) * t115 + (-t41 + t46) * t71 + (-mrSges(5,2) + mrSges(6,3)) * t61 + (-t36 + t114) * t62 + (-mrSges(7,1) + t116) * t33 - t94;
t4 = m(4) * t101 - qJDD(3) * mrSges(4,2) + t66 * mrSges(4,3) - qJD(3) * t69 - t64 * t109 + t123 * t9 - t77 * t7;
t68 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t109;
t5 = m(4) * t96 + qJDD(3) * mrSges(4,1) - t67 * mrSges(4,3) + qJD(3) * t68 - t64 * t110 - t127;
t103 = t80 * t4 - t78 * t5;
t92 = -m(3) * (-qJDD(1) * pkin(1) + t91) - t78 * t4 - t80 * t5;
t90 = -t34 * mrSges(7,2) - t63 * t43 + t105;
t89 = m(4) * t85 - t66 * mrSges(4,1) + t67 * mrSges(4,2) + t68 * t109 + t69 * t110 + t123 * t7 + t77 * t9;
t86 = -m(3) * (t83 * pkin(1) + t129) + t89;
t2 = m(2) * t98 + t119 * qJDD(1) - (t121 * t83) + t86;
t1 = m(2) * t102 + t121 * qJDD(1) + t119 * t83 + t92;
t3 = [-m(1) * g(1) - t79 * t1 + t81 * t2, t2, -m(3) * g(3) + t103, t4, t9, -t33 * mrSges(6,2) - t62 * t44 + t128 + t90, t90; -m(1) * g(2) + t81 * t1 + t79 * t2, t1, -(t83 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t86, t5, t7, -t61 * mrSges(6,3) - t71 * t46 + (t36 + t39) * t62 + (mrSges(6,1) + mrSges(7,1)) * t33 + t94, -t61 * mrSges(7,3) - t71 * t45 + t106; (-m(1) + t125) * g(3) + t103, t125 * g(3) + t103, qJDD(1) * mrSges(3,2) - t83 * mrSges(3,3) - t92, t89, t127, t34 * mrSges(6,1) + t111 * t71 + t117 * t61 + t63 * t39 + t95, -t33 * mrSges(7,1) - t62 * t36 + t104;];
f_new  = t3;
