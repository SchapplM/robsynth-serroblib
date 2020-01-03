% Calculate vector of cutting forces with Newton-Euler
% S5RRPPP1
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
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:24
% EndTime: 2019-12-31 19:23:27
% DurationCPUTime: 1.32s
% Computational Cost: add. (11033->182), mult. (27747->224), div. (0->0), fcn. (19247->8), ass. (0->88)
t80 = cos(qJ(2));
t107 = qJD(1) * t80;
t78 = sin(qJ(2));
t108 = qJD(1) * t78;
t76 = sin(pkin(8));
t109 = cos(pkin(8));
t110 = cos(pkin(5));
t90 = t110 * t109;
t77 = sin(pkin(5));
t97 = t77 * t109;
t52 = -qJD(2) * t97 - t107 * t90 + t108 * t76;
t106 = qJD(2) * t77;
t98 = t76 * t110;
t53 = t76 * t106 + (t109 * t78 + t80 * t98) * qJD(1);
t32 = pkin(3) * t52 - qJ(4) * t53;
t126 = (2 * qJD(3)) + t32;
t63 = -qJD(2) * t110 + t107 * t77;
t41 = mrSges(5,1) * t53 - mrSges(5,2) * t63;
t105 = qJD(1) * qJD(2);
t67 = qJDD(1) * t78 + t105 * t80;
t68 = qJDD(1) * t80 - t105 * t78;
t45 = -qJDD(2) * t97 + t67 * t76 - t68 * t90;
t88 = qJDD(2) * t77 + t110 * t68;
t46 = t109 * t67 + t76 * t88;
t117 = t52 * t63;
t121 = -2 * qJD(4);
t112 = qJ(3) * t67;
t96 = qJD(1) * t110;
t58 = (t80 * t96 + t106) * qJ(3);
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t100 = t79 * g(1) - t81 * g(2);
t82 = qJD(1) ^ 2;
t59 = -qJDD(1) * pkin(1) - t82 * pkin(7) - t100;
t111 = qJ(3) * t78;
t65 = qJD(2) * pkin(2) - t111 * t96;
t25 = -t77 * t112 - t68 * pkin(2) + (-t58 * t80 + t65 * t78) * qJD(1) + t59;
t118 = t80 * g(3);
t95 = -g(1) * t81 - g(2) * t79;
t60 = -pkin(1) * t82 + qJDD(1) * pkin(7) + t95;
t61 = (-pkin(2) * t80 - t111 * t77) * qJD(1);
t26 = -t110 * t112 + qJDD(2) * pkin(2) - t118 + qJD(2) * t58 + (-qJD(1) * t61 - t60) * t78;
t94 = t110 * t25 - t77 * t26 + qJDD(3);
t84 = (-t46 - t117) * qJ(4) + t94 + (-pkin(3) * t63 + t121) * t53;
t125 = m(5) * (t45 * pkin(3) + t84) - t53 * t41 - t46 * mrSges(5,3);
t39 = mrSges(5,1) * t52 + mrSges(5,3) * t63;
t40 = -mrSges(6,1) * t52 - mrSges(6,2) * t63;
t113 = -t39 + t40;
t34 = -mrSges(5,2) * t52 - mrSges(5,3) * t53;
t114 = -mrSges(4,1) * t52 - mrSges(4,2) * t53 - t34;
t115 = -mrSges(4,3) - mrSges(5,1);
t116 = mrSges(5,2) - mrSges(6,3);
t123 = -2 * qJD(3);
t35 = mrSges(4,2) * t63 - mrSges(4,3) * t52;
t55 = qJDD(2) * t110 - t68 * t77;
t99 = -t78 * g(3) + t80 * t60;
t27 = qJ(3) * t88 - qJD(2) * t65 + t107 * t61 + t99;
t83 = t25 * t97 + t26 * t90 - t27 * t76;
t120 = 2 * qJD(5);
t62 = t63 ^ 2;
t16 = -t55 * pkin(3) - t62 * qJ(4) + t126 * t53 + qJDD(4) - t83;
t31 = -mrSges(6,2) * t53 + mrSges(6,3) * t52;
t104 = m(6) * (t63 * t120 + (t52 * t53 - t55) * qJ(5) + (t46 - t117) * pkin(4) + t16) + t53 * t31 + t46 * mrSges(6,1);
t92 = m(5) * t16 + t104;
t7 = m(4) * t83 + (m(4) * t123 + t114) * t53 + t115 * t46 + (-t35 - t113) * t63 + (mrSges(4,1) - t116) * t55 - t92;
t70 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t108;
t71 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t107;
t102 = t25 * t76 * t77 + t109 * t27 + t26 * t98;
t36 = -mrSges(4,1) * t63 - mrSges(4,3) * t53;
t49 = t52 * t123;
t37 = pkin(4) * t53 + qJ(5) * t63;
t38 = mrSges(6,1) * t53 + mrSges(6,3) * t63;
t51 = t52 ^ 2;
t86 = t62 * pkin(3) - t55 * qJ(4) - t102;
t101 = t63 * t38 - m(6) * (-t45 * pkin(4) - t51 * qJ(5) - t52 * t32 + qJDD(5) + t49 + (t121 - t37) * t63 - t86) - t55 * mrSges(6,2);
t91 = m(5) * (0.2e1 * qJD(4) * t63 + t126 * t52 + t86) + t101;
t8 = m(4) * (t49 + t102) + (t36 - t41) * t63 + (-mrSges(4,2) + mrSges(5,3)) * t55 + (-t31 + t114) * t52 + (-mrSges(6,1) + t115) * t45 - t91;
t103 = m(6) * (-t51 * pkin(4) + t52 * t120 - t53 * t37 + (pkin(3) + qJ(5)) * t45 + t84) + t52 * t40 + t45 * mrSges(6,3);
t9 = m(4) * t94 + (t36 - t38) * t53 + (t35 - t39) * t52 + (mrSges(4,2) - mrSges(6,2)) * t46 + (mrSges(4,1) - mrSges(5,2)) * t45 + t103 + t125;
t124 = m(3) * t59 - t68 * mrSges(3,1) + t67 * mrSges(3,2) + (t70 * t78 - t71 * t80) * qJD(1) + t110 * t9 + t77 * (t109 * t7 + t76 * t8);
t66 = (-mrSges(3,1) * t80 + mrSges(3,2) * t78) * qJD(1);
t4 = m(3) * (-t60 * t78 - t118) - t67 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t66 * t108 + qJD(2) * t71 + t8 * t98 + t7 * t90 - t77 * t9;
t6 = m(3) * t99 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t68 - qJD(2) * t70 + t107 * t66 + t109 * t8 - t7 * t76;
t119 = t4 * t80 + t6 * t78;
t87 = -t46 * mrSges(6,2) - t53 * t38 + t103;
t2 = m(2) * t100 + qJDD(1) * mrSges(2,1) - t82 * mrSges(2,2) - t124;
t1 = m(2) * t95 - mrSges(2,1) * t82 - qJDD(1) * mrSges(2,2) - t4 * t78 + t6 * t80;
t3 = [-m(1) * g(1) + t1 * t81 - t2 * t79, t1, t6, t8, -t45 * mrSges(5,2) - t52 * t39 + t125 + t87, t87; -m(1) * g(2) + t1 * t79 + t2 * t81, t2, t4, t7, -t55 * mrSges(5,3) + t63 * t41 + (t31 + t34) * t52 + (mrSges(5,1) + mrSges(6,1)) * t45 + t91, -t55 * mrSges(6,3) + t63 * t40 + t104; (-m(1) - m(2)) * g(3) + t119, -m(2) * g(3) + t119, t124, t9, t46 * mrSges(5,1) + t113 * t63 + t116 * t55 + t53 * t34 + t92, -t45 * mrSges(6,1) - t52 * t31 - t101;];
f_new = t3;
