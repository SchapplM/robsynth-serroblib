% Calculate vector of cutting forces with Newton-Euler
% S6RPRRPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRPP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:15:18
% EndTime: 2019-05-05 21:15:22
% DurationCPUTime: 1.97s
% Computational Cost: add. (25999->173), mult. (51397->221), div. (0->0), fcn. (32954->10), ass. (0->87)
t87 = sin(qJ(1));
t90 = cos(qJ(1));
t106 = t87 * g(1) - t90 * g(2);
t68 = qJDD(1) * pkin(1) + t106;
t100 = -t90 * g(1) - t87 * g(2);
t92 = qJD(1) ^ 2;
t70 = -t92 * pkin(1) + t100;
t83 = sin(pkin(9));
t84 = cos(pkin(9));
t101 = t84 * t68 - t83 * t70;
t46 = -qJDD(1) * pkin(2) - t92 * pkin(7) - t101;
t113 = cos(pkin(10));
t86 = sin(qJ(3));
t112 = qJD(1) * t86;
t85 = sin(qJ(4));
t88 = cos(qJ(4));
t66 = t88 * qJD(3) - t85 * t112;
t67 = t85 * qJD(3) + t88 * t112;
t82 = sin(pkin(10));
t52 = -t113 * t66 + t82 * t67;
t53 = t113 * t67 + t82 * t66;
t36 = t52 * mrSges(7,1) - t53 * mrSges(7,3);
t116 = -t52 * mrSges(6,1) - t53 * mrSges(6,2) - t36;
t118 = -mrSges(6,3) - mrSges(7,2);
t35 = t52 * pkin(5) - t53 * qJ(6);
t110 = qJD(1) * qJD(3);
t105 = t86 * t110;
t89 = cos(qJ(3));
t73 = t89 * qJDD(1) - t105;
t65 = qJDD(4) - t73;
t111 = t89 * qJD(1);
t77 = qJD(4) - t111;
t76 = t77 ^ 2;
t104 = t89 * t110;
t72 = t86 * qJDD(1) + t104;
t27 = (-t72 - t104) * pkin(8) + (-t73 + t105) * pkin(3) + t46;
t114 = t83 * t68 + t84 * t70;
t47 = -t92 * pkin(2) + qJDD(1) * pkin(7) + t114;
t81 = -g(3) + qJDD(2);
t115 = t89 * t47 + t86 * t81;
t71 = (-pkin(3) * t89 - pkin(8) * t86) * qJD(1);
t91 = qJD(3) ^ 2;
t34 = -t91 * pkin(3) + qJDD(3) * pkin(8) + t71 * t111 + t115;
t103 = t88 * t27 - t85 * t34;
t51 = t66 * qJD(4) + t85 * qJDD(3) + t88 * t72;
t18 = (t66 * t77 - t51) * qJ(5) + (t66 * t67 + t65) * pkin(4) + t103;
t117 = t85 * t27 + t88 * t34;
t50 = -t67 * qJD(4) + t88 * qJDD(3) - t85 * t72;
t56 = t77 * pkin(4) - t67 * qJ(5);
t64 = t66 ^ 2;
t20 = -t64 * pkin(4) + t50 * qJ(5) - t77 * t56 + t117;
t98 = t113 * t18 - t82 * t20;
t119 = m(7) * (-t65 * pkin(5) - t76 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t35) * t53 - t98);
t120 = -2 * qJD(5);
t29 = t113 * t51 + t82 * t50;
t39 = -t77 * mrSges(6,2) - t52 * mrSges(6,3);
t42 = -t52 * mrSges(7,2) + t77 * mrSges(7,3);
t10 = m(6) * t98 - t119 + (t39 + t42) * t77 + (mrSges(6,1) + mrSges(7,1)) * t65 + (m(6) * t120 + t116) * t53 + t118 * t29;
t54 = -t66 * mrSges(5,1) + t67 * mrSges(5,2);
t55 = -t77 * mrSges(5,2) + t66 * mrSges(5,3);
t107 = t113 * t20 + t52 * t120 + t82 * t18;
t41 = -t77 * mrSges(7,1) + t53 * mrSges(7,2);
t108 = m(7) * (-t76 * pkin(5) + t65 * qJ(6) + 0.2e1 * qJD(6) * t77 - t52 * t35 + t107) + t77 * t41 + t65 * mrSges(7,3);
t28 = -t113 * t50 + t82 * t51;
t40 = t77 * mrSges(6,1) - t53 * mrSges(6,3);
t9 = m(6) * t107 - t65 * mrSges(6,2) + t116 * t52 + t118 * t28 - t77 * t40 + t108;
t7 = m(5) * t103 + t65 * mrSges(5,1) - t51 * mrSges(5,3) + t113 * t10 - t67 * t54 + t77 * t55 + t82 * t9;
t74 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t112;
t75 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t111;
t57 = t77 * mrSges(5,1) - t67 * mrSges(5,3);
t8 = m(5) * t117 - t65 * mrSges(5,2) + t50 * mrSges(5,3) - t82 * t10 + t113 * t9 + t66 * t54 - t77 * t57;
t121 = m(4) * t46 - t73 * mrSges(4,1) + t72 * mrSges(4,2) + t88 * t7 + t85 * t8 + (t86 * t74 - t89 * t75) * qJD(1);
t102 = -t86 * t47 + t89 * t81;
t69 = (-mrSges(4,1) * t89 + mrSges(4,2) * t86) * qJD(1);
t33 = -qJDD(3) * pkin(3) - t91 * pkin(8) + t71 * t112 - t102;
t94 = -t50 * pkin(4) - t64 * qJ(5) + t67 * t56 + qJDD(5) + t33;
t97 = t29 * mrSges(7,3) + t53 * t41 - m(7) * (-0.2e1 * qJD(6) * t53 + (t52 * t77 - t29) * qJ(6) + (t53 * t77 + t28) * pkin(5) + t94) - t28 * mrSges(7,1) - t52 * t42;
t95 = m(6) * t94 + t28 * mrSges(6,1) + t29 * mrSges(6,2) + t52 * t39 + t53 * t40 - t97;
t93 = m(5) * t33 - t50 * mrSges(5,1) + t51 * mrSges(5,2) - t66 * t55 + t67 * t57 + t95;
t12 = m(4) * t102 + qJDD(3) * mrSges(4,1) - t72 * mrSges(4,3) + qJD(3) * t75 - t69 * t112 - t93;
t6 = m(4) * t115 - qJDD(3) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(3) * t74 + t69 * t111 - t85 * t7 + t88 * t8;
t109 = m(3) * t81 + t89 * t12 + t86 * t6;
t4 = m(3) * t101 + qJDD(1) * mrSges(3,1) - t92 * mrSges(3,2) - t121;
t3 = m(3) * t114 - t92 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t86 * t12 + t89 * t6;
t2 = m(2) * t100 - t92 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t84 * t3 - t83 * t4;
t1 = m(2) * t106 + qJDD(1) * mrSges(2,1) - t92 * mrSges(2,2) + t83 * t3 + t84 * t4;
t5 = [-m(1) * g(1) - t87 * t1 + t90 * t2, t2, t3, t6, t8, t9, -t28 * mrSges(7,2) - t52 * t36 + t108; -m(1) * g(2) + t90 * t1 + t87 * t2, t1, t4, t12, t7, t10, -t97; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t109, t121, t93, t95, -t65 * mrSges(7,1) + t29 * mrSges(7,2) + t53 * t36 - t77 * t42 + t119;];
f_new  = t5;
