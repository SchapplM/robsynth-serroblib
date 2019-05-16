% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:17:38
% EndTime: 2019-05-06 01:17:43
% DurationCPUTime: 1.97s
% Computational Cost: add. (27090->173), mult. (51969->219), div. (0->0), fcn. (33610->10), ass. (0->88)
t89 = sin(qJ(1));
t92 = cos(qJ(1));
t107 = t89 * g(1) - t92 * g(2);
t68 = qJDD(1) * pkin(1) + t107;
t102 = -t92 * g(1) - t89 * g(2);
t94 = qJD(1) ^ 2;
t70 = -t94 * pkin(1) + t102;
t84 = sin(pkin(10));
t85 = cos(pkin(10));
t103 = t85 * t68 - t84 * t70;
t46 = -qJDD(1) * pkin(2) - t94 * pkin(7) - t103;
t110 = qJD(1) * qJD(3);
t91 = cos(qJ(3));
t106 = t91 * t110;
t88 = sin(qJ(3));
t72 = t88 * qJDD(1) + t106;
t80 = t88 * t110;
t73 = t91 * qJDD(1) - t80;
t29 = (-t72 - t106) * pkin(8) + (-t73 + t80) * pkin(3) + t46;
t111 = t91 * qJD(1);
t113 = t84 * t68 + t85 * t70;
t47 = -t94 * pkin(2) + qJDD(1) * pkin(7) + t113;
t83 = -g(3) + qJDD(2);
t114 = t91 * t47 + t88 * t83;
t71 = (-pkin(3) * t91 - pkin(8) * t88) * qJD(1);
t93 = qJD(3) ^ 2;
t34 = -t93 * pkin(3) + qJDD(3) * pkin(8) + t71 * t111 + t114;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t105 = t90 * t29 - t87 * t34;
t119 = cos(qJ(5));
t112 = qJD(1) * t88;
t66 = t90 * qJD(3) - t87 * t112;
t50 = t66 * qJD(4) + t87 * qJDD(3) + t90 * t72;
t65 = qJDD(4) - t73;
t67 = t87 * qJD(3) + t90 * t112;
t78 = qJD(4) - t111;
t18 = (t66 * t78 - t50) * pkin(9) + (t66 * t67 + t65) * pkin(4) + t105;
t116 = t87 * t29 + t90 * t34;
t49 = -t67 * qJD(4) + t90 * qJDD(3) - t87 * t72;
t56 = t78 * pkin(4) - t67 * pkin(9);
t64 = t66 ^ 2;
t20 = -t64 * pkin(4) + t49 * pkin(9) - t78 * t56 + t116;
t86 = sin(qJ(5));
t117 = t119 * t20 + t86 * t18;
t51 = -t119 * t66 + t86 * t67;
t52 = t119 * t67 + t86 * t66;
t35 = t51 * pkin(5) - t52 * qJ(6);
t77 = qJD(5) + t78;
t42 = -t77 * mrSges(7,1) + t52 * mrSges(7,2);
t63 = qJDD(5) + t65;
t76 = t77 ^ 2;
t108 = m(7) * (-t76 * pkin(5) + t63 * qJ(6) + 0.2e1 * qJD(6) * t77 - t51 * t35 + t117) + t77 * t42 + t63 * mrSges(7,3);
t36 = t51 * mrSges(7,1) - t52 * mrSges(7,3);
t115 = -t51 * mrSges(6,1) - t52 * mrSges(6,2) - t36;
t118 = -mrSges(6,3) - mrSges(7,2);
t25 = t52 * qJD(5) - t119 * t49 + t86 * t50;
t41 = t77 * mrSges(6,1) - t52 * mrSges(6,3);
t11 = m(6) * t117 - t63 * mrSges(6,2) + t115 * t51 + t118 * t25 - t77 * t41 + t108;
t100 = t119 * t18 - t86 * t20;
t120 = m(7) * (-t63 * pkin(5) - t76 * qJ(6) + t52 * t35 + qJDD(6) - t100);
t26 = -t51 * qJD(5) + t119 * t50 + t86 * t49;
t39 = -t51 * mrSges(7,2) + t77 * mrSges(7,3);
t40 = -t77 * mrSges(6,2) - t51 * mrSges(6,3);
t12 = m(6) * t100 - t120 + (t40 + t39) * t77 + (mrSges(6,1) + mrSges(7,1)) * t63 + t115 * t52 + t118 * t26;
t53 = -t66 * mrSges(5,1) + t67 * mrSges(5,2);
t54 = -t78 * mrSges(5,2) + t66 * mrSges(5,3);
t7 = m(5) * t105 + t65 * mrSges(5,1) - t50 * mrSges(5,3) + t86 * t11 + t119 * t12 - t67 * t53 + t78 * t54;
t74 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t112;
t75 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t111;
t55 = t78 * mrSges(5,1) - t67 * mrSges(5,3);
t8 = m(5) * t116 - t65 * mrSges(5,2) + t49 * mrSges(5,3) + t119 * t11 - t86 * t12 + t66 * t53 - t78 * t55;
t121 = m(4) * t46 - t73 * mrSges(4,1) + t72 * mrSges(4,2) + t90 * t7 + t87 * t8 + (t88 * t74 - t91 * t75) * qJD(1);
t104 = -t88 * t47 + t91 * t83;
t69 = (-mrSges(4,1) * t91 + mrSges(4,2) * t88) * qJD(1);
t33 = -qJDD(3) * pkin(3) - t93 * pkin(8) + t71 * t112 - t104;
t97 = -t49 * pkin(4) - t64 * pkin(9) + t67 * t56 + t33;
t99 = t26 * mrSges(7,3) + t52 * t42 - m(7) * (-0.2e1 * qJD(6) * t52 + (t51 * t77 - t26) * qJ(6) + (t52 * t77 + t25) * pkin(5) + t97) - t25 * mrSges(7,1) - t51 * t39;
t96 = m(6) * t97 + t25 * mrSges(6,1) + t26 * mrSges(6,2) + t51 * t40 + t52 * t41 - t99;
t95 = m(5) * t33 - t49 * mrSges(5,1) + t50 * mrSges(5,2) - t66 * t54 + t67 * t55 + t96;
t10 = m(4) * t104 + qJDD(3) * mrSges(4,1) - t72 * mrSges(4,3) + qJD(3) * t75 - t69 * t112 - t95;
t6 = m(4) * t114 - qJDD(3) * mrSges(4,2) + t73 * mrSges(4,3) - qJD(3) * t74 + t69 * t111 - t87 * t7 + t90 * t8;
t109 = m(3) * t83 + t91 * t10 + t88 * t6;
t4 = m(3) * t103 + qJDD(1) * mrSges(3,1) - t94 * mrSges(3,2) - t121;
t3 = m(3) * t113 - t94 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t88 * t10 + t91 * t6;
t2 = m(2) * t102 - t94 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t85 * t3 - t84 * t4;
t1 = m(2) * t107 + qJDD(1) * mrSges(2,1) - t94 * mrSges(2,2) + t84 * t3 + t85 * t4;
t5 = [-m(1) * g(1) - t89 * t1 + t92 * t2, t2, t3, t6, t8, t11, -t25 * mrSges(7,2) - t51 * t36 + t108; -m(1) * g(2) + t92 * t1 + t89 * t2, t1, t4, t10, t7, t12, -t99; (-m(1) - m(2)) * g(3) + t109, -m(2) * g(3) + t109, t109, t121, t95, t96, -t63 * mrSges(7,1) + t26 * mrSges(7,2) + t52 * t36 - t77 * t39 + t120;];
f_new  = t5;
