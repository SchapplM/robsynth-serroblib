% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:42:34
% EndTime: 2019-05-05 09:42:41
% DurationCPUTime: 2.84s
% Computational Cost: add. (38725->173), mult. (75179->224), div. (0->0), fcn. (53395->12), ass. (0->92)
t100 = qJD(3) ^ 2;
t101 = qJD(2) ^ 2;
t89 = sin(pkin(6));
t95 = sin(qJ(2));
t127 = t89 * t95;
t88 = sin(pkin(11));
t90 = cos(pkin(11));
t78 = t88 * g(1) - t90 * g(2);
t91 = cos(pkin(6));
t128 = t78 * t91;
t79 = -t90 * g(1) - t88 * g(2);
t87 = -g(3) + qJDD(1);
t99 = cos(qJ(2));
t115 = t87 * t127 + t95 * t128 + t99 * t79;
t46 = -t101 * pkin(2) + qJDD(2) * pkin(8) + t115;
t63 = -t89 * t78 + t91 * t87;
t94 = sin(qJ(3));
t98 = cos(qJ(3));
t108 = -t94 * t46 + t98 * t63;
t121 = qJD(2) * t94;
t75 = (-pkin(3) * t98 - pkin(9) * t94) * qJD(2);
t33 = -qJDD(3) * pkin(3) - t100 * pkin(9) + t75 * t121 - t108;
t93 = sin(qJ(4));
t97 = cos(qJ(4));
t73 = t93 * qJD(3) + t97 * t121;
t119 = qJD(2) * qJD(3);
t111 = t98 * t119;
t76 = t94 * qJDD(2) + t111;
t53 = -t73 * qJD(4) + t97 * qJDD(3) - t93 * t76;
t120 = t98 * qJD(2);
t84 = qJD(4) - t120;
t62 = t84 * pkin(4) - t73 * pkin(10);
t72 = t97 * qJD(3) - t93 * t121;
t68 = t72 ^ 2;
t103 = -t53 * pkin(4) - t68 * pkin(10) + t73 * t62 + t33;
t54 = t72 * qJD(4) + t93 * qJDD(3) + t97 * t76;
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t57 = t92 * t72 + t96 * t73;
t29 = -t57 * qJD(5) + t96 * t53 - t92 * t54;
t56 = t96 * t72 - t92 * t73;
t30 = t56 * qJD(5) + t92 * t53 + t96 * t54;
t83 = qJD(5) + t84;
t49 = t83 * pkin(5) - t57 * qJ(6);
t50 = t83 * mrSges(7,1) - t57 * mrSges(7,3);
t55 = t56 ^ 2;
t116 = m(7) * (-t29 * pkin(5) - t55 * qJ(6) + t57 * t49 + qJDD(6) + t103) + t30 * mrSges(7,2) + t57 * t50;
t47 = -t83 * mrSges(7,2) + t56 * mrSges(7,3);
t48 = -t83 * mrSges(6,2) + t56 * mrSges(6,3);
t51 = t83 * mrSges(6,1) - t57 * mrSges(6,3);
t134 = m(6) * t103 + t30 * mrSges(6,2) + t57 * t51 + t116 - (t48 + t47) * t56 - (mrSges(6,1) + mrSges(7,1)) * t29;
t60 = -t84 * mrSges(5,2) + t72 * mrSges(5,3);
t61 = t84 * mrSges(5,1) - t73 * mrSges(5,3);
t133 = m(5) * t33 - t53 * mrSges(5,1) + t54 * mrSges(5,2) - t72 * t60 + t73 * t61 + t134;
t131 = (t87 * t89 + t128) * t99 - t95 * t79;
t123 = t98 * t46 + t94 * t63;
t34 = -t100 * pkin(3) + qJDD(3) * pkin(9) + t75 * t120 + t123;
t45 = -qJDD(2) * pkin(2) - t101 * pkin(8) - t131;
t85 = t94 * t119;
t77 = t98 * qJDD(2) - t85;
t37 = (-t76 - t111) * pkin(9) + (-t77 + t85) * pkin(3) + t45;
t109 = -t93 * t34 + t97 * t37;
t69 = qJDD(4) - t77;
t21 = (t72 * t84 - t54) * pkin(10) + (t72 * t73 + t69) * pkin(4) + t109;
t124 = t97 * t34 + t93 * t37;
t23 = -t68 * pkin(4) + t53 * pkin(10) - t84 * t62 + t124;
t110 = t96 * t21 - t92 * t23;
t67 = qJDD(5) + t69;
t118 = m(7) * (-0.2e1 * qJD(6) * t57 + (t56 * t83 - t30) * qJ(6) + (t56 * t57 + t67) * pkin(5) + t110) + t83 * t47 + t67 * mrSges(7,1);
t40 = -t56 * mrSges(7,1) + t57 * mrSges(7,2);
t41 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t12 = m(6) * t110 + t67 * mrSges(6,1) + t83 * t48 + (-t41 - t40) * t57 + (-mrSges(6,3) - mrSges(7,3)) * t30 + t118;
t125 = t92 * t21 + t96 * t23;
t117 = m(7) * (-t55 * pkin(5) + t29 * qJ(6) + 0.2e1 * qJD(6) * t56 - t83 * t49 + t125) + t29 * mrSges(7,3) + t56 * t40;
t13 = m(6) * t125 + t29 * mrSges(6,3) + t56 * t41 + (-t51 - t50) * t83 + (-mrSges(6,2) - mrSges(7,2)) * t67 + t117;
t58 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t10 = m(5) * t109 + t69 * mrSges(5,1) - t54 * mrSges(5,3) + t96 * t12 + t92 * t13 - t73 * t58 + t84 * t60;
t11 = m(5) * t124 - t69 * mrSges(5,2) + t53 * mrSges(5,3) - t92 * t12 + t96 * t13 + t72 * t58 - t84 * t61;
t80 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t121;
t81 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t120;
t130 = m(4) * t45 - t77 * mrSges(4,1) + t76 * mrSges(4,2) + t97 * t10 + t93 * t11 + (t94 * t80 - t98 * t81) * qJD(2);
t8 = m(3) * t131 + qJDD(2) * mrSges(3,1) - t101 * mrSges(3,2) - t130;
t129 = t8 * t99;
t74 = (-mrSges(4,1) * t98 + mrSges(4,2) * t94) * qJD(2);
t14 = m(4) * t108 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t81 - t74 * t121 - t133;
t9 = m(4) * t123 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t80 - t93 * t10 + t97 * t11 + t74 * t120;
t4 = m(3) * t115 - t101 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t94 * t14 + t98 * t9;
t6 = m(3) * t63 + t98 * t14 + t94 * t9;
t114 = m(2) * t87 + t4 * t127 + t89 * t129 + t91 * t6;
t2 = m(2) * t79 + t99 * t4 - t95 * t8;
t1 = m(2) * t78 - t89 * t6 + (t4 * t95 + t129) * t91;
t3 = [-m(1) * g(1) - t88 * t1 + t90 * t2, t2, t4, t9, t11, t13, -t67 * mrSges(7,2) - t83 * t50 + t117; -m(1) * g(2) + t90 * t1 + t88 * t2, t1, t8, t14, t10, t12, -t30 * mrSges(7,3) - t57 * t40 + t118; -m(1) * g(3) + t114, t114, t6, t130, t133, t134, -t29 * mrSges(7,1) - t56 * t47 + t116;];
f_new  = t3;
