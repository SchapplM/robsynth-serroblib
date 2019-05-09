% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRP6
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:11:06
% EndTime: 2019-05-05 10:11:15
% DurationCPUTime: 4.81s
% Computational Cost: add. (68909->180), mult. (141136->244), div. (0->0), fcn. (111754->14), ass. (0->100)
t103 = qJD(2) ^ 2;
t102 = cos(qJ(2));
t92 = sin(pkin(6));
t121 = t102 * t92;
t90 = sin(pkin(12));
t93 = cos(pkin(12));
t81 = t90 * g(1) - t93 * g(2);
t95 = cos(pkin(6));
t129 = t81 * t95;
t82 = -t93 * g(1) - t90 * g(2);
t89 = -g(3) + qJDD(1);
t99 = sin(qJ(2));
t110 = t102 * t129 + t89 * t121 - t99 * t82;
t91 = sin(pkin(7));
t133 = pkin(9) * t91;
t51 = qJDD(2) * pkin(2) + t103 * t133 + t110;
t67 = -t92 * t81 + t95 * t89;
t94 = cos(pkin(7));
t137 = t51 * t94 + t67 * t91;
t101 = cos(qJ(3));
t119 = qJD(2) * qJD(3);
t98 = sin(qJ(3));
t74 = (-qJDD(2) * t101 + t98 * t119) * t91;
t128 = t92 * t99;
t115 = t102 * t82 + t89 * t128 + t99 * t129;
t52 = -t103 * pkin(2) + qJDD(2) * t133 + t115;
t136 = t137 * t101 - t98 * t52;
t100 = cos(qJ(4));
t120 = qJD(2) * t91;
t113 = t101 * t120;
t116 = t101 * t52 + t137 * t98;
t72 = (-pkin(3) * t101 - pkin(10) * t98) * t120;
t87 = t94 * qJD(2) + qJD(3);
t85 = t87 ^ 2;
t86 = t94 * qJDD(2) + qJDD(3);
t28 = -t85 * pkin(3) + t86 * pkin(10) + t72 * t113 + t116;
t61 = t94 * t67;
t73 = (qJDD(2) * t98 + t101 * t119) * t91;
t30 = t74 * pkin(3) - t73 * pkin(10) + t61 + (-t51 + (pkin(3) * t98 - pkin(10) * t101) * t87 * qJD(2)) * t91;
t97 = sin(qJ(4));
t111 = t100 * t30 - t97 * t28;
t114 = t98 * t120;
t65 = t100 * t87 - t97 * t114;
t66 = t100 * t114 + t97 * t87;
t54 = -t65 * pkin(4) - t66 * pkin(11);
t68 = qJDD(4) + t74;
t80 = qJD(4) - t113;
t79 = t80 ^ 2;
t21 = -t68 * pkin(4) - t79 * pkin(11) + t66 * t54 - t111;
t132 = cos(qJ(5));
t47 = t65 * qJD(4) + t100 * t73 + t97 * t86;
t96 = sin(qJ(5));
t56 = t132 * t66 + t96 * t80;
t32 = t56 * qJD(5) - t132 * t68 + t96 * t47;
t55 = -t132 * t80 + t96 * t66;
t33 = -t55 * qJD(5) + t132 * t47 + t96 * t68;
t63 = qJD(5) - t65;
t40 = -t55 * mrSges(7,2) + t63 * mrSges(7,3);
t117 = m(7) * (-0.2e1 * qJD(6) * t56 + (t55 * t63 - t33) * qJ(6) + (t56 * t63 + t32) * pkin(5) + t21) + t32 * mrSges(7,1) + t55 * t40;
t41 = -t63 * mrSges(6,2) - t55 * mrSges(6,3);
t42 = t63 * mrSges(6,1) - t56 * mrSges(6,3);
t43 = -t63 * mrSges(7,1) + t56 * mrSges(7,2);
t135 = m(6) * t21 + t32 * mrSges(6,1) + (t42 - t43) * t56 + (mrSges(6,2) - mrSges(7,3)) * t33 + t55 * t41 + t117;
t124 = t100 * t28 + t97 * t30;
t22 = -t79 * pkin(4) + t68 * pkin(11) + t65 * t54 + t124;
t27 = -t86 * pkin(3) - t85 * pkin(10) + t72 * t114 - t136;
t46 = -t66 * qJD(4) + t100 * t86 - t97 * t73;
t24 = (-t65 * t80 - t47) * pkin(11) + (t66 * t80 - t46) * pkin(4) + t27;
t107 = t132 * t24 - t96 * t22;
t36 = t55 * pkin(5) - t56 * qJ(6);
t45 = qJDD(5) - t46;
t62 = t63 ^ 2;
t134 = m(7) * (-t45 * pkin(5) - t62 * qJ(6) + t56 * t36 + qJDD(6) - t107);
t126 = -mrSges(6,3) - mrSges(7,2);
t125 = t132 * t22 + t96 * t24;
t37 = t55 * mrSges(7,1) - t56 * mrSges(7,3);
t123 = -t55 * mrSges(6,1) - t56 * mrSges(6,2) - t37;
t118 = m(7) * (-t62 * pkin(5) + t45 * qJ(6) + 0.2e1 * qJD(6) * t63 - t55 * t36 + t125) + t63 * t43 + t45 * mrSges(7,3);
t14 = m(6) * t125 - t45 * mrSges(6,2) + t123 * t55 + t126 * t32 - t63 * t42 + t118;
t15 = m(6) * t107 - t134 + (t41 + t40) * t63 + t123 * t56 + (mrSges(6,1) + mrSges(7,1)) * t45 + t126 * t33;
t53 = -t65 * mrSges(5,1) + t66 * mrSges(5,2);
t58 = t80 * mrSges(5,1) - t66 * mrSges(5,3);
t12 = m(5) * t124 - t68 * mrSges(5,2) + t46 * mrSges(5,3) + t132 * t14 - t96 * t15 + t65 * t53 - t80 * t58;
t57 = -t80 * mrSges(5,2) + t65 * mrSges(5,3);
t13 = m(5) * t111 + t68 * mrSges(5,1) - t47 * mrSges(5,3) - t66 * t53 + t80 * t57 - t135;
t69 = t87 * mrSges(4,1) - mrSges(4,3) * t114;
t70 = -t87 * mrSges(4,2) + mrSges(4,3) * t113;
t10 = m(4) * (-t91 * t51 + t61) + t73 * mrSges(4,2) + t74 * mrSges(4,1) + t97 * t12 + t100 * t13 + (-t101 * t70 + t69 * t98) * t120;
t104 = m(5) * t27 - t46 * mrSges(5,1) + t47 * mrSges(5,2) + t132 * t15 + t96 * t14 - t65 * t57 + t66 * t58;
t71 = (-mrSges(4,1) * t101 + mrSges(4,2) * t98) * t120;
t11 = m(4) * t136 + t86 * mrSges(4,1) - t73 * mrSges(4,3) - t71 * t114 + t87 * t70 - t104;
t9 = m(4) * t116 - t86 * mrSges(4,2) - t74 * mrSges(4,3) + t100 * t12 + t71 * t113 - t97 * t13 - t87 * t69;
t109 = t101 * t11 + t98 * t9;
t4 = m(3) * t110 + qJDD(2) * mrSges(3,1) - t103 * mrSges(3,2) - t91 * t10 + t109 * t94;
t6 = m(3) * t67 + t94 * t10 + t109 * t91;
t8 = m(3) * t115 - t103 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t101 * t9 - t98 * t11;
t112 = m(2) * t89 + t4 * t121 + t8 * t128 + t95 * t6;
t2 = m(2) * t82 + t102 * t8 - t99 * t4;
t1 = m(2) * t81 - t92 * t6 + (t102 * t4 + t8 * t99) * t95;
t3 = [-m(1) * g(1) - t90 * t1 + t93 * t2, t2, t8, t9, t12, t14, -t32 * mrSges(7,2) - t55 * t37 + t118; -m(1) * g(2) + t93 * t1 + t90 * t2, t1, t4, t11, t13, t15, -t33 * mrSges(7,3) - t56 * t43 + t117; -m(1) * g(3) + t112, t112, t6, t10, t104, t135, -t45 * mrSges(7,1) + t33 * mrSges(7,2) + t56 * t37 - t63 * t40 + t134;];
f_new  = t3;
