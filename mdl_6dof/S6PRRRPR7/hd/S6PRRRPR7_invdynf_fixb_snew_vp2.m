% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR7
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 08:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 08:42:28
% EndTime: 2019-05-05 08:42:41
% DurationCPUTime: 7.67s
% Computational Cost: add. (145128->183), mult. (304580->260), div. (0->0), fcn. (245242->16), ass. (0->107)
t103 = cos(pkin(7));
t112 = qJD(2) ^ 2;
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t100 = sin(pkin(6));
t129 = t100 * t111;
t104 = cos(pkin(6));
t102 = cos(pkin(12));
t98 = sin(pkin(12));
t88 = t98 * g(1) - t102 * g(2);
t132 = t104 * t88;
t89 = -t102 * g(1) - t98 * g(2);
t96 = -g(3) + qJDD(1);
t120 = -t108 * t89 + t111 * t132 + t96 * t129;
t99 = sin(pkin(7));
t137 = pkin(9) * t99;
t56 = qJDD(2) * pkin(2) + t112 * t137 + t120;
t73 = -t100 * t88 + t104 * t96;
t139 = t103 * t56 + t73 * t99;
t107 = sin(qJ(3));
t110 = cos(qJ(3));
t128 = qJD(2) * qJD(3);
t80 = (-qJDD(2) * t110 + t107 * t128) * t99;
t130 = t100 * t108;
t125 = t108 * t132 + t111 * t89 + t96 * t130;
t57 = -t112 * pkin(2) + qJDD(2) * t137 + t125;
t138 = -t107 * t57 + t139 * t110;
t136 = cos(qJ(4));
t106 = sin(qJ(4));
t131 = qJD(2) * t99;
t123 = t110 * t131;
t126 = t139 * t107 + t110 * t57;
t78 = (-pkin(3) * t110 - pkin(10) * t107) * t131;
t94 = t103 * qJD(2) + qJD(3);
t92 = t94 ^ 2;
t93 = t103 * qJDD(2) + qJDD(3);
t35 = -t92 * pkin(3) + t93 * pkin(10) + t78 * t123 + t126;
t69 = t103 * t73;
t79 = (qJDD(2) * t107 + t110 * t128) * t99;
t38 = t80 * pkin(3) - t79 * pkin(10) + t69 + (-t56 + (pkin(3) * t107 - pkin(10) * t110) * t94 * qJD(2)) * t99;
t134 = t106 * t38 + t136 * t35;
t101 = cos(pkin(13));
t124 = t107 * t131;
t71 = t106 * t124 - t136 * t94;
t72 = t106 * t94 + t136 * t124;
t58 = t71 * pkin(4) - t72 * qJ(5);
t74 = qJDD(4) + t80;
t87 = qJD(4) - t123;
t86 = t87 ^ 2;
t24 = -t86 * pkin(4) + t74 * qJ(5) - t71 * t58 + t134;
t34 = -t93 * pkin(3) - t92 * pkin(10) + t78 * t124 - t138;
t51 = t72 * qJD(4) + t106 * t79 - t136 * t93;
t52 = -t71 * qJD(4) + t106 * t93 + t136 * t79;
t27 = (t71 * t87 - t52) * qJ(5) + (t72 * t87 + t51) * pkin(4) + t34;
t97 = sin(pkin(13));
t63 = t101 * t87 - t97 * t72;
t127 = 0.2e1 * qJD(5) * t63 + t101 * t24 + t97 * t27;
t105 = sin(qJ(6));
t109 = cos(qJ(6));
t64 = t101 * t72 + t97 * t87;
t121 = -0.2e1 * qJD(5) * t64 + t101 * t27 - t97 * t24;
t44 = t101 * t52 + t97 * t74;
t18 = (t63 * t71 - t44) * pkin(11) + (t63 * t64 + t51) * pkin(5) + t121;
t43 = t101 * t74 - t97 * t52;
t49 = t71 * pkin(5) - t64 * pkin(11);
t62 = t63 ^ 2;
t19 = -t62 * pkin(5) + t43 * pkin(11) - t71 * t49 + t127;
t41 = -t105 * t64 + t109 * t63;
t30 = t41 * qJD(6) + t105 * t43 + t109 * t44;
t42 = t105 * t63 + t109 * t64;
t36 = -t41 * mrSges(7,1) + t42 * mrSges(7,2);
t70 = qJD(6) + t71;
t39 = -t70 * mrSges(7,2) + t41 * mrSges(7,3);
t50 = qJDD(6) + t51;
t16 = m(7) * (-t105 * t19 + t109 * t18) - t30 * mrSges(7,3) + t50 * mrSges(7,1) - t42 * t36 + t70 * t39;
t29 = -t42 * qJD(6) - t105 * t44 + t109 * t43;
t40 = t70 * mrSges(7,1) - t42 * mrSges(7,3);
t17 = m(7) * (t105 * t18 + t109 * t19) + t29 * mrSges(7,3) - t50 * mrSges(7,2) + t41 * t36 - t70 * t40;
t45 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t47 = -t71 * mrSges(6,2) + t63 * mrSges(6,3);
t13 = m(6) * t121 + t51 * mrSges(6,1) - t44 * mrSges(6,3) + t105 * t17 + t109 * t16 - t64 * t45 + t71 * t47;
t48 = t71 * mrSges(6,1) - t64 * mrSges(6,3);
t14 = m(6) * t127 - t51 * mrSges(6,2) + t43 * mrSges(6,3) - t105 * t16 + t109 * t17 + t63 * t45 - t71 * t48;
t59 = t71 * mrSges(5,1) + t72 * mrSges(5,2);
t66 = t87 * mrSges(5,1) - t72 * mrSges(5,3);
t12 = m(5) * t134 - t74 * mrSges(5,2) - t51 * mrSges(5,3) + t101 * t14 - t97 * t13 - t71 * t59 - t87 * t66;
t119 = -t106 * t35 + t136 * t38;
t23 = -t74 * pkin(4) - t86 * qJ(5) + t72 * t58 + qJDD(5) - t119;
t115 = t29 * mrSges(7,1) + t41 * t39 - m(7) * (-t43 * pkin(5) - t62 * pkin(11) + t64 * t49 + t23) - t30 * mrSges(7,2) - t42 * t40;
t113 = m(6) * t23 - t43 * mrSges(6,1) + t44 * mrSges(6,2) - t63 * t47 + t64 * t48 - t115;
t65 = -t87 * mrSges(5,2) - t71 * mrSges(5,3);
t15 = m(5) * t119 + t74 * mrSges(5,1) - t52 * mrSges(5,3) - t72 * t59 + t87 * t65 - t113;
t75 = t94 * mrSges(4,1) - mrSges(4,3) * t124;
t76 = -t94 * mrSges(4,2) + mrSges(4,3) * t123;
t10 = m(4) * (-t99 * t56 + t69) + t79 * mrSges(4,2) + t80 * mrSges(4,1) + t106 * t12 + t136 * t15 + (t107 * t75 - t110 * t76) * t131;
t114 = m(5) * t34 + t51 * mrSges(5,1) + t52 * mrSges(5,2) + t101 * t13 + t97 * t14 + t71 * t65 + t72 * t66;
t77 = (-mrSges(4,1) * t110 + mrSges(4,2) * t107) * t131;
t11 = m(4) * t138 + t93 * mrSges(4,1) - t79 * mrSges(4,3) - t77 * t124 + t94 * t76 - t114;
t9 = m(4) * t126 - t93 * mrSges(4,2) - t80 * mrSges(4,3) - t106 * t15 + t136 * t12 + t77 * t123 - t94 * t75;
t117 = t107 * t9 + t110 * t11;
t4 = m(3) * t120 + qJDD(2) * mrSges(3,1) - t112 * mrSges(3,2) - t99 * t10 + t117 * t103;
t6 = m(3) * t73 + t103 * t10 + t117 * t99;
t8 = m(3) * t125 - t112 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t107 * t11 + t110 * t9;
t122 = m(2) * t96 + t104 * t6 + t4 * t129 + t8 * t130;
t2 = m(2) * t89 - t108 * t4 + t111 * t8;
t1 = m(2) * t88 - t100 * t6 + (t108 * t8 + t111 * t4) * t104;
t3 = [-m(1) * g(1) - t98 * t1 + t102 * t2, t2, t8, t9, t12, t14, t17; -m(1) * g(2) + t102 * t1 + t98 * t2, t1, t4, t11, t15, t13, t16; -m(1) * g(3) + t122, t122, t6, t10, t114, t113, -t115;];
f_new  = t3;
