% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:31:54
% EndTime: 2019-05-05 12:32:22
% DurationCPUTime: 17.23s
% Computational Cost: add. (328119->197), mult. (771990->286), div. (0->0), fcn. (652045->18), ass. (0->121)
t101 = cos(pkin(7));
t113 = qJD(2) ^ 2;
t107 = sin(qJ(2));
t112 = cos(qJ(2));
t98 = sin(pkin(6));
t135 = t112 * t98;
t102 = cos(pkin(6));
t95 = sin(pkin(14));
t99 = cos(pkin(14));
t87 = t95 * g(1) - t99 * g(2);
t138 = t102 * t87;
t88 = -t99 * g(1) - t95 * g(2);
t94 = -g(3) + qJDD(1);
t123 = -t107 * t88 + t112 * t138 + t94 * t135;
t97 = sin(pkin(7));
t147 = pkin(10) * t97;
t63 = qJDD(2) * pkin(2) + t113 * t147 + t123;
t77 = t102 * t94 - t98 * t87;
t149 = t101 * t63 + t77 * t97;
t105 = sin(qJ(4));
t110 = cos(qJ(4));
t100 = cos(pkin(8));
t106 = sin(qJ(3));
t134 = qJD(2) * t97;
t111 = cos(qJ(3));
t141 = t149 * t111;
t145 = pkin(11) * t100;
t136 = t107 * t98;
t129 = t107 * t138 + t112 * t88 + t94 * t136;
t64 = -t113 * pkin(2) + qJDD(2) * t147 + t129;
t126 = t111 * t134;
t92 = t101 * qJD(2) + qJD(3);
t96 = sin(pkin(8));
t143 = t92 * t96;
t74 = (t100 * t126 + t143) * pkin(11);
t146 = pkin(11) * t96;
t78 = (-pkin(3) * t111 - t106 * t146) * t134;
t132 = qJD(2) * qJD(3);
t82 = (qJDD(2) * t106 + t111 * t132) * t97;
t91 = t101 * qJDD(2) + qJDD(3);
t33 = -t82 * t145 + t91 * pkin(3) + t92 * t74 + (-t78 * t134 - t64) * t106 + t141;
t140 = t100 * t33;
t83 = (qJDD(2) * t111 - t106 * t132) * t97;
t121 = t100 * t83 + t91 * t96;
t130 = t149 * t106 + t111 * t64;
t127 = t106 * t134;
t76 = t92 * pkin(3) - t127 * t145;
t34 = t121 * pkin(11) + t78 * t126 - t92 * t76 + t130;
t72 = t101 * t77;
t42 = -t82 * t146 - t83 * pkin(3) + t72 + (-t63 + (t106 * t76 - t111 * t74) * qJD(2)) * t97;
t148 = -t105 * t34 + (t42 * t96 + t140) * t110;
t133 = t100 * t111;
t137 = t105 * t96;
t68 = t92 * t137 + (t105 * t133 + t106 * t110) * t134;
t51 = -t68 * qJD(4) - t105 * t82 + t121 * t110;
t67 = (-t105 * t106 + t110 * t133) * t134 + t110 * t143;
t104 = sin(qJ(5));
t109 = cos(qJ(5));
t131 = t105 * t140 + t110 * t34 + t42 * t137;
t54 = -t67 * pkin(4) - t68 * pkin(12);
t69 = t100 * t91 - t96 * t83 + qJDD(4);
t75 = t100 * t92 - t96 * t126 + qJD(4);
t73 = t75 ^ 2;
t24 = -t73 * pkin(4) + t69 * pkin(12) + t67 * t54 + t131;
t124 = t100 * t42 - t96 * t33;
t52 = t67 * qJD(4) + t121 * t105 + t110 * t82;
t26 = (-t67 * t75 - t52) * pkin(12) + (t68 * t75 - t51) * pkin(4) + t124;
t142 = t104 * t26 + t109 * t24;
t103 = sin(qJ(6));
t108 = cos(qJ(6));
t56 = -t104 * t68 + t109 * t75;
t57 = t104 * t75 + t109 * t68;
t44 = -t56 * pkin(5) - t57 * pkin(13);
t50 = qJDD(5) - t51;
t66 = qJD(5) - t67;
t65 = t66 ^ 2;
t20 = -t65 * pkin(5) + t50 * pkin(13) + t56 * t44 + t142;
t23 = -t69 * pkin(4) - t73 * pkin(12) + t68 * t54 - t148;
t36 = -t57 * qJD(5) - t104 * t52 + t109 * t69;
t37 = t56 * qJD(5) + t104 * t69 + t109 * t52;
t21 = (-t56 * t66 - t37) * pkin(13) + (t57 * t66 - t36) * pkin(5) + t23;
t46 = -t103 * t57 + t108 * t66;
t28 = t46 * qJD(6) + t103 * t50 + t108 * t37;
t47 = t103 * t66 + t108 * t57;
t29 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t35 = qJDD(6) - t36;
t55 = qJD(6) - t56;
t38 = -t55 * mrSges(7,2) + t46 * mrSges(7,3);
t17 = m(7) * (-t103 * t20 + t108 * t21) - t28 * mrSges(7,3) + t35 * mrSges(7,1) - t47 * t29 + t55 * t38;
t27 = -t47 * qJD(6) - t103 * t37 + t108 * t50;
t39 = t55 * mrSges(7,1) - t47 * mrSges(7,3);
t18 = m(7) * (t103 * t21 + t108 * t20) + t27 * mrSges(7,3) - t35 * mrSges(7,2) + t46 * t29 - t55 * t39;
t43 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t49 = t66 * mrSges(6,1) - t57 * mrSges(6,3);
t15 = m(6) * t142 - t50 * mrSges(6,2) + t36 * mrSges(6,3) - t103 * t17 + t108 * t18 + t56 * t43 - t66 * t49;
t119 = -t104 * t24 + t109 * t26;
t115 = m(7) * (-t50 * pkin(5) - t65 * pkin(13) + t57 * t44 - t119) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t46 * t38 + t47 * t39;
t48 = -t66 * mrSges(6,2) + t56 * mrSges(6,3);
t16 = m(6) * t119 + t50 * mrSges(6,1) - t37 * mrSges(6,3) - t57 * t43 + t66 * t48 - t115;
t53 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t59 = t75 * mrSges(5,1) - t68 * mrSges(5,3);
t12 = m(5) * t131 - t69 * mrSges(5,2) + t51 * mrSges(5,3) - t104 * t16 + t109 * t15 + t67 * t53 - t75 * t59;
t114 = m(6) * t23 - t36 * mrSges(6,1) + t37 * mrSges(6,2) + t103 * t18 + t108 * t17 - t56 * t48 + t57 * t49;
t58 = -t75 * mrSges(5,2) + t67 * mrSges(5,3);
t14 = m(5) * t148 + t69 * mrSges(5,1) - t52 * mrSges(5,3) - t68 * t53 + t75 * t58 - t114;
t118 = t105 * t12 + t110 * t14;
t13 = m(5) * t124 - t51 * mrSges(5,1) + t52 * mrSges(5,2) + t104 * t15 + t109 * t16 - t67 * t58 + t68 * t59;
t79 = t92 * mrSges(4,1) - mrSges(4,3) * t127;
t80 = -t92 * mrSges(4,2) + mrSges(4,3) * t126;
t10 = m(4) * (-t97 * t63 + t72) + t82 * mrSges(4,2) - t83 * mrSges(4,1) + t100 * t13 + t118 * t96 + (t106 * t79 - t111 * t80) * t134;
t81 = (-mrSges(4,1) * t111 + mrSges(4,2) * t106) * t134;
t11 = m(4) * t130 - t91 * mrSges(4,2) + t83 * mrSges(4,3) - t105 * t14 + t110 * t12 + t81 * t126 - t92 * t79;
t9 = m(4) * (-t106 * t64 + t141) - t82 * mrSges(4,3) + t91 * mrSges(4,1) - t81 * t127 + t92 * t80 - t96 * t13 + t118 * t100;
t120 = t106 * t11 + t111 * t9;
t4 = m(3) * t123 + qJDD(2) * mrSges(3,1) - t113 * mrSges(3,2) - t97 * t10 + t120 * t101;
t6 = m(3) * t77 + t101 * t10 + t120 * t97;
t8 = m(3) * t129 - t113 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t106 * t9 + t111 * t11;
t125 = m(2) * t94 + t102 * t6 + t4 * t135 + t8 * t136;
t2 = m(2) * t88 - t107 * t4 + t112 * t8;
t1 = m(2) * t87 - t98 * t6 + (t107 * t8 + t112 * t4) * t102;
t3 = [-m(1) * g(1) - t95 * t1 + t99 * t2, t2, t8, t11, t12, t15, t18; -m(1) * g(2) + t99 * t1 + t95 * t2, t1, t4, t9, t14, t16, t17; -m(1) * g(3) + t125, t125, t6, t10, t13, t114, t115;];
f_new  = t3;
