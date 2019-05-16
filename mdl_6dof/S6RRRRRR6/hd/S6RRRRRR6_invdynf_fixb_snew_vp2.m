% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 11:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 11:07:11
% EndTime: 2019-05-08 11:07:39
% DurationCPUTime: 10.28s
% Computational Cost: add. (199669->216), mult. (422222->293), div. (0->0), fcn. (346508->14), ass. (0->117)
t112 = cos(pkin(6));
t150 = t112 * g(3);
t114 = sin(qJ(5));
t120 = cos(qJ(5));
t111 = sin(pkin(6));
t123 = cos(qJ(2));
t141 = qJD(1) * t123;
t137 = t111 * t141;
t104 = qJD(3) - t137;
t103 = qJD(4) + t104;
t102 = t103 ^ 2;
t115 = sin(qJ(4));
t121 = cos(qJ(4));
t116 = sin(qJ(3));
t122 = cos(qJ(3));
t108 = t112 * qJD(1) + qJD(2);
t106 = t108 ^ 2;
t107 = t112 * qJDD(1) + qJDD(2);
t117 = sin(qJ(2));
t125 = qJD(1) ^ 2;
t118 = sin(qJ(1));
t124 = cos(qJ(1));
t136 = t118 * g(1) - t124 * g(2);
t96 = t125 * t111 * pkin(8) + qJDD(1) * pkin(1) + t136;
t145 = t112 * t96;
t132 = -t124 * g(1) - t118 * g(2);
t139 = qJDD(1) * t111;
t97 = -t125 * pkin(1) + pkin(8) * t139 + t132;
t146 = t117 * t145 + t123 * t97;
t142 = qJD(1) * t111;
t99 = (-pkin(2) * t123 - pkin(9) * t117) * t142;
t64 = -t106 * pkin(2) + t107 * pkin(9) + (-g(3) * t117 + t141 * t99) * t111 + t146;
t100 = (qJD(2) * t141 + qJDD(1) * t117) * t111;
t138 = t117 * t142;
t101 = -qJD(2) * t138 + t123 * t139;
t65 = -t101 * pkin(2) - t100 * pkin(9) - t150 + (-t96 + (pkin(2) * t117 - pkin(9) * t123) * t108 * qJD(1)) * t111;
t134 = -t116 * t64 + t122 * t65;
t88 = t122 * t108 - t116 * t138;
t74 = t88 * qJD(3) + t122 * t100 + t116 * t107;
t89 = t116 * t108 + t122 * t138;
t93 = qJDD(3) - t101;
t33 = (t104 * t88 - t74) * pkin(10) + (t88 * t89 + t93) * pkin(3) + t134;
t147 = t116 * t65 + t122 * t64;
t73 = -t89 * qJD(3) - t116 * t100 + t122 * t107;
t83 = t104 * pkin(3) - t89 * pkin(10);
t87 = t88 ^ 2;
t37 = -t87 * pkin(3) + t73 * pkin(10) - t104 * t83 + t147;
t148 = t115 * t33 + t121 * t37;
t78 = -t115 * t89 + t121 * t88;
t79 = t115 * t88 + t121 * t89;
t59 = -t78 * pkin(4) - t79 * pkin(11);
t92 = qJDD(4) + t93;
t25 = -t102 * pkin(4) + t92 * pkin(11) + t78 * t59 + t148;
t143 = t111 * t123;
t131 = -g(3) * t143 - t117 * t97 + t123 * t145;
t63 = -t107 * pkin(2) - t106 * pkin(9) + t99 * t138 - t131;
t128 = -t73 * pkin(3) - t87 * pkin(10) + t89 * t83 + t63;
t48 = -t79 * qJD(4) - t115 * t74 + t121 * t73;
t49 = t78 * qJD(4) + t115 * t73 + t121 * t74;
t28 = (-t103 * t78 - t49) * pkin(11) + (t103 * t79 - t48) * pkin(4) + t128;
t149 = t114 * t28 + t120 * t25;
t144 = t111 * t117;
t113 = sin(qJ(6));
t119 = cos(qJ(6));
t135 = -t114 * t25 + t120 * t28;
t67 = t120 * t103 - t114 * t79;
t40 = t67 * qJD(5) + t114 * t92 + t120 * t49;
t47 = qJDD(5) - t48;
t68 = t114 * t103 + t120 * t79;
t77 = qJD(5) - t78;
t19 = (t67 * t77 - t40) * pkin(12) + (t67 * t68 + t47) * pkin(5) + t135;
t39 = -t68 * qJD(5) - t114 * t49 + t120 * t92;
t56 = t77 * pkin(5) - t68 * pkin(12);
t66 = t67 ^ 2;
t20 = -t66 * pkin(5) + t39 * pkin(12) - t77 * t56 + t149;
t51 = -t113 * t68 + t119 * t67;
t31 = t51 * qJD(6) + t113 * t39 + t119 * t40;
t52 = t113 * t67 + t119 * t68;
t38 = -t51 * mrSges(7,1) + t52 * mrSges(7,2);
t75 = qJD(6) + t77;
t43 = -t75 * mrSges(7,2) + t51 * mrSges(7,3);
t45 = qJDD(6) + t47;
t17 = m(7) * (-t113 * t20 + t119 * t19) - t31 * mrSges(7,3) + t45 * mrSges(7,1) - t52 * t38 + t75 * t43;
t30 = -t52 * qJD(6) - t113 * t40 + t119 * t39;
t44 = t75 * mrSges(7,1) - t52 * mrSges(7,3);
t18 = m(7) * (t113 * t19 + t119 * t20) + t30 * mrSges(7,3) - t45 * mrSges(7,2) + t51 * t38 - t75 * t44;
t53 = -t67 * mrSges(6,1) + t68 * mrSges(6,2);
t54 = -t77 * mrSges(6,2) + t67 * mrSges(6,3);
t14 = m(6) * t135 + t47 * mrSges(6,1) - t40 * mrSges(6,3) + t113 * t18 + t119 * t17 - t68 * t53 + t77 * t54;
t55 = t77 * mrSges(6,1) - t68 * mrSges(6,3);
t15 = m(6) * t149 - t47 * mrSges(6,2) + t39 * mrSges(6,3) - t113 * t17 + t119 * t18 + t67 * t53 - t77 * t55;
t69 = -t103 * mrSges(5,2) + t78 * mrSges(5,3);
t70 = t103 * mrSges(5,1) - t79 * mrSges(5,3);
t129 = -m(5) * t128 + t48 * mrSges(5,1) - t49 * mrSges(5,2) - t114 * t15 - t120 * t14 + t78 * t69 - t79 * t70;
t81 = -t104 * mrSges(4,2) + t88 * mrSges(4,3);
t82 = t104 * mrSges(4,1) - t89 * mrSges(4,3);
t126 = m(4) * t63 - t73 * mrSges(4,1) + t74 * mrSges(4,2) - t88 * t81 + t89 * t82 - t129;
t95 = -t108 * mrSges(3,2) + mrSges(3,3) * t137;
t98 = (-mrSges(3,1) * t123 + mrSges(3,2) * t117) * t142;
t10 = m(3) * t131 + t107 * mrSges(3,1) - t100 * mrSges(3,3) + t108 * t95 - t138 * t98 - t126;
t58 = -t78 * mrSges(5,1) + t79 * mrSges(5,2);
t11 = m(5) * t148 - t92 * mrSges(5,2) + t48 * mrSges(5,3) - t103 * t70 - t114 * t14 + t120 * t15 + t78 * t58;
t133 = -t115 * t37 + t121 * t33;
t24 = -t92 * pkin(4) - t102 * pkin(11) + t79 * t59 - t133;
t130 = t30 * mrSges(7,1) + t51 * t43 - m(7) * (-t39 * pkin(5) - t66 * pkin(12) + t68 * t56 + t24) - t31 * mrSges(7,2) - t52 * t44;
t127 = m(6) * t24 - t39 * mrSges(6,1) + t40 * mrSges(6,2) - t67 * t54 + t68 * t55 - t130;
t16 = m(5) * t133 + t92 * mrSges(5,1) - t49 * mrSges(5,3) + t103 * t69 - t79 * t58 - t127;
t80 = -t88 * mrSges(4,1) + t89 * mrSges(4,2);
t7 = m(4) * t134 + t93 * mrSges(4,1) - t74 * mrSges(4,3) + t104 * t81 + t115 * t11 + t121 * t16 - t89 * t80;
t8 = m(4) * t147 - t93 * mrSges(4,2) + t73 * mrSges(4,3) - t104 * t82 + t121 * t11 - t115 * t16 + t88 * t80;
t94 = t108 * mrSges(3,1) - mrSges(3,3) * t138;
t4 = m(3) * (-g(3) * t144 + t146) + t101 * mrSges(3,3) - t107 * mrSges(3,2) + t98 * t137 - t108 * t94 + t122 * t8 - t116 * t7;
t6 = m(3) * (-t111 * t96 - t150) + t100 * mrSges(3,2) - t101 * mrSges(3,1) + t116 * t8 + t122 * t7 + (t117 * t94 - t123 * t95) * t142;
t140 = t10 * t143 + t112 * t6 + t4 * t144;
t2 = m(2) * t132 - t125 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t117 * t10 + t123 * t4;
t1 = m(2) * t136 + qJDD(1) * mrSges(2,1) - t125 * mrSges(2,2) - t111 * t6 + (t123 * t10 + t117 * t4) * t112;
t3 = [-m(1) * g(1) - t118 * t1 + t124 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t124 * t1 + t118 * t2, t1, t10, t7, t16, t14, t17; (-m(1) - m(2)) * g(3) + t140, -m(2) * g(3) + t140, t6, t126, -t129, t127, -t130;];
f_new  = t3;
