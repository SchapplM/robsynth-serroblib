% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 11:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:27:24
% EndTime: 2019-05-07 11:28:07
% DurationCPUTime: 11.16s
% Computational Cost: add. (201959->215), mult. (447969->295), div. (0->0), fcn. (365701->14), ass. (0->115)
t110 = cos(pkin(6));
t147 = t110 * g(3);
t112 = sin(qJ(5));
t117 = cos(qJ(5));
t108 = sin(pkin(6));
t119 = cos(qJ(2));
t139 = qJD(1) * t119;
t134 = t108 * t139;
t100 = qJD(3) - t134;
t107 = sin(pkin(12));
t109 = cos(pkin(12));
t113 = sin(qJ(3));
t118 = cos(qJ(3));
t104 = t110 * qJD(1) + qJD(2);
t102 = t104 ^ 2;
t103 = t110 * qJDD(1) + qJDD(2);
t114 = sin(qJ(2));
t121 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t120 = cos(qJ(1));
t133 = t115 * g(1) - t120 * g(2);
t92 = t121 * t108 * pkin(8) + qJDD(1) * pkin(1) + t133;
t143 = t110 * t92;
t130 = -t120 * g(1) - t115 * g(2);
t138 = qJDD(1) * t108;
t93 = -t121 * pkin(1) + pkin(8) * t138 + t130;
t144 = t114 * t143 + t119 * t93;
t140 = qJD(1) * t108;
t95 = (-pkin(2) * t119 - pkin(9) * t114) * t140;
t64 = -t102 * pkin(2) + t103 * pkin(9) + (-g(3) * t114 + t139 * t95) * t108 + t144;
t96 = (qJD(2) * t139 + qJDD(1) * t114) * t108;
t135 = t114 * t140;
t97 = -qJD(2) * t135 + t119 * t138;
t65 = -t97 * pkin(2) - t96 * pkin(9) - t147 + (-t92 + (pkin(2) * t114 - pkin(9) * t119) * t104 * qJD(1)) * t108;
t132 = -t113 * t64 + t118 * t65;
t84 = t118 * t104 - t113 * t135;
t71 = t84 * qJD(3) + t113 * t103 + t118 * t96;
t85 = t113 * t104 + t118 * t135;
t89 = qJDD(3) - t97;
t36 = (t100 * t84 - t71) * qJ(4) + (t84 * t85 + t89) * pkin(3) + t132;
t145 = t113 * t65 + t118 * t64;
t70 = -t85 * qJD(3) + t118 * t103 - t113 * t96;
t79 = t100 * pkin(3) - t85 * qJ(4);
t83 = t84 ^ 2;
t38 = -t83 * pkin(3) + t70 * qJ(4) - t100 * t79 + t145;
t76 = t107 * t84 + t109 * t85;
t131 = -0.2e1 * qJD(4) * t76 - t107 * t38 + t109 * t36;
t53 = t107 * t70 + t109 * t71;
t75 = -t107 * t85 + t109 * t84;
t23 = (t100 * t75 - t53) * pkin(10) + (t75 * t76 + t89) * pkin(4) + t131;
t136 = 0.2e1 * qJD(4) * t75 + t107 * t36 + t109 * t38;
t52 = -t107 * t71 + t109 * t70;
t68 = t100 * pkin(4) - t76 * pkin(10);
t74 = t75 ^ 2;
t25 = -t74 * pkin(4) + t52 * pkin(10) - t100 * t68 + t136;
t146 = t112 * t23 + t117 * t25;
t142 = t108 * t114;
t141 = t108 * t119;
t128 = -g(3) * t141 - t114 * t93 + t119 * t143;
t63 = -t103 * pkin(2) - t102 * pkin(9) + t95 * t135 - t128;
t124 = -t70 * pkin(3) - t83 * qJ(4) + t85 * t79 + qJDD(4) + t63;
t111 = sin(qJ(6));
t116 = cos(qJ(6));
t123 = -t52 * pkin(4) - t74 * pkin(10) + t76 * t68 + t124;
t57 = -t112 * t76 + t117 * t75;
t58 = t112 * t75 + t117 * t76;
t46 = -t57 * pkin(5) - t58 * pkin(11);
t88 = qJDD(5) + t89;
t99 = qJD(5) + t100;
t98 = t99 ^ 2;
t20 = -t98 * pkin(5) + t88 * pkin(11) + t57 * t46 + t146;
t34 = -t58 * qJD(5) - t112 * t53 + t117 * t52;
t35 = t57 * qJD(5) + t112 * t52 + t117 * t53;
t21 = t123 + (-t57 * t99 - t35) * pkin(11) + (t58 * t99 - t34) * pkin(5);
t47 = -t111 * t58 + t116 * t99;
t29 = t47 * qJD(6) + t111 * t88 + t116 * t35;
t33 = qJDD(6) - t34;
t48 = t111 * t99 + t116 * t58;
t39 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t56 = qJD(6) - t57;
t40 = -t56 * mrSges(7,2) + t47 * mrSges(7,3);
t17 = m(7) * (-t111 * t20 + t116 * t21) - t29 * mrSges(7,3) + t33 * mrSges(7,1) - t48 * t39 + t56 * t40;
t28 = -t48 * qJD(6) - t111 * t35 + t116 * t88;
t41 = t56 * mrSges(7,1) - t48 * mrSges(7,3);
t18 = m(7) * (t111 * t21 + t116 * t20) + t28 * mrSges(7,3) - t33 * mrSges(7,2) + t47 * t39 - t56 * t41;
t49 = -t99 * mrSges(6,2) + t57 * mrSges(6,3);
t50 = t99 * mrSges(6,1) - t58 * mrSges(6,3);
t127 = -m(6) * t123 + t34 * mrSges(6,1) - t35 * mrSges(6,2) - t111 * t18 - t116 * t17 + t57 * t49 - t58 * t50;
t66 = -t100 * mrSges(5,2) + t75 * mrSges(5,3);
t67 = t100 * mrSges(5,1) - t76 * mrSges(5,3);
t125 = -m(5) * t124 + t52 * mrSges(5,1) - t53 * mrSges(5,2) + t75 * t66 - t76 * t67 + t127;
t78 = -t100 * mrSges(4,2) + t84 * mrSges(4,3);
t80 = t100 * mrSges(4,1) - t85 * mrSges(4,3);
t122 = m(4) * t63 - t70 * mrSges(4,1) + t71 * mrSges(4,2) - t84 * t78 + t85 * t80 - t125;
t91 = -t104 * mrSges(3,2) + mrSges(3,3) * t134;
t94 = (-mrSges(3,1) * t119 + mrSges(3,2) * t114) * t140;
t13 = m(3) * t128 + t103 * mrSges(3,1) - t96 * mrSges(3,3) + t104 * t91 - t135 * t94 - t122;
t45 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t11 = m(6) * t146 - t88 * mrSges(6,2) + t34 * mrSges(6,3) - t111 * t17 + t116 * t18 + t57 * t45 - t99 * t50;
t129 = -t112 * t25 + t117 * t23;
t126 = m(7) * (-t88 * pkin(5) - t98 * pkin(11) + t58 * t46 - t129) - t28 * mrSges(7,1) + t29 * mrSges(7,2) - t47 * t40 + t48 * t41;
t14 = m(6) * t129 + t88 * mrSges(6,1) - t35 * mrSges(6,3) - t58 * t45 + t99 * t49 - t126;
t59 = -t75 * mrSges(5,1) + t76 * mrSges(5,2);
t10 = m(5) * t136 - t89 * mrSges(5,2) + t52 * mrSges(5,3) - t100 * t67 + t117 * t11 - t112 * t14 + t75 * t59;
t77 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t9 = m(5) * t131 + t89 * mrSges(5,1) - t53 * mrSges(5,3) + t100 * t66 + t112 * t11 + t117 * t14 - t76 * t59;
t7 = m(4) * t132 + t89 * mrSges(4,1) - t71 * mrSges(4,3) + t107 * t10 + t100 * t78 + t109 * t9 - t85 * t77;
t8 = m(4) * t145 - t89 * mrSges(4,2) + t70 * mrSges(4,3) + t109 * t10 - t100 * t80 - t107 * t9 + t84 * t77;
t90 = t104 * mrSges(3,1) - mrSges(3,3) * t135;
t4 = m(3) * (-g(3) * t142 + t144) + t97 * mrSges(3,3) - t103 * mrSges(3,2) + t94 * t134 - t104 * t90 + t118 * t8 - t113 * t7;
t6 = m(3) * (-t108 * t92 - t147) + t96 * mrSges(3,2) - t97 * mrSges(3,1) + t113 * t8 + t118 * t7 + (t114 * t90 - t119 * t91) * t140;
t137 = t110 * t6 + t13 * t141 + t4 * t142;
t2 = m(2) * t130 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t114 * t13 + t119 * t4;
t1 = m(2) * t133 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t108 * t6 + (t114 * t4 + t119 * t13) * t110;
t3 = [-m(1) * g(1) - t115 * t1 + t120 * t2, t2, t4, t8, t10, t11, t18; -m(1) * g(2) + t120 * t1 + t115 * t2, t1, t13, t7, t9, t14, t17; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t6, t122, -t125, -t127, t126;];
f_new  = t3;
