% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:00:47
% EndTime: 2019-05-06 23:01:09
% DurationCPUTime: 10.73s
% Computational Cost: add. (193161->215), mult. (440155->295), div. (0->0), fcn. (360910->14), ass. (0->115)
t110 = cos(pkin(6));
t147 = t110 * g(3);
t112 = sin(qJ(5));
t117 = cos(qJ(5));
t108 = sin(pkin(6));
t119 = cos(qJ(2));
t139 = qJD(1) * t119;
t134 = t108 * t139;
t100 = qJD(4) - t134;
t113 = sin(qJ(4));
t118 = cos(qJ(4));
t107 = sin(pkin(12));
t109 = cos(pkin(12));
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
t94 = (-pkin(2) * t119 - qJ(3) * t114) * t140;
t64 = -t102 * pkin(2) + t103 * qJ(3) + (-g(3) * t114 + t94 * t139) * t108 + t144;
t96 = (qJD(2) * t139 + qJDD(1) * t114) * t108;
t135 = t114 * t140;
t97 = -qJD(2) * t135 + t119 * t138;
t65 = -t97 * pkin(2) - t147 - t96 * qJ(3) + (-t92 + (pkin(2) * t114 - qJ(3) * t119) * t104 * qJD(1)) * t108;
t85 = t107 * t104 + t109 * t135;
t131 = -0.2e1 * qJD(3) * t85 - t107 * t64 + t109 * t65;
t77 = t107 * t103 + t109 * t96;
t84 = t109 * t104 - t107 * t135;
t36 = (-t84 * t134 - t77) * pkin(9) + (t84 * t85 - t97) * pkin(3) + t131;
t136 = 0.2e1 * qJD(3) * t84 + t107 * t65 + t109 * t64;
t76 = t109 * t103 - t107 * t96;
t78 = -pkin(3) * t134 - t85 * pkin(9);
t83 = t84 ^ 2;
t38 = -t83 * pkin(3) + t76 * pkin(9) + t78 * t134 + t136;
t132 = -t113 * t38 + t118 * t36;
t71 = -t113 * t85 + t118 * t84;
t51 = t71 * qJD(4) + t113 * t76 + t118 * t77;
t72 = t113 * t84 + t118 * t85;
t89 = qJDD(4) - t97;
t23 = (t100 * t71 - t51) * pkin(10) + (t71 * t72 + t89) * pkin(4) + t132;
t145 = t113 * t36 + t118 * t38;
t50 = -t72 * qJD(4) - t113 * t77 + t118 * t76;
t68 = t100 * pkin(4) - t72 * pkin(10);
t70 = t71 ^ 2;
t25 = -t70 * pkin(4) + t50 * pkin(10) - t100 * t68 + t145;
t146 = t112 * t23 + t117 * t25;
t142 = t108 * t114;
t141 = t108 * t119;
t128 = -g(3) * t141 - t114 * t93 + t119 * t143;
t63 = -t103 * pkin(2) - t102 * qJ(3) + t94 * t135 + qJDD(3) - t128;
t124 = -t76 * pkin(3) - t83 * pkin(9) + t85 * t78 + t63;
t111 = sin(qJ(6));
t116 = cos(qJ(6));
t123 = -t50 * pkin(4) - t70 * pkin(10) + t72 * t68 + t124;
t57 = -t112 * t72 + t117 * t71;
t58 = t112 * t71 + t117 * t72;
t44 = -t57 * pkin(5) - t58 * pkin(11);
t88 = qJDD(5) + t89;
t99 = qJD(5) + t100;
t98 = t99 ^ 2;
t20 = -t98 * pkin(5) + t88 * pkin(11) + t57 * t44 + t146;
t32 = -t58 * qJD(5) - t112 * t51 + t117 * t50;
t33 = t57 * qJD(5) + t112 * t50 + t117 * t51;
t21 = t123 + (-t57 * t99 - t33) * pkin(11) + (t58 * t99 - t32) * pkin(5);
t47 = -t111 * t58 + t116 * t99;
t27 = t47 * qJD(6) + t111 * t88 + t116 * t33;
t31 = qJDD(6) - t32;
t48 = t111 * t99 + t116 * t58;
t39 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t56 = qJD(6) - t57;
t40 = -t56 * mrSges(7,2) + t47 * mrSges(7,3);
t17 = m(7) * (-t111 * t20 + t116 * t21) - t27 * mrSges(7,3) + t31 * mrSges(7,1) - t48 * t39 + t56 * t40;
t26 = -t48 * qJD(6) - t111 * t33 + t116 * t88;
t41 = t56 * mrSges(7,1) - t48 * mrSges(7,3);
t18 = m(7) * (t111 * t21 + t116 * t20) + t26 * mrSges(7,3) - t31 * mrSges(7,2) + t47 * t39 - t56 * t41;
t52 = -t99 * mrSges(6,2) + t57 * mrSges(6,3);
t53 = t99 * mrSges(6,1) - t58 * mrSges(6,3);
t127 = -m(6) * t123 + t32 * mrSges(6,1) - t33 * mrSges(6,2) - t111 * t18 - t116 * t17 + t57 * t52 - t58 * t53;
t66 = -t100 * mrSges(5,2) + t71 * mrSges(5,3);
t67 = t100 * mrSges(5,1) - t72 * mrSges(5,3);
t125 = -m(5) * t124 + t50 * mrSges(5,1) - t51 * mrSges(5,2) + t71 * t66 - t72 * t67 + t127;
t74 = mrSges(4,2) * t134 + t84 * mrSges(4,3);
t75 = -mrSges(4,1) * t134 - t85 * mrSges(4,3);
t122 = m(4) * t63 - t76 * mrSges(4,1) + t77 * mrSges(4,2) - t84 * t74 + t85 * t75 - t125;
t91 = -t104 * mrSges(3,2) + mrSges(3,3) * t134;
t95 = (-mrSges(3,1) * t119 + mrSges(3,2) * t114) * t140;
t12 = m(3) * t128 + t103 * mrSges(3,1) - t96 * mrSges(3,3) + t104 * t91 - t95 * t135 - t122;
t43 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t13 = m(6) * t146 - t88 * mrSges(6,2) + t32 * mrSges(6,3) - t111 * t17 + t116 * t18 + t57 * t43 - t99 * t53;
t129 = -t112 * t25 + t117 * t23;
t126 = m(7) * (-t88 * pkin(5) - t98 * pkin(11) + t58 * t44 - t129) - t26 * mrSges(7,1) + t27 * mrSges(7,2) - t47 * t40 + t48 * t41;
t14 = m(6) * t129 + t88 * mrSges(6,1) - t33 * mrSges(6,3) - t58 * t43 + t99 * t52 - t126;
t59 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t10 = m(5) * t145 - t89 * mrSges(5,2) + t50 * mrSges(5,3) - t100 * t67 - t112 * t14 + t117 * t13 + t71 * t59;
t73 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t9 = m(5) * t132 + t89 * mrSges(5,1) - t51 * mrSges(5,3) + t100 * t66 + t112 * t13 + t117 * t14 - t72 * t59;
t7 = m(4) * t131 - t97 * mrSges(4,1) - t77 * mrSges(4,3) + t113 * t10 + t118 * t9 - t74 * t134 - t85 * t73;
t8 = m(4) * t136 + t97 * mrSges(4,2) + t76 * mrSges(4,3) + t118 * t10 - t113 * t9 + t75 * t134 + t84 * t73;
t90 = t104 * mrSges(3,1) - mrSges(3,3) * t135;
t4 = m(3) * (-g(3) * t142 + t144) + t97 * mrSges(3,3) - t103 * mrSges(3,2) + t95 * t134 - t104 * t90 + t109 * t8 - t107 * t7;
t6 = m(3) * (-t108 * t92 - t147) + t96 * mrSges(3,2) - t97 * mrSges(3,1) + t107 * t8 + t109 * t7 + (t114 * t90 - t119 * t91) * t140;
t137 = t110 * t6 + t12 * t141 + t4 * t142;
t2 = m(2) * t130 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t114 * t12 + t119 * t4;
t1 = m(2) * t133 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t108 * t6 + (t114 * t4 + t119 * t12) * t110;
t3 = [-m(1) * g(1) - t115 * t1 + t120 * t2, t2, t4, t8, t10, t13, t18; -m(1) * g(2) + t120 * t1 + t115 * t2, t1, t12, t7, t9, t14, t17; (-m(1) - m(2)) * g(3) + t137, -m(2) * g(3) + t137, t6, t122, -t125, -t127, t126;];
f_new  = t3;
