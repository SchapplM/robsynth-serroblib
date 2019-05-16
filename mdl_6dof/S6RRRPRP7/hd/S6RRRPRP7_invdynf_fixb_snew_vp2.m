% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:12:53
% EndTime: 2019-05-07 08:13:08
% DurationCPUTime: 4.91s
% Computational Cost: add. (85542->211), mult. (188647->282), div. (0->0), fcn. (150376->12), ass. (0->107)
t103 = sin(pkin(6));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t130 = qJD(1) * qJD(2);
t93 = (-qJDD(1) * t111 + t108 * t130) * t103;
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t107 = sin(qJ(3));
t110 = cos(qJ(3));
t132 = qJD(1) * t111;
t105 = cos(pkin(6));
t113 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t124 = t109 * g(1) - t112 * g(2);
t146 = pkin(8) * t103;
t88 = qJDD(1) * pkin(1) + t113 * t146 + t124;
t137 = t105 * t88;
t121 = -t112 * g(1) - t109 * g(2);
t89 = -t113 * pkin(1) + qJDD(1) * t146 + t121;
t138 = t108 * t137 + t111 * t89;
t133 = qJD(1) * t103;
t91 = (-pkin(2) * t111 - pkin(9) * t108) * t133;
t99 = t105 * qJD(1) + qJD(2);
t97 = t99 ^ 2;
t98 = t105 * qJDD(1) + qJDD(2);
t58 = -t97 * pkin(2) + t98 * pkin(9) + (-g(3) * t108 + t91 * t132) * t103 + t138;
t145 = t105 * g(3);
t92 = (qJDD(1) * t108 + t111 * t130) * t103;
t59 = t93 * pkin(2) - t92 * pkin(9) - t145 + (-t88 + (pkin(2) * t108 - pkin(9) * t111) * t99 * qJD(1)) * t103;
t122 = -t107 * t58 + t110 * t59;
t126 = t108 * t133;
t81 = -t107 * t126 + t110 * t99;
t66 = t81 * qJD(3) + t107 * t98 + t110 * t92;
t82 = t107 * t99 + t110 * t126;
t85 = qJDD(3) + t93;
t125 = t103 * t132;
t96 = qJD(3) - t125;
t27 = (t81 * t96 - t66) * qJ(4) + (t81 * t82 + t85) * pkin(3) + t122;
t139 = t107 * t59 + t110 * t58;
t65 = -t82 * qJD(3) - t107 * t92 + t110 * t98;
t75 = t96 * pkin(3) - t82 * qJ(4);
t80 = t81 ^ 2;
t30 = -t80 * pkin(3) + t65 * qJ(4) - t96 * t75 + t139;
t72 = t102 * t81 + t104 * t82;
t150 = -0.2e1 * qJD(4) * t72 - t102 * t30 + t104 * t27;
t71 = -t102 * t82 + t104 * t81;
t53 = -t71 * pkin(4) - t72 * pkin(10);
t95 = t96 ^ 2;
t22 = -t85 * pkin(4) - t95 * pkin(10) + t72 * t53 - t150;
t106 = sin(qJ(5));
t147 = cos(qJ(5));
t50 = t102 * t65 + t104 * t66;
t61 = t106 * t96 + t147 * t72;
t34 = t61 * qJD(5) + t106 * t50 - t147 * t85;
t60 = t106 * t72 - t147 * t96;
t35 = -t60 * qJD(5) + t106 * t85 + t147 * t50;
t70 = qJD(5) - t71;
t42 = -t60 * mrSges(7,2) + t70 * mrSges(7,3);
t128 = m(7) * (-0.2e1 * qJD(6) * t61 + (t60 * t70 - t35) * qJ(6) + (t61 * t70 + t34) * pkin(5) + t22) + t34 * mrSges(7,1) + t60 * t42;
t43 = -t70 * mrSges(6,2) - t60 * mrSges(6,3);
t44 = t70 * mrSges(6,1) - t61 * mrSges(6,3);
t45 = -t70 * mrSges(7,1) + t61 * mrSges(7,2);
t149 = m(6) * t22 + t34 * mrSges(6,1) + (t44 - t45) * t61 + (mrSges(6,2) - mrSges(7,3)) * t35 + t60 * t43 + t128;
t127 = 0.2e1 * qJD(4) * t71 + t102 * t27 + t104 * t30;
t23 = -t95 * pkin(4) + t85 * pkin(10) + t71 * t53 + t127;
t134 = t103 * t111;
t120 = -g(3) * t134 - t108 * t89 + t111 * t137;
t57 = -t98 * pkin(2) - t97 * pkin(9) + t91 * t126 - t120;
t115 = -t65 * pkin(3) - t80 * qJ(4) + t82 * t75 + qJDD(4) + t57;
t49 = -t102 * t66 + t104 * t65;
t25 = (-t71 * t96 - t50) * pkin(10) + (t72 * t96 - t49) * pkin(4) + t115;
t119 = -t106 * t23 + t147 * t25;
t39 = t60 * pkin(5) - t61 * qJ(6);
t48 = qJDD(5) - t49;
t69 = t70 ^ 2;
t148 = m(7) * (-t48 * pkin(5) - t69 * qJ(6) + t61 * t39 + qJDD(6) - t119);
t143 = -mrSges(6,3) - mrSges(7,2);
t142 = t106 * t25 + t147 * t23;
t40 = t60 * mrSges(7,1) - t61 * mrSges(7,3);
t141 = -t60 * mrSges(6,1) - t61 * mrSges(6,2) - t40;
t135 = t103 * t108;
t129 = m(7) * (-t69 * pkin(5) + t48 * qJ(6) + 0.2e1 * qJD(6) * t70 - t60 * t39 + t142) + t70 * t45 + t48 * mrSges(7,3);
t14 = m(6) * t142 - t48 * mrSges(6,2) + t141 * t60 + t143 * t34 - t70 * t44 + t129;
t16 = m(6) * t119 - t148 + (t43 + t42) * t70 + t141 * t61 + (mrSges(6,1) + mrSges(7,1)) * t48 + t143 * t35;
t62 = -t96 * mrSges(5,2) + t71 * mrSges(5,3);
t63 = t96 * mrSges(5,1) - t72 * mrSges(5,3);
t117 = -m(5) * t115 + t49 * mrSges(5,1) - t50 * mrSges(5,2) - t106 * t14 - t147 * t16 + t71 * t62 - t72 * t63;
t74 = -t96 * mrSges(4,2) + t81 * mrSges(4,3);
t76 = t96 * mrSges(4,1) - t82 * mrSges(4,3);
t114 = m(4) * t57 - t65 * mrSges(4,1) + t66 * mrSges(4,2) - t81 * t74 + t82 * t76 - t117;
t87 = -t99 * mrSges(3,2) + mrSges(3,3) * t125;
t90 = (-mrSges(3,1) * t111 + mrSges(3,2) * t108) * t133;
t10 = m(3) * t120 + t98 * mrSges(3,1) - t92 * mrSges(3,3) - t90 * t126 + t99 * t87 - t114;
t52 = -t71 * mrSges(5,1) + t72 * mrSges(5,2);
t11 = m(5) * t127 - t85 * mrSges(5,2) + t49 * mrSges(5,3) - t106 * t16 + t147 * t14 + t71 * t52 - t96 * t63;
t12 = m(5) * t150 + t85 * mrSges(5,1) - t50 * mrSges(5,3) - t72 * t52 + t96 * t62 - t149;
t73 = -t81 * mrSges(4,1) + t82 * mrSges(4,2);
t7 = m(4) * t122 + t85 * mrSges(4,1) - t66 * mrSges(4,3) + t102 * t11 + t104 * t12 - t82 * t73 + t96 * t74;
t8 = m(4) * t139 - t85 * mrSges(4,2) + t65 * mrSges(4,3) - t102 * t12 + t104 * t11 + t81 * t73 - t96 * t76;
t86 = t99 * mrSges(3,1) - mrSges(3,3) * t126;
t4 = m(3) * (-g(3) * t135 + t138) - t93 * mrSges(3,3) - t98 * mrSges(3,2) + t90 * t125 - t99 * t86 + t110 * t8 - t107 * t7;
t6 = m(3) * (-t103 * t88 - t145) + t92 * mrSges(3,2) + t93 * mrSges(3,1) + t107 * t8 + t110 * t7 + (t108 * t86 - t111 * t87) * t133;
t131 = t10 * t134 + t105 * t6 + t4 * t135;
t2 = m(2) * t121 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t108 * t10 + t111 * t4;
t1 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) - t103 * t6 + (t111 * t10 + t108 * t4) * t105;
t3 = [-m(1) * g(1) - t109 * t1 + t112 * t2, t2, t4, t8, t11, t14, -t34 * mrSges(7,2) - t60 * t40 + t129; -m(1) * g(2) + t112 * t1 + t109 * t2, t1, t10, t7, t12, t16, -t35 * mrSges(7,3) - t61 * t45 + t128; (-m(1) - m(2)) * g(3) + t131, -m(2) * g(3) + t131, t6, t114, -t117, t149, -t48 * mrSges(7,1) + t35 * mrSges(7,2) + t61 * t40 - t70 * t42 + t148;];
f_new  = t3;
