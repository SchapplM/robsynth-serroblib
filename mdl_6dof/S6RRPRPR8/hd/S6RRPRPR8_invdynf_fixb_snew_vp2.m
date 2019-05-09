% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:54:37
% EndTime: 2019-05-06 14:54:46
% DurationCPUTime: 3.18s
% Computational Cost: add. (39827->202), mult. (86002->256), div. (0->0), fcn. (59810->10), ass. (0->100)
t111 = qJD(2) ^ 2;
t106 = sin(qJ(2));
t134 = qJD(1) * t106;
t109 = cos(qJ(2));
t112 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t110 = cos(qJ(1));
t123 = -g(1) * t110 - g(2) * t107;
t81 = -pkin(1) * t112 + qJDD(1) * pkin(7) + t123;
t135 = -t109 * g(3) - t106 * t81;
t89 = (-pkin(2) * t109 - qJ(3) * t106) * qJD(1);
t118 = qJDD(2) * pkin(2) + qJ(3) * t111 - t89 * t134 - qJDD(3) + t135;
t102 = sin(pkin(10));
t103 = cos(pkin(10));
t132 = qJD(1) * qJD(2);
t128 = t109 * t132;
t91 = qJDD(1) * t106 + t128;
t71 = qJDD(2) * t103 - t102 * t91;
t133 = qJD(1) * t109;
t87 = qJD(2) * t102 + t103 * t134;
t73 = -pkin(3) * t133 - pkin(8) * t87;
t86 = qJD(2) * t103 - t102 * t134;
t85 = t86 ^ 2;
t115 = pkin(3) * t71 + pkin(8) * t85 - t87 * t73 + t118;
t105 = sin(qJ(4));
t140 = cos(qJ(4));
t65 = t105 * t87 - t140 * t86;
t98 = qJD(4) - t133;
t139 = t65 * t98;
t72 = qJDD(2) * t102 + t103 * t91;
t43 = -t65 * qJD(4) + t105 * t71 + t140 * t72;
t144 = (-t43 + t139) * qJ(5) - t115;
t127 = g(1) * t107 - t110 * g(2);
t80 = -qJDD(1) * pkin(1) - pkin(7) * t112 - t127;
t99 = t106 * t132;
t92 = qJDD(1) * t109 - t99;
t53 = (-t91 - t128) * qJ(3) + (-t92 + t99) * pkin(2) + t80;
t130 = -g(3) * t106 + t109 * t81;
t57 = -pkin(2) * t111 + qJDD(2) * qJ(3) + t89 * t133 + t130;
t126 = -0.2e1 * qJD(3) * t87 - t102 * t57 + t103 * t53;
t67 = -mrSges(4,1) * t86 + mrSges(4,2) * t87;
t69 = mrSges(4,2) * t133 + mrSges(4,3) * t86;
t104 = sin(qJ(6));
t108 = cos(qJ(6));
t27 = (-t86 * t133 - t72) * pkin(8) + (t86 * t87 - t92) * pkin(3) + t126;
t131 = 0.2e1 * qJD(3) * t86 + t102 * t53 + t103 * t57;
t30 = -pkin(3) * t85 + pkin(8) * t71 + t73 * t133 + t131;
t124 = -t105 * t30 + t140 * t27;
t66 = t105 * t86 + t140 * t87;
t48 = t65 * pkin(4) - qJ(5) * t66;
t88 = qJDD(4) - t92;
t97 = t98 ^ 2;
t19 = -t88 * pkin(4) - t97 * qJ(5) + t66 * t48 + qJDD(5) - t124;
t14 = (-t43 - t139) * pkin(9) + (t65 * t66 - t88) * pkin(5) + t19;
t137 = t105 * t27 + t140 * t30;
t142 = 2 * qJD(5);
t120 = -pkin(4) * t97 + t88 * qJ(5) + t98 * t142 - t65 * t48 + t137;
t42 = qJD(4) * t66 + t105 * t72 - t140 * t71;
t62 = -pkin(5) * t98 - pkin(9) * t66;
t64 = t65 ^ 2;
t15 = -t64 * pkin(5) + t42 * pkin(9) + t62 * t98 + t120;
t46 = -t104 * t66 + t108 * t65;
t25 = t46 * qJD(6) + t104 * t42 + t108 * t43;
t47 = t104 * t65 + t108 * t66;
t33 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t96 = qJD(6) - t98;
t38 = -mrSges(7,2) * t96 + t46 * mrSges(7,3);
t84 = qJDD(6) - t88;
t12 = m(7) * (-t104 * t15 + t108 * t14) - t25 * mrSges(7,3) + t84 * mrSges(7,1) - t47 * t33 + t96 * t38;
t24 = -t47 * qJD(6) - t104 * t43 + t108 * t42;
t39 = mrSges(7,1) * t96 - t47 * mrSges(7,3);
t13 = m(7) * (t104 * t14 + t108 * t15) + t24 * mrSges(7,3) - t84 * mrSges(7,2) + t46 * t33 - t96 * t39;
t60 = -mrSges(6,1) * t98 + mrSges(6,2) * t66;
t121 = m(6) * t120 + t88 * mrSges(6,3) - t104 * t12 + t108 * t13 + t98 * t60;
t49 = t65 * mrSges(6,1) - mrSges(6,3) * t66;
t136 = -t65 * mrSges(5,1) - mrSges(5,2) * t66 - t49;
t138 = -mrSges(5,3) - mrSges(6,2);
t59 = mrSges(5,1) * t98 - mrSges(5,3) * t66;
t7 = m(5) * t137 - t88 * mrSges(5,2) + t136 * t65 + t138 * t42 - t98 * t59 + t121;
t119 = -m(6) * t19 - t104 * t13 - t108 * t12;
t58 = -mrSges(5,2) * t98 - t65 * mrSges(5,3);
t61 = -t65 * mrSges(6,2) + mrSges(6,3) * t98;
t8 = m(5) * t124 + (t58 + t61) * t98 + (mrSges(5,1) + mrSges(6,1)) * t88 + t136 * t66 + t138 * t43 + t119;
t5 = m(4) * t126 - t92 * mrSges(4,1) - t72 * mrSges(4,3) + t105 * t7 - t69 * t133 + t140 * t8 - t87 * t67;
t70 = -mrSges(4,1) * t133 - mrSges(4,3) * t87;
t6 = m(4) * t131 + t92 * mrSges(4,2) + t71 * mrSges(4,3) - t105 * t8 + t70 * t133 + t140 * t7 + t86 * t67;
t93 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t143 = m(3) * t80 - t92 * mrSges(3,1) + t91 * mrSges(3,2) + t102 * t6 + t103 * t5 + (t106 * t93 - t109 * t94) * qJD(1);
t125 = m(7) * (-t64 * pkin(9) + (-pkin(4) - pkin(5)) * t42 + (-pkin(4) * t98 + t142 + t62) * t66 - t144) + t25 * mrSges(7,2) - t24 * mrSges(7,1) + t47 * t39 - t46 * t38;
t116 = t43 * mrSges(6,3) + t66 * t60 + t125 - m(6) * (-0.2e1 * qJD(5) * t66 + (t66 * t98 + t42) * pkin(4) + t144) - t42 * mrSges(6,1) - t65 * t61;
t114 = -m(5) * t115 + t42 * mrSges(5,1) + t43 * mrSges(5,2) + t65 * t58 + t66 * t59 - t116;
t113 = -m(4) * t118 - t71 * mrSges(4,1) + t72 * mrSges(4,2) - t86 * t69 + t87 * t70 + t114;
t90 = (-mrSges(3,1) * t109 + mrSges(3,2) * t106) * qJD(1);
t10 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t91 * mrSges(3,3) + qJD(2) * t94 - t90 * t134 - t113;
t4 = m(3) * t130 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t93 - t102 * t5 + t103 * t6 + t90 * t133;
t141 = t109 * t10 + t106 * t4;
t2 = m(2) * t127 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t143;
t1 = m(2) * t123 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t10 + t109 * t4;
t3 = [-m(1) * g(1) + t1 * t110 - t107 * t2, t1, t4, t6, t7, -t42 * mrSges(6,2) - t65 * t49 + t121, t13; -m(1) * g(2) + t1 * t107 + t110 * t2, t2, t10, t5, t8, -t116, t12; (-m(1) - m(2)) * g(3) + t141, -m(2) * g(3) + t141, t143, t113, t114, -t88 * mrSges(6,1) + t43 * mrSges(6,2) + t66 * t49 - t98 * t61 - t119, t125;];
f_new  = t3;
