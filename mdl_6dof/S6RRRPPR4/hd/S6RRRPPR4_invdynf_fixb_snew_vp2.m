% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-05-07 04:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:40:40
% EndTime: 2019-05-07 04:40:52
% DurationCPUTime: 3.29s
% Computational Cost: add. (42329->203), mult. (87856->256), div. (0->0), fcn. (61037->10), ass. (0->101)
t111 = qJD(2) ^ 2;
t105 = sin(qJ(2));
t135 = qJD(1) * t105;
t109 = cos(qJ(2));
t112 = qJD(1) ^ 2;
t106 = sin(qJ(1));
t110 = cos(qJ(1));
t124 = -t110 * g(1) - t106 * g(2);
t81 = -t112 * pkin(1) + qJDD(1) * pkin(7) + t124;
t137 = -t109 * g(3) - t105 * t81;
t90 = (-pkin(2) * t109 - pkin(8) * t105) * qJD(1);
t118 = qJDD(2) * pkin(2) + t111 * pkin(8) - t90 * t135 + t137;
t104 = sin(qJ(3));
t108 = cos(qJ(3));
t88 = t104 * qJD(2) + t108 * t135;
t133 = qJD(1) * qJD(2);
t128 = t109 * t133;
t91 = t105 * qJDD(1) + t128;
t67 = -t88 * qJD(3) + t108 * qJDD(2) - t104 * t91;
t134 = t109 * qJD(1);
t98 = qJD(3) - t134;
t74 = t98 * pkin(3) - t88 * qJ(4);
t87 = t108 * qJD(2) - t104 * t135;
t85 = t87 ^ 2;
t115 = t67 * pkin(3) + t85 * qJ(4) - t88 * t74 - qJDD(4) + t118;
t102 = sin(pkin(10));
t136 = cos(pkin(10));
t70 = t102 * t88 - t136 * t87;
t141 = t70 * t98;
t68 = t87 * qJD(3) + t104 * qJDD(2) + t108 * t91;
t43 = t102 * t67 + t136 * t68;
t146 = (-t43 + t141) * qJ(5) - t115;
t129 = t105 * t133;
t127 = t106 * g(1) - t110 * g(2);
t80 = -qJDD(1) * pkin(1) - t112 * pkin(7) - t127;
t92 = t109 * qJDD(1) - t129;
t53 = (-t91 - t128) * pkin(8) + (-t92 + t129) * pkin(2) + t80;
t131 = -t105 * g(3) + t109 * t81;
t57 = -t111 * pkin(2) + qJDD(2) * pkin(8) + t134 * t90 + t131;
t126 = -t104 * t57 + t108 * t53;
t103 = sin(qJ(6));
t107 = cos(qJ(6));
t144 = -2 * qJD(4);
t86 = qJDD(3) - t92;
t27 = (t87 * t98 - t68) * qJ(4) + (t87 * t88 + t86) * pkin(3) + t126;
t138 = t104 * t53 + t108 * t57;
t30 = -t85 * pkin(3) + t67 * qJ(4) - t98 * t74 + t138;
t132 = t102 * t27 + t136 * t30 + t144 * t70;
t143 = 2 * qJD(5);
t71 = t102 * t87 + t136 * t88;
t48 = t70 * pkin(4) - t71 * qJ(5);
t97 = t98 ^ 2;
t119 = -t97 * pkin(4) + t86 * qJ(5) + t98 * t143 - t70 * t48 + t132;
t123 = -t102 * t30 + t136 * t27;
t19 = -t86 * pkin(4) - t97 * qJ(5) + qJDD(5) - t123 + ((2 * qJD(4)) + t48) * t71;
t14 = (-t43 - t141) * pkin(9) + (t70 * t71 - t86) * pkin(5) + t19;
t42 = t102 * t68 - t136 * t67;
t62 = -t98 * pkin(5) - t71 * pkin(9);
t69 = t70 ^ 2;
t15 = -t69 * pkin(5) + t42 * pkin(9) + t98 * t62 + t119;
t46 = -t103 * t71 + t107 * t70;
t25 = t46 * qJD(6) + t103 * t42 + t107 * t43;
t47 = t103 * t70 + t107 * t71;
t33 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t96 = qJD(6) - t98;
t38 = -t96 * mrSges(7,2) + t46 * mrSges(7,3);
t84 = qJDD(6) - t86;
t12 = m(7) * (-t103 * t15 + t107 * t14) - t25 * mrSges(7,3) + t84 * mrSges(7,1) - t47 * t33 + t96 * t38;
t24 = -t47 * qJD(6) - t103 * t43 + t107 * t42;
t39 = t96 * mrSges(7,1) - t47 * mrSges(7,3);
t13 = m(7) * (t103 * t14 + t107 * t15) + t24 * mrSges(7,3) - t84 * mrSges(7,2) + t46 * t33 - t96 * t39;
t60 = -t98 * mrSges(6,1) + t71 * mrSges(6,2);
t121 = m(6) * t119 + t86 * mrSges(6,3) - t103 * t12 + t107 * t13 + t98 * t60;
t49 = t70 * mrSges(6,1) - t71 * mrSges(6,3);
t139 = -t70 * mrSges(5,1) - t71 * mrSges(5,2) - t49;
t140 = -mrSges(5,3) - mrSges(6,2);
t59 = t98 * mrSges(5,1) - t71 * mrSges(5,3);
t7 = m(5) * t132 - t86 * mrSges(5,2) + t139 * t70 + t140 * t42 - t98 * t59 + t121;
t72 = -t87 * mrSges(4,1) + t88 * mrSges(4,2);
t73 = -t98 * mrSges(4,2) + t87 * mrSges(4,3);
t120 = -m(6) * t19 - t103 * t13 - t107 * t12;
t58 = -t98 * mrSges(5,2) - t70 * mrSges(5,3);
t61 = -t70 * mrSges(6,2) + t98 * mrSges(6,3);
t8 = m(5) * t123 + (t58 + t61) * t98 + (mrSges(5,1) + mrSges(6,1)) * t86 + (m(5) * t144 + t139) * t71 + t140 * t43 + t120;
t5 = m(4) * t126 + t86 * mrSges(4,1) - t68 * mrSges(4,3) + t102 * t7 + t136 * t8 - t88 * t72 + t98 * t73;
t75 = t98 * mrSges(4,1) - t88 * mrSges(4,3);
t6 = m(4) * t138 - t86 * mrSges(4,2) + t67 * mrSges(4,3) - t102 * t8 + t136 * t7 + t87 * t72 - t98 * t75;
t93 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t135;
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t134;
t145 = m(3) * t80 - t92 * mrSges(3,1) + t91 * mrSges(3,2) + t104 * t6 + t108 * t5 + (t105 * t93 - t109 * t94) * qJD(1);
t125 = m(7) * (-t69 * pkin(9) + (-pkin(4) - pkin(5)) * t42 + (-pkin(4) * t98 + t143 + t62) * t71 - t146) + t25 * mrSges(7,2) - t24 * mrSges(7,1) + t47 * t39 - t46 * t38;
t116 = t43 * mrSges(6,3) + t71 * t60 + t125 - m(6) * (-0.2e1 * qJD(5) * t71 + (t71 * t98 + t42) * pkin(4) + t146) - t42 * mrSges(6,1) - t70 * t61;
t114 = -m(5) * t115 + t42 * mrSges(5,1) + t43 * mrSges(5,2) + t70 * t58 + t71 * t59 - t116;
t113 = -m(4) * t118 - t67 * mrSges(4,1) + t68 * mrSges(4,2) - t87 * t73 + t88 * t75 + t114;
t89 = (-mrSges(3,1) * t109 + mrSges(3,2) * t105) * qJD(1);
t10 = m(3) * t137 + qJDD(2) * mrSges(3,1) - t91 * mrSges(3,3) + qJD(2) * t94 - t135 * t89 - t113;
t4 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t92 * mrSges(3,3) - qJD(2) * t93 - t104 * t5 + t108 * t6 + t134 * t89;
t142 = t109 * t10 + t105 * t4;
t2 = m(2) * t127 + qJDD(1) * mrSges(2,1) - t112 * mrSges(2,2) - t145;
t1 = m(2) * t124 - t112 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t105 * t10 + t109 * t4;
t3 = [-m(1) * g(1) + t110 * t1 - t106 * t2, t1, t4, t6, t7, -t42 * mrSges(6,2) - t70 * t49 + t121, t13; -m(1) * g(2) + t106 * t1 + t110 * t2, t2, t10, t5, t8, -t116, t12; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t145, t113, t114, -t86 * mrSges(6,1) + t43 * mrSges(6,2) + t71 * t49 - t98 * t61 - t120, t125;];
f_new  = t3;
