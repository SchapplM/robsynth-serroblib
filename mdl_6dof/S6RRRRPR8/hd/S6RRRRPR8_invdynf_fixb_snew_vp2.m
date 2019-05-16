% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 21:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:46:22
% EndTime: 2019-05-07 21:46:33
% DurationCPUTime: 3.24s
% Computational Cost: add. (45254->203), mult. (90428->254), div. (0->0), fcn. (63581->10), ass. (0->102)
t112 = qJD(2) ^ 2;
t106 = sin(qJ(2));
t134 = qJD(1) * t106;
t110 = cos(qJ(2));
t113 = qJD(1) ^ 2;
t107 = sin(qJ(1));
t111 = cos(qJ(1));
t124 = -t111 * g(1) - t107 * g(2);
t82 = -t113 * pkin(1) + qJDD(1) * pkin(7) + t124;
t135 = -t110 * g(3) - t106 * t82;
t89 = (-pkin(2) * t110 - pkin(8) * t106) * qJD(1);
t119 = qJDD(2) * pkin(2) + t112 * pkin(8) - t89 * t134 + t135;
t105 = sin(qJ(3));
t109 = cos(qJ(3));
t87 = t105 * qJD(2) + t109 * t134;
t132 = qJD(1) * qJD(2);
t129 = t110 * t132;
t90 = t106 * qJDD(1) + t129;
t65 = -t87 * qJD(3) + t109 * qJDD(2) - t105 * t90;
t133 = t110 * qJD(1);
t98 = qJD(3) - t133;
t73 = t98 * pkin(3) - t87 * pkin(9);
t86 = t109 * qJD(2) - t105 * t134;
t84 = t86 ^ 2;
t116 = t65 * pkin(3) + t84 * pkin(9) - t87 * t73 + t119;
t104 = sin(qJ(4));
t141 = cos(qJ(4));
t68 = t104 * t87 - t141 * t86;
t97 = qJD(4) + t98;
t140 = t68 * t97;
t66 = t86 * qJD(3) + t105 * qJDD(2) + t109 * t90;
t39 = -t68 * qJD(4) + t104 * t65 + t141 * t66;
t145 = (-t39 + t140) * qJ(5) - t116;
t128 = t107 * g(1) - t111 * g(2);
t81 = -qJDD(1) * pkin(1) - t113 * pkin(7) - t128;
t99 = t106 * t132;
t91 = t110 * qJDD(1) - t99;
t53 = (-t90 - t129) * pkin(8) + (-t91 + t99) * pkin(2) + t81;
t131 = -t106 * g(3) + t110 * t82;
t57 = -t112 * pkin(2) + qJDD(2) * pkin(8) + t89 * t133 + t131;
t127 = -t105 * t57 + t109 * t53;
t103 = sin(qJ(6));
t108 = cos(qJ(6));
t85 = qJDD(3) - t91;
t27 = (t86 * t98 - t66) * pkin(9) + (t86 * t87 + t85) * pkin(3) + t127;
t136 = t105 * t53 + t109 * t57;
t30 = -t84 * pkin(3) + t65 * pkin(9) - t98 * t73 + t136;
t125 = -t104 * t30 + t141 * t27;
t69 = t104 * t86 + t141 * t87;
t48 = t68 * pkin(4) - t69 * qJ(5);
t83 = qJDD(4) + t85;
t96 = t97 ^ 2;
t19 = -t83 * pkin(4) - t96 * qJ(5) + t69 * t48 + qJDD(5) - t125;
t14 = (-t39 - t140) * pkin(10) + (t68 * t69 - t83) * pkin(5) + t19;
t138 = t104 * t27 + t141 * t30;
t143 = 2 * qJD(5);
t121 = -t96 * pkin(4) + t83 * qJ(5) + t97 * t143 - t68 * t48 + t138;
t38 = t69 * qJD(4) + t104 * t66 - t141 * t65;
t62 = -t97 * pkin(5) - t69 * pkin(10);
t67 = t68 ^ 2;
t15 = -t67 * pkin(5) + t38 * pkin(10) + t97 * t62 + t121;
t46 = -t103 * t69 + t108 * t68;
t25 = t46 * qJD(6) + t103 * t38 + t108 * t39;
t47 = t103 * t68 + t108 * t69;
t33 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t93 = qJD(6) - t97;
t42 = -t93 * mrSges(7,2) + t46 * mrSges(7,3);
t79 = qJDD(6) - t83;
t12 = m(7) * (-t103 * t15 + t108 * t14) - t25 * mrSges(7,3) + t79 * mrSges(7,1) - t47 * t33 + t93 * t42;
t24 = -t47 * qJD(6) - t103 * t39 + t108 * t38;
t43 = t93 * mrSges(7,1) - t47 * mrSges(7,3);
t13 = m(7) * (t103 * t14 + t108 * t15) + t24 * mrSges(7,3) - t79 * mrSges(7,2) + t46 * t33 - t93 * t43;
t61 = -t97 * mrSges(6,1) + t69 * mrSges(6,2);
t122 = m(6) * t121 + t83 * mrSges(6,3) - t103 * t12 + t108 * t13 + t97 * t61;
t49 = t68 * mrSges(6,1) - t69 * mrSges(6,3);
t137 = -t68 * mrSges(5,1) - t69 * mrSges(5,2) - t49;
t139 = -mrSges(5,3) - mrSges(6,2);
t60 = t97 * mrSges(5,1) - t69 * mrSges(5,3);
t7 = m(5) * t138 - t83 * mrSges(5,2) + t137 * t68 + t139 * t38 - t97 * t60 + t122;
t70 = -t86 * mrSges(4,1) + t87 * mrSges(4,2);
t71 = -t98 * mrSges(4,2) + t86 * mrSges(4,3);
t120 = -m(6) * t19 - t103 * t13 - t108 * t12;
t58 = -t68 * mrSges(6,2) + t97 * mrSges(6,3);
t59 = -t97 * mrSges(5,2) - t68 * mrSges(5,3);
t8 = m(5) * t125 + (t59 + t58) * t97 + (mrSges(5,1) + mrSges(6,1)) * t83 + t137 * t69 + t139 * t39 + t120;
t5 = m(4) * t127 + t85 * mrSges(4,1) - t66 * mrSges(4,3) + t104 * t7 + t141 * t8 - t87 * t70 + t98 * t71;
t72 = t98 * mrSges(4,1) - t87 * mrSges(4,3);
t6 = m(4) * t136 - t85 * mrSges(4,2) + t65 * mrSges(4,3) - t104 * t8 + t141 * t7 + t86 * t70 - t98 * t72;
t94 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t134;
t95 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133;
t144 = m(3) * t81 - t91 * mrSges(3,1) + t90 * mrSges(3,2) + t105 * t6 + t109 * t5 + (t106 * t94 - t110 * t95) * qJD(1);
t126 = m(7) * (-t67 * pkin(10) + (-pkin(4) - pkin(5)) * t38 + (-pkin(4) * t97 + t143 + t62) * t69 - t145) + t25 * mrSges(7,2) - t24 * mrSges(7,1) + t47 * t43 - t46 * t42;
t117 = t39 * mrSges(6,3) + t69 * t61 + t126 - m(6) * (-0.2e1 * qJD(5) * t69 + (t69 * t97 + t38) * pkin(4) + t145) - t38 * mrSges(6,1) - t68 * t58;
t115 = -m(5) * t116 + t38 * mrSges(5,1) + t39 * mrSges(5,2) + t68 * t59 + t69 * t60 - t117;
t114 = -m(4) * t119 - t65 * mrSges(4,1) + t66 * mrSges(4,2) - t86 * t71 + t87 * t72 + t115;
t88 = (-mrSges(3,1) * t110 + mrSges(3,2) * t106) * qJD(1);
t10 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t90 * mrSges(3,3) + qJD(2) * t95 - t88 * t134 - t114;
t4 = m(3) * t131 - qJDD(2) * mrSges(3,2) + t91 * mrSges(3,3) - qJD(2) * t94 - t105 * t5 + t109 * t6 + t88 * t133;
t142 = t110 * t10 + t106 * t4;
t2 = m(2) * t128 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) - t144;
t1 = m(2) * t124 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t10 + t110 * t4;
t3 = [-m(1) * g(1) + t111 * t1 - t107 * t2, t1, t4, t6, t7, -t38 * mrSges(6,2) - t68 * t49 + t122, t13; -m(1) * g(2) + t107 * t1 + t111 * t2, t2, t10, t5, t8, -t117, t12; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t144, t114, t115, -t83 * mrSges(6,1) + t39 * mrSges(6,2) + t69 * t49 - t97 * t58 - t120, t126;];
f_new  = t3;
