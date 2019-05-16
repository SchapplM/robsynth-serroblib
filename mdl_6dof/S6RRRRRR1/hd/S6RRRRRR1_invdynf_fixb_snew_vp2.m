% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 08:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:57:33
% EndTime: 2019-05-08 07:57:53
% DurationCPUTime: 7.43s
% Computational Cost: add. (129849->209), mult. (296961->274), div. (0->0), fcn. (228472->12), ass. (0->108)
t106 = sin(qJ(2));
t112 = cos(qJ(2));
t114 = qJD(1) ^ 2;
t101 = t112 ^ 2;
t107 = sin(qJ(1));
t113 = cos(qJ(1));
t130 = t107 * g(1) - t113 * g(2);
t123 = -qJDD(1) * pkin(1) - t130;
t133 = qJD(1) * t106;
t131 = qJD(1) * qJD(2);
t89 = t112 * qJDD(1) - t106 * t131;
t92 = qJD(2) * pkin(2) - pkin(8) * t133;
t120 = -t89 * pkin(2) + t92 * t133 + (-pkin(8) * t101 - pkin(7)) * t114 + t123;
t105 = sin(qJ(3));
t111 = cos(qJ(3));
t83 = (t105 * t112 + t106 * t111) * qJD(1);
t88 = t106 * qJDD(1) + t112 * t131;
t65 = -t83 * qJD(3) - t105 * t88 + t111 * t89;
t100 = qJD(2) + qJD(3);
t79 = t100 * pkin(3) - t83 * pkin(9);
t82 = (-t105 * t106 + t111 * t112) * qJD(1);
t81 = t82 ^ 2;
t118 = -t65 * pkin(3) - t81 * pkin(9) + t83 * t79 + t120;
t102 = sin(qJ(6));
t108 = cos(qJ(6));
t104 = sin(qJ(4));
t110 = cos(qJ(4));
t66 = t82 * qJD(3) + t105 * t89 + t111 * t88;
t75 = t104 * t82 + t110 * t83;
t46 = -t75 * qJD(4) - t104 * t66 + t110 * t65;
t97 = qJD(4) + t100;
t70 = t97 * pkin(4) - t75 * pkin(10);
t74 = -t104 * t83 + t110 * t82;
t73 = t74 ^ 2;
t116 = -t46 * pkin(4) - t73 * pkin(10) + t75 * t70 + t118;
t103 = sin(qJ(5));
t109 = cos(qJ(5));
t126 = -t113 * g(1) - t107 * g(2);
t85 = -t114 * pkin(1) + qJDD(1) * pkin(7) + t126;
t134 = t106 * t85;
t138 = pkin(2) * t114;
t61 = qJDD(2) * pkin(2) - t88 * pkin(8) - t134 + (pkin(8) * t131 + t106 * t138 - g(3)) * t112;
t129 = -t106 * g(3) + t112 * t85;
t62 = t89 * pkin(8) - qJD(2) * t92 - t101 * t138 + t129;
t127 = -t105 * t62 + t111 * t61;
t99 = qJDD(2) + qJDD(3);
t34 = (t100 * t82 - t66) * pkin(9) + (t82 * t83 + t99) * pkin(3) + t127;
t135 = t105 * t61 + t111 * t62;
t36 = -t81 * pkin(3) + t65 * pkin(9) - t100 * t79 + t135;
t128 = -t104 * t36 + t110 * t34;
t47 = t74 * qJD(4) + t104 * t65 + t110 * t66;
t96 = qJDD(4) + t99;
t21 = (t74 * t97 - t47) * pkin(10) + (t74 * t75 + t96) * pkin(4) + t128;
t136 = t104 * t34 + t110 * t36;
t23 = -t73 * pkin(4) + t46 * pkin(10) - t97 * t70 + t136;
t137 = t103 * t21 + t109 * t23;
t55 = -t103 * t75 + t109 * t74;
t56 = t103 * t74 + t109 * t75;
t42 = -t55 * pkin(5) - t56 * pkin(11);
t93 = qJDD(5) + t96;
t95 = qJD(5) + t97;
t94 = t95 ^ 2;
t18 = -t94 * pkin(5) + t93 * pkin(11) + t55 * t42 + t137;
t30 = -t56 * qJD(5) - t103 * t47 + t109 * t46;
t31 = t55 * qJD(5) + t103 * t46 + t109 * t47;
t19 = (-t55 * t95 - t31) * pkin(11) + (t56 * t95 - t30) * pkin(5) + t116;
t48 = -t102 * t56 + t108 * t95;
t25 = t48 * qJD(6) + t102 * t93 + t108 * t31;
t29 = qJDD(6) - t30;
t49 = t102 * t95 + t108 * t56;
t37 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t54 = qJD(6) - t55;
t38 = -t54 * mrSges(7,2) + t48 * mrSges(7,3);
t15 = m(7) * (-t102 * t18 + t108 * t19) - t25 * mrSges(7,3) + t29 * mrSges(7,1) - t49 * t37 + t54 * t38;
t24 = -t49 * qJD(6) - t102 * t31 + t108 * t93;
t39 = t54 * mrSges(7,1) - t49 * mrSges(7,3);
t16 = m(7) * (t102 * t19 + t108 * t18) + t24 * mrSges(7,3) - t29 * mrSges(7,2) + t48 * t37 - t54 * t39;
t50 = -t95 * mrSges(6,2) + t55 * mrSges(6,3);
t51 = t95 * mrSges(6,1) - t56 * mrSges(6,3);
t122 = -m(6) * t116 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t102 * t16 - t108 * t15 + t55 * t50 - t56 * t51;
t68 = -t97 * mrSges(5,2) + t74 * mrSges(5,3);
t69 = t97 * mrSges(5,1) - t75 * mrSges(5,3);
t119 = -m(5) * t118 + t46 * mrSges(5,1) - t47 * mrSges(5,2) + t74 * t68 - t75 * t69 + t122;
t77 = -t100 * mrSges(4,2) + t82 * mrSges(4,3);
t78 = t100 * mrSges(4,1) - t83 * mrSges(4,3);
t117 = -m(4) * t120 + t65 * mrSges(4,1) - t66 * mrSges(4,2) + t82 * t77 - t83 * t78 + t119;
t90 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t133;
t132 = qJD(1) * t112;
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t132;
t140 = (t106 * t90 - t112 * t91) * qJD(1) + m(3) * (-t114 * pkin(7) + t123) - t89 * mrSges(3,1) + t88 * mrSges(3,2) - t117;
t76 = -t82 * mrSges(4,1) + t83 * mrSges(4,2);
t41 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t11 = m(6) * t137 - t93 * mrSges(6,2) + t30 * mrSges(6,3) - t102 * t15 + t108 * t16 + t55 * t41 - t95 * t51;
t125 = -t103 * t23 + t109 * t21;
t121 = m(7) * (-t93 * pkin(5) - t94 * pkin(11) + t56 * t42 - t125) - t24 * mrSges(7,1) + t25 * mrSges(7,2) - t48 * t38 + t49 * t39;
t12 = m(6) * t125 + t93 * mrSges(6,1) - t31 * mrSges(6,3) - t56 * t41 + t95 * t50 - t121;
t57 = -t74 * mrSges(5,1) + t75 * mrSges(5,2);
t8 = m(5) * t128 + t96 * mrSges(5,1) - t47 * mrSges(5,3) + t103 * t11 + t109 * t12 - t75 * t57 + t97 * t68;
t9 = m(5) * t136 - t96 * mrSges(5,2) + t46 * mrSges(5,3) - t103 * t12 + t109 * t11 + t74 * t57 - t97 * t69;
t6 = m(4) * t127 + t99 * mrSges(4,1) - t66 * mrSges(4,3) + t100 * t77 + t104 * t9 + t110 * t8 - t83 * t76;
t7 = m(4) * t135 - t99 * mrSges(4,2) + t65 * mrSges(4,3) - t100 * t78 - t104 * t8 + t110 * t9 + t82 * t76;
t87 = (-mrSges(3,1) * t112 + mrSges(3,2) * t106) * qJD(1);
t4 = m(3) * (-t112 * g(3) - t134) - t88 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t87 * t133 + qJD(2) * t91 + t105 * t7 + t111 * t6;
t5 = m(3) * t129 - qJDD(2) * mrSges(3,2) + t89 * mrSges(3,3) - qJD(2) * t90 - t105 * t6 + t111 * t7 + t87 * t132;
t139 = t106 * t5 + t112 * t4;
t10 = m(2) * t130 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t140;
t1 = m(2) * t126 - t114 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t112 * t5;
t2 = [-m(1) * g(1) + t113 * t1 - t107 * t10, t1, t5, t7, t9, t11, t16; -m(1) * g(2) + t107 * t1 + t113 * t10, t10, t4, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t139, -m(2) * g(3) + t139, t140, -t117, -t119, -t122, t121;];
f_new  = t2;
