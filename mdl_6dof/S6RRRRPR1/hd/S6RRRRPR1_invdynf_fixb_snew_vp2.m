% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 19:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:26:14
% EndTime: 2019-05-07 19:26:30
% DurationCPUTime: 7.05s
% Computational Cost: add. (124234->209), mult. (286949->276), div. (0->0), fcn. (218616->12), ass. (0->107)
t106 = sin(qJ(2));
t111 = cos(qJ(2));
t113 = qJD(1) ^ 2;
t100 = t111 ^ 2;
t107 = sin(qJ(1));
t112 = cos(qJ(1));
t129 = g(1) * t107 - t112 * g(2);
t122 = -qJDD(1) * pkin(1) - t129;
t133 = qJD(1) * t106;
t131 = qJD(1) * qJD(2);
t90 = qJDD(1) * t111 - t106 * t131;
t93 = qJD(2) * pkin(2) - pkin(8) * t133;
t119 = -pkin(2) * t90 + t93 * t133 + (-pkin(8) * t100 - pkin(7)) * t113 + t122;
t105 = sin(qJ(3));
t110 = cos(qJ(3));
t84 = (t105 * t111 + t106 * t110) * qJD(1);
t89 = qJDD(1) * t106 + t111 * t131;
t66 = -qJD(3) * t84 - t105 * t89 + t110 * t90;
t99 = qJD(2) + qJD(3);
t80 = pkin(3) * t99 - pkin(9) * t84;
t83 = (-t105 * t106 + t110 * t111) * qJD(1);
t82 = t83 ^ 2;
t117 = -pkin(3) * t66 - pkin(9) * t82 + t84 * t80 + t119;
t103 = sin(qJ(6));
t108 = cos(qJ(6));
t104 = sin(qJ(4));
t109 = cos(qJ(4));
t67 = qJD(3) * t83 + t105 * t90 + t110 * t89;
t76 = t104 * t83 + t109 * t84;
t46 = -qJD(4) * t76 - t104 * t67 + t109 * t66;
t96 = qJD(4) + t99;
t70 = pkin(4) * t96 - qJ(5) * t76;
t75 = -t104 * t84 + t109 * t83;
t74 = t75 ^ 2;
t115 = -pkin(4) * t46 - qJ(5) * t74 + t76 * t70 + qJDD(5) + t117;
t101 = sin(pkin(11));
t102 = cos(pkin(11));
t139 = 2 * qJD(5);
t125 = -g(1) * t112 - g(2) * t107;
t86 = -pkin(1) * t113 + qJDD(1) * pkin(7) + t125;
t134 = t106 * t86;
t137 = pkin(2) * t113;
t62 = qJDD(2) * pkin(2) - t89 * pkin(8) - t134 + (pkin(8) * t131 + t106 * t137 - g(3)) * t111;
t128 = -g(3) * t106 + t111 * t86;
t63 = pkin(8) * t90 - qJD(2) * t93 - t100 * t137 + t128;
t126 = -t105 * t63 + t110 * t62;
t98 = qJDD(2) + qJDD(3);
t34 = (t83 * t99 - t67) * pkin(9) + (t83 * t84 + t98) * pkin(3) + t126;
t135 = t105 * t62 + t110 * t63;
t36 = -pkin(3) * t82 + pkin(9) * t66 - t80 * t99 + t135;
t127 = -t104 * t36 + t109 * t34;
t47 = qJD(4) * t75 + t104 * t66 + t109 * t67;
t95 = qJDD(4) + t98;
t21 = (t75 * t96 - t47) * qJ(5) + (t75 * t76 + t95) * pkin(4) + t127;
t136 = t104 * t34 + t109 * t36;
t23 = -pkin(4) * t74 + qJ(5) * t46 - t70 * t96 + t136;
t56 = -t101 * t76 + t102 * t75;
t130 = t101 * t21 + t102 * t23 + t56 * t139;
t57 = t101 * t75 + t102 * t76;
t42 = -pkin(5) * t56 - pkin(10) * t57;
t94 = t96 ^ 2;
t18 = -pkin(5) * t94 + pkin(10) * t95 + t42 * t56 + t130;
t30 = -t101 * t47 + t102 * t46;
t31 = t101 * t46 + t102 * t47;
t19 = (-t56 * t96 - t31) * pkin(10) + (t57 * t96 - t30) * pkin(5) + t115;
t48 = -t103 * t57 + t108 * t96;
t27 = qJD(6) * t48 + t103 * t95 + t108 * t31;
t29 = qJDD(6) - t30;
t49 = t103 * t96 + t108 * t57;
t37 = -mrSges(7,1) * t48 + mrSges(7,2) * t49;
t55 = qJD(6) - t56;
t38 = -mrSges(7,2) * t55 + mrSges(7,3) * t48;
t15 = m(7) * (-t103 * t18 + t108 * t19) - t27 * mrSges(7,3) + t29 * mrSges(7,1) - t49 * t37 + t55 * t38;
t26 = -qJD(6) * t49 - t103 * t31 + t108 * t95;
t39 = mrSges(7,1) * t55 - mrSges(7,3) * t49;
t16 = m(7) * (t103 * t19 + t108 * t18) + t26 * mrSges(7,3) - t29 * mrSges(7,2) + t48 * t37 - t55 * t39;
t50 = -mrSges(6,2) * t96 + mrSges(6,3) * t56;
t51 = mrSges(6,1) * t96 - mrSges(6,3) * t57;
t121 = -m(6) * t115 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t103 * t16 - t108 * t15 + t56 * t50 - t57 * t51;
t69 = -mrSges(5,2) * t96 + mrSges(5,3) * t75;
t71 = mrSges(5,1) * t96 - mrSges(5,3) * t76;
t118 = -m(5) * t117 + t46 * mrSges(5,1) - t47 * mrSges(5,2) + t75 * t69 - t76 * t71 + t121;
t78 = -mrSges(4,2) * t99 + mrSges(4,3) * t83;
t79 = mrSges(4,1) * t99 - mrSges(4,3) * t84;
t116 = -m(4) * t119 + t66 * mrSges(4,1) - t67 * mrSges(4,2) + t83 * t78 - t84 * t79 + t118;
t91 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t133;
t132 = qJD(1) * t111;
t92 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t132;
t140 = (t106 * t91 - t111 * t92) * qJD(1) + m(3) * (-pkin(7) * t113 + t122) - t90 * mrSges(3,1) + t89 * mrSges(3,2) - t116;
t77 = -mrSges(4,1) * t83 + mrSges(4,2) * t84;
t41 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t11 = m(6) * t130 - t95 * mrSges(6,2) + t30 * mrSges(6,3) - t103 * t15 + t108 * t16 + t56 * t41 - t96 * t51;
t124 = t101 * t23 - t102 * t21;
t120 = m(7) * (-t95 * pkin(5) - t94 * pkin(10) + (t139 + t42) * t57 + t124) - t26 * mrSges(7,1) + t27 * mrSges(7,2) - t48 * t38 + t49 * t39;
t12 = m(6) * (-0.2e1 * qJD(5) * t57 - t124) - t31 * mrSges(6,3) + t95 * mrSges(6,1) - t57 * t41 + t96 * t50 - t120;
t58 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t8 = m(5) * t127 + t95 * mrSges(5,1) - t47 * mrSges(5,3) + t101 * t11 + t102 * t12 - t76 * t58 + t96 * t69;
t9 = m(5) * t136 - t95 * mrSges(5,2) + t46 * mrSges(5,3) - t101 * t12 + t102 * t11 + t75 * t58 - t96 * t71;
t6 = m(4) * t126 + t98 * mrSges(4,1) - t67 * mrSges(4,3) + t104 * t9 + t109 * t8 - t84 * t77 + t99 * t78;
t7 = m(4) * t135 - t98 * mrSges(4,2) + t66 * mrSges(4,3) - t104 * t8 + t109 * t9 + t83 * t77 - t99 * t79;
t88 = (-mrSges(3,1) * t111 + mrSges(3,2) * t106) * qJD(1);
t4 = m(3) * (-t111 * g(3) - t134) - t89 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t88 * t133 + qJD(2) * t92 + t105 * t7 + t110 * t6;
t5 = m(3) * t128 - qJDD(2) * mrSges(3,2) + t90 * mrSges(3,3) - qJD(2) * t91 - t105 * t6 + t110 * t7 + t88 * t132;
t138 = t106 * t5 + t111 * t4;
t10 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t113 * mrSges(2,2) - t140;
t1 = m(2) * t125 - t113 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t4 + t111 * t5;
t2 = [-m(1) * g(1) + t1 * t112 - t10 * t107, t1, t5, t7, t9, t11, t16; -m(1) * g(2) + t1 * t107 + t10 * t112, t10, t4, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t140, -t116, -t118, -t121, t120;];
f_new  = t2;
