% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:14:21
% EndTime: 2019-05-06 13:14:34
% DurationCPUTime: 5.98s
% Computational Cost: add. (101289->207), mult. (235550->277), div. (0->0), fcn. (171718->12), ass. (0->104)
t105 = sin(pkin(10));
t107 = cos(pkin(10));
t110 = sin(qJ(2));
t114 = cos(qJ(2));
t135 = qJD(1) * qJD(2);
t117 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t127 = -g(1) * t115 - g(2) * t111;
t93 = -pkin(1) * t117 + qJDD(1) * pkin(7) + t127;
t139 = t110 * t93;
t141 = pkin(2) * t117;
t96 = qJDD(1) * t110 + t114 * t135;
t62 = qJDD(2) * pkin(2) - t96 * qJ(3) - t139 + (qJ(3) * t135 + t110 * t141 - g(3)) * t114;
t103 = t114 ^ 2;
t132 = -g(3) * t110 + t114 * t93;
t97 = qJDD(1) * t114 - t110 * t135;
t137 = qJD(1) * t110;
t98 = qJD(2) * pkin(2) - qJ(3) * t137;
t63 = qJ(3) * t97 - qJD(2) * t98 - t103 * t141 + t132;
t90 = (t105 * t114 + t107 * t110) * qJD(1);
t144 = -0.2e1 * qJD(3) * t90 - t105 * t63 + t107 * t62;
t136 = qJD(1) * t114;
t100 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t136;
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t116 = qJD(2) ^ 2;
t89 = -t105 * t137 + t107 * t136;
t133 = 0.2e1 * qJD(3) * t89 + t105 * t62 + t107 * t63;
t73 = -pkin(3) * t89 - pkin(8) * t90;
t35 = -pkin(3) * t116 + qJDD(2) * pkin(8) + t73 * t89 + t133;
t131 = t111 * g(1) - t115 * g(2);
t125 = -qJDD(1) * pkin(1) - t131;
t120 = -t97 * pkin(2) + qJDD(3) + t98 * t137 + (-qJ(3) * t103 - pkin(7)) * t117 + t125;
t77 = -t105 * t96 + t107 * t97;
t78 = t105 * t97 + t107 * t96;
t38 = (-qJD(2) * t89 - t78) * pkin(8) + (qJD(2) * t90 - t77) * pkin(3) + t120;
t129 = -t109 * t35 + t113 * t38;
t80 = qJD(2) * t113 - t109 * t90;
t53 = qJD(4) * t80 + qJDD(2) * t109 + t113 * t78;
t76 = qJDD(4) - t77;
t81 = qJD(2) * t109 + t113 * t90;
t88 = qJD(4) - t89;
t23 = (t80 * t88 - t53) * qJ(5) + (t80 * t81 + t76) * pkin(4) + t129;
t140 = t109 * t38 + t113 * t35;
t52 = -qJD(4) * t81 + qJDD(2) * t113 - t109 * t78;
t69 = pkin(4) * t88 - qJ(5) * t81;
t79 = t80 ^ 2;
t25 = -pkin(4) * t79 + qJ(5) * t52 - t69 * t88 + t140;
t61 = t104 * t80 + t106 * t81;
t128 = -0.2e1 * qJD(5) * t61 - t104 * t25 + t106 * t23;
t43 = t104 * t52 + t106 * t53;
t60 = -t104 * t81 + t106 * t80;
t17 = (t60 * t88 - t43) * pkin(9) + (t60 * t61 + t76) * pkin(5) + t128;
t134 = 0.2e1 * qJD(5) * t60 + t104 * t23 + t106 * t25;
t42 = -t104 * t53 + t106 * t52;
t50 = pkin(5) * t88 - pkin(9) * t61;
t57 = t60 ^ 2;
t18 = -pkin(5) * t57 + pkin(9) * t42 - t50 * t88 + t134;
t45 = -t108 * t61 + t112 * t60;
t28 = qJD(6) * t45 + t108 * t42 + t112 * t43;
t46 = t108 * t60 + t112 * t61;
t32 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t84 = qJD(6) + t88;
t39 = -mrSges(7,2) * t84 + mrSges(7,3) * t45;
t74 = qJDD(6) + t76;
t14 = m(7) * (-t108 * t18 + t112 * t17) - t28 * mrSges(7,3) + t74 * mrSges(7,1) - t46 * t32 + t84 * t39;
t27 = -qJD(6) * t46 - t108 * t43 + t112 * t42;
t40 = mrSges(7,1) * t84 - mrSges(7,3) * t46;
t15 = m(7) * (t108 * t17 + t112 * t18) + t27 * mrSges(7,3) - t74 * mrSges(7,2) + t45 * t32 - t84 * t40;
t47 = -mrSges(6,1) * t60 + mrSges(6,2) * t61;
t48 = -mrSges(6,2) * t88 + mrSges(6,3) * t60;
t12 = m(6) * t128 + t76 * mrSges(6,1) - t43 * mrSges(6,3) + t108 * t15 + t112 * t14 - t61 * t47 + t88 * t48;
t49 = mrSges(6,1) * t88 - mrSges(6,3) * t61;
t13 = m(6) * t134 - t76 * mrSges(6,2) + t42 * mrSges(6,3) - t108 * t14 + t112 * t15 + t60 * t47 - t88 * t49;
t64 = -mrSges(5,1) * t80 + mrSges(5,2) * t81;
t68 = -mrSges(5,2) * t88 + mrSges(5,3) * t80;
t10 = m(5) * t129 + t76 * mrSges(5,1) - t53 * mrSges(5,3) + t104 * t13 + t106 * t12 - t81 * t64 + t88 * t68;
t70 = mrSges(5,1) * t88 - mrSges(5,3) * t81;
t11 = m(5) * t140 - t76 * mrSges(5,2) + t52 * mrSges(5,3) - t104 * t12 + t106 * t13 + t80 * t64 - t88 * t70;
t82 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t89;
t83 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t90;
t123 = -m(4) * t120 + t77 * mrSges(4,1) - t78 * mrSges(4,2) - t113 * t10 - t109 * t11 + t89 * t82 - t90 * t83;
t99 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t137;
t143 = -(t100 * t114 - t110 * t99) * qJD(1) + m(3) * (-pkin(7) * t117 + t125) - t97 * mrSges(3,1) + t96 * mrSges(3,2) - t123;
t34 = -qJDD(2) * pkin(3) - pkin(8) * t116 + t90 * t73 - t144;
t119 = -pkin(4) * t52 - qJ(5) * t79 + t81 * t69 + qJDD(5) + t34;
t124 = t27 * mrSges(7,1) + t45 * t39 - m(7) * (-pkin(5) * t42 - pkin(9) * t57 + t50 * t61 + t119) - t28 * mrSges(7,2) - t46 * t40;
t122 = -m(6) * t119 + t42 * mrSges(6,1) - t43 * mrSges(6,2) + t60 * t48 - t61 * t49 + t124;
t118 = m(5) * t34 - t52 * mrSges(5,1) + t53 * mrSges(5,2) - t80 * t68 + t81 * t70 - t122;
t72 = -mrSges(4,1) * t89 + mrSges(4,2) * t90;
t16 = m(4) * t144 + qJDD(2) * mrSges(4,1) - t78 * mrSges(4,3) + qJD(2) * t82 - t90 * t72 - t118;
t7 = m(4) * t133 - qJDD(2) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(2) * t83 - t109 * t10 + t113 * t11 + t89 * t72;
t95 = (-mrSges(3,1) * t114 + mrSges(3,2) * t110) * qJD(1);
t4 = m(3) * (-t114 * g(3) - t139) - t96 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t95 * t137 + qJD(2) * t100 + t105 * t7 + t107 * t16;
t5 = m(3) * t132 - qJDD(2) * mrSges(3,2) + t97 * mrSges(3,3) - qJD(2) * t99 - t105 * t16 + t107 * t7 + t95 * t136;
t142 = t110 * t5 + t114 * t4;
t6 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t143;
t1 = m(2) * t127 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t4 + t114 * t5;
t2 = [-m(1) * g(1) + t1 * t115 - t111 * t6, t1, t5, t7, t11, t13, t15; -m(1) * g(2) + t1 * t111 + t115 * t6, t6, t4, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, -t123, t118, -t122, -t124;];
f_new  = t2;
