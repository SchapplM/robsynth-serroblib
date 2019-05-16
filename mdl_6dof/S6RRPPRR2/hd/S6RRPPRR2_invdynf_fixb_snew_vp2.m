% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:43:47
% EndTime: 2019-05-06 09:44:05
% DurationCPUTime: 6.03s
% Computational Cost: add. (97013->206), mult. (230479->276), div. (0->0), fcn. (167488->12), ass. (0->105)
t145 = -2 * qJD(3);
t105 = sin(pkin(10));
t107 = cos(pkin(10));
t110 = sin(qJ(2));
t114 = cos(qJ(2));
t135 = qJD(1) * qJD(2);
t117 = qJD(1) ^ 2;
t111 = sin(qJ(1));
t115 = cos(qJ(1));
t127 = -t115 * g(1) - t111 * g(2);
t93 = -t117 * pkin(1) + qJDD(1) * pkin(7) + t127;
t139 = t110 * t93;
t141 = pkin(2) * t117;
t96 = t110 * qJDD(1) + t114 * t135;
t58 = qJDD(2) * pkin(2) - t96 * qJ(3) - t139 + (qJ(3) * t135 + t110 * t141 - g(3)) * t114;
t103 = t114 ^ 2;
t132 = -t110 * g(3) + t114 * t93;
t97 = t114 * qJDD(1) - t110 * t135;
t137 = qJD(1) * t110;
t98 = qJD(2) * pkin(2) - qJ(3) * t137;
t60 = t97 * qJ(3) - qJD(2) * t98 - t103 * t141 + t132;
t90 = (t105 * t114 + t107 * t110) * qJD(1);
t144 = -t105 * t60 + t107 * t58 + t90 * t145;
t136 = qJD(1) * t114;
t100 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t136;
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t104 = sin(pkin(11));
t106 = cos(pkin(11));
t116 = qJD(2) ^ 2;
t89 = t105 * t137 - t107 * t136;
t133 = t105 * t58 + t107 * t60 + t89 * t145;
t70 = t89 * pkin(3) - t90 * qJ(4);
t35 = -t116 * pkin(3) + qJDD(2) * qJ(4) - t89 * t70 + t133;
t131 = t111 * g(1) - t115 * g(2);
t125 = -qJDD(1) * pkin(1) - t131;
t120 = -t97 * pkin(2) + qJDD(3) + t98 * t137 + (-qJ(3) * t103 - pkin(7)) * t117 + t125;
t75 = t105 * t96 - t107 * t97;
t76 = t105 * t97 + t107 * t96;
t38 = (qJD(2) * t89 - t76) * qJ(4) + (qJD(2) * t90 + t75) * pkin(3) + t120;
t81 = t104 * qJD(2) + t106 * t90;
t128 = -0.2e1 * qJD(4) * t81 - t104 * t35 + t106 * t38;
t68 = t104 * qJDD(2) + t106 * t76;
t80 = t106 * qJD(2) - t104 * t90;
t23 = (t80 * t89 - t68) * pkin(8) + (t80 * t81 + t75) * pkin(4) + t128;
t134 = 0.2e1 * qJD(4) * t80 + t104 * t38 + t106 * t35;
t66 = t89 * pkin(4) - t81 * pkin(8);
t67 = t106 * qJDD(2) - t104 * t76;
t79 = t80 ^ 2;
t25 = -t79 * pkin(4) + t67 * pkin(8) - t89 * t66 + t134;
t129 = -t109 * t25 + t113 * t23;
t54 = -t109 * t81 + t113 * t80;
t41 = t54 * qJD(5) + t109 * t67 + t113 * t68;
t55 = t109 * t80 + t113 * t81;
t74 = qJDD(5) + t75;
t88 = qJD(5) + t89;
t17 = (t54 * t88 - t41) * pkin(9) + (t54 * t55 + t74) * pkin(5) + t129;
t140 = t109 * t23 + t113 * t25;
t40 = -t55 * qJD(5) - t109 * t68 + t113 * t67;
t50 = t88 * pkin(5) - t55 * pkin(9);
t53 = t54 ^ 2;
t18 = -t53 * pkin(5) + t40 * pkin(9) - t88 * t50 + t140;
t45 = -t108 * t55 + t112 * t54;
t28 = t45 * qJD(6) + t108 * t40 + t112 * t41;
t46 = t108 * t54 + t112 * t55;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t84 = qJD(6) + t88;
t42 = -t84 * mrSges(7,2) + t45 * mrSges(7,3);
t72 = qJDD(6) + t74;
t14 = m(7) * (-t108 * t18 + t112 * t17) - t28 * mrSges(7,3) + t72 * mrSges(7,1) - t46 * t32 + t84 * t42;
t27 = -t46 * qJD(6) - t108 * t41 + t112 * t40;
t43 = t84 * mrSges(7,1) - t46 * mrSges(7,3);
t15 = m(7) * (t108 * t17 + t112 * t18) + t27 * mrSges(7,3) - t72 * mrSges(7,2) + t45 * t32 - t84 * t43;
t47 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t48 = -t88 * mrSges(6,2) + t54 * mrSges(6,3);
t12 = m(6) * t129 + t74 * mrSges(6,1) - t41 * mrSges(6,3) + t108 * t15 + t112 * t14 - t55 * t47 + t88 * t48;
t49 = t88 * mrSges(6,1) - t55 * mrSges(6,3);
t13 = m(6) * t140 - t74 * mrSges(6,2) + t40 * mrSges(6,3) - t108 * t14 + t112 * t15 + t54 * t47 - t88 * t49;
t59 = -t80 * mrSges(5,1) + t81 * mrSges(5,2);
t64 = -t89 * mrSges(5,2) + t80 * mrSges(5,3);
t10 = m(5) * t128 + t75 * mrSges(5,1) - t68 * mrSges(5,3) + t109 * t13 + t113 * t12 - t81 * t59 + t89 * t64;
t65 = t89 * mrSges(5,1) - t81 * mrSges(5,3);
t11 = m(5) * t134 - t75 * mrSges(5,2) + t67 * mrSges(5,3) - t109 * t12 + t113 * t13 + t80 * t59 - t89 * t65;
t82 = -qJD(2) * mrSges(4,2) - t89 * mrSges(4,3);
t83 = qJD(2) * mrSges(4,1) - t90 * mrSges(4,3);
t123 = m(4) * t120 + t75 * mrSges(4,1) + t76 * mrSges(4,2) + t106 * t10 + t104 * t11 + t89 * t82 + t90 * t83;
t99 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t137;
t143 = -(t114 * t100 - t110 * t99) * qJD(1) + m(3) * (-t117 * pkin(7) + t125) - t97 * mrSges(3,1) + t96 * mrSges(3,2) + t123;
t34 = -qJDD(2) * pkin(3) - t116 * qJ(4) + t90 * t70 + qJDD(4) - t144;
t119 = -t67 * pkin(4) - t79 * pkin(8) + t81 * t66 + t34;
t124 = t27 * mrSges(7,1) + t45 * t42 - m(7) * (-t40 * pkin(5) - t53 * pkin(9) + t55 * t50 + t119) - t28 * mrSges(7,2) - t46 * t43;
t122 = -m(6) * t119 + t40 * mrSges(6,1) - t41 * mrSges(6,2) + t54 * t48 - t55 * t49 + t124;
t118 = m(5) * t34 - t67 * mrSges(5,1) + t68 * mrSges(5,2) - t80 * t64 + t81 * t65 - t122;
t71 = t89 * mrSges(4,1) + t90 * mrSges(4,2);
t16 = m(4) * t144 + qJDD(2) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(2) * t82 - t90 * t71 - t118;
t7 = m(4) * t133 - qJDD(2) * mrSges(4,2) - t75 * mrSges(4,3) - qJD(2) * t83 - t104 * t10 + t106 * t11 - t89 * t71;
t95 = (-mrSges(3,1) * t114 + mrSges(3,2) * t110) * qJD(1);
t4 = m(3) * (-t114 * g(3) - t139) - t96 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t95 * t137 + qJD(2) * t100 + t105 * t7 + t107 * t16;
t5 = m(3) * t132 - qJDD(2) * mrSges(3,2) + t97 * mrSges(3,3) - qJD(2) * t99 - t105 * t16 + t107 * t7 + t95 * t136;
t142 = t110 * t5 + t114 * t4;
t6 = m(2) * t131 + qJDD(1) * mrSges(2,1) - t117 * mrSges(2,2) - t143;
t1 = m(2) * t127 - t117 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t110 * t4 + t114 * t5;
t2 = [-m(1) * g(1) + t115 * t1 - t111 * t6, t1, t5, t7, t11, t13, t15; -m(1) * g(2) + t111 * t1 + t115 * t6, t6, t4, t16, t10, t12, t14; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t143, t123, t118, -t122, -t124;];
f_new  = t2;
