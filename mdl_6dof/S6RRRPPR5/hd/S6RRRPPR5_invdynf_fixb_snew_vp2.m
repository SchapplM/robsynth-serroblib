% Calculate vector of cutting forces with Newton-Euler
% S6RRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 05:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:59:47
% EndTime: 2019-05-07 05:00:07
% DurationCPUTime: 9.55s
% Computational Cost: add. (174051->213), mult. (387960->297), div. (0->0), fcn. (313732->14), ass. (0->115)
t153 = -2 * qJD(4);
t111 = sin(pkin(6));
t116 = sin(qJ(2));
t120 = cos(qJ(2));
t139 = qJD(1) * qJD(2);
t101 = (-qJDD(1) * t120 + t116 * t139) * t111;
t110 = sin(pkin(11));
t146 = cos(pkin(11));
t141 = qJD(1) * t120;
t135 = t111 * t141;
t103 = qJD(3) - t135;
t115 = sin(qJ(3));
t119 = cos(qJ(3));
t113 = cos(pkin(6));
t106 = t113 * qJD(1) + qJD(2);
t104 = t106 ^ 2;
t105 = t113 * qJDD(1) + qJDD(2);
t122 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t121 = cos(qJ(1));
t134 = t117 * g(1) - t121 * g(2);
t151 = pkin(8) * t111;
t96 = qJDD(1) * pkin(1) + t122 * t151 + t134;
t147 = t113 * t96;
t131 = -t121 * g(1) - t117 * g(2);
t97 = -t122 * pkin(1) + qJDD(1) * t151 + t131;
t148 = t116 * t147 + t120 * t97;
t142 = qJD(1) * t111;
t99 = (-pkin(2) * t120 - pkin(9) * t116) * t142;
t63 = -t104 * pkin(2) + t105 * pkin(9) + (-g(3) * t116 + t141 * t99) * t111 + t148;
t100 = (qJDD(1) * t116 + t120 * t139) * t111;
t150 = t113 * g(3);
t64 = t101 * pkin(2) - t100 * pkin(9) - t150 + (-t96 + (pkin(2) * t116 - pkin(9) * t120) * t106 * qJD(1)) * t111;
t133 = -t115 * t63 + t119 * t64;
t136 = t116 * t142;
t89 = t119 * t106 - t115 * t136;
t75 = t89 * qJD(3) + t119 * t100 + t115 * t105;
t90 = t115 * t106 + t119 * t136;
t93 = qJDD(3) + t101;
t33 = (t103 * t89 - t75) * qJ(4) + (t89 * t90 + t93) * pkin(3) + t133;
t149 = t115 * t64 + t119 * t63;
t74 = -t90 * qJD(3) - t115 * t100 + t119 * t105;
t83 = t103 * pkin(3) - t90 * qJ(4);
t88 = t89 ^ 2;
t37 = -t88 * pkin(3) + t74 * qJ(4) - t103 * t83 + t149;
t80 = t110 * t89 + t146 * t90;
t152 = -t110 * t37 + t146 * t33 + t153 * t80;
t144 = t111 * t116;
t143 = t111 * t120;
t109 = sin(pkin(12));
t112 = cos(pkin(12));
t129 = -g(3) * t143 - t116 * t97 + t120 * t147;
t62 = -t105 * pkin(2) - t104 * pkin(9) + t99 * t136 - t129;
t125 = -t74 * pkin(3) - t88 * qJ(4) + t90 * t83 + qJDD(4) + t62;
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t102 = t103 ^ 2;
t79 = t110 * t90 - t146 * t89;
t137 = t110 * t33 + t146 * t37 + t153 * t79;
t57 = t79 * pkin(4) - t80 * qJ(5);
t25 = -t102 * pkin(4) + t93 * qJ(5) - t79 * t57 + t137;
t54 = t110 * t75 - t146 * t74;
t55 = t110 * t74 + t146 * t75;
t28 = (t103 * t79 - t55) * qJ(5) + (t103 * t80 + t54) * pkin(4) + t125;
t69 = t109 * t103 + t112 * t80;
t132 = -0.2e1 * qJD(5) * t69 - t109 * t25 + t112 * t28;
t45 = t109 * t93 + t112 * t55;
t68 = t112 * t103 - t109 * t80;
t19 = (t68 * t79 - t45) * pkin(10) + (t68 * t69 + t54) * pkin(5) + t132;
t138 = 0.2e1 * qJD(5) * t68 + t109 * t28 + t112 * t25;
t44 = -t109 * t55 + t112 * t93;
t51 = t79 * pkin(5) - t69 * pkin(10);
t67 = t68 ^ 2;
t20 = -t67 * pkin(5) + t44 * pkin(10) - t79 * t51 + t138;
t46 = -t114 * t69 + t118 * t68;
t31 = t46 * qJD(6) + t114 * t44 + t118 * t45;
t47 = t114 * t68 + t118 * t69;
t38 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t78 = qJD(6) + t79;
t41 = -t78 * mrSges(7,2) + t46 * mrSges(7,3);
t53 = qJDD(6) + t54;
t17 = m(7) * (-t114 * t20 + t118 * t19) - t31 * mrSges(7,3) + t53 * mrSges(7,1) - t47 * t38 + t78 * t41;
t30 = -t47 * qJD(6) - t114 * t45 + t118 * t44;
t42 = t78 * mrSges(7,1) - t47 * mrSges(7,3);
t18 = m(7) * (t114 * t19 + t118 * t20) + t30 * mrSges(7,3) - t53 * mrSges(7,2) + t46 * t38 - t78 * t42;
t48 = -t68 * mrSges(6,1) + t69 * mrSges(6,2);
t49 = -t79 * mrSges(6,2) + t68 * mrSges(6,3);
t14 = m(6) * t132 + t54 * mrSges(6,1) - t45 * mrSges(6,3) + t114 * t18 + t118 * t17 - t69 * t48 + t79 * t49;
t50 = t79 * mrSges(6,1) - t69 * mrSges(6,3);
t15 = m(6) * t138 - t54 * mrSges(6,2) + t44 * mrSges(6,3) - t114 * t17 + t118 * t18 + t68 * t48 - t79 * t50;
t70 = -t103 * mrSges(5,2) - t79 * mrSges(5,3);
t71 = t103 * mrSges(5,1) - t80 * mrSges(5,3);
t126 = m(5) * t125 + t54 * mrSges(5,1) + t55 * mrSges(5,2) + t109 * t15 + t112 * t14 + t79 * t70 + t80 * t71;
t82 = -t103 * mrSges(4,2) + t89 * mrSges(4,3);
t84 = t103 * mrSges(4,1) - t90 * mrSges(4,3);
t123 = m(4) * t62 - t74 * mrSges(4,1) + t75 * mrSges(4,2) - t89 * t82 + t90 * t84 + t126;
t95 = -t106 * mrSges(3,2) + mrSges(3,3) * t135;
t98 = (-mrSges(3,1) * t120 + mrSges(3,2) * t116) * t142;
t10 = m(3) * t129 + t105 * mrSges(3,1) - t100 * mrSges(3,3) + t106 * t95 - t136 * t98 - t123;
t58 = t79 * mrSges(5,1) + t80 * mrSges(5,2);
t11 = m(5) * t137 - t93 * mrSges(5,2) - t54 * mrSges(5,3) - t103 * t71 - t109 * t14 + t112 * t15 - t79 * t58;
t24 = -t93 * pkin(4) - t102 * qJ(5) + t80 * t57 + qJDD(5) - t152;
t127 = t30 * mrSges(7,1) + t46 * t41 - m(7) * (-t44 * pkin(5) - t67 * pkin(10) + t69 * t51 + t24) - t31 * mrSges(7,2) - t47 * t42;
t124 = m(6) * t24 - t44 * mrSges(6,1) + t45 * mrSges(6,2) - t68 * t49 + t69 * t50 - t127;
t16 = m(5) * t152 + t93 * mrSges(5,1) - t55 * mrSges(5,3) + t103 * t70 - t80 * t58 - t124;
t81 = -t89 * mrSges(4,1) + t90 * mrSges(4,2);
t7 = m(4) * t133 + t93 * mrSges(4,1) - t75 * mrSges(4,3) + t103 * t82 + t110 * t11 + t146 * t16 - t90 * t81;
t8 = m(4) * t149 - t93 * mrSges(4,2) + t74 * mrSges(4,3) - t103 * t84 + t11 * t146 - t110 * t16 + t89 * t81;
t94 = t106 * mrSges(3,1) - mrSges(3,3) * t136;
t4 = m(3) * (-g(3) * t144 + t148) - t101 * mrSges(3,3) - t105 * mrSges(3,2) + t98 * t135 - t106 * t94 + t119 * t8 - t115 * t7;
t6 = m(3) * (-t111 * t96 - t150) + t100 * mrSges(3,2) + t101 * mrSges(3,1) + t115 * t8 + t119 * t7 + (t116 * t94 - t120 * t95) * t142;
t140 = t10 * t143 + t113 * t6 + t4 * t144;
t2 = m(2) * t131 - t122 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t116 * t10 + t120 * t4;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t122 * mrSges(2,2) - t111 * t6 + (t120 * t10 + t116 * t4) * t113;
t3 = [-m(1) * g(1) - t117 * t1 + t121 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t121 * t1 + t117 * t2, t1, t10, t7, t16, t14, t17; (-m(1) - m(2)) * g(3) + t140, -m(2) * g(3) + t140, t6, t123, t126, t124, -t127;];
f_new  = t3;
