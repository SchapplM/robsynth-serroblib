% Calculate vector of cutting forces with Newton-Euler
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 14:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRPR5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:55:41
% EndTime: 2019-05-06 13:56:00
% DurationCPUTime: 9.11s
% Computational Cost: add. (149556->212), mult. (396340->297), div. (0->0), fcn. (315561->14), ass. (0->115)
t153 = -2 * qJD(3);
t109 = sin(pkin(11));
t112 = cos(pkin(11));
t116 = sin(qJ(2));
t119 = cos(qJ(2));
t110 = sin(pkin(6));
t143 = qJD(1) * t110;
t89 = (t109 * t116 - t112 * t119) * t143;
t113 = cos(pkin(6));
t103 = t113 * qJDD(1) + qJDD(2);
t104 = t113 * qJD(1) + qJD(2);
t121 = qJD(1) ^ 2;
t117 = sin(qJ(1));
t120 = cos(qJ(1));
t135 = t117 * g(1) - t120 * g(2);
t150 = pkin(8) * t110;
t95 = qJDD(1) * pkin(1) + t121 * t150 + t135;
t148 = t113 * t95;
t130 = -t120 * g(1) - t117 * g(2);
t96 = -t121 * pkin(1) + qJDD(1) * t150 + t130;
t133 = -t116 * t96 + t119 * t148;
t146 = t110 ^ 2 * t121;
t141 = qJD(1) * qJD(2);
t98 = (qJDD(1) * t116 + t119 * t141) * t110;
t51 = t103 * pkin(2) - t98 * qJ(3) + (pkin(2) * t116 * t146 + (qJ(3) * qJD(1) * t104 - g(3)) * t110) * t119 + t133;
t145 = t110 * t116;
t128 = -g(3) * t145 + t116 * t148 + t119 * t96;
t138 = t119 ^ 2 * t146;
t137 = t116 * t143;
t92 = t104 * pkin(2) - qJ(3) * t137;
t99 = (qJDD(1) * t119 - t116 * t141) * t110;
t55 = -pkin(2) * t138 + t99 * qJ(3) - t104 * t92 + t128;
t90 = (t109 * t119 + t112 * t116) * t143;
t152 = -t109 * t55 + t112 * t51 + t90 * t153;
t151 = cos(qJ(4));
t115 = sin(qJ(4));
t102 = t104 ^ 2;
t139 = t109 * t51 + t112 * t55 + t89 * t153;
t72 = t89 * pkin(3) - t90 * pkin(9);
t35 = -t102 * pkin(3) + t103 * pkin(9) - t89 * t72 + t139;
t129 = -t113 * g(3) - t110 * t95;
t124 = -t99 * pkin(2) - qJ(3) * t138 + t92 * t137 + qJDD(3) + t129;
t75 = -t109 * t98 + t112 * t99;
t76 = t109 * t99 + t112 * t98;
t39 = (t104 * t89 - t76) * pkin(9) + (t104 * t90 - t75) * pkin(3) + t124;
t149 = t115 * t39 + t151 * t35;
t144 = t110 * t119;
t108 = sin(pkin(12));
t111 = cos(pkin(12));
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t78 = -t151 * t104 + t115 * t90;
t79 = t115 * t104 + t151 * t90;
t60 = t78 * pkin(4) - t79 * qJ(5);
t74 = qJDD(4) - t75;
t88 = qJD(4) + t89;
t87 = t88 ^ 2;
t25 = -t87 * pkin(4) + t74 * qJ(5) - t78 * t60 + t149;
t34 = -t103 * pkin(3) - t102 * pkin(9) + t90 * t72 - t152;
t57 = t79 * qJD(4) - t151 * t103 + t115 * t76;
t58 = -t78 * qJD(4) + t115 * t103 + t151 * t76;
t28 = (t78 * t88 - t58) * qJ(5) + (t79 * t88 + t57) * pkin(4) + t34;
t67 = t108 * t88 + t111 * t79;
t132 = -0.2e1 * qJD(5) * t67 - t108 * t25 + t111 * t28;
t43 = t108 * t74 + t111 * t58;
t66 = -t108 * t79 + t111 * t88;
t19 = (t66 * t78 - t43) * pkin(10) + (t66 * t67 + t57) * pkin(5) + t132;
t140 = 0.2e1 * qJD(5) * t66 + t108 * t28 + t111 * t25;
t42 = -t108 * t58 + t111 * t74;
t54 = t78 * pkin(5) - t67 * pkin(10);
t65 = t66 ^ 2;
t20 = -t65 * pkin(5) + t42 * pkin(10) - t78 * t54 + t140;
t44 = -t114 * t67 + t118 * t66;
t31 = t44 * qJD(6) + t114 * t42 + t118 * t43;
t45 = t114 * t66 + t118 * t67;
t37 = -t44 * mrSges(7,1) + t45 * mrSges(7,2);
t77 = qJD(6) + t78;
t40 = -t77 * mrSges(7,2) + t44 * mrSges(7,3);
t56 = qJDD(6) + t57;
t17 = m(7) * (-t114 * t20 + t118 * t19) - t31 * mrSges(7,3) + t56 * mrSges(7,1) - t45 * t37 + t77 * t40;
t30 = -t45 * qJD(6) - t114 * t43 + t118 * t42;
t41 = t77 * mrSges(7,1) - t45 * mrSges(7,3);
t18 = m(7) * (t114 * t19 + t118 * t20) + t30 * mrSges(7,3) - t56 * mrSges(7,2) + t44 * t37 - t77 * t41;
t46 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t52 = -t78 * mrSges(6,2) + t66 * mrSges(6,3);
t13 = m(6) * t132 + t57 * mrSges(6,1) - t43 * mrSges(6,3) + t114 * t18 + t118 * t17 - t67 * t46 + t78 * t52;
t53 = t78 * mrSges(6,1) - t67 * mrSges(6,3);
t14 = m(6) * t140 - t57 * mrSges(6,2) + t42 * mrSges(6,3) - t114 * t17 + t118 * t18 + t66 * t46 - t78 * t53;
t68 = -t88 * mrSges(5,2) - t78 * mrSges(5,3);
t69 = t88 * mrSges(5,1) - t79 * mrSges(5,3);
t123 = m(5) * t34 + t57 * mrSges(5,1) + t58 * mrSges(5,2) + t108 * t14 + t111 * t13 + t78 * t68 + t79 * t69;
t71 = t89 * mrSges(4,1) + t90 * mrSges(4,2);
t80 = -t104 * mrSges(4,2) - t89 * mrSges(4,3);
t10 = m(4) * t152 + t103 * mrSges(4,1) - t76 * mrSges(4,3) + t104 * t80 - t90 * t71 - t123;
t61 = t78 * mrSges(5,1) + t79 * mrSges(5,2);
t12 = m(5) * t149 - t74 * mrSges(5,2) - t57 * mrSges(5,3) - t108 * t13 + t111 * t14 - t78 * t61 - t88 * t69;
t131 = -t115 * t35 + t151 * t39;
t24 = -t74 * pkin(4) - t87 * qJ(5) + t79 * t60 + qJDD(5) - t131;
t126 = t30 * mrSges(7,1) + t44 * t40 - m(7) * (-t42 * pkin(5) - t65 * pkin(10) + t67 * t54 + t24) - t31 * mrSges(7,2) - t45 * t41;
t122 = m(6) * t24 - t42 * mrSges(6,1) + t43 * mrSges(6,2) - t66 * t52 + t67 * t53 - t126;
t16 = m(5) * t131 + t74 * mrSges(5,1) - t58 * mrSges(5,3) - t79 * t61 + t88 * t68 - t122;
t81 = t104 * mrSges(4,1) - t90 * mrSges(4,3);
t7 = m(4) * t139 - t103 * mrSges(4,2) + t75 * mrSges(4,3) - t104 * t81 - t115 * t16 + t151 * t12 - t89 * t71;
t136 = t119 * t143;
t94 = -t104 * mrSges(3,2) + mrSges(3,3) * t136;
t97 = (-mrSges(3,1) * t119 + mrSges(3,2) * t116) * t143;
t5 = m(3) * (-g(3) * t144 + t133) - t98 * mrSges(3,3) + t103 * mrSges(3,1) - t97 * t137 + t104 * t94 + t109 * t7 + t112 * t10;
t93 = t104 * mrSges(3,1) - mrSges(3,3) * t137;
t6 = m(3) * t128 - t103 * mrSges(3,2) + t99 * mrSges(3,3) - t109 * t10 - t104 * t93 + t112 * t7 + t97 * t136;
t125 = m(4) * t124 - t75 * mrSges(4,1) + t76 * mrSges(4,2) + t115 * t12 + t151 * t16 + t89 * t80 + t90 * t81;
t9 = m(3) * t129 + t98 * mrSges(3,2) - t99 * mrSges(3,1) + (t116 * t93 - t119 * t94) * t143 + t125;
t142 = t113 * t9 + t5 * t144 + t6 * t145;
t2 = m(2) * t130 - t121 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t116 * t5 + t119 * t6;
t1 = m(2) * t135 + qJDD(1) * mrSges(2,1) - t121 * mrSges(2,2) - t110 * t9 + (t116 * t6 + t119 * t5) * t113;
t3 = [-m(1) * g(1) - t117 * t1 + t120 * t2, t2, t6, t7, t12, t14, t18; -m(1) * g(2) + t120 * t1 + t117 * t2, t1, t5, t10, t16, t13, t17; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t9, t125, t123, t122, -t126;];
f_new  = t3;
