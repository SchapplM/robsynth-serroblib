% Calculate vector of cutting forces with Newton-Euler
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRRR10_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 18:03:35
% EndTime: 2019-05-08 18:06:19
% DurationCPUTime: 56.31s
% Computational Cost: add. (1004595->238), mult. (2575649->346), div. (0->0), fcn. (2232196->18), ass. (0->141)
t104 = sin(pkin(8));
t107 = cos(pkin(8));
t113 = sin(qJ(3));
t119 = cos(qJ(3));
t105 = sin(pkin(7));
t151 = t105 * t119;
t108 = cos(pkin(7));
t109 = cos(pkin(6));
t101 = t109 * qJDD(1) + qJDD(2);
t102 = t109 * qJD(1) + qJD(2);
t106 = sin(pkin(6));
t120 = cos(qJ(2));
t114 = sin(qJ(2));
t122 = qJD(1) ^ 2;
t115 = sin(qJ(1));
t121 = cos(qJ(1));
t138 = t115 * g(1) - t121 * g(2);
t165 = pkin(10) * t106;
t96 = qJDD(1) * pkin(1) + t122 * t165 + t138;
t154 = t109 * t96;
t133 = -t121 * g(1) - t115 * g(2);
t97 = -t122 * pkin(1) + qJDD(1) * t165 + t133;
t135 = -t114 * t97 + t120 * t154;
t146 = qJD(1) * t114;
t163 = pkin(11) * t108;
t145 = qJD(1) * t120;
t139 = t106 * t145;
t89 = (t102 * t105 + t108 * t139) * pkin(11);
t147 = qJD(1) * t106;
t164 = pkin(11) * t105;
t93 = (-pkin(2) * t120 - t114 * t164) * t147;
t143 = qJD(1) * qJD(2);
t99 = (qJDD(1) * t114 + t120 * t143) * t106;
t63 = -t99 * t163 + t101 * pkin(2) + t102 * t89 + (-g(3) * t120 - t93 * t146) * t106 + t135;
t155 = t108 * t63;
t100 = (qJDD(1) * t120 - t114 * t143) * t106;
t125 = t100 * t108 + t101 * t105;
t158 = t114 * t154 + t120 * t97;
t140 = t106 * t146;
t92 = t102 * pkin(2) - t140 * t163;
t64 = -t102 * t92 + (-g(3) * t114 + t93 * t145) * t106 + t125 * pkin(11) + t158;
t160 = t109 * g(3);
t69 = -t99 * t164 - t100 * pkin(2) - t160 + (-t96 + (t114 * t92 - t120 * t89) * qJD(1)) * t106;
t134 = -t113 * t64 + t119 * t155 + t69 * t151;
t161 = pkin(12) * t107;
t148 = t108 * t120;
t84 = t102 * t151 + (-t113 * t114 + t119 * t148) * t147;
t75 = t84 * qJD(3) + t125 * t113 + t119 * t99;
t162 = pkin(12) * t104;
t152 = t105 * t113;
t85 = t102 * t152 + (t113 * t148 + t114 * t119) * t147;
t76 = -t84 * pkin(3) - t85 * t162;
t90 = t108 * t102 - t105 * t139 + qJD(3);
t129 = t104 * t90 + t107 * t84;
t79 = t129 * pkin(12);
t86 = -t105 * t100 + t108 * t101 + qJDD(3);
t32 = t86 * pkin(3) - t75 * t161 - t85 * t76 + t90 * t79 + t134;
t136 = -t105 * t63 + t108 * t69;
t74 = -t85 * qJD(3) - t113 * t99 + t125 * t119;
t81 = t90 * pkin(3) - t85 * t161;
t39 = -t74 * pkin(3) - t75 * t162 - t84 * t79 + t85 * t81 + t136;
t167 = t104 * t39 + t107 * t32;
t112 = sin(qJ(4));
t118 = cos(qJ(4));
t130 = t104 * t86 + t107 * t74;
t141 = t113 * t155 + t119 * t64 + t69 * t152;
t33 = t130 * pkin(12) + t84 * t76 - t90 * t81 + t141;
t166 = -t112 * t33 + t167 * t118;
t73 = t129 * t112 + t118 * t85;
t46 = -t73 * qJD(4) - t112 * t75 + t130 * t118;
t72 = -t112 * t85 + t129 * t118;
t111 = sin(qJ(5));
t117 = cos(qJ(5));
t142 = t167 * t112 + t118 * t33;
t54 = -t72 * pkin(4) - t73 * pkin(13);
t65 = -t104 * t74 + t107 * t86 + qJDD(4);
t80 = -t104 * t84 + t107 * t90 + qJD(4);
t78 = t80 ^ 2;
t24 = -t78 * pkin(4) + t65 * pkin(13) + t72 * t54 + t142;
t137 = -t104 * t32 + t107 * t39;
t47 = t72 * qJD(4) + t130 * t112 + t118 * t75;
t26 = (-t72 * t80 - t47) * pkin(13) + (t73 * t80 - t46) * pkin(4) + t137;
t159 = t111 * t26 + t117 * t24;
t150 = t106 * t114;
t149 = t106 * t120;
t110 = sin(qJ(6));
t116 = cos(qJ(6));
t56 = -t111 * t73 + t117 * t80;
t57 = t111 * t80 + t117 * t73;
t44 = -t56 * pkin(5) - t57 * pkin(14);
t45 = qJDD(5) - t46;
t71 = qJD(5) - t72;
t70 = t71 ^ 2;
t20 = -t70 * pkin(5) + t45 * pkin(14) + t56 * t44 + t159;
t23 = -t65 * pkin(4) - t78 * pkin(13) + t73 * t54 - t166;
t35 = -t57 * qJD(5) - t111 * t47 + t117 * t65;
t36 = t56 * qJD(5) + t111 * t65 + t117 * t47;
t21 = (-t56 * t71 - t36) * pkin(14) + (t57 * t71 - t35) * pkin(5) + t23;
t49 = -t110 * t57 + t116 * t71;
t28 = t49 * qJD(6) + t110 * t45 + t116 * t36;
t34 = qJDD(6) - t35;
t50 = t110 * t71 + t116 * t57;
t40 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t55 = qJD(6) - t56;
t41 = -t55 * mrSges(7,2) + t49 * mrSges(7,3);
t17 = m(7) * (-t110 * t20 + t116 * t21) - t28 * mrSges(7,3) + t34 * mrSges(7,1) - t50 * t40 + t55 * t41;
t27 = -t50 * qJD(6) - t110 * t36 + t116 * t45;
t42 = t55 * mrSges(7,1) - t50 * mrSges(7,3);
t18 = m(7) * (t110 * t21 + t116 * t20) + t27 * mrSges(7,3) - t34 * mrSges(7,2) + t49 * t40 - t55 * t42;
t43 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t52 = t71 * mrSges(6,1) - t57 * mrSges(6,3);
t15 = m(6) * t159 - t45 * mrSges(6,2) + t35 * mrSges(6,3) - t110 * t17 + t116 * t18 + t56 * t43 - t71 * t52;
t128 = -t111 * t24 + t117 * t26;
t124 = m(7) * (-t45 * pkin(5) - t70 * pkin(14) + t57 * t44 - t128) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t49 * t41 + t50 * t42;
t51 = -t71 * mrSges(6,2) + t56 * mrSges(6,3);
t16 = m(6) * t128 + t45 * mrSges(6,1) - t36 * mrSges(6,3) - t57 * t43 + t71 * t51 - t124;
t53 = -t72 * mrSges(5,1) + t73 * mrSges(5,2);
t59 = t80 * mrSges(5,1) - t73 * mrSges(5,3);
t12 = m(5) * t142 - t65 * mrSges(5,2) + t46 * mrSges(5,3) - t111 * t16 + t117 * t15 + t72 * t53 - t80 * t59;
t123 = m(6) * t23 - t35 * mrSges(6,1) + t36 * mrSges(6,2) + t110 * t18 + t116 * t17 - t56 * t51 + t57 * t52;
t58 = -t80 * mrSges(5,2) + t72 * mrSges(5,3);
t14 = m(5) * t166 + t65 * mrSges(5,1) - t47 * mrSges(5,3) - t73 * t53 + t80 * t58 - t123;
t127 = t112 * t12 + t118 * t14;
t13 = m(5) * t137 - t46 * mrSges(5,1) + t47 * mrSges(5,2) + t111 * t15 + t117 * t16 - t72 * t58 + t73 * t59;
t82 = -t90 * mrSges(4,2) + t84 * mrSges(4,3);
t83 = t90 * mrSges(4,1) - t85 * mrSges(4,3);
t10 = m(4) * t136 - t74 * mrSges(4,1) + t75 * mrSges(4,2) + t127 * t104 + t107 * t13 - t84 * t82 + t85 * t83;
t77 = -t84 * mrSges(4,1) + t85 * mrSges(4,2);
t11 = m(4) * t141 - t86 * mrSges(4,2) + t74 * mrSges(4,3) - t112 * t14 + t118 * t12 + t84 * t77 - t90 * t83;
t9 = m(4) * t134 + t86 * mrSges(4,1) - t75 * mrSges(4,3) - t104 * t13 + t127 * t107 - t85 * t77 + t90 * t82;
t132 = t113 * t11 + t119 * t9;
t95 = -t102 * mrSges(3,2) + mrSges(3,3) * t139;
t98 = (-mrSges(3,1) * t120 + mrSges(3,2) * t114) * t147;
t4 = m(3) * (-g(3) * t149 + t135) - t99 * mrSges(3,3) + t101 * mrSges(3,1) - t98 * t140 + t102 * t95 - t105 * t10 + t132 * t108;
t94 = t102 * mrSges(3,1) - mrSges(3,3) * t140;
t6 = m(3) * (-t106 * t96 - t160) + t99 * mrSges(3,2) - t100 * mrSges(3,1) + t108 * t10 + t132 * t105 + (t114 * t94 - t120 * t95) * t147;
t8 = m(3) * (-g(3) * t150 + t158) + t100 * mrSges(3,3) - t101 * mrSges(3,2) + t98 * t139 - t102 * t94 + t119 * t11 - t113 * t9;
t144 = t109 * t6 + t4 * t149 + t8 * t150;
t2 = m(2) * t133 - t122 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t114 * t4 + t120 * t8;
t1 = m(2) * t138 + qJDD(1) * mrSges(2,1) - t122 * mrSges(2,2) - t106 * t6 + (t114 * t8 + t120 * t4) * t109;
t3 = [-m(1) * g(1) - t115 * t1 + t121 * t2, t2, t8, t11, t12, t15, t18; -m(1) * g(2) + t121 * t1 + t115 * t2, t1, t4, t9, t14, t16, t17; (-m(1) - m(2)) * g(3) + t144, -m(2) * g(3) + t144, t6, t10, t13, t123, t124;];
f_new  = t3;
