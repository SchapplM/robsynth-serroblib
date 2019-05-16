% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 07:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 06:23:16
% EndTime: 2019-05-06 06:24:05
% DurationCPUTime: 44.74s
% Computational Cost: add. (734798->229), mult. (2396616->341), div. (0->0), fcn. (2107126->18), ass. (0->140)
t103 = sin(pkin(8));
t107 = cos(pkin(8));
t113 = sin(qJ(3));
t108 = cos(pkin(7));
t118 = cos(qJ(3));
t151 = t108 * t118;
t104 = sin(pkin(7));
t155 = t104 * t118;
t102 = sin(pkin(14));
t105 = sin(pkin(6));
t109 = cos(pkin(6));
t106 = cos(pkin(14));
t167 = pkin(10) * t102;
t129 = -pkin(2) * t106 - t104 * t167;
t150 = qJD(1) * t105;
t159 = pkin(10) * qJDD(1);
t126 = qJD(1) * t129 * t150 + t108 * t159;
t145 = qJD(2) * t150;
t154 = t105 * t106;
t120 = qJD(1) ^ 2;
t114 = sin(qJ(1));
t119 = cos(qJ(1));
t144 = t114 * g(1) - g(2) * t119;
t158 = qJ(2) * t105;
t94 = qJDD(1) * pkin(1) + t120 * t158 + t144;
t161 = t109 * t94;
t130 = -g(3) * t154 - 0.2e1 * t102 * t145 + t106 * t161;
t153 = t105 * t108;
t91 = (t104 * t109 + t106 * t153) * qJD(1) * pkin(10);
t139 = -g(1) * t119 - g(2) * t114;
t95 = -pkin(1) * t120 + qJDD(1) * t158 + t139;
t63 = (pkin(2) * qJDD(1) + qJD(1) * t91) * t109 + (-t126 * t105 - t95) * t102 + t130;
t146 = t102 * t161 + (0.2e1 * t145 + t95) * t106;
t96 = (pkin(2) * t109 - t153 * t167) * qJD(1);
t64 = (-qJD(1) * t96 + t104 * t159) * t109 + (-g(3) * t102 + t126 * t106) * t105 + t146;
t141 = -g(3) * t109 + qJDD(2);
t72 = (-t94 + t129 * qJDD(1) + (t102 * t96 - t106 * t91) * qJD(1)) * t105 + t141;
t140 = -t113 * t64 + t151 * t63 + t155 * t72;
t165 = pkin(11) * t107;
t166 = pkin(11) * t103;
t122 = t109 * t155 + (-t102 * t113 + t106 * t151) * t105;
t84 = t122 * qJD(1);
t152 = t108 * t113;
t156 = t104 * t113;
t123 = t109 * t156 + (t102 * t118 + t106 * t152) * t105;
t85 = t123 * qJD(1);
t74 = -pkin(3) * t84 - t166 * t85;
t77 = qJD(3) * t84 + qJDD(1) * t123;
t125 = -t104 * t154 + t108 * t109;
t92 = qJD(1) * t125 + qJD(3);
t134 = t103 * t92 + t107 * t84;
t79 = t134 * pkin(11);
t89 = qJDD(1) * t125 + qJDD(3);
t32 = pkin(3) * t89 - t165 * t77 - t74 * t85 + t79 * t92 + t140;
t142 = -t104 * t63 + t108 * t72;
t76 = -qJD(3) * t85 + qJDD(1) * t122;
t81 = pkin(3) * t92 - t165 * t85;
t36 = -pkin(3) * t76 - t166 * t77 - t79 * t84 + t81 * t85 + t142;
t169 = t103 * t36 + t107 * t32;
t112 = sin(qJ(4));
t117 = cos(qJ(4));
t135 = t103 * t89 + t107 * t76;
t147 = t118 * t64 + t152 * t63 + t156 * t72;
t33 = pkin(11) * t135 + t74 * t84 - t81 * t92 + t147;
t168 = -t112 * t33 + t117 * t169;
t71 = t112 * t134 + t117 * t85;
t49 = -qJD(4) * t71 - t112 * t77 + t117 * t135;
t70 = -t112 * t85 + t117 * t134;
t111 = sin(qJ(5));
t116 = cos(qJ(5));
t148 = t112 * t169 + t117 * t33;
t54 = -pkin(4) * t70 - pkin(12) * t71;
t73 = -t103 * t76 + t107 * t89 + qJDD(4);
t80 = -t103 * t84 + t107 * t92 + qJD(4);
t78 = t80 ^ 2;
t24 = -pkin(4) * t78 + pkin(12) * t73 + t54 * t70 + t148;
t143 = -t103 * t32 + t107 * t36;
t50 = qJD(4) * t70 + t112 * t135 + t117 * t77;
t26 = (-t70 * t80 - t50) * pkin(12) + (t71 * t80 - t49) * pkin(4) + t143;
t164 = t111 * t26 + t116 * t24;
t157 = t102 * t105;
t110 = sin(qJ(6));
t115 = cos(qJ(6));
t56 = -t111 * t71 + t116 * t80;
t57 = t111 * t80 + t116 * t71;
t44 = -pkin(5) * t56 - pkin(13) * t57;
t48 = qJDD(5) - t49;
t68 = qJD(5) - t70;
t67 = t68 ^ 2;
t20 = -pkin(5) * t67 + pkin(13) * t48 + t44 * t56 + t164;
t23 = -pkin(4) * t73 - pkin(12) * t78 + t54 * t71 - t168;
t39 = -qJD(5) * t57 - t111 * t50 + t116 * t73;
t40 = qJD(5) * t56 + t111 * t73 + t116 * t50;
t21 = (-t56 * t68 - t40) * pkin(13) + (t57 * t68 - t39) * pkin(5) + t23;
t46 = -t110 * t57 + t115 * t68;
t28 = qJD(6) * t46 + t110 * t48 + t115 * t40;
t47 = t110 * t68 + t115 * t57;
t37 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t38 = qJDD(6) - t39;
t55 = qJD(6) - t56;
t41 = -mrSges(7,2) * t55 + mrSges(7,3) * t46;
t17 = m(7) * (-t110 * t20 + t115 * t21) - t28 * mrSges(7,3) + t38 * mrSges(7,1) - t47 * t37 + t55 * t41;
t27 = -qJD(6) * t47 - t110 * t40 + t115 * t48;
t42 = mrSges(7,1) * t55 - mrSges(7,3) * t47;
t18 = m(7) * (t110 * t21 + t115 * t20) + t27 * mrSges(7,3) - t38 * mrSges(7,2) + t46 * t37 - t55 * t42;
t43 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t52 = mrSges(6,1) * t68 - mrSges(6,3) * t57;
t15 = m(6) * t164 - mrSges(6,2) * t48 + mrSges(6,3) * t39 - t110 * t17 + t115 * t18 + t43 * t56 - t52 * t68;
t133 = -t111 * t24 + t116 * t26;
t124 = m(7) * (-pkin(5) * t48 - pkin(13) * t67 + t44 * t57 - t133) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t46 * t41 + t47 * t42;
t51 = -mrSges(6,2) * t68 + mrSges(6,3) * t56;
t16 = m(6) * t133 + mrSges(6,1) * t48 - mrSges(6,3) * t40 - t43 * t57 + t51 * t68 - t124;
t58 = -mrSges(5,2) * t80 + mrSges(5,3) * t70;
t59 = mrSges(5,1) * t80 - mrSges(5,3) * t71;
t13 = m(5) * t143 - mrSges(5,1) * t49 + mrSges(5,2) * t50 + t111 * t15 + t116 * t16 - t58 * t70 + t59 * t71;
t53 = -mrSges(5,1) * t70 + mrSges(5,2) * t71;
t12 = m(5) * t148 - mrSges(5,2) * t73 + mrSges(5,3) * t49 - t111 * t16 + t116 * t15 + t53 * t70 - t59 * t80;
t121 = m(6) * t23 - t39 * mrSges(6,1) + t40 * mrSges(6,2) + t110 * t18 + t115 * t17 - t56 * t51 + t57 * t52;
t14 = m(5) * t168 + t73 * mrSges(5,1) - t50 * mrSges(5,3) - t71 * t53 + t80 * t58 - t121;
t132 = t112 * t12 + t117 * t14;
t82 = -mrSges(4,2) * t92 + mrSges(4,3) * t84;
t83 = mrSges(4,1) * t92 - mrSges(4,3) * t85;
t10 = m(4) * t142 - t76 * mrSges(4,1) + t77 * mrSges(4,2) + t103 * t132 + t107 * t13 - t84 * t82 + t85 * t83;
t128 = mrSges(3,1) * t109 - mrSges(3,3) * t157;
t75 = -mrSges(4,1) * t84 + mrSges(4,2) * t85;
t11 = m(4) * t147 - mrSges(4,2) * t89 + mrSges(4,3) * t76 - t112 * t14 + t117 * t12 + t75 * t84 - t83 * t92;
t9 = m(4) * t140 + t89 * mrSges(4,1) - t77 * mrSges(4,3) - t103 * t13 + t107 * t132 - t85 * t75 + t92 * t82;
t138 = t11 * t113 + t118 * t9;
t137 = -mrSges(3,1) * t106 + mrSges(3,2) * t102;
t93 = t137 * t150;
t127 = -mrSges(3,2) * t109 + mrSges(3,3) * t154;
t98 = t127 * qJD(1);
t4 = m(3) * (-t102 * t95 + t130) - t104 * t10 + t138 * t108 + t128 * qJDD(1) + (t109 * t98 - t157 * t93) * qJD(1);
t97 = t128 * qJD(1);
t6 = m(3) * t141 + t108 * t10 + t138 * t104 + (-m(3) * t94 + t137 * qJDD(1) + (t102 * t97 - t106 * t98) * qJD(1)) * t105;
t8 = m(3) * (-g(3) * t157 + t146) + t118 * t11 - t113 * t9 + t127 * qJDD(1) + (-t109 * t97 + t154 * t93) * qJD(1);
t149 = t109 * t6 + t154 * t4 + t157 * t8;
t2 = m(2) * t139 - mrSges(2,1) * t120 - qJDD(1) * mrSges(2,2) - t102 * t4 + t106 * t8;
t1 = m(2) * t144 + qJDD(1) * mrSges(2,1) - t120 * mrSges(2,2) - t105 * t6 + (t102 * t8 + t106 * t4) * t109;
t3 = [-m(1) * g(1) - t1 * t114 + t119 * t2, t2, t8, t11, t12, t15, t18; -m(1) * g(2) + t1 * t119 + t114 * t2, t1, t4, t9, t14, t16, t17; (-m(1) - m(2)) * g(3) + t149, -m(2) * g(3) + t149, t6, t10, t13, t121, t124;];
f_new  = t3;
