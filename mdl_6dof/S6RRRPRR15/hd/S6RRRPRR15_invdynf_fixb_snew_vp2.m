% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 17:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRR15_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:05:38
% EndTime: 2019-05-07 17:05:57
% DurationCPUTime: 8.09s
% Computational Cost: add. (140377->227), mult. (351516->311), div. (0->0), fcn. (291611->14), ass. (0->126)
t107 = sin(pkin(7));
t111 = sin(qJ(3));
t155 = cos(pkin(7));
t138 = t111 * t155;
t166 = cos(qJ(3));
t156 = cos(pkin(6));
t105 = qJDD(1) * t156 + qJDD(2);
t108 = sin(pkin(6));
t116 = cos(qJ(2));
t131 = qJD(1) * t156 + qJD(2);
t112 = sin(qJ(2));
t136 = t116 * t156;
t118 = qJD(1) ^ 2;
t113 = sin(qJ(1));
t117 = cos(qJ(1));
t141 = t113 * g(1) - t117 * g(2);
t165 = pkin(9) * t108;
t96 = qJDD(1) * pkin(1) + t118 * t165 + t141;
t133 = -t117 * g(1) - t113 * g(2);
t97 = -t118 * pkin(1) + qJDD(1) * t165 + t133;
t139 = -t112 * t97 + t96 * t136;
t150 = qJD(1) * t112;
t147 = qJD(1) * qJD(2);
t99 = (qJDD(1) * t112 + t116 * t147) * t108;
t167 = pkin(10) * t99;
t126 = t131 * t107;
t151 = qJD(1) * t108;
t134 = t155 * t151;
t168 = t116 * t134 + t126;
t87 = t168 * pkin(10);
t164 = pkin(10) * t112;
t93 = (-pkin(2) * t116 - t107 * t164) * t151;
t44 = -t155 * t167 + t105 * pkin(2) + t131 * t87 + (-t116 * g(3) - t150 * t93) * t108 + t139;
t100 = (qJDD(1) * t116 - t112 * t147) * t108;
t127 = t100 * t155 + t105 * t107;
t149 = qJD(1) * t116;
t137 = t112 * t156;
t157 = t116 * t97 + t96 * t137;
t92 = pkin(2) * t131 - t134 * t164;
t45 = -t131 * t92 + (-t112 * g(3) + t149 * t93) * t108 + t127 * pkin(10) + t157;
t144 = t156 * g(3);
t50 = -t107 * t167 - t144 - t100 * pkin(2) + (-t96 + (t112 * t92 - t116 * t87) * qJD(1)) * t108;
t146 = t107 * t111 * t50 + t44 * t138 + t166 * t45;
t143 = t108 * t150;
t79 = t111 * t143 - t168 * t166;
t80 = t111 * t126 + (t112 * t166 + t116 * t138) * t151;
t63 = t79 * pkin(3) - t80 * qJ(4);
t82 = -t107 * t100 + t105 * t155 + qJDD(3);
t142 = t108 * t149;
t89 = t107 * t142 - t131 * t155 - qJD(3);
t86 = t89 ^ 2;
t169 = t86 * pkin(3) - t82 * qJ(4) + 0.2e1 * qJD(4) * t89 + t79 * t63 - t146;
t163 = t79 * t89;
t162 = mrSges(4,1) - mrSges(5,2);
t161 = -mrSges(4,3) - mrSges(5,1);
t110 = sin(qJ(5));
t115 = cos(qJ(5));
t135 = t155 * t166;
t145 = t107 * t166;
t120 = -t111 * t45 + t135 * t44 + t145 * t50;
t30 = -t82 * pkin(3) - t86 * qJ(4) + t80 * t63 + qJDD(4) - t120;
t61 = -t79 * qJD(3) + t111 * t127 + t166 * t99;
t24 = (t79 * t80 - t82) * pkin(11) + (t61 - t163) * pkin(4) + t30;
t140 = -t107 * t44 + t155 * t50;
t123 = (-t61 - t163) * qJ(4) + t140 + (-t89 * pkin(3) - 0.2e1 * qJD(4)) * t80;
t60 = t80 * qJD(3) - t100 * t135 - t105 * t145 + t111 * t99;
t73 = t80 * pkin(4) + t89 * pkin(11);
t78 = t79 ^ 2;
t25 = -t78 * pkin(4) - t80 * t73 + (pkin(3) + pkin(11)) * t60 + t123;
t160 = t110 * t24 + t115 * t25;
t65 = -t79 * mrSges(5,2) - t80 * mrSges(5,3);
t159 = -t79 * mrSges(4,1) - t80 * mrSges(4,2) - t65;
t70 = t79 * mrSges(5,1) + t89 * mrSges(5,3);
t158 = t89 * mrSges(4,2) - t79 * mrSges(4,3) - t70;
t153 = t108 * t112;
t152 = t108 * t116;
t109 = sin(qJ(6));
t114 = cos(qJ(6));
t67 = t110 * t89 + t115 * t79;
t68 = t110 * t79 - t115 * t89;
t47 = -t67 * pkin(5) - t68 * pkin(12);
t59 = qJDD(5) + t61;
t77 = qJD(5) + t80;
t76 = t77 ^ 2;
t20 = -t76 * pkin(5) + t59 * pkin(12) + t67 * t47 + t160;
t119 = -t60 * pkin(4) - t78 * pkin(11) - t89 * t73 - t169;
t36 = -t68 * qJD(5) - t110 * t82 + t115 * t60;
t37 = t67 * qJD(5) + t110 * t60 + t115 * t82;
t21 = (-t67 * t77 - t37) * pkin(12) + (t68 * t77 - t36) * pkin(5) + t119;
t53 = -t109 * t68 + t114 * t77;
t32 = t53 * qJD(6) + t109 * t59 + t114 * t37;
t54 = t109 * t77 + t114 * t68;
t33 = -mrSges(7,1) * t53 + mrSges(7,2) * t54;
t35 = qJDD(6) - t36;
t66 = qJD(6) - t67;
t38 = -t66 * mrSges(7,2) + t53 * mrSges(7,3);
t17 = m(7) * (-t109 * t20 + t114 * t21) - t32 * mrSges(7,3) + t35 * mrSges(7,1) - t54 * t33 + t66 * t38;
t31 = -t54 * qJD(6) - t109 * t37 + t114 * t59;
t39 = t66 * mrSges(7,1) - t54 * mrSges(7,3);
t18 = m(7) * (t109 * t21 + t114 * t20) + t31 * mrSges(7,3) - t35 * mrSges(7,2) + t53 * t33 - t66 * t39;
t46 = -t67 * mrSges(6,1) + t68 * mrSges(6,2);
t56 = t77 * mrSges(6,1) - t68 * mrSges(6,3);
t13 = m(6) * t160 - t59 * mrSges(6,2) + t36 * mrSges(6,3) - t109 * t17 + t114 * t18 + t67 * t46 - t77 * t56;
t132 = -t110 * t25 + t115 * t24;
t121 = m(7) * (-t59 * pkin(5) - t76 * pkin(12) + t68 * t47 - t132) - t31 * mrSges(7,1) + t32 * mrSges(7,2) - t53 * t38 + t54 * t39;
t55 = -t77 * mrSges(6,2) + t67 * mrSges(6,3);
t14 = m(6) * t132 + t59 * mrSges(6,1) - t37 * mrSges(6,3) - t68 * t46 + t77 * t55 - t121;
t125 = -m(5) * t30 - t110 * t13 - t115 * t14;
t10 = m(4) * t120 - t158 * t89 + t159 * t80 + t161 * t61 + t162 * t82 + t125;
t124 = m(6) * t119 - t36 * mrSges(6,1) + t37 * mrSges(6,2) + t109 * t18 + t114 * t17 - t67 * t55 + t68 * t56;
t122 = -m(5) * t169 + t124;
t71 = t80 * mrSges(5,1) - t89 * mrSges(5,2);
t72 = -t89 * mrSges(4,1) - t80 * mrSges(4,3);
t11 = m(4) * t146 + (t72 - t71) * t89 + (-mrSges(4,2) + mrSges(5,3)) * t82 + t159 * t79 + t161 * t60 + t122;
t129 = -t110 * t14 + t115 * t13 + m(5) * (t60 * pkin(3) + t123) - t80 * t71 - t61 * mrSges(5,3);
t9 = m(4) * t140 + t61 * mrSges(4,2) + t158 * t79 + t162 * t60 + t80 * t72 + t129;
t95 = -mrSges(3,2) * t131 + mrSges(3,3) * t142;
t98 = (-mrSges(3,1) * t116 + mrSges(3,2) * t112) * t151;
t4 = m(3) * (-g(3) * t152 + t139) - t99 * mrSges(3,3) + t105 * mrSges(3,1) - t98 * t143 + t131 * t95 + t11 * t138 + t10 * t135 - t107 * t9;
t94 = mrSges(3,1) * t131 - mrSges(3,3) * t143;
t6 = m(3) * (-t108 * t96 - t144) + t99 * mrSges(3,2) - t100 * mrSges(3,1) + t155 * t9 + (t10 * t166 + t111 * t11) * t107 + (t112 * t94 - t116 * t95) * t151;
t8 = m(3) * (-g(3) * t153 + t157) + t100 * mrSges(3,3) - t105 * mrSges(3,2) + t98 * t142 - t131 * t94 + t166 * t11 - t111 * t10;
t148 = t4 * t152 + t8 * t153 + t156 * t6;
t2 = m(2) * t133 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t112 * t4 + t116 * t8;
t1 = m(2) * t141 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t108 * t6 + t136 * t4 + t137 * t8;
t3 = [-m(1) * g(1) - t113 * t1 + t117 * t2, t2, t8, t11, -t60 * mrSges(5,2) - t79 * t70 + t129, t13, t18; -m(1) * g(2) + t117 * t1 + t113 * t2, t1, t4, t10, t60 * mrSges(5,1) - t82 * mrSges(5,3) + t79 * t65 + t89 * t71 - t122, t14, t17; (-m(1) - m(2)) * g(3) + t148, -m(2) * g(3) + t148, t6, t9, t61 * mrSges(5,1) + t82 * mrSges(5,2) + t80 * t65 - t89 * t70 - t125, t124, t121;];
f_new  = t3;
