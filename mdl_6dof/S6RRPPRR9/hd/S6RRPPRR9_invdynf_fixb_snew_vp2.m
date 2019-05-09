% Calculate vector of cutting forces with Newton-Euler
% S6RRPPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 11:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPPRR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:31:43
% EndTime: 2019-05-06 11:31:47
% DurationCPUTime: 2.00s
% Computational Cost: add. (20434->224), mult. (47684->277), div. (0->0), fcn. (31862->10), ass. (0->113)
t174 = -2 * qJD(3);
t173 = -2 * qJD(4);
t114 = cos(pkin(6));
t106 = qJDD(1) * t114 + qJDD(2);
t107 = qJD(1) * t114 + qJD(2);
t117 = sin(qJ(2));
t113 = sin(pkin(6));
t121 = cos(qJ(2));
t152 = qJD(1) * t121;
t142 = t113 * t152;
t140 = t107 * t142;
t123 = qJD(1) ^ 2;
t156 = t113 ^ 2 * t123;
t144 = t121 * t156;
t105 = t107 ^ 2;
t153 = qJD(1) * t113;
t143 = t117 * t153;
t154 = t113 * t121;
t118 = sin(qJ(1));
t122 = cos(qJ(1));
t141 = g(1) * t118 - g(2) * t122;
t81 = pkin(8) * t113 * t123 + qJDD(1) * pkin(1) + t141;
t160 = t114 * t81;
t139 = -g(1) * t122 - g(2) * t118;
t148 = qJDD(1) * t113;
t82 = -pkin(1) * t123 + pkin(8) * t148 + t139;
t146 = -g(3) * t154 - t117 * t82 + t121 * t160;
t83 = (-pkin(2) * t121 - qJ(3) * t117) * t153;
t36 = -pkin(2) * t106 - qJ(3) * t105 + t143 * t83 + qJDD(3) - t146;
t88 = (qJD(2) * t152 + qJDD(1) * t117) * t113;
t172 = t36 - (t117 * t144 + t106) * qJ(4) - (t140 - t88) * pkin(3) + t107 * t173;
t165 = t117 * t160 + t121 * t82;
t171 = t105 * pkin(2) - qJ(3) * t106 + t107 * t174 - t142 * t83 - t165;
t170 = t114 * g(3);
t169 = t88 * mrSges(5,2);
t168 = mrSges(4,2) - mrSges(5,3);
t167 = -mrSges(3,3) - mrSges(4,1);
t116 = sin(qJ(5));
t120 = cos(qJ(5));
t111 = t117 ^ 2;
t112 = t121 ^ 2;
t89 = -qJD(2) * t143 + t121 * t148;
t133 = t143 * t174 - t170 + (t107 * t143 - t89) * pkin(2);
t130 = -qJ(4) * t89 + t142 * t173 + t133;
t157 = qJ(3) * t107;
t76 = pkin(3) * t143 - qJ(4) * t107;
t159 = t117 * t76;
t87 = pkin(4) * t142 - pkin(9) * t107;
t22 = (pkin(9) - qJ(3)) * t88 + (-pkin(3) * t112 - pkin(4) * t111) * t156 + (-t81 + (-t159 + (-t87 - t157) * t121) * qJD(1)) * t113 + t130;
t145 = t112 * t156;
t125 = -qJ(4) * t145 + t107 * t76 + qJDD(4) - t171;
t26 = -t106 * pkin(9) + (pkin(3) + pkin(4)) * t89 + (pkin(9) * t144 + (pkin(4) * qJD(1) * t107 - g(3)) * t113) * t117 + t125;
t166 = t116 * t26 + t120 * t22;
t77 = mrSges(5,1) * t143 - mrSges(5,3) * t107;
t80 = mrSges(4,1) * t143 + mrSges(4,2) * t107;
t164 = -t77 - t80;
t78 = -mrSges(4,1) * t142 - mrSges(4,3) * t107;
t79 = mrSges(5,1) * t142 + mrSges(5,2) * t107;
t163 = t78 - t79;
t84 = (mrSges(4,2) * t121 - mrSges(4,3) * t117) * t153;
t162 = t84 + (-mrSges(3,1) * t121 + mrSges(3,2) * t117) * t153;
t161 = t113 * t81;
t158 = t121 * t79;
t155 = t113 * t117;
t115 = sin(qJ(6));
t119 = cos(qJ(6));
t64 = -t107 * t116 + t120 * t143;
t65 = t107 * t120 + t116 * t143;
t48 = -pkin(5) * t64 - pkin(10) * t65;
t72 = qJDD(5) + t89;
t96 = qJD(5) + t142;
t93 = t96 ^ 2;
t19 = -pkin(5) * t93 + pkin(10) * t72 + t48 * t64 + t166;
t124 = -t111 * pkin(9) * t156 - t88 * pkin(4) + t107 * t87 - t172;
t45 = -qJD(5) * t65 - t106 * t116 + t120 * t88;
t46 = qJD(5) * t64 + t106 * t120 + t116 * t88;
t20 = (-t64 * t96 - t46) * pkin(10) + (t65 * t96 - t45) * pkin(5) + t124;
t49 = -t115 * t65 + t119 * t96;
t33 = qJD(6) * t49 + t115 * t72 + t119 * t46;
t50 = t115 * t96 + t119 * t65;
t37 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t63 = qJD(6) - t64;
t38 = -mrSges(7,2) * t63 + mrSges(7,3) * t49;
t42 = qJDD(6) - t45;
t16 = m(7) * (-t115 * t19 + t119 * t20) - t33 * mrSges(7,3) + t42 * mrSges(7,1) - t50 * t37 + t63 * t38;
t32 = -qJD(6) * t50 - t115 * t46 + t119 * t72;
t39 = mrSges(7,1) * t63 - mrSges(7,3) * t50;
t17 = m(7) * (t115 * t20 + t119 * t19) + t32 * mrSges(7,3) - t42 * mrSges(7,2) + t49 * t37 - t63 * t39;
t47 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t52 = mrSges(6,1) * t96 - mrSges(6,3) * t65;
t11 = m(6) * t166 - mrSges(6,2) * t72 + mrSges(6,3) * t45 - t115 * t16 + t119 * t17 + t47 * t64 - t52 * t96;
t138 = -t116 * t22 + t120 * t26;
t126 = m(7) * (-pkin(5) * t72 - pkin(10) * t93 + t48 * t65 - t138) - t32 * mrSges(7,1) + t33 * mrSges(7,2) - t49 * t38 + t50 * t39;
t51 = -mrSges(6,2) * t96 + mrSges(6,3) * t64;
t13 = m(6) * t138 + mrSges(6,1) * t72 - mrSges(6,3) * t46 - t47 * t65 + t51 * t96 - t126;
t137 = t120 * t11 - t116 * t13 + m(5) * (-pkin(3) * t145 - t88 * qJ(3) + (-t81 + (-t121 * t157 - t159) * qJD(1)) * t113 + t130) - t89 * mrSges(5,3);
t131 = t137 + m(4) * (-t161 + (-t88 - t140) * qJ(3) + t133) + t78 * t142 - t88 * mrSges(4,3);
t74 = mrSges(3,1) * t107 - mrSges(3,3) * t143;
t75 = -mrSges(3,2) * t107 + mrSges(3,3) * t142;
t4 = m(3) * (-t161 - t170) + (-mrSges(3,1) + mrSges(4,2)) * t89 + (mrSges(3,2) - mrSges(5,2)) * t88 + ((-t75 - t79) * t121 + (t74 + t164) * t117) * t153 + t131;
t147 = g(3) * t155;
t86 = (-mrSges(5,2) * t117 - mrSges(5,3) * t121) * t153;
t136 = t116 * t11 + t120 * t13 + t106 * mrSges(5,2) + m(5) * (t89 * pkin(3) + t125 - t147) + t107 * t77 + t86 * t142;
t132 = m(4) * (t147 + t171) - t136;
t6 = m(3) * (-t147 + t165) + (-t74 + t80) * t107 + (-mrSges(3,2) + mrSges(4,3)) * t106 + t162 * t142 + (mrSges(5,1) - t167) * t89 - t132;
t134 = m(6) * t124 - mrSges(6,1) * t45 + mrSges(6,2) * t46 + t115 * t17 + t119 * t16 - t51 * t64 + t52 * t65;
t129 = m(5) * t172 + mrSges(5,1) * t88 + t143 * t86 - t134;
t127 = m(4) * t36 + t129;
t8 = -t162 * t143 - t127 + t167 * t88 + (t75 - t163) * t107 + (mrSges(3,1) - t168) * t106 + m(3) * t146;
t149 = t114 * t4 + t154 * t8 + t155 * t6;
t2 = m(2) * t139 - mrSges(2,1) * t123 - qJDD(1) * mrSges(2,2) - t117 * t8 + t121 * t6;
t1 = m(2) * t141 + qJDD(1) * mrSges(2,1) - t123 * mrSges(2,2) - t113 * t4 + (t117 * t6 + t121 * t8) * t114;
t3 = [-m(1) * g(1) - t1 * t118 + t122 * t2, t2, t6, t89 * mrSges(4,2) - t169 + (t117 * t164 - t158) * t153 + t131, -t169 + (-t117 * t77 - t158) * t153 + t137, t11, t17; -m(1) * g(2) + t1 * t122 + t118 * t2, t1, t8, -t84 * t142 - t106 * mrSges(4,3) - t107 * t80 + (-mrSges(4,1) - mrSges(5,1)) * t89 + t132, -t106 * mrSges(5,3) - t107 * t79 + t129, t13, t16; (-m(1) - m(2)) * g(3) + t149, -m(2) * g(3) + t149, t4, t88 * mrSges(4,1) + t106 * t168 + t107 * t163 + t143 * t84 + t127, t89 * mrSges(5,1) + t136, t134, t126;];
f_new  = t3;
