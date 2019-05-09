% Calculate vector of cutting forces with Newton-Euler
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 00:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRPRRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:37:11
% EndTime: 2019-05-07 00:37:22
% DurationCPUTime: 4.80s
% Computational Cost: add. (62474->216), mult. (139621->283), div. (0->0), fcn. (104490->12), ass. (0->114)
t158 = -2 * qJD(3);
t107 = cos(pkin(6));
t101 = qJD(1) * t107 + qJD(2);
t111 = sin(qJ(2));
t106 = sin(pkin(6));
t142 = qJD(1) * t106;
t99 = t111 * t142;
t157 = (pkin(2) * t101 + t158) * t99;
t100 = qJDD(1) * t107 + qJDD(2);
t116 = cos(qJ(2));
t144 = t106 * t111;
t118 = qJD(1) ^ 2;
t112 = sin(qJ(1));
t117 = cos(qJ(1));
t134 = t112 * g(1) - g(2) * t117;
t82 = pkin(8) * t106 * t118 + qJDD(1) * pkin(1) + t134;
t147 = t107 * t82;
t132 = -g(1) * t117 - g(2) * t112;
t139 = qJDD(1) * t106;
t83 = -pkin(1) * t118 + pkin(8) * t139 + t132;
t129 = -g(3) * t144 + t111 * t147 + t116 * t83;
t141 = qJD(1) * t116;
t135 = t106 * t141;
t84 = (-pkin(2) * t116 - qJ(3) * t111) * t142;
t98 = t101 ^ 2;
t156 = pkin(2) * t98 - t100 * qJ(3) + t101 * t158 - t84 * t135 - t129;
t155 = -pkin(2) - pkin(9);
t154 = t107 * g(3);
t153 = mrSges(3,1) - mrSges(4,2);
t152 = mrSges(3,3) + mrSges(4,1);
t109 = sin(qJ(5));
t114 = cos(qJ(5));
t110 = sin(qJ(4));
t115 = cos(qJ(4));
t145 = t106 ^ 2 * t118;
t136 = t116 ^ 2 * t145;
t87 = pkin(3) * t99 - pkin(9) * t101;
t88 = (qJD(2) * t141 + qJDD(1) * t111) * t106;
t89 = -qJD(2) * t99 + t116 * t139;
t36 = -pkin(3) * t136 - t154 - t88 * qJ(3) + t155 * t89 + (-t82 + (-qJ(3) * t101 * t116 - t111 * t87) * qJD(1)) * t106 + t157;
t143 = t106 * t116;
t146 = g(3) * t143 + t111 * t83;
t127 = -qJ(3) * t98 + t84 * t99 + qJDD(3) + t146;
t39 = pkin(3) * t88 + t155 * t100 + (-pkin(3) * t101 * t142 - pkin(9) * t111 * t145 - t147) * t116 + t127;
t133 = -t110 * t36 + t115 * t39;
t71 = -t101 * t110 - t115 * t135;
t57 = qJD(4) * t71 + t100 * t115 - t110 * t89;
t72 = t101 * t115 - t110 * t135;
t77 = qJDD(4) + t88;
t94 = t99 + qJD(4);
t22 = (t71 * t94 - t57) * pkin(10) + (t71 * t72 + t77) * pkin(4) + t133;
t150 = t110 * t39 + t115 * t36;
t56 = -qJD(4) * t72 - t100 * t110 - t115 * t89;
t64 = pkin(4) * t94 - pkin(10) * t72;
t70 = t71 ^ 2;
t24 = -pkin(4) * t70 + t56 * pkin(10) - t64 * t94 + t150;
t151 = t109 * t22 + t114 * t24;
t81 = mrSges(4,1) * t99 + mrSges(4,2) * t101;
t149 = mrSges(3,1) * t101 - mrSges(3,3) * t99 - t81;
t85 = (mrSges(4,2) * t116 - mrSges(4,3) * t111) * t142;
t148 = t85 + (-mrSges(3,1) * t116 + mrSges(3,2) * t111) * t142;
t121 = pkin(3) * t89 - pkin(9) * t136 + t101 * t87 - t156;
t108 = sin(qJ(6));
t113 = cos(qJ(6));
t119 = -t56 * pkin(4) - pkin(10) * t70 + t72 * t64 + t121;
t59 = -t109 * t72 + t114 * t71;
t60 = t109 * t71 + t114 * t72;
t47 = -pkin(5) * t59 - pkin(11) * t60;
t75 = qJDD(5) + t77;
t92 = qJD(5) + t94;
t91 = t92 ^ 2;
t19 = -pkin(5) * t91 + pkin(11) * t75 + t59 * t47 + t151;
t31 = -t60 * qJD(5) - t109 * t57 + t114 * t56;
t32 = t59 * qJD(5) + t109 * t56 + t114 * t57;
t20 = t119 + (-t59 * t92 - t32) * pkin(11) + (t60 * t92 - t31) * pkin(5);
t49 = -t108 * t60 + t113 * t92;
t28 = t49 * qJD(6) + t108 * t75 + t113 * t32;
t30 = qJDD(6) - t31;
t50 = t108 * t92 + t113 * t60;
t40 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t58 = qJD(6) - t59;
t41 = -mrSges(7,2) * t58 + mrSges(7,3) * t49;
t16 = m(7) * (-t108 * t19 + t113 * t20) - t28 * mrSges(7,3) + t30 * mrSges(7,1) - t50 * t40 + t58 * t41;
t27 = -t50 * qJD(6) - t108 * t32 + t113 * t75;
t42 = mrSges(7,1) * t58 - mrSges(7,3) * t50;
t17 = m(7) * (t108 * t20 + t113 * t19) + t27 * mrSges(7,3) - t30 * mrSges(7,2) + t49 * t40 - t58 * t42;
t51 = -mrSges(6,2) * t92 + t59 * mrSges(6,3);
t52 = mrSges(6,1) * t92 - t60 * mrSges(6,3);
t124 = m(6) * t119 - t31 * mrSges(6,1) + t32 * mrSges(6,2) + t108 * t17 + t113 * t16 - t59 * t51 + t60 * t52;
t62 = -mrSges(5,2) * t94 + mrSges(5,3) * t71;
t63 = mrSges(5,1) * t94 - mrSges(5,3) * t72;
t122 = m(5) * t121 - t56 * mrSges(5,1) + t57 * mrSges(5,2) - t71 * t62 + t72 * t63 + t124;
t120 = -m(4) * t156 + t122;
t11 = m(3) * t129 + t120 + t148 * t135 + t152 * t89 - t149 * t101 + (-mrSges(3,2) + mrSges(4,3)) * t100;
t131 = -t106 * t82 - t154;
t46 = -mrSges(6,1) * t59 + mrSges(6,2) * t60;
t12 = m(6) * t151 - t75 * mrSges(6,2) + t31 * mrSges(6,3) - t108 * t16 + t113 * t17 + t59 * t46 - t92 * t52;
t130 = -t109 * t24 + t114 * t22;
t123 = m(7) * (-pkin(5) * t75 - pkin(11) * t91 + t60 * t47 - t130) - t27 * mrSges(7,1) + t28 * mrSges(7,2) - t49 * t41 + t50 * t42;
t13 = m(6) * t130 + t75 * mrSges(6,1) - t32 * mrSges(6,3) - t60 * t46 + t92 * t51 - t123;
t61 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t8 = m(5) * t133 + t77 * mrSges(5,1) - t57 * mrSges(5,3) + t109 * t12 + t114 * t13 - t72 * t61 + t94 * t62;
t80 = -mrSges(4,1) * t135 - mrSges(4,3) * t101;
t9 = m(5) * t150 - t77 * mrSges(5,2) + t56 * mrSges(5,3) - t109 * t13 + t114 * t12 + t71 * t61 - t94 * t63;
t128 = -t110 * t8 + t115 * t9 + m(4) * (-t89 * pkin(2) + (-t101 * t135 - t88) * qJ(3) + t131 + t157) + t80 * t135 - t88 * mrSges(4,3);
t79 = -mrSges(3,2) * t101 + mrSges(3,3) * t135;
t5 = m(3) * t131 + t88 * mrSges(3,2) - t153 * t89 + (t149 * t111 - t116 * t79) * t142 + t128;
t137 = t116 * t147;
t126 = -m(4) * (-pkin(2) * t100 + t127 - t137) - t110 * t9 - t115 * t8;
t6 = m(3) * (t137 - t146) - t152 * t88 + (t79 - t80) * t101 + t153 * t100 - t148 * t99 + t126;
t138 = t107 * t5 + t11 * t144 + t6 * t143;
t2 = m(2) * t132 - t118 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t116 * t11 - t111 * t6;
t1 = m(2) * t134 + qJDD(1) * mrSges(2,1) - t118 * mrSges(2,2) - t106 * t5 + (t11 * t111 + t116 * t6) * t107;
t3 = [-m(1) * g(1) - t1 * t112 + t117 * t2, t2, t11, t89 * mrSges(4,2) - t81 * t99 + t128, t9, t12, t17; -m(1) * g(2) + t1 * t117 + t112 * t2, t1, t6, -t89 * mrSges(4,1) - t100 * mrSges(4,3) - t101 * t81 - t85 * t135 - t120, t8, t13, t16; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t5, t88 * mrSges(4,1) + t100 * mrSges(4,2) + t101 * t80 + t85 * t99 - t126, t122, t124, t123;];
f_new  = t3;
