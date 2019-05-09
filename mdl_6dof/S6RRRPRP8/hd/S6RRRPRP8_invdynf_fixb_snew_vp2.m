% Calculate vector of cutting forces with Newton-Euler
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 08:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRPRP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:29:16
% EndTime: 2019-05-07 08:29:21
% DurationCPUTime: 1.80s
% Computational Cost: add. (19166->198), mult. (37773->239), div. (0->0), fcn. (24742->8), ass. (0->93)
t107 = sin(qJ(3));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t146 = cos(qJ(3));
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t137 = qJD(1) * t108;
t87 = -qJD(2) * t146 + t107 * t137;
t136 = t111 * qJD(1);
t98 = qJD(3) - t136;
t145 = t87 * t98;
t135 = qJD(1) * qJD(2);
t130 = t111 * t135;
t131 = t108 * t135;
t114 = qJD(1) ^ 2;
t109 = sin(qJ(1));
t112 = cos(qJ(1));
t129 = t109 * g(1) - g(2) * t112;
t80 = -qJDD(1) * pkin(1) - t114 * pkin(7) - t129;
t91 = qJDD(1) * t108 + t130;
t92 = qJDD(1) * t111 - t131;
t42 = (-t91 - t130) * pkin(8) + (-t92 + t131) * pkin(2) + t80;
t113 = qJD(2) ^ 2;
t125 = -g(1) * t112 - g(2) * t109;
t81 = -pkin(1) * t114 + qJDD(1) * pkin(7) + t125;
t132 = -g(3) * t108 + t111 * t81;
t90 = (-pkin(2) * t111 - pkin(8) * t108) * qJD(1);
t47 = -pkin(2) * t113 + qJDD(2) * pkin(8) + t136 * t90 + t132;
t126 = -t107 * t47 + t146 * t42;
t88 = qJD(2) * t107 + t137 * t146;
t66 = pkin(3) * t87 - qJ(4) * t88;
t86 = qJDD(3) - t92;
t97 = t98 ^ 2;
t26 = -t86 * pkin(3) - t97 * qJ(4) + t66 * t88 + qJDD(4) - t126;
t61 = -qJD(3) * t87 + qJDD(2) * t107 + t146 * t91;
t18 = (-t61 - t145) * pkin(9) + (t87 * t88 - t86) * pkin(4) + t26;
t141 = t107 * t42 + t146 * t47;
t149 = 2 * qJD(4);
t122 = -pkin(3) * t97 + qJ(4) * t86 + t149 * t98 - t66 * t87 + t141;
t60 = qJD(3) * t88 - qJDD(2) * t146 + t107 * t91;
t73 = -pkin(4) * t98 - pkin(9) * t88;
t85 = t87 ^ 2;
t21 = -pkin(4) * t85 + t60 * pkin(9) + t73 * t98 + t122;
t142 = t106 * t18 + t110 * t21;
t64 = t106 * t87 + t110 * t88;
t32 = -t64 * qJD(5) - t106 * t61 + t110 * t60;
t63 = -t106 * t88 + t110 * t87;
t39 = -mrSges(7,1) * t63 + mrSges(7,2) * t64;
t96 = qJD(5) - t98;
t50 = pkin(5) * t96 - t64 * qJ(6);
t62 = t63 ^ 2;
t133 = m(7) * (-t62 * pkin(5) + t32 * qJ(6) + 0.2e1 * qJD(6) * t63 - t50 * t96 + t142) + t32 * mrSges(7,3) + t63 * t39;
t40 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t51 = mrSges(7,1) * t96 - t64 * mrSges(7,3);
t52 = mrSges(6,1) * t96 - t64 * mrSges(6,3);
t84 = qJDD(5) - t86;
t11 = m(6) * t142 + t32 * mrSges(6,3) + t63 * t40 + (-t52 - t51) * t96 + (-mrSges(6,2) - mrSges(7,2)) * t84 + t133;
t71 = -mrSges(5,1) * t98 + mrSges(5,2) * t88;
t128 = -t106 * t21 + t110 * t18;
t33 = t63 * qJD(5) + t106 * t60 + t110 * t61;
t48 = -mrSges(7,2) * t96 + t63 * mrSges(7,3);
t134 = m(7) * (-0.2e1 * qJD(6) * t64 + (t63 * t96 - t33) * qJ(6) + (t63 * t64 + t84) * pkin(5) + t128) + t96 * t48 + t84 * mrSges(7,1);
t49 = -mrSges(6,2) * t96 + t63 * mrSges(6,3);
t9 = m(6) * t128 + t84 * mrSges(6,1) + t96 * t49 + (-t40 - t39) * t64 + (-mrSges(6,3) - mrSges(7,3)) * t33 + t134;
t123 = m(5) * t122 + mrSges(5,3) * t86 - t106 * t9 + t11 * t110 + t71 * t98;
t67 = mrSges(5,1) * t87 - mrSges(5,3) * t88;
t140 = -mrSges(4,1) * t87 - mrSges(4,2) * t88 - t67;
t143 = -mrSges(4,3) - mrSges(5,2);
t70 = mrSges(4,1) * t98 - mrSges(4,3) * t88;
t5 = m(4) * t141 - t86 * mrSges(4,2) + t140 * t87 + t143 * t60 - t98 * t70 + t123;
t121 = -m(5) * t26 - t106 * t11 - t110 * t9;
t69 = -mrSges(4,2) * t98 - mrSges(4,3) * t87;
t72 = -mrSges(5,2) * t87 + mrSges(5,3) * t98;
t6 = m(4) * t126 + (t69 + t72) * t98 + t140 * t88 + (mrSges(4,1) + mrSges(5,1)) * t86 + t143 * t61 + t121;
t93 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t137;
t94 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t136;
t151 = m(3) * t80 - t92 * mrSges(3,1) + t91 * mrSges(3,2) + (t108 * t93 - t111 * t94) * qJD(1) + t107 * t5 + t146 * t6;
t138 = -g(3) * t111 - t108 * t81;
t46 = -qJDD(2) * pkin(2) - pkin(8) * t113 + t137 * t90 - t138;
t120 = t60 * pkin(3) + t46 + (t145 - t61) * qJ(4);
t147 = pkin(3) * t98;
t116 = -t60 * pkin(4) - t85 * pkin(9) - t120 + (-t147 + t149 + t73) * t88;
t127 = m(7) * (-t32 * pkin(5) - t62 * qJ(6) + t64 * t50 + qJDD(6) + t116) + t33 * mrSges(7,2) - t32 * mrSges(7,1) + t64 * t51 - t63 * t48;
t119 = m(6) * t116 - t32 * mrSges(6,1) + t33 * mrSges(6,2) - t63 * t49 + t64 * t52 + t127;
t118 = m(5) * ((-(2 * qJD(4)) + t147) * t88 + t120) + t60 * mrSges(5,1) + t87 * t72 - t119;
t150 = m(4) * t46 + t60 * mrSges(4,1) + (t70 - t71) * t88 + (mrSges(4,2) - mrSges(5,3)) * t61 + t87 * t69 + t118;
t89 = (-mrSges(3,1) * t111 + mrSges(3,2) * t108) * qJD(1);
t4 = m(3) * t132 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t92 - qJD(2) * t93 - t107 * t6 + t136 * t89 + t146 * t5;
t8 = m(3) * t138 + qJDD(2) * mrSges(3,1) - t91 * mrSges(3,3) + qJD(2) * t94 - t137 * t89 - t150;
t148 = t108 * t4 + t111 * t8;
t2 = m(2) * t129 + qJDD(1) * mrSges(2,1) - t114 * mrSges(2,2) - t151;
t1 = m(2) * t125 - mrSges(2,1) * t114 - qJDD(1) * mrSges(2,2) - t108 * t8 + t111 * t4;
t3 = [-m(1) * g(1) + t1 * t112 - t109 * t2, t1, t4, t5, -t60 * mrSges(5,2) - t87 * t67 + t123, t11, -t84 * mrSges(7,2) - t96 * t51 + t133; -m(1) * g(2) + t1 * t109 + t112 * t2, t2, t8, t6, -t61 * mrSges(5,3) - t88 * t71 + t118, t9, -t33 * mrSges(7,3) - t64 * t39 + t134; (-m(1) - m(2)) * g(3) + t148, -m(2) * g(3) + t148, t151, t150, -t86 * mrSges(5,1) + t61 * mrSges(5,2) + t88 * t67 - t98 * t72 - t121, t119, t127;];
f_new  = t3;
