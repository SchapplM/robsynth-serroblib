% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 18:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:01:14
% EndTime: 2019-05-07 18:01:20
% DurationCPUTime: 1.81s
% Computational Cost: add. (20951->205), mult. (42068->244), div. (0->0), fcn. (28896->8), ass. (0->93)
t101 = sin(qJ(4));
t140 = cos(qJ(4));
t105 = cos(qJ(2));
t100 = t105 ^ 2;
t107 = qJD(1) ^ 2;
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t124 = g(1) * t104 - t106 * g(2);
t116 = -qJDD(1) * pkin(1) - t124;
t103 = sin(qJ(2));
t129 = qJD(1) * t103;
t127 = qJD(1) * qJD(2);
t89 = qJDD(1) * t105 - t103 * t127;
t92 = qJD(2) * pkin(2) - pkin(8) * t129;
t109 = -pkin(2) * t89 + t92 * t129 + (-pkin(8) * t100 - pkin(7)) * t107 + t116;
t102 = sin(qJ(3));
t141 = cos(qJ(3));
t82 = (t102 * t105 + t141 * t103) * qJD(1);
t88 = qJDD(1) * t103 + t105 * t127;
t59 = -qJD(3) * t82 - t102 * t88 + t141 * t89;
t128 = qJD(1) * t105;
t81 = -t102 * t129 + t141 * t128;
t60 = t81 * qJD(3) + t102 * t89 + t141 * t88;
t99 = qJD(2) + qJD(3);
t23 = (-t81 * t99 - t60) * pkin(9) + (t82 * t99 - t59) * pkin(3) + t109;
t119 = -g(1) * t106 - g(2) * t104;
t84 = -pkin(1) * t107 + qJDD(1) * pkin(7) + t119;
t130 = t103 * t84;
t139 = pkin(2) * t107;
t50 = qJDD(2) * pkin(2) - pkin(8) * t88 - t130 + (pkin(8) * t127 + t103 * t139 - g(3)) * t105;
t123 = -g(3) * t103 + t105 * t84;
t51 = pkin(8) * t89 - qJD(2) * t92 - t100 * t139 + t123;
t133 = t102 * t50 + t141 * t51;
t71 = -pkin(3) * t81 - pkin(9) * t82;
t97 = t99 ^ 2;
t98 = qJDD(2) + qJDD(3);
t27 = -pkin(3) * t97 + pkin(9) * t98 + t71 * t81 + t133;
t135 = t101 * t23 + t140 * t27;
t143 = 2 * qJD(5);
t73 = t101 * t82 - t140 * t99;
t74 = t101 * t99 + t140 * t82;
t46 = pkin(4) * t73 - qJ(5) * t74;
t57 = qJDD(4) - t59;
t80 = qJD(4) - t81;
t79 = t80 ^ 2;
t115 = -pkin(4) * t79 + t57 * qJ(5) + t80 * t143 - t73 * t46 + t135;
t68 = -mrSges(6,1) * t80 + mrSges(6,2) * t74;
t148 = m(6) * t115 + t57 * mrSges(6,3) + t80 * t68;
t134 = -t102 * t51 + t141 * t50;
t113 = pkin(3) * t98 + pkin(9) * t97 - t82 * t71 + t134;
t138 = t73 * t80;
t33 = -t73 * qJD(4) + t101 * t98 + t140 * t60;
t147 = (-t33 + t138) * qJ(5) - t113;
t32 = qJD(4) * t74 + t101 * t60 - t140 * t98;
t48 = -mrSges(7,1) * t73 + mrSges(7,2) * t74;
t65 = -pkin(5) * t80 - qJ(6) * t74;
t72 = t73 ^ 2;
t126 = m(7) * (-pkin(5) * t72 + t32 * qJ(6) + 0.2e1 * qJD(6) * t73 + t65 * t80 + t115) + t32 * mrSges(7,3) + t73 * t48;
t47 = mrSges(6,1) * t73 - mrSges(6,3) * t74;
t132 = -mrSges(5,1) * t73 - mrSges(5,2) * t74 - t47;
t136 = -mrSges(5,3) - mrSges(6,2);
t66 = -mrSges(7,1) * t80 - mrSges(7,3) * t74;
t67 = mrSges(5,1) * t80 - mrSges(5,3) * t74;
t12 = m(5) * t135 + (-t67 + t66) * t80 + t132 * t73 + (-mrSges(5,2) + mrSges(7,2)) * t57 + t136 * t32 + t126 + t148;
t75 = -mrSges(4,2) * t99 + mrSges(4,3) * t81;
t76 = mrSges(4,1) * t99 - mrSges(4,3) * t82;
t144 = -0.2e1 * t74;
t120 = -t101 * t27 + t140 * t23;
t19 = -t57 * pkin(4) - t79 * qJ(5) + t74 * t46 + qJDD(5) - t120;
t63 = mrSges(7,2) * t80 + mrSges(7,3) * t73;
t125 = t80 * t63 + t57 * mrSges(7,1) - m(7) * (qJD(6) * t144 + (-t33 - t138) * qJ(6) + (t73 * t74 - t57) * pkin(5) + t19);
t118 = m(6) * t19 - t125;
t62 = -mrSges(6,2) * t73 + mrSges(6,3) * t80;
t64 = -mrSges(5,2) * t80 - mrSges(5,3) * t73;
t9 = m(5) * t120 + (t64 + t62) * t80 + (mrSges(5,1) + mrSges(6,1)) * t57 + (t48 + t132) * t74 + (mrSges(7,3) + t136) * t33 - t118;
t111 = -m(4) * t109 + t59 * mrSges(4,1) - t60 * mrSges(4,2) - t101 * t12 - t140 * t9 + t81 * t75 - t82 * t76;
t90 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t129;
t91 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t128;
t146 = (t103 * t90 - t105 * t91) * qJD(1) + m(3) * (-pkin(7) * t107 + t116) - t89 * mrSges(3,1) + t88 * mrSges(3,2) - t111;
t121 = m(7) * (-qJ(6) * t72 + qJDD(6) + (-pkin(4) - pkin(5)) * t32 + (-pkin(4) * t80 + t143 + t65) * t74 - t147) + t33 * mrSges(7,2) - t32 * mrSges(7,1) + t74 * t66 - t73 * t63;
t114 = m(6) * (qJD(5) * t144 + (t74 * t80 + t32) * pkin(4) + t147) + t32 * mrSges(6,1) + t73 * t62 - t121;
t145 = -m(5) * t113 + t32 * mrSges(5,1) + (t67 - t68) * t74 + (mrSges(5,2) - mrSges(6,3)) * t33 + t73 * t64 + t114;
t70 = -mrSges(4,1) * t81 + mrSges(4,2) * t82;
t10 = m(4) * t134 + t98 * mrSges(4,1) - t60 * mrSges(4,3) - t82 * t70 + t99 * t75 - t145;
t7 = m(4) * t133 - t98 * mrSges(4,2) + t59 * mrSges(4,3) - t101 * t9 + t140 * t12 + t81 * t70 - t99 * t76;
t87 = (-mrSges(3,1) * t105 + mrSges(3,2) * t103) * qJD(1);
t4 = m(3) * (-g(3) * t105 - t130) - t88 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t87 * t129 + qJD(2) * t91 + t102 * t7 + t141 * t10;
t5 = m(3) * t123 - qJDD(2) * mrSges(3,2) + t89 * mrSges(3,3) - qJD(2) * t90 - t102 * t10 + t87 * t128 + t141 * t7;
t142 = t103 * t5 + t105 * t4;
t112 = t57 * mrSges(7,2) + t80 * t66 + t126;
t6 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t107 * mrSges(2,2) - t146;
t1 = m(2) * t119 - t107 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t103 * t4 + t105 * t5;
t2 = [-m(1) * g(1) + t1 * t106 - t104 * t6, t1, t5, t7, t12, -t32 * mrSges(6,2) - t73 * t47 + t112 + t148, t112; -m(1) * g(2) + t1 * t104 + t106 * t6, t6, t4, t10, t9, -t33 * mrSges(6,3) - t74 * t68 + t114, -t33 * mrSges(7,3) - t74 * t48 - t125; (-m(1) - m(2)) * g(3) + t142, -m(2) * g(3) + t142, t146, -t111, t145, -t57 * mrSges(6,1) - t80 * t62 + (t47 - t48) * t74 + (mrSges(6,2) - mrSges(7,3)) * t33 + t118, t121;];
f_new  = t2;
