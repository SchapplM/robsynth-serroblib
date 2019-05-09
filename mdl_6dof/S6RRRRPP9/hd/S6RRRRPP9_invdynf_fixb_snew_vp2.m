% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-05-07 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:11:56
% EndTime: 2019-05-07 19:12:04
% DurationCPUTime: 2.94s
% Computational Cost: add. (37199->215), mult. (78965->265), div. (0->0), fcn. (60858->10), ass. (0->103)
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t125 = qJD(1) * qJD(2);
t97 = sin(pkin(6));
t87 = (-qJDD(1) * t104 + t101 * t125) * t97;
t143 = cos(qJ(4));
t100 = sin(qJ(3));
t103 = cos(qJ(3));
t128 = qJD(1) * t97;
t121 = t101 * t128;
t98 = cos(pkin(6));
t94 = t98 * qJD(1) + qJD(2);
t76 = t100 * t94 + t103 * t121;
t127 = qJD(1) * t104;
t120 = t97 * t127;
t91 = qJD(3) - t120;
t99 = sin(qJ(4));
t66 = -t143 * t91 + t99 * t76;
t75 = -t100 * t121 + t103 * t94;
t73 = qJD(4) - t75;
t142 = t66 * t73;
t146 = -2 * qJD(5);
t106 = qJD(1) ^ 2;
t102 = sin(qJ(1));
t105 = cos(qJ(1));
t119 = t102 * g(1) - t105 * g(2);
t145 = pkin(8) * t97;
t82 = qJDD(1) * pkin(1) + t106 * t145 + t119;
t141 = t82 * t98;
t116 = -t105 * g(1) - t102 * g(2);
t83 = -t106 * pkin(1) + qJDD(1) * t145 + t116;
t131 = t101 * t141 + t104 * t83;
t85 = (-pkin(2) * t104 - pkin(9) * t101) * t128;
t92 = t94 ^ 2;
t93 = t98 * qJDD(1) + qJDD(2);
t42 = -t92 * pkin(2) + t93 * pkin(9) + (-g(3) * t101 + t127 * t85) * t97 + t131;
t144 = t98 * g(3);
t86 = (qJDD(1) * t101 + t104 * t125) * t97;
t43 = t87 * pkin(2) - t86 * pkin(9) - t144 + (-t82 + (pkin(2) * t101 - pkin(9) * t104) * t94 * qJD(1)) * t97;
t118 = -t100 * t42 + t103 * t43;
t63 = -t75 * pkin(3) - t76 * pkin(10);
t79 = qJDD(3) + t87;
t90 = t91 ^ 2;
t25 = -t79 * pkin(3) - t90 * pkin(10) + t76 * t63 - t118;
t61 = t75 * qJD(3) + t100 * t93 + t103 * t86;
t33 = -t66 * qJD(4) + t143 * t61 + t99 * t79;
t67 = t143 * t76 + t99 * t91;
t108 = (-t33 + t142) * qJ(5) + t25 + (t73 * pkin(4) + t146) * t67;
t32 = t67 * qJD(4) - t143 * t79 + t99 * t61;
t54 = t67 * mrSges(6,1) + t73 * mrSges(6,2);
t148 = m(6) * (t32 * pkin(4) + t108) - t33 * mrSges(6,3) - t67 * t54;
t50 = t67 * pkin(5) - t73 * qJ(6);
t53 = -t66 * mrSges(7,1) + t73 * mrSges(7,2);
t65 = t66 ^ 2;
t122 = m(7) * (-t65 * pkin(5) + 0.2e1 * qJD(6) * t66 - t67 * t50 + (pkin(4) + qJ(6)) * t32 + t108) + t32 * mrSges(7,3) + t66 * t53;
t52 = t66 * mrSges(6,1) - t73 * mrSges(6,3);
t133 = t73 * mrSges(5,2) + t66 * mrSges(5,3) + t52;
t51 = t67 * mrSges(7,1) - t73 * mrSges(7,3);
t56 = t73 * mrSges(5,1) - t67 * mrSges(5,3);
t147 = m(5) * t25 - (-t56 + t51) * t67 - t133 * t66 + (mrSges(5,2) - mrSges(7,2)) * t33 + (mrSges(5,1) - mrSges(6,2)) * t32 + t122 + t148;
t138 = mrSges(6,2) - mrSges(7,3);
t137 = -mrSges(5,3) - mrSges(6,1);
t135 = t100 * t43 + t103 * t42;
t26 = -t90 * pkin(3) + t79 * pkin(10) + t75 * t63 + t135;
t129 = t104 * t97;
t113 = -g(3) * t129 - t101 * t83 + t104 * t141;
t41 = -t93 * pkin(2) - t92 * pkin(9) + t85 * t121 - t113;
t60 = -t76 * qJD(3) - t100 * t86 + t103 * t93;
t28 = (-t75 * t91 - t61) * pkin(10) + (t76 * t91 - t60) * pkin(3) + t41;
t136 = t143 * t26 + t99 * t28;
t48 = -t66 * mrSges(6,2) - t67 * mrSges(6,3);
t134 = -t66 * mrSges(5,1) - t67 * mrSges(5,2) - t48;
t130 = t101 * t97;
t62 = -t75 * mrSges(4,1) + t76 * mrSges(4,2);
t68 = -t91 * mrSges(4,2) + t75 * mrSges(4,3);
t10 = m(4) * t118 + t79 * mrSges(4,1) - t61 * mrSges(4,3) - t76 * t62 + t91 * t68 - t147;
t80 = t94 * mrSges(3,1) - mrSges(3,3) * t121;
t84 = (-mrSges(3,1) * t104 + mrSges(3,2) * t101) * t128;
t117 = t143 * t28 - t99 * t26;
t46 = t66 * pkin(4) - t67 * qJ(5);
t59 = qJDD(4) - t60;
t72 = t73 ^ 2;
t20 = -t59 * pkin(4) - t72 * qJ(5) + t67 * t46 + qJDD(5) - t117;
t45 = -t67 * mrSges(7,2) + t66 * mrSges(7,3);
t124 = m(7) * (-0.2e1 * qJD(6) * t73 + (t66 * t67 - t59) * qJ(6) + (t33 + t142) * pkin(5) + t20) + t33 * mrSges(7,1) + t67 * t45;
t115 = m(6) * t20 + t124;
t11 = m(5) * t117 + t134 * t67 + t137 * t33 + (t53 - t133) * t73 + (mrSges(5,1) - t138) * t59 - t115;
t109 = -t72 * pkin(4) + t59 * qJ(5) - t66 * t46 + t136;
t123 = m(7) * (-t32 * pkin(5) - t65 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t50) * t73 + t109) + t73 * t51 + t59 * mrSges(7,2);
t114 = m(6) * (t73 * t146 - t109) - t123;
t12 = m(5) * t136 + (-t56 + t54) * t73 + (-mrSges(5,2) + mrSges(6,3)) * t59 + (-t45 + t134) * t66 + (-mrSges(7,1) + t137) * t32 - t114;
t69 = t91 * mrSges(4,1) - t76 * mrSges(4,3);
t9 = m(4) * t135 - t79 * mrSges(4,2) + t60 * mrSges(4,3) - t99 * t11 + t12 * t143 + t75 * t62 - t91 * t69;
t4 = m(3) * (-g(3) * t130 + t131) - t87 * mrSges(3,3) - t93 * mrSges(3,2) + t84 * t120 - t94 * t80 + t103 * t9 - t100 * t10;
t81 = -t94 * mrSges(3,2) + mrSges(3,3) * t120;
t6 = m(3) * (-t97 * t82 - t144) + t86 * mrSges(3,2) + t87 * mrSges(3,1) + t100 * t9 + t103 * t10 + (t101 * t80 - t104 * t81) * t128;
t107 = m(4) * t41 - t60 * mrSges(4,1) + t61 * mrSges(4,2) + t143 * t11 + t99 * t12 - t75 * t68 + t76 * t69;
t8 = m(3) * t113 + t93 * mrSges(3,1) - t86 * mrSges(3,3) - t121 * t84 + t94 * t81 - t107;
t126 = t8 * t129 + t4 * t130 + t98 * t6;
t111 = -t33 * mrSges(7,2) - t67 * t51 + t122;
t2 = m(2) * t116 - t106 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t101 * t8 + t104 * t4;
t1 = m(2) * t119 + qJDD(1) * mrSges(2,1) - t106 * mrSges(2,2) - t97 * t6 + (t101 * t4 + t104 * t8) * t98;
t3 = [-m(1) * g(1) - t102 * t1 + t105 * t2, t2, t4, t9, t12, -t32 * mrSges(6,2) - t66 * t52 + t111 + t148, t111; -m(1) * g(2) + t105 * t1 + t102 * t2, t1, t8, t10, t11, -t59 * mrSges(6,3) - t73 * t54 + (t45 + t48) * t66 + (mrSges(6,1) + mrSges(7,1)) * t32 + t114, -t59 * mrSges(7,3) - t73 * t53 + t124; (-m(1) - m(2)) * g(3) + t126, -m(2) * g(3) + t126, t6, t107, t147, t33 * mrSges(6,1) + t67 * t48 + (t52 - t53) * t73 + t138 * t59 + t115, -t32 * mrSges(7,1) - t66 * t45 + t123;];
f_new  = t3;
