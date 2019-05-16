% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP8
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
% Datum: 2019-05-07 19:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:57:21
% EndTime: 2019-05-07 18:57:28
% DurationCPUTime: 2.97s
% Computational Cost: add. (37268->212), mult. (79191->265), div. (0->0), fcn. (61064->10), ass. (0->103)
t103 = sin(pkin(6));
t107 = sin(qJ(2));
t109 = cos(qJ(2));
t130 = qJD(1) * qJD(2);
t92 = (-qJDD(1) * t109 + t107 * t130) * t103;
t105 = sin(qJ(4));
t148 = cos(qJ(4));
t106 = sin(qJ(3));
t149 = cos(qJ(3));
t132 = qJD(1) * t109;
t104 = cos(pkin(6));
t111 = qJD(1) ^ 2;
t108 = sin(qJ(1));
t110 = cos(qJ(1));
t124 = t108 * g(1) - t110 * g(2);
t147 = pkin(8) * t103;
t87 = qJDD(1) * pkin(1) + t111 * t147 + t124;
t136 = t104 * t87;
t121 = -t110 * g(1) - t108 * g(2);
t88 = -t111 * pkin(1) + qJDD(1) * t147 + t121;
t137 = t107 * t136 + t109 * t88;
t133 = qJD(1) * t103;
t90 = (-pkin(2) * t109 - pkin(9) * t107) * t133;
t100 = t104 * qJD(1) + qJD(2);
t98 = t100 ^ 2;
t99 = t104 * qJDD(1) + qJDD(2);
t45 = -t98 * pkin(2) + t99 * pkin(9) + (-g(3) * t107 + t90 * t132) * t103 + t137;
t146 = t104 * g(3);
t91 = (qJDD(1) * t107 + t109 * t130) * t103;
t46 = t92 * pkin(2) - t91 * pkin(9) - t146 + (-t87 + (pkin(2) * t107 - pkin(9) * t109) * t100 * qJD(1)) * t103;
t140 = t106 * t46 + t149 * t45;
t127 = t107 * t133;
t80 = t149 * t100 - t106 * t127;
t81 = t106 * t100 + t149 * t127;
t68 = -t80 * pkin(3) - t81 * pkin(10);
t84 = qJDD(3) + t92;
t126 = t103 * t132;
t97 = qJD(3) - t126;
t96 = t97 ^ 2;
t25 = -t96 * pkin(3) + t84 * pkin(10) + t80 * t68 + t140;
t134 = t103 * t109;
t119 = -g(3) * t134 - t107 * t88 + t109 * t136;
t44 = -t99 * pkin(2) - t98 * pkin(9) + t90 * t127 - t119;
t65 = -t81 * qJD(3) - t106 * t91 + t149 * t99;
t66 = t80 * qJD(3) + t106 * t99 + t149 * t91;
t27 = (-t80 * t97 - t66) * pkin(10) + (t81 * t97 - t65) * pkin(3) + t44;
t142 = t105 * t27 + t148 * t25;
t150 = 2 * qJD(5);
t70 = t105 * t81 - t148 * t97;
t71 = t105 * t97 + t148 * t81;
t49 = t70 * pkin(4) - t71 * qJ(5);
t64 = qJDD(4) - t65;
t78 = qJD(4) - t80;
t77 = t78 ^ 2;
t118 = -t77 * pkin(4) + t64 * qJ(5) + t78 * t150 - t70 * t49 + t142;
t60 = -t78 * mrSges(6,1) + t71 * mrSges(6,2);
t154 = m(6) * t118 + t64 * mrSges(6,3) + t78 * t60;
t141 = -t106 * t45 + t149 * t46;
t115 = t84 * pkin(3) + t96 * pkin(10) - t81 * t68 + t141;
t145 = t70 * t78;
t33 = -t70 * qJD(4) + t105 * t84 + t148 * t66;
t153 = (-t33 + t145) * qJ(5) - t115;
t32 = t71 * qJD(4) + t105 * t66 - t148 * t84;
t55 = t78 * mrSges(7,2) + t70 * mrSges(7,3);
t57 = -t78 * pkin(5) - t71 * qJ(6);
t58 = -t78 * mrSges(7,1) - t71 * mrSges(7,3);
t69 = t70 ^ 2;
t123 = m(7) * (-t69 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t32 + (-pkin(4) * t78 + t150 + t57) * t71 - t153) + t33 * mrSges(7,2) - t32 * mrSges(7,1) + t71 * t58 - t70 * t55;
t151 = -0.2e1 * t71;
t54 = -t70 * mrSges(6,2) + t78 * mrSges(6,3);
t116 = m(6) * (qJD(5) * t151 + (t71 * t78 + t32) * pkin(4) + t153) + t32 * mrSges(6,1) + t70 * t54 - t123;
t56 = -t78 * mrSges(5,2) - t70 * mrSges(5,3);
t59 = t78 * mrSges(5,1) - t71 * mrSges(5,3);
t152 = -m(5) * t115 + t32 * mrSges(5,1) + (t59 - t60) * t71 + (mrSges(5,2) - mrSges(6,3)) * t33 + t70 * t56 + t116;
t143 = -mrSges(5,3) - mrSges(6,2);
t50 = t70 * mrSges(6,1) - t71 * mrSges(6,3);
t139 = -t70 * mrSges(5,1) - t71 * mrSges(5,2) - t50;
t135 = t103 * t107;
t67 = -t80 * mrSges(4,1) + t81 * mrSges(4,2);
t72 = -t97 * mrSges(4,2) + t80 * mrSges(4,3);
t10 = m(4) * t141 + t84 * mrSges(4,1) - t66 * mrSges(4,3) - t81 * t67 + t97 * t72 - t152;
t85 = t100 * mrSges(3,1) - mrSges(3,3) * t127;
t89 = (-mrSges(3,1) * t109 + mrSges(3,2) * t107) * t133;
t122 = -t105 * t25 + t148 * t27;
t19 = -t64 * pkin(4) - t77 * qJ(5) + t71 * t49 + qJDD(5) - t122;
t128 = t78 * t55 + t64 * mrSges(7,1) - m(7) * (qJD(6) * t151 + (-t33 - t145) * qJ(6) + (t70 * t71 - t64) * pkin(5) + t19);
t120 = m(6) * t19 - t128;
t51 = -t70 * mrSges(7,1) + t71 * mrSges(7,2);
t11 = m(5) * t122 + (t56 + t54) * t78 + (mrSges(5,1) + mrSges(6,1)) * t64 + (t51 + t139) * t71 + (mrSges(7,3) + t143) * t33 - t120;
t129 = m(7) * (-t69 * pkin(5) + t32 * qJ(6) + 0.2e1 * qJD(6) * t70 + t78 * t57 + t118) + t32 * mrSges(7,3) + t70 * t51;
t12 = m(5) * t142 + (-t59 + t58) * t78 + t139 * t70 + (-mrSges(5,2) + mrSges(7,2)) * t64 + t143 * t32 + t129 + t154;
t73 = t97 * mrSges(4,1) - t81 * mrSges(4,3);
t9 = m(4) * t140 - t84 * mrSges(4,2) + t65 * mrSges(4,3) - t105 * t11 + t148 * t12 + t80 * t67 - t97 * t73;
t4 = m(3) * (-g(3) * t135 + t137) - t92 * mrSges(3,3) - t99 * mrSges(3,2) + t89 * t126 - t100 * t85 + t149 * t9 - t106 * t10;
t86 = -t100 * mrSges(3,2) + mrSges(3,3) * t126;
t6 = m(3) * (-t103 * t87 - t146) + t91 * mrSges(3,2) + t92 * mrSges(3,1) + t106 * t9 + t149 * t10 + (t107 * t85 - t109 * t86) * t133;
t112 = m(4) * t44 - t65 * mrSges(4,1) + t66 * mrSges(4,2) + t105 * t12 + t148 * t11 - t80 * t72 + t81 * t73;
t8 = m(3) * t119 + t99 * mrSges(3,1) - t91 * mrSges(3,3) + t100 * t86 - t89 * t127 - t112;
t131 = t104 * t6 + t8 * t134 + t4 * t135;
t114 = t64 * mrSges(7,2) + t78 * t58 + t129;
t2 = m(2) * t121 - t111 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t107 * t8 + t109 * t4;
t1 = m(2) * t124 + qJDD(1) * mrSges(2,1) - t111 * mrSges(2,2) - t103 * t6 + (t107 * t4 + t109 * t8) * t104;
t3 = [-m(1) * g(1) - t108 * t1 + t110 * t2, t2, t4, t9, t12, -t32 * mrSges(6,2) - t70 * t50 + t114 + t154, t114; -m(1) * g(2) + t110 * t1 + t108 * t2, t1, t8, t10, t11, -t33 * mrSges(6,3) - t71 * t60 + t116, -t33 * mrSges(7,3) - t71 * t51 - t128; (-m(1) - m(2)) * g(3) + t131, -m(2) * g(3) + t131, t6, t112, t152, -t64 * mrSges(6,1) - t78 * t54 + (t50 - t51) * t71 + (mrSges(6,2) - mrSges(7,3)) * t33 + t120, t123;];
f_new  = t3;
