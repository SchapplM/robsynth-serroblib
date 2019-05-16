% Calculate vector of cutting forces with Newton-Euler
% S6RRRRPP3
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
% Datum: 2019-05-07 18:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RRRRPP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:07:50
% EndTime: 2019-05-07 18:07:54
% DurationCPUTime: 1.78s
% Computational Cost: add. (20915->208), mult. (41986->244), div. (0->0), fcn. (28824->8), ass. (0->93)
t137 = cos(qJ(4));
t100 = cos(qJ(2));
t96 = sin(qJ(3));
t97 = sin(qJ(2));
t99 = cos(qJ(3));
t77 = (t100 * t96 + t97 * t99) * qJD(1);
t93 = qJD(2) + qJD(3);
t95 = sin(qJ(4));
t69 = -t137 * t93 + t77 * t95;
t123 = qJD(1) * t100;
t124 = qJD(1) * t97;
t76 = t99 * t123 - t96 * t124;
t75 = qJD(4) - t76;
t135 = t69 * t75;
t139 = -2 * qJD(5);
t122 = qJD(1) * qJD(2);
t102 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t98 = sin(qJ(1));
t114 = -g(1) * t101 - g(2) * t98;
t79 = -pkin(1) * t102 + qJDD(1) * pkin(7) + t114;
t134 = t79 * t97;
t136 = pkin(2) * t102;
t82 = qJDD(1) * t97 + t100 * t122;
t46 = qJDD(2) * pkin(2) - pkin(8) * t82 - t134 + (pkin(8) * t122 + t97 * t136 - g(3)) * t100;
t117 = -g(3) * t97 + t100 * t79;
t83 = qJDD(1) * t100 - t97 * t122;
t86 = qJD(2) * pkin(2) - pkin(8) * t124;
t94 = t100 ^ 2;
t47 = pkin(8) * t83 - qJD(2) * t86 - t94 * t136 + t117;
t116 = t46 * t99 - t96 * t47;
t66 = -pkin(3) * t76 - pkin(9) * t77;
t91 = t93 ^ 2;
t92 = qJDD(2) + qJDD(3);
t27 = -pkin(3) * t92 - pkin(9) * t91 + t77 * t66 - t116;
t54 = qJD(3) * t76 + t82 * t99 + t83 * t96;
t33 = -t69 * qJD(4) + t137 * t54 + t95 * t92;
t70 = t137 * t77 + t95 * t93;
t104 = (-t33 + t135) * qJ(5) + t27 + (t75 * pkin(4) + t139) * t70;
t32 = qJD(4) * t70 - t137 * t92 + t54 * t95;
t60 = mrSges(6,1) * t70 + mrSges(6,2) * t75;
t142 = m(6) * (t32 * pkin(4) + t104) - t33 * mrSges(6,3) - t70 * t60;
t118 = g(1) * t98 - t101 * g(2);
t110 = -qJDD(1) * pkin(1) - t118;
t105 = -pkin(2) * t83 + t86 * t124 + (-pkin(8) * t94 - pkin(7)) * t102 + t110;
t53 = -qJD(3) * t77 - t82 * t96 + t83 * t99;
t24 = (-t76 * t93 - t54) * pkin(9) + (t77 * t93 - t53) * pkin(3) + t105;
t128 = t96 * t46 + t99 * t47;
t28 = -pkin(3) * t91 + pkin(9) * t92 + t66 * t76 + t128;
t115 = t137 * t24 - t95 * t28;
t43 = pkin(4) * t69 - qJ(5) * t70;
t51 = qJDD(4) - t53;
t74 = t75 ^ 2;
t20 = -t51 * pkin(4) - t74 * qJ(5) + t70 * t43 + qJDD(5) - t115;
t42 = -mrSges(7,2) * t70 + t69 * mrSges(7,3);
t121 = m(7) * (-0.2e1 * qJD(6) * t75 + (t69 * t70 - t51) * qJ(6) + (t33 + t135) * pkin(5) + t20) + t33 * mrSges(7,1) + t70 * t42;
t112 = m(6) * t20 + t121;
t58 = t69 * mrSges(6,1) - mrSges(6,3) * t75;
t126 = mrSges(5,2) * t75 + mrSges(5,3) * t69 + t58;
t45 = -t69 * mrSges(6,2) - mrSges(6,3) * t70;
t127 = -mrSges(5,1) * t69 - mrSges(5,2) * t70 - t45;
t130 = -mrSges(5,3) - mrSges(6,1);
t131 = mrSges(6,2) - mrSges(7,3);
t59 = -t69 * mrSges(7,1) + mrSges(7,2) * t75;
t10 = m(5) * t115 + t127 * t70 + t130 * t33 + (t59 - t126) * t75 + (mrSges(5,1) - t131) * t51 - t112;
t129 = t137 * t28 + t95 * t24;
t107 = -pkin(4) * t74 + t51 * qJ(5) - t43 * t69 + t129;
t56 = pkin(5) * t70 - qJ(6) * t75;
t57 = mrSges(7,1) * t70 - mrSges(7,3) * t75;
t68 = t69 ^ 2;
t120 = m(7) * (-t32 * pkin(5) - t68 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t56) * t75 + t107) + t75 * t57 + t51 * mrSges(7,2);
t111 = m(6) * (t75 * t139 - t107) - t120;
t62 = mrSges(5,1) * t75 - mrSges(5,3) * t70;
t12 = m(5) * t129 + (-t62 + t60) * t75 + (-mrSges(5,2) + mrSges(6,3)) * t51 + (-t42 + t127) * t69 + (-mrSges(7,1) + t130) * t32 - t111;
t71 = -mrSges(4,2) * t93 + mrSges(4,3) * t76;
t72 = mrSges(4,1) * t93 - mrSges(4,3) * t77;
t106 = -m(4) * t105 + t53 * mrSges(4,1) - t54 * mrSges(4,2) - t137 * t10 - t95 * t12 + t76 * t71 - t77 * t72;
t84 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t124;
t85 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t123;
t141 = -(t100 * t85 - t84 * t97) * qJD(1) + m(3) * (-pkin(7) * t102 + t110) - t83 * mrSges(3,1) + t82 * mrSges(3,2) - t106;
t119 = m(7) * (-t68 * pkin(5) + 0.2e1 * qJD(6) * t69 - t56 * t70 + (pkin(4) + qJ(6)) * t32 + t104) + t32 * mrSges(7,3) + t69 * t59;
t140 = m(5) * t27 - (-t62 + t57) * t70 - t126 * t69 + (mrSges(5,2) - mrSges(7,2)) * t33 + (mrSges(5,1) - mrSges(6,2)) * t32 + t119 + t142;
t65 = -mrSges(4,1) * t76 + mrSges(4,2) * t77;
t7 = m(4) * t128 - t92 * mrSges(4,2) + t53 * mrSges(4,3) - t95 * t10 + t137 * t12 + t76 * t65 - t93 * t72;
t8 = m(4) * t116 + t92 * mrSges(4,1) - t54 * mrSges(4,3) - t77 * t65 + t93 * t71 - t140;
t81 = (-mrSges(3,1) * t100 + mrSges(3,2) * t97) * qJD(1);
t4 = m(3) * (-g(3) * t100 - t134) - t82 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t81 * t124 + qJD(2) * t85 + t96 * t7 + t99 * t8;
t5 = m(3) * t117 - qJDD(2) * mrSges(3,2) + t83 * mrSges(3,3) - qJD(2) * t84 + t81 * t123 + t99 * t7 - t96 * t8;
t138 = t100 * t4 + t97 * t5;
t109 = -t33 * mrSges(7,2) - t70 * t57 + t119;
t6 = m(2) * t118 + qJDD(1) * mrSges(2,1) - t102 * mrSges(2,2) - t141;
t1 = m(2) * t114 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t100 * t5 - t97 * t4;
t2 = [-m(1) * g(1) + t1 * t101 - t6 * t98, t1, t5, t7, t12, -t32 * mrSges(6,2) - t69 * t58 + t109 + t142, t109; -m(1) * g(2) + t1 * t98 + t101 * t6, t6, t4, t8, t10, -t51 * mrSges(6,3) - t75 * t60 + (t42 + t45) * t69 + (mrSges(6,1) + mrSges(7,1)) * t32 + t111, -t51 * mrSges(7,3) - t75 * t59 + t121; (-m(1) - m(2)) * g(3) + t138, -m(2) * g(3) + t138, t141, -t106, t140, t33 * mrSges(6,1) + t70 * t45 + (t58 - t59) * t75 + t131 * t51 + t112, -t32 * mrSges(7,1) - t69 * t42 + t120;];
f_new  = t2;
