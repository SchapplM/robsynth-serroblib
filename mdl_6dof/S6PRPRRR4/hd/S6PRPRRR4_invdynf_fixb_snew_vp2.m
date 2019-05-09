% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 01:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:57:16
% EndTime: 2019-05-05 00:57:24
% DurationCPUTime: 4.45s
% Computational Cost: add. (59354->168), mult. (132240->224), div. (0->0), fcn. (101206->14), ass. (0->100)
t106 = qJD(2) ^ 2;
t94 = cos(pkin(12));
t89 = t94 ^ 2;
t91 = sin(pkin(12));
t127 = t91 ^ 2 + t89;
t135 = t127 * mrSges(4,3);
t100 = sin(qJ(2));
t104 = cos(qJ(2));
t92 = sin(pkin(11));
t95 = cos(pkin(11));
t80 = g(1) * t92 - g(2) * t95;
t96 = cos(pkin(6));
t132 = t80 * t96;
t81 = -g(1) * t95 - g(2) * t92;
t90 = -g(3) + qJDD(1);
t93 = sin(pkin(6));
t134 = -t100 * t81 + (t90 * t93 + t132) * t104;
t133 = pkin(3) * t106;
t99 = sin(qJ(4));
t131 = t91 * t99;
t102 = cos(qJ(5));
t105 = qJD(4) ^ 2;
t103 = cos(qJ(4));
t122 = qJD(2) * qJD(3);
t70 = -t80 * t93 + t90 * t96;
t128 = -0.2e1 * t122 * t91 + t70 * t94;
t125 = t100 * t93;
t120 = t100 * t132 + t104 * t81 + t125 * t90;
t54 = -pkin(2) * t106 + qJDD(2) * qJ(3) + t120;
t37 = (-pkin(8) * qJDD(2) + t133 * t94 - t54) * t91 + t128;
t121 = t91 * t70 + (0.2e1 * t122 + t54) * t94;
t123 = qJDD(2) * t94;
t38 = pkin(8) * t123 - t133 * t89 + t121;
t129 = t103 * t38 + t37 * t99;
t74 = (t103 * t94 - t131) * qJD(2);
t114 = t103 * t91 + t94 * t99;
t75 = t114 * qJD(2);
t60 = -pkin(4) * t74 - pkin(9) * t75;
t25 = -pkin(4) * t105 + qJDD(4) * pkin(9) + t60 * t74 + t129;
t111 = qJDD(3) - t134;
t107 = (-pkin(3) * t94 - pkin(2)) * qJDD(2) + (-pkin(8) * t127 - qJ(3)) * t106 + t111;
t124 = t74 * qJD(4);
t72 = t75 * qJD(4);
t61 = -qJDD(2) * t131 + t103 * t123 - t72;
t62 = qJDD(2) * t114 + t124;
t31 = (-t62 - t124) * pkin(9) + (-t61 + t72) * pkin(4) + t107;
t98 = sin(qJ(5));
t130 = t102 * t25 + t31 * t98;
t101 = cos(qJ(6));
t118 = t102 * t31 - t98 * t25;
t64 = qJD(4) * t102 - t75 * t98;
t42 = qJD(5) * t64 + qJDD(4) * t98 + t102 * t62;
t59 = qJDD(5) - t61;
t65 = qJD(4) * t98 + t102 * t75;
t73 = qJD(5) - t74;
t19 = (t64 * t73 - t42) * pkin(10) + (t64 * t65 + t59) * pkin(5) + t118;
t41 = -qJD(5) * t65 + qJDD(4) * t102 - t62 * t98;
t53 = pkin(5) * t73 - pkin(10) * t65;
t63 = t64 ^ 2;
t20 = -pkin(5) * t63 + pkin(10) * t41 - t53 * t73 + t130;
t97 = sin(qJ(6));
t45 = t101 * t64 - t65 * t97;
t28 = qJD(6) * t45 + t101 * t42 + t41 * t97;
t46 = t101 * t65 + t64 * t97;
t33 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t71 = qJD(6) + t73;
t39 = -mrSges(7,2) * t71 + mrSges(7,3) * t45;
t56 = qJDD(6) + t59;
t17 = m(7) * (t101 * t19 - t20 * t97) - t28 * mrSges(7,3) + t56 * mrSges(7,1) - t46 * t33 + t71 * t39;
t27 = -qJD(6) * t46 + t101 * t41 - t42 * t97;
t40 = mrSges(7,1) * t71 - mrSges(7,3) * t46;
t18 = m(7) * (t101 * t20 + t19 * t97) + t27 * mrSges(7,3) - t56 * mrSges(7,2) + t45 * t33 - t71 * t40;
t47 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t51 = -mrSges(6,2) * t73 + mrSges(6,3) * t64;
t14 = m(6) * t118 + mrSges(6,1) * t59 - mrSges(6,3) * t42 + t101 * t17 + t18 * t97 - t47 * t65 + t51 * t73;
t52 = mrSges(6,1) * t73 - mrSges(6,3) * t65;
t15 = m(6) * t130 - mrSges(6,2) * t59 + mrSges(6,3) * t41 + t101 * t18 - t17 * t97 + t47 * t64 - t52 * t73;
t68 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t74;
t69 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t75;
t110 = -m(5) * t107 + t61 * mrSges(5,1) - mrSges(5,2) * t62 - t102 * t14 - t15 * t98 + t74 * t68 - t69 * t75;
t109 = m(4) * (-qJDD(2) * pkin(2) - t106 * qJ(3) + t111) - t110;
t116 = -mrSges(4,1) * t94 + mrSges(4,2) * t91;
t10 = m(3) * t134 + (-mrSges(3,2) + t135) * t106 + (mrSges(3,1) - t116) * qJDD(2) - t109;
t126 = t10 * t104;
t57 = -mrSges(5,1) * t74 + mrSges(5,2) * t75;
t11 = m(5) * t129 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t61 - qJD(4) * t69 + t102 * t15 - t14 * t98 + t57 * t74;
t113 = mrSges(4,3) * qJDD(2) + t106 * t116;
t117 = t103 * t37 - t38 * t99;
t24 = -qJDD(4) * pkin(4) - pkin(9) * t105 + t60 * t75 - t117;
t112 = t27 * mrSges(7,1) + t45 * t39 - m(7) * (-pkin(5) * t41 - pkin(10) * t63 + t53 * t65 + t24) - t28 * mrSges(7,2) - t46 * t40;
t108 = m(6) * t24 - mrSges(6,1) * t41 + mrSges(6,2) * t42 - t51 * t64 + t52 * t65 - t112;
t16 = m(5) * t117 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t62 + qJD(4) * t68 - t57 * t75 - t108;
t7 = m(4) * t128 + t99 * t11 + t103 * t16 + (-m(4) * t54 - t113) * t91;
t8 = m(4) * t121 + t103 * t11 + t113 * t94 - t99 * t16;
t4 = m(3) * t120 - mrSges(3,1) * t106 - qJDD(2) * mrSges(3,2) - t7 * t91 + t8 * t94;
t6 = m(3) * t70 + t7 * t94 + t8 * t91;
t119 = m(2) * t90 + t125 * t4 + t126 * t93 + t6 * t96;
t2 = m(2) * t81 - t10 * t100 + t104 * t4;
t1 = m(2) * t80 - t93 * t6 + (t100 * t4 + t126) * t96;
t3 = [-m(1) * g(1) - t1 * t92 + t2 * t95, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t1 * t95 + t2 * t92, t1, t10, t7, t16, t14, t17; -m(1) * g(3) + t119, t119, t6, qJDD(2) * t116 - t106 * t135 + t109, -t110, t108, -t112;];
f_new  = t3;
