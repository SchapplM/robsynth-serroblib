% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:28:19
% EndTime: 2019-05-05 04:28:28
% DurationCPUTime: 4.81s
% Computational Cost: add. (63944->178), mult. (138380->244), div. (0->0), fcn. (100720->14), ass. (0->99)
t103 = sin(qJ(3));
t107 = cos(qJ(3));
t110 = qJD(2) ^ 2;
t104 = sin(qJ(2));
t108 = cos(qJ(2));
t97 = sin(pkin(6));
t130 = t104 * t97;
t100 = cos(pkin(6));
t96 = sin(pkin(11));
t99 = cos(pkin(11));
t85 = g(1) * t96 - g(2) * t99;
t131 = t100 * t85;
t86 = -g(1) * t99 - g(2) * t96;
t94 = -g(3) + qJDD(1);
t124 = t104 * t131 + t108 * t86 + t130 * t94;
t54 = -pkin(2) * t110 + qJDD(2) * pkin(8) + t124;
t70 = t100 * t94 - t85 * t97;
t119 = -t103 * t54 + t107 * t70;
t126 = qJD(2) * qJD(3);
t121 = t107 * t126;
t83 = qJDD(2) * t103 + t121;
t37 = (-t83 + t121) * qJ(4) + (t103 * t107 * t110 + qJDD(3)) * pkin(3) + t119;
t133 = t103 * t70 + t107 * t54;
t84 = qJDD(2) * t107 - t103 * t126;
t128 = qJD(2) * t103;
t87 = qJD(3) * pkin(3) - qJ(4) * t128;
t93 = t107 ^ 2;
t38 = -pkin(3) * t110 * t93 + qJ(4) * t84 - qJD(3) * t87 + t133;
t95 = sin(pkin(12));
t98 = cos(pkin(12));
t76 = (t103 * t98 + t107 * t95) * qJD(2);
t137 = -0.2e1 * qJD(4) * t76 + t37 * t98 - t38 * t95;
t136 = -t104 * t86 + (t94 * t97 + t131) * t108;
t114 = -qJDD(2) * pkin(2) - t136;
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t111 = -t84 * pkin(3) + qJDD(4) + t87 * t128 + (-qJ(4) * t93 - pkin(8)) * t110 + t114;
t101 = sin(qJ(6));
t105 = cos(qJ(6));
t109 = qJD(3) ^ 2;
t127 = qJD(2) * t107;
t75 = t127 * t98 - t128 * t95;
t125 = 0.2e1 * qJD(4) * t75 + t37 * t95 + t38 * t98;
t57 = -pkin(4) * t75 - pkin(9) * t76;
t25 = -pkin(4) * t109 + qJDD(3) * pkin(9) + t57 * t75 + t125;
t61 = -t83 * t95 + t84 * t98;
t62 = t83 * t98 + t84 * t95;
t31 = (-qJD(3) * t75 - t62) * pkin(9) + (qJD(3) * t76 - t61) * pkin(4) + t111;
t120 = -t102 * t25 + t106 * t31;
t64 = qJD(3) * t106 - t102 * t76;
t44 = qJD(5) * t64 + qJDD(3) * t102 + t106 * t62;
t60 = qJDD(5) - t61;
t65 = qJD(3) * t102 + t106 * t76;
t74 = qJD(5) - t75;
t19 = (t64 * t74 - t44) * pkin(10) + (t64 * t65 + t60) * pkin(5) + t120;
t134 = t102 * t31 + t106 * t25;
t43 = -qJD(5) * t65 + qJDD(3) * t106 - t102 * t62;
t52 = pkin(5) * t74 - pkin(10) * t65;
t63 = t64 ^ 2;
t20 = -pkin(5) * t63 + pkin(10) * t43 - t52 * t74 + t134;
t45 = -t101 * t65 + t105 * t64;
t28 = qJD(6) * t45 + t101 * t43 + t105 * t44;
t46 = t101 * t64 + t105 * t65;
t33 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t71 = qJD(6) + t74;
t41 = -mrSges(7,2) * t71 + mrSges(7,3) * t45;
t58 = qJDD(6) + t60;
t17 = m(7) * (-t101 * t20 + t105 * t19) - t28 * mrSges(7,3) + t58 * mrSges(7,1) - t46 * t33 + t71 * t41;
t27 = -qJD(6) * t46 - t101 * t44 + t105 * t43;
t42 = mrSges(7,1) * t71 - mrSges(7,3) * t46;
t18 = m(7) * (t101 * t19 + t105 * t20) + t27 * mrSges(7,3) - t58 * mrSges(7,2) + t45 * t33 - t71 * t42;
t47 = -mrSges(6,1) * t64 + mrSges(6,2) * t65;
t50 = -mrSges(6,2) * t74 + mrSges(6,3) * t64;
t14 = m(6) * t120 + mrSges(6,1) * t60 - mrSges(6,3) * t44 + t101 * t18 + t105 * t17 - t47 * t65 + t50 * t74;
t51 = mrSges(6,1) * t74 - mrSges(6,3) * t65;
t15 = m(6) * t134 - mrSges(6,2) * t60 + mrSges(6,3) * t43 - t101 * t17 + t105 * t18 + t47 * t64 - t51 * t74;
t68 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t75;
t69 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t76;
t115 = -m(5) * t111 + t61 * mrSges(5,1) - mrSges(5,2) * t62 - t102 * t15 - t106 * t14 + t75 * t68 - t69 * t76;
t88 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t128;
t89 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t127;
t135 = (t103 * t88 - t107 * t89) * qJD(2) + m(4) * (-t110 * pkin(8) + t114) - t84 * mrSges(4,1) + t83 * mrSges(4,2) - t115;
t10 = m(3) * t136 + qJDD(2) * mrSges(3,1) - t110 * mrSges(3,2) - t135;
t132 = t10 * t108;
t56 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t11 = m(5) * t125 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t61 - qJD(3) * t69 - t102 * t14 + t106 * t15 + t56 * t75;
t24 = -qJDD(3) * pkin(4) - pkin(9) * t109 + t57 * t76 - t137;
t116 = t27 * mrSges(7,1) + t45 * t41 - m(7) * (-pkin(5) * t43 - pkin(10) * t63 + t52 * t65 + t24) - t28 * mrSges(7,2) - t46 * t42;
t112 = m(6) * t24 - t43 * mrSges(6,1) + t44 * mrSges(6,2) - t64 * t50 + t65 * t51 - t116;
t16 = m(5) * t137 + qJDD(3) * mrSges(5,1) - t62 * mrSges(5,3) + qJD(3) * t68 - t76 * t56 - t112;
t82 = (-mrSges(4,1) * t107 + mrSges(4,2) * t103) * qJD(2);
t7 = m(4) * t119 + qJDD(3) * mrSges(4,1) - t83 * mrSges(4,3) + qJD(3) * t89 + t95 * t11 - t128 * t82 + t98 * t16;
t8 = m(4) * t133 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t84 - qJD(3) * t88 + t11 * t98 + t127 * t82 - t16 * t95;
t4 = m(3) * t124 - mrSges(3,1) * t110 - qJDD(2) * mrSges(3,2) - t103 * t7 + t107 * t8;
t6 = m(3) * t70 + t103 * t8 + t107 * t7;
t123 = m(2) * t94 + t100 * t6 + t130 * t4 + t132 * t97;
t2 = m(2) * t86 - t10 * t104 + t108 * t4;
t1 = m(2) * t85 - t97 * t6 + (t104 * t4 + t132) * t100;
t3 = [-m(1) * g(1) - t1 * t96 + t2 * t99, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t1 * t99 + t2 * t96, t1, t10, t7, t16, t14, t17; -m(1) * g(3) + t123, t123, t6, t135, -t115, t112, -t116;];
f_new  = t3;
