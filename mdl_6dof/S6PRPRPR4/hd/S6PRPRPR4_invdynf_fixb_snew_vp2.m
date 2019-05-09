% Calculate vector of cutting forces with Newton-Euler
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:45:38
% EndTime: 2019-05-04 22:45:47
% DurationCPUTime: 4.31s
% Computational Cost: add. (55885->165), mult. (128445->224), div. (0->0), fcn. (97564->14), ass. (0->98)
t134 = cos(qJ(4));
t91 = sin(pkin(11));
t95 = cos(pkin(11));
t99 = sin(qJ(4));
t137 = -t95 * t134 + t91 * t99;
t104 = qJD(2) ^ 2;
t88 = t95 ^ 2;
t128 = t91 ^ 2 + t88;
t136 = t128 * mrSges(4,3);
t100 = sin(qJ(2));
t102 = cos(qJ(2));
t92 = sin(pkin(10));
t96 = cos(pkin(10));
t79 = g(1) * t92 - g(2) * t96;
t97 = cos(pkin(6));
t132 = t79 * t97;
t80 = -g(1) * t96 - g(2) * t92;
t89 = -g(3) + qJDD(1);
t93 = sin(pkin(6));
t135 = -t100 * t80 + (t89 * t93 + t132) * t102;
t133 = pkin(3) * t104;
t125 = pkin(8) * qJDD(2);
t122 = qJD(2) * qJD(3);
t71 = -t79 * t93 + t89 * t97;
t129 = -0.2e1 * t91 * t122 + t95 * t71;
t126 = t100 * t93;
t119 = t100 * t132 + t102 * t80 + t89 * t126;
t49 = -pkin(2) * t104 + qJDD(2) * qJ(3) + t119;
t37 = (t95 * t133 - t125 - t49) * t91 + t129;
t120 = t91 * t71 + (0.2e1 * t122 + t49) * t95;
t38 = t95 * t125 - t88 * t133 + t120;
t130 = t134 * t38 + t99 * t37;
t109 = qJDD(3) - t135;
t105 = (-pkin(3) * t95 - pkin(2)) * qJDD(2) + (-t128 * pkin(8) - qJ(3)) * t104 + t109;
t101 = cos(qJ(6));
t103 = qJD(4) ^ 2;
t73 = t137 * qJD(2);
t111 = t134 * t91 + t95 * t99;
t74 = t111 * qJD(2);
t56 = pkin(4) * t73 - qJ(5) * t74;
t25 = -pkin(4) * t103 + qJDD(4) * qJ(5) - t56 * t73 + t130;
t123 = t73 * qJD(4);
t124 = qJD(4) * t74;
t60 = qJDD(2) * t137 + t124;
t61 = t111 * qJDD(2) - t123;
t28 = (-t61 + t123) * qJ(5) + (t60 + t124) * pkin(4) + t105;
t90 = sin(pkin(12));
t94 = cos(pkin(12));
t66 = qJD(4) * t90 + t74 * t94;
t116 = -0.2e1 * qJD(5) * t66 - t90 * t25 + t94 * t28;
t54 = qJDD(4) * t90 + t61 * t94;
t65 = qJD(4) * t94 - t74 * t90;
t19 = (t65 * t73 - t54) * pkin(9) + (t65 * t66 + t60) * pkin(5) + t116;
t121 = 0.2e1 * qJD(5) * t65 + t94 * t25 + t90 * t28;
t52 = pkin(5) * t73 - t66 * pkin(9);
t53 = qJDD(4) * t94 - t61 * t90;
t64 = t65 ^ 2;
t20 = -t64 * pkin(5) + t53 * pkin(9) - t52 * t73 + t121;
t98 = sin(qJ(6));
t43 = t101 * t65 - t66 * t98;
t31 = t43 * qJD(6) + t101 * t54 + t53 * t98;
t44 = t101 * t66 + t65 * t98;
t33 = -mrSges(7,1) * t43 + mrSges(7,2) * t44;
t72 = qJD(6) + t73;
t39 = -mrSges(7,2) * t72 + t43 * mrSges(7,3);
t59 = qJDD(6) + t60;
t17 = m(7) * (t101 * t19 - t20 * t98) - t31 * mrSges(7,3) + t59 * mrSges(7,1) - t44 * t33 + t72 * t39;
t30 = -t44 * qJD(6) + t101 * t53 - t54 * t98;
t40 = mrSges(7,1) * t72 - t44 * mrSges(7,3);
t18 = m(7) * (t101 * t20 + t19 * t98) + t30 * mrSges(7,3) - t59 * mrSges(7,2) + t43 * t33 - t72 * t40;
t45 = -mrSges(6,1) * t65 + mrSges(6,2) * t66;
t50 = -mrSges(6,2) * t73 + t65 * mrSges(6,3);
t14 = m(6) * t116 + t60 * mrSges(6,1) - t54 * mrSges(6,3) + t101 * t17 + t98 * t18 - t66 * t45 + t73 * t50;
t51 = mrSges(6,1) * t73 - t66 * mrSges(6,3);
t15 = m(6) * t121 - t60 * mrSges(6,2) + t53 * mrSges(6,3) + t101 * t18 - t98 * t17 + t65 * t45 - t73 * t51;
t69 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t73;
t70 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t74;
t108 = m(5) * t105 + t60 * mrSges(5,1) + t61 * mrSges(5,2) + t94 * t14 + t90 * t15 + t73 * t69 + t74 * t70;
t107 = m(4) * (-qJDD(2) * pkin(2) - t104 * qJ(3) + t109) + t108;
t114 = -mrSges(4,1) * t95 + mrSges(4,2) * t91;
t10 = m(3) * t135 + (-mrSges(3,2) + t136) * t104 + (mrSges(3,1) - t114) * qJDD(2) - t107;
t127 = t10 * t102;
t57 = mrSges(5,1) * t73 + mrSges(5,2) * t74;
t11 = m(5) * t130 - qJDD(4) * mrSges(5,2) - t60 * mrSges(5,3) - qJD(4) * t70 - t90 * t14 + t94 * t15 - t73 * t57;
t112 = mrSges(4,3) * qJDD(2) + t104 * t114;
t115 = t134 * t37 - t99 * t38;
t24 = -qJDD(4) * pkin(4) - t103 * qJ(5) + t74 * t56 + qJDD(5) - t115;
t110 = t30 * mrSges(7,1) + t43 * t39 - m(7) * (-t53 * pkin(5) - t64 * pkin(9) + t66 * t52 + t24) - t31 * mrSges(7,2) - t44 * t40;
t106 = m(6) * t24 - t53 * mrSges(6,1) + t54 * mrSges(6,2) - t65 * t50 + t66 * t51 - t110;
t16 = m(5) * t115 + qJDD(4) * mrSges(5,1) - t61 * mrSges(5,3) + qJD(4) * t69 - t74 * t57 - t106;
t7 = m(4) * t129 + t99 * t11 + t134 * t16 + (-m(4) * t49 - t112) * t91;
t8 = m(4) * t120 + t134 * t11 + t112 * t95 - t99 * t16;
t4 = m(3) * t119 - t104 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t91 * t7 + t95 * t8;
t6 = m(3) * t71 + t7 * t95 + t8 * t91;
t117 = m(2) * t89 + t4 * t126 + t93 * t127 + t97 * t6;
t2 = m(2) * t80 - t10 * t100 + t102 * t4;
t1 = m(2) * t79 - t93 * t6 + (t100 * t4 + t127) * t97;
t3 = [-m(1) * g(1) - t1 * t92 + t2 * t96, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t1 * t96 + t2 * t92, t1, t10, t7, t16, t14, t17; -m(1) * g(3) + t117, t117, t6, t114 * qJDD(2) - t104 * t136 + t107, t108, t106, -t110;];
f_new  = t3;
