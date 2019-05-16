% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 10:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:41:40
% EndTime: 2019-05-05 10:41:50
% DurationCPUTime: 5.06s
% Computational Cost: add. (75370->179), mult. (147401->242), div. (0->0), fcn. (108408->14), ass. (0->101)
t105 = sin(qJ(2));
t110 = cos(qJ(2));
t100 = cos(pkin(6));
t97 = sin(pkin(12));
t99 = cos(pkin(12));
t84 = t97 * g(1) - t99 * g(2);
t130 = t100 * t84;
t85 = -t99 * g(1) - t97 * g(2);
t96 = -g(3) + qJDD(1);
t98 = sin(pkin(6));
t136 = -t105 * t85 + (t96 * t98 + t130) * t110;
t104 = sin(qJ(3));
t109 = cos(qJ(3));
t111 = qJD(2) ^ 2;
t115 = -qJDD(2) * pkin(2) - t136;
t102 = sin(qJ(5));
t107 = cos(qJ(5));
t128 = qJD(2) * t104;
t126 = qJD(2) * qJD(3);
t83 = t109 * qJDD(2) - t104 * t126;
t89 = qJD(3) * pkin(3) - pkin(9) * t128;
t95 = t109 ^ 2;
t113 = -t83 * pkin(3) + t89 * t128 + (-pkin(9) * t95 - pkin(8)) * t111 + t115;
t101 = sin(qJ(6));
t106 = cos(qJ(6));
t103 = sin(qJ(4));
t108 = cos(qJ(4));
t129 = t105 * t98;
t125 = t105 * t130 + t110 * t85 + t96 * t129;
t59 = -t111 * pkin(2) + qJDD(2) * pkin(8) + t125;
t70 = t100 * t96 - t98 * t84;
t121 = -t104 * t59 + t109 * t70;
t123 = t109 * t126;
t82 = t104 * qJDD(2) + t123;
t37 = (-t82 + t123) * pkin(9) + (t104 * t109 * t111 + qJDD(3)) * pkin(3) + t121;
t132 = t104 * t70 + t109 * t59;
t38 = -t95 * t111 * pkin(3) + t83 * pkin(9) - qJD(3) * t89 + t132;
t133 = t103 * t37 + t108 * t38;
t127 = qJD(2) * t109;
t75 = -t103 * t128 + t108 * t127;
t76 = (t103 * t109 + t104 * t108) * qJD(2);
t62 = -t75 * pkin(4) - t76 * pkin(10);
t94 = qJD(3) + qJD(4);
t92 = t94 ^ 2;
t93 = qJDD(3) + qJDD(4);
t25 = -t92 * pkin(4) + t93 * pkin(10) + t75 * t62 + t133;
t53 = -t76 * qJD(4) - t103 * t82 + t108 * t83;
t54 = t75 * qJD(4) + t103 * t83 + t108 * t82;
t31 = (-t75 * t94 - t54) * pkin(10) + (t76 * t94 - t53) * pkin(4) + t113;
t122 = -t102 * t25 + t107 * t31;
t64 = -t102 * t76 + t107 * t94;
t40 = t64 * qJD(5) + t102 * t93 + t107 * t54;
t51 = qJDD(5) - t53;
t65 = t102 * t94 + t107 * t76;
t73 = qJD(5) - t75;
t19 = (t64 * t73 - t40) * pkin(11) + (t64 * t65 + t51) * pkin(5) + t122;
t134 = t102 * t31 + t107 * t25;
t39 = -t65 * qJD(5) - t102 * t54 + t107 * t93;
t57 = t73 * pkin(5) - t65 * pkin(11);
t63 = t64 ^ 2;
t20 = -t63 * pkin(5) + t39 * pkin(11) - t73 * t57 + t134;
t45 = -t101 * t65 + t106 * t64;
t28 = t45 * qJD(6) + t101 * t39 + t106 * t40;
t46 = t101 * t64 + t106 * t65;
t33 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t71 = qJD(6) + t73;
t41 = -t71 * mrSges(7,2) + t45 * mrSges(7,3);
t48 = qJDD(6) + t51;
t17 = m(7) * (-t101 * t20 + t106 * t19) - t28 * mrSges(7,3) + t48 * mrSges(7,1) - t46 * t33 + t71 * t41;
t27 = -t46 * qJD(6) - t101 * t40 + t106 * t39;
t42 = t71 * mrSges(7,1) - t46 * mrSges(7,3);
t18 = m(7) * (t101 * t19 + t106 * t20) + t27 * mrSges(7,3) - t48 * mrSges(7,2) + t45 * t33 - t71 * t42;
t47 = -t64 * mrSges(6,1) + t65 * mrSges(6,2);
t55 = -t73 * mrSges(6,2) + t64 * mrSges(6,3);
t14 = m(6) * t122 + t51 * mrSges(6,1) - t40 * mrSges(6,3) + t101 * t18 + t106 * t17 - t65 * t47 + t73 * t55;
t56 = t73 * mrSges(6,1) - t65 * mrSges(6,3);
t15 = m(6) * t134 - t51 * mrSges(6,2) + t39 * mrSges(6,3) - t101 * t17 + t106 * t18 + t64 * t47 - t73 * t56;
t68 = -t94 * mrSges(5,2) + t75 * mrSges(5,3);
t69 = t94 * mrSges(5,1) - t76 * mrSges(5,3);
t116 = -m(5) * t113 + t53 * mrSges(5,1) - t54 * mrSges(5,2) - t102 * t15 - t107 * t14 + t75 * t68 - t76 * t69;
t86 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t128;
t87 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t127;
t135 = (t104 * t86 - t109 * t87) * qJD(2) + m(4) * (-t111 * pkin(8) + t115) - t83 * mrSges(4,1) + t82 * mrSges(4,2) - t116;
t10 = m(3) * t136 + qJDD(2) * mrSges(3,1) - t111 * mrSges(3,2) - t135;
t131 = t10 * t110;
t61 = -t75 * mrSges(5,1) + t76 * mrSges(5,2);
t11 = m(5) * t133 - t93 * mrSges(5,2) + t53 * mrSges(5,3) - t102 * t14 + t107 * t15 + t75 * t61 - t94 * t69;
t120 = -t103 * t38 + t108 * t37;
t24 = -t93 * pkin(4) - t92 * pkin(10) + t76 * t62 - t120;
t117 = t27 * mrSges(7,1) + t45 * t41 - m(7) * (-t39 * pkin(5) - t63 * pkin(11) + t65 * t57 + t24) - t28 * mrSges(7,2) - t46 * t42;
t112 = m(6) * t24 - t39 * mrSges(6,1) + t40 * mrSges(6,2) - t64 * t55 + t65 * t56 - t117;
t16 = m(5) * t120 + t93 * mrSges(5,1) - t54 * mrSges(5,3) - t76 * t61 + t94 * t68 - t112;
t81 = (-mrSges(4,1) * t109 + mrSges(4,2) * t104) * qJD(2);
t7 = m(4) * t121 + qJDD(3) * mrSges(4,1) - t82 * mrSges(4,3) + qJD(3) * t87 + t103 * t11 + t108 * t16 - t81 * t128;
t8 = m(4) * t132 - qJDD(3) * mrSges(4,2) + t83 * mrSges(4,3) - qJD(3) * t86 - t103 * t16 + t108 * t11 + t81 * t127;
t4 = m(3) * t125 - t111 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t104 * t7 + t109 * t8;
t6 = m(3) * t70 + t104 * t8 + t109 * t7;
t124 = m(2) * t96 + t100 * t6 + t4 * t129 + t98 * t131;
t2 = m(2) * t85 - t105 * t10 + t110 * t4;
t1 = m(2) * t84 - t98 * t6 + (t105 * t4 + t131) * t100;
t3 = [-m(1) * g(1) - t97 * t1 + t99 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t99 * t1 + t97 * t2, t1, t10, t7, t16, t14, t17; -m(1) * g(3) + t124, t124, t6, t135, -t116, t112, -t117;];
f_new  = t3;
