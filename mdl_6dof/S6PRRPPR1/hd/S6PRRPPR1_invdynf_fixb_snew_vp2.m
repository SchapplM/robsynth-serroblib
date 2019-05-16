% Calculate vector of cutting forces with Newton-Euler
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-05-05 02:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:30:16
% EndTime: 2019-05-05 02:30:25
% DurationCPUTime: 4.55s
% Computational Cost: add. (60475->176), mult. (134585->245), div. (0->0), fcn. (97078->14), ass. (0->98)
t138 = -2 * qJD(4);
t130 = cos(pkin(11));
t103 = sin(qJ(3));
t106 = cos(qJ(3));
t109 = qJD(2) ^ 2;
t104 = sin(qJ(2));
t107 = cos(qJ(2));
t98 = sin(pkin(6));
t131 = t104 * t98;
t101 = cos(pkin(6));
t100 = cos(pkin(10));
t97 = sin(pkin(10));
t85 = g(1) * t97 - g(2) * t100;
t132 = t101 * t85;
t86 = -g(1) * t100 - g(2) * t97;
t94 = -g(3) + qJDD(1);
t123 = t104 * t132 + t107 * t86 + t94 * t131;
t51 = -pkin(2) * t109 + qJDD(2) * pkin(8) + t123;
t71 = t101 * t94 - t85 * t98;
t120 = -t103 * t51 + t106 * t71;
t126 = qJD(2) * qJD(3);
t121 = t106 * t126;
t83 = qJDD(2) * t103 + t121;
t37 = (-t83 + t121) * qJ(4) + (t103 * t106 * t109 + qJDD(3)) * pkin(3) + t120;
t134 = t103 * t71 + t106 * t51;
t84 = qJDD(2) * t106 - t103 * t126;
t128 = qJD(2) * t103;
t87 = qJD(3) * pkin(3) - qJ(4) * t128;
t93 = t106 ^ 2;
t38 = -pkin(3) * t109 * t93 + qJ(4) * t84 - qJD(3) * t87 + t134;
t96 = sin(pkin(11));
t76 = (t130 * t103 + t106 * t96) * qJD(2);
t137 = t130 * t37 + t76 * t138 - t96 * t38;
t136 = -t104 * t86 + (t94 * t98 + t132) * t107;
t113 = -qJDD(2) * pkin(2) - t136;
t110 = -t84 * pkin(3) + qJDD(4) + t87 * t128 + (-qJ(4) * t93 - pkin(8)) * t109 + t113;
t102 = sin(qJ(6));
t105 = cos(qJ(6));
t108 = qJD(3) ^ 2;
t127 = qJD(2) * t106;
t75 = -t130 * t127 + t96 * t128;
t124 = t130 * t38 + t75 * t138 + t96 * t37;
t56 = pkin(4) * t75 - qJ(5) * t76;
t25 = -pkin(4) * t108 + qJDD(3) * qJ(5) - t56 * t75 + t124;
t60 = -t130 * t84 + t83 * t96;
t61 = t130 * t83 + t96 * t84;
t28 = (qJD(3) * t75 - t61) * qJ(5) + (qJD(3) * t76 + t60) * pkin(4) + t110;
t95 = sin(pkin(12));
t99 = cos(pkin(12));
t66 = qJD(3) * t95 + t76 * t99;
t119 = -0.2e1 * qJD(5) * t66 - t95 * t25 + t99 * t28;
t54 = qJDD(3) * t95 + t61 * t99;
t65 = qJD(3) * t99 - t76 * t95;
t19 = (t65 * t75 - t54) * pkin(9) + (t65 * t66 + t60) * pkin(5) + t119;
t125 = 0.2e1 * qJD(5) * t65 + t99 * t25 + t95 * t28;
t52 = pkin(5) * t75 - t66 * pkin(9);
t53 = qJDD(3) * t99 - t61 * t95;
t64 = t65 ^ 2;
t20 = -t64 * pkin(5) + t53 * pkin(9) - t52 * t75 + t125;
t43 = -t102 * t66 + t105 * t65;
t31 = t43 * qJD(6) + t102 * t53 + t105 * t54;
t44 = t102 * t65 + t105 * t66;
t33 = -mrSges(7,1) * t43 + mrSges(7,2) * t44;
t74 = qJD(6) + t75;
t41 = -mrSges(7,2) * t74 + t43 * mrSges(7,3);
t59 = qJDD(6) + t60;
t17 = m(7) * (-t102 * t20 + t105 * t19) - t31 * mrSges(7,3) + t59 * mrSges(7,1) - t44 * t33 + t74 * t41;
t30 = -t44 * qJD(6) - t102 * t54 + t105 * t53;
t42 = mrSges(7,1) * t74 - t44 * mrSges(7,3);
t18 = m(7) * (t102 * t19 + t105 * t20) + t30 * mrSges(7,3) - t59 * mrSges(7,2) + t43 * t33 - t74 * t42;
t45 = -mrSges(6,1) * t65 + mrSges(6,2) * t66;
t48 = -mrSges(6,2) * t75 + t65 * mrSges(6,3);
t14 = m(6) * t119 + t60 * mrSges(6,1) - t54 * mrSges(6,3) + t102 * t18 + t105 * t17 - t66 * t45 + t75 * t48;
t49 = mrSges(6,1) * t75 - t66 * mrSges(6,3);
t15 = m(6) * t125 - t60 * mrSges(6,2) + t53 * mrSges(6,3) - t102 * t17 + t105 * t18 + t65 * t45 - t75 * t49;
t69 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t75;
t70 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t76;
t114 = m(5) * t110 + t60 * mrSges(5,1) + t61 * mrSges(5,2) + t99 * t14 + t95 * t15 + t75 * t69 + t76 * t70;
t88 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t128;
t89 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t127;
t135 = (t103 * t88 - t106 * t89) * qJD(2) + m(4) * (-t109 * pkin(8) + t113) - t84 * mrSges(4,1) + t83 * mrSges(4,2) + t114;
t10 = m(3) * t136 + qJDD(2) * mrSges(3,1) - t109 * mrSges(3,2) - t135;
t133 = t10 * t107;
t57 = mrSges(5,1) * t75 + mrSges(5,2) * t76;
t11 = m(5) * t124 - qJDD(3) * mrSges(5,2) - t60 * mrSges(5,3) - qJD(3) * t70 - t95 * t14 + t99 * t15 - t75 * t57;
t24 = -qJDD(3) * pkin(4) - t108 * qJ(5) + t76 * t56 + qJDD(5) - t137;
t115 = t30 * mrSges(7,1) + t43 * t41 - m(7) * (-t53 * pkin(5) - t64 * pkin(9) + t66 * t52 + t24) - t31 * mrSges(7,2) - t44 * t42;
t111 = m(6) * t24 - t53 * mrSges(6,1) + t54 * mrSges(6,2) - t65 * t48 + t66 * t49 - t115;
t16 = m(5) * t137 + qJDD(3) * mrSges(5,1) - t61 * mrSges(5,3) + qJD(3) * t69 - t76 * t57 - t111;
t82 = (-mrSges(4,1) * t106 + mrSges(4,2) * t103) * qJD(2);
t7 = m(4) * t120 + qJDD(3) * mrSges(4,1) - t83 * mrSges(4,3) + qJD(3) * t89 + t96 * t11 - t82 * t128 + t130 * t16;
t8 = m(4) * t134 - qJDD(3) * mrSges(4,2) + t84 * mrSges(4,3) - qJD(3) * t88 + t130 * t11 + t82 * t127 - t96 * t16;
t4 = m(3) * t123 - t109 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t103 * t7 + t106 * t8;
t6 = m(3) * t71 + t103 * t8 + t106 * t7;
t122 = m(2) * t94 + t101 * t6 + t4 * t131 + t98 * t133;
t2 = m(2) * t86 - t10 * t104 + t107 * t4;
t1 = m(2) * t85 - t98 * t6 + (t104 * t4 + t133) * t101;
t3 = [-m(1) * g(1) - t1 * t97 + t100 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t1 * t100 + t2 * t97, t1, t10, t7, t16, t14, t17; -m(1) * g(3) + t122, t122, t6, t135, t114, t111, -t115;];
f_new  = t3;
