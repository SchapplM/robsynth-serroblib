% Calculate vector of cutting forces with Newton-Euler
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-05 07:18
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:12:46
% EndTime: 2019-05-05 07:12:56
% DurationCPUTime: 4.95s
% Computational Cost: add. (71294->177), mult. (143606->244), div. (0->0), fcn. (104766->14), ass. (0->99)
t105 = sin(qJ(2));
t108 = cos(qJ(2));
t101 = cos(pkin(6));
t100 = cos(pkin(11));
t97 = sin(pkin(11));
t83 = t97 * g(1) - t100 * g(2);
t129 = t101 * t83;
t84 = -t100 * g(1) - t97 * g(2);
t95 = -g(3) + qJDD(1);
t98 = sin(pkin(6));
t135 = -t105 * t84 + t108 * (t95 * t98 + t129);
t104 = sin(qJ(3));
t107 = cos(qJ(3));
t109 = qJD(2) ^ 2;
t113 = -qJDD(2) * pkin(2) - t135;
t127 = qJD(2) * t104;
t125 = qJD(2) * qJD(3);
t82 = t107 * qJDD(2) - t104 * t125;
t88 = qJD(3) * pkin(3) - pkin(9) * t127;
t94 = t107 ^ 2;
t111 = -t82 * pkin(3) + t88 * t127 + (-pkin(9) * t94 - pkin(8)) * t109 + t113;
t102 = sin(qJ(6));
t106 = cos(qJ(6));
t103 = sin(qJ(4));
t133 = cos(qJ(4));
t128 = t105 * t98;
t123 = t105 * t129 + t108 * t84 + t95 * t128;
t58 = -t109 * pkin(2) + qJDD(2) * pkin(8) + t123;
t71 = t101 * t95 - t98 * t83;
t120 = -t104 * t58 + t107 * t71;
t121 = t107 * t125;
t81 = t104 * qJDD(2) + t121;
t37 = (-t81 + t121) * pkin(9) + (t104 * t107 * t109 + qJDD(3)) * pkin(3) + t120;
t131 = t104 * t71 + t107 * t58;
t38 = -t94 * t109 * pkin(3) + t82 * pkin(9) - qJD(3) * t88 + t131;
t132 = t103 * t37 + t133 * t38;
t126 = qJD(2) * t107;
t74 = t103 * t127 - t133 * t126;
t75 = (t103 * t107 + t104 * t133) * qJD(2);
t60 = t74 * pkin(4) - t75 * qJ(5);
t93 = qJD(3) + qJD(4);
t91 = t93 ^ 2;
t92 = qJDD(3) + qJDD(4);
t25 = -t91 * pkin(4) + t92 * qJ(5) - t74 * t60 + t132;
t52 = t75 * qJD(4) + t103 * t81 - t133 * t82;
t53 = -t74 * qJD(4) + t103 * t82 + t133 * t81;
t28 = (t74 * t93 - t53) * qJ(5) + (t75 * t93 + t52) * pkin(4) + t111;
t96 = sin(pkin(12));
t99 = cos(pkin(12));
t66 = t99 * t75 + t96 * t93;
t119 = -0.2e1 * qJD(5) * t66 - t96 * t25 + t99 * t28;
t46 = t99 * t53 + t96 * t92;
t65 = -t96 * t75 + t99 * t93;
t19 = (t65 * t74 - t46) * pkin(10) + (t65 * t66 + t52) * pkin(5) + t119;
t124 = 0.2e1 * qJD(5) * t65 + t99 * t25 + t96 * t28;
t45 = -t96 * t53 + t99 * t92;
t56 = t74 * pkin(5) - t66 * pkin(10);
t64 = t65 ^ 2;
t20 = -t64 * pkin(5) + t45 * pkin(10) - t74 * t56 + t124;
t43 = -t102 * t66 + t106 * t65;
t31 = t43 * qJD(6) + t102 * t45 + t106 * t46;
t44 = t102 * t65 + t106 * t66;
t33 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t72 = qJD(6) + t74;
t39 = -t72 * mrSges(7,2) + t43 * mrSges(7,3);
t50 = qJDD(6) + t52;
t17 = m(7) * (-t102 * t20 + t106 * t19) - t31 * mrSges(7,3) + t50 * mrSges(7,1) - t44 * t33 + t72 * t39;
t30 = -t44 * qJD(6) - t102 * t46 + t106 * t45;
t40 = t72 * mrSges(7,1) - t44 * mrSges(7,3);
t18 = m(7) * (t102 * t19 + t106 * t20) + t30 * mrSges(7,3) - t50 * mrSges(7,2) + t43 * t33 - t72 * t40;
t47 = -t65 * mrSges(6,1) + t66 * mrSges(6,2);
t54 = -t74 * mrSges(6,2) + t65 * mrSges(6,3);
t14 = m(6) * t119 + t52 * mrSges(6,1) - t46 * mrSges(6,3) + t102 * t18 + t106 * t17 - t66 * t47 + t74 * t54;
t55 = t74 * mrSges(6,1) - t66 * mrSges(6,3);
t15 = m(6) * t124 - t52 * mrSges(6,2) + t45 * mrSges(6,3) - t102 * t17 + t106 * t18 + t65 * t47 - t74 * t55;
t69 = -t93 * mrSges(5,2) - t74 * mrSges(5,3);
t70 = t93 * mrSges(5,1) - t75 * mrSges(5,3);
t114 = m(5) * t111 + t52 * mrSges(5,1) + t53 * mrSges(5,2) + t99 * t14 + t96 * t15 + t74 * t69 + t75 * t70;
t85 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t127;
t86 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t126;
t134 = (t104 * t85 - t107 * t86) * qJD(2) + m(4) * (-t109 * pkin(8) + t113) - t82 * mrSges(4,1) + t81 * mrSges(4,2) + t114;
t10 = m(3) * t135 + qJDD(2) * mrSges(3,1) - t109 * mrSges(3,2) - t134;
t130 = t10 * t108;
t61 = t74 * mrSges(5,1) + t75 * mrSges(5,2);
t11 = m(5) * t132 - t92 * mrSges(5,2) - t52 * mrSges(5,3) - t96 * t14 + t99 * t15 - t74 * t61 - t93 * t70;
t118 = -t103 * t38 + t133 * t37;
t24 = -t92 * pkin(4) - t91 * qJ(5) + t75 * t60 + qJDD(5) - t118;
t115 = t30 * mrSges(7,1) + t43 * t39 - m(7) * (-t45 * pkin(5) - t64 * pkin(10) + t66 * t56 + t24) - t31 * mrSges(7,2) - t44 * t40;
t110 = m(6) * t24 - t45 * mrSges(6,1) + t46 * mrSges(6,2) - t65 * t54 + t66 * t55 - t115;
t16 = m(5) * t118 + t92 * mrSges(5,1) - t53 * mrSges(5,3) - t75 * t61 + t93 * t69 - t110;
t80 = (-mrSges(4,1) * t107 + mrSges(4,2) * t104) * qJD(2);
t7 = m(4) * t120 + qJDD(3) * mrSges(4,1) - t81 * mrSges(4,3) + qJD(3) * t86 + t103 * t11 - t80 * t127 + t133 * t16;
t8 = m(4) * t131 - qJDD(3) * mrSges(4,2) + t82 * mrSges(4,3) - qJD(3) * t85 - t103 * t16 + t11 * t133 + t80 * t126;
t4 = m(3) * t123 - t109 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t104 * t7 + t107 * t8;
t6 = m(3) * t71 + t104 * t8 + t107 * t7;
t122 = m(2) * t95 + t101 * t6 + t4 * t128 + t98 * t130;
t2 = m(2) * t84 - t105 * t10 + t108 * t4;
t1 = m(2) * t83 - t98 * t6 + (t105 * t4 + t130) * t101;
t3 = [-m(1) * g(1) - t97 * t1 + t100 * t2, t2, t4, t8, t11, t15, t18; -m(1) * g(2) + t100 * t1 + t97 * t2, t1, t10, t7, t16, t14, t17; -m(1) * g(3) + t122, t122, t6, t134, t114, t110, -t115;];
f_new  = t3;
