% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 05:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:06:30
% EndTime: 2019-05-05 05:06:34
% DurationCPUTime: 1.88s
% Computational Cost: add. (22819->178), mult. (45072->233), div. (0->0), fcn. (29033->12), ass. (0->95)
t104 = sin(qJ(3));
t108 = cos(qJ(3));
t129 = qJD(2) * qJD(3);
t124 = t108 * t129;
t111 = qJD(2) ^ 2;
t105 = sin(qJ(2));
t109 = cos(qJ(2));
t98 = sin(pkin(6));
t133 = t109 * t98;
t100 = cos(pkin(6));
t97 = sin(pkin(11));
t99 = cos(pkin(11));
t75 = t97 * g(1) - t99 * g(2);
t136 = t100 * t75;
t76 = -t99 * g(1) - t97 * g(2);
t95 = -g(3) + qJDD(1);
t128 = -t105 * t76 + t109 * t136 + t95 * t133;
t43 = -qJDD(2) * pkin(2) - t111 * pkin(8) - t128;
t72 = t104 * qJDD(2) + t124;
t125 = t104 * t129;
t73 = t108 * qJDD(2) - t125;
t117 = -t73 * pkin(3) + t43 + (-t124 - t72) * qJ(4);
t102 = sin(qJ(6));
t106 = cos(qJ(6));
t131 = qJD(2) * t104;
t132 = t111 * t108 ^ 2;
t143 = 2 * qJD(4);
t83 = -qJD(3) * pkin(4) - pkin(9) * t131;
t112 = -pkin(3) * t125 + t73 * pkin(4) - pkin(9) * t132 - t117 + (t143 + t83) * t131;
t103 = sin(qJ(5));
t107 = cos(qJ(5));
t110 = qJD(3) ^ 2;
t130 = qJD(2) * t108;
t135 = t105 * t98;
t127 = t105 * t136 + t109 * t76 + t95 * t135;
t44 = -t111 * pkin(2) + qJDD(2) * pkin(8) + t127;
t54 = t100 * t95 - t98 * t75;
t139 = t104 * t54 + t108 * t44;
t69 = (-pkin(3) * t108 - qJ(4) * t104) * qJD(2);
t119 = -t110 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t143 + t69 * t130 + t139;
t23 = -pkin(4) * t132 - t73 * pkin(9) + qJD(3) * t83 + t119;
t37 = t104 * t44;
t121 = -t110 * qJ(4) + t69 * t131 + qJDD(4) + t37;
t24 = -t72 * pkin(9) + (-pkin(3) - pkin(4)) * qJDD(3) + (-pkin(4) * t104 * t111 + pkin(9) * t129 - t54) * t108 + t121;
t140 = t103 * t24 + t107 * t23;
t61 = (t103 * t104 + t107 * t108) * qJD(2);
t62 = (-t103 * t108 + t104 * t107) * qJD(2);
t48 = t61 * pkin(5) - t62 * pkin(10);
t91 = -qJD(3) + qJD(5);
t89 = t91 ^ 2;
t90 = -qJDD(3) + qJDD(5);
t19 = -t89 * pkin(5) + t90 * pkin(10) - t61 * t48 + t140;
t39 = -t62 * qJD(5) - t103 * t72 - t107 * t73;
t40 = -t61 * qJD(5) - t103 * t73 + t107 * t72;
t20 = (t62 * t91 - t39) * pkin(5) + (t61 * t91 - t40) * pkin(10) + t112;
t49 = -t102 * t62 + t106 * t91;
t32 = t49 * qJD(6) + t102 * t90 + t106 * t40;
t50 = t102 * t91 + t106 * t62;
t33 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t36 = qJDD(6) - t39;
t58 = qJD(6) + t61;
t41 = -t58 * mrSges(7,2) + t49 * mrSges(7,3);
t16 = m(7) * (-t102 * t19 + t106 * t20) - t32 * mrSges(7,3) + t36 * mrSges(7,1) - t50 * t33 + t58 * t41;
t31 = -t50 * qJD(6) - t102 * t40 + t106 * t90;
t42 = t58 * mrSges(7,1) - t50 * mrSges(7,3);
t17 = m(7) * (t102 * t20 + t106 * t19) + t31 * mrSges(7,3) - t36 * mrSges(7,2) + t49 * t33 - t58 * t42;
t52 = -t91 * mrSges(6,2) - t61 * mrSges(6,3);
t53 = t91 * mrSges(6,1) - t62 * mrSges(6,3);
t120 = m(6) * t112 - t39 * mrSges(6,1) + t40 * mrSges(6,2) + t102 * t17 + t106 * t16 + t61 * t52 + t62 * t53;
t115 = m(5) * ((pkin(3) * qJD(3) - (2 * qJD(4))) * t131 + t117) - t73 * mrSges(5,1) - t120;
t80 = mrSges(5,2) * t130 + qJD(3) * mrSges(5,3);
t137 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t130 + t80;
t77 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t131;
t78 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t131;
t146 = (-t137 * t108 + (t77 - t78) * t104) * qJD(2) + m(4) * t43 - t73 * mrSges(4,1) + (mrSges(4,2) - mrSges(5,3)) * t72 + t115;
t141 = mrSges(4,3) + mrSges(5,2);
t134 = t108 * t54;
t10 = m(3) * t128 + qJDD(2) * mrSges(3,1) - t111 * mrSges(3,2) - t146;
t47 = t61 * mrSges(6,1) + t62 * mrSges(6,2);
t12 = m(6) * t140 - t90 * mrSges(6,2) + t39 * mrSges(6,3) - t102 * t16 + t106 * t17 - t61 * t47 - t91 * t53;
t123 = -t103 * t23 + t107 * t24;
t114 = m(7) * (-t90 * pkin(5) - t89 * pkin(10) + t62 * t48 - t123) - t31 * mrSges(7,1) + t32 * mrSges(7,2) - t49 * t41 + t50 * t42;
t13 = m(6) * t123 + t90 * mrSges(6,1) - t40 * mrSges(6,3) - t62 * t47 + t91 * t52 - t114;
t70 = (-mrSges(5,1) * t108 - mrSges(5,3) * t104) * qJD(2);
t118 = m(5) * t119 + qJDD(3) * mrSges(5,3) + qJD(3) * t78 - t103 * t13 + t107 * t12 + t70 * t130;
t71 = (-mrSges(4,1) * t108 + mrSges(4,2) * t104) * qJD(2);
t7 = m(4) * t139 - qJDD(3) * mrSges(4,2) - qJD(3) * t77 + t71 * t130 + t141 * t73 + t118;
t116 = -m(5) * (-qJDD(3) * pkin(3) + t121 - t134) - t103 * t12 - t107 * t13;
t8 = m(4) * (-t37 + t134) - t141 * t72 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + t137 * qJD(3) + (-t70 - t71) * t131 + t116;
t4 = m(3) * t127 - t111 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t104 * t8 + t108 * t7;
t6 = m(3) * t54 + t104 * t7 + t108 * t8;
t126 = m(2) * t95 + t10 * t133 + t100 * t6 + t4 * t135;
t2 = m(2) * t76 - t105 * t10 + t109 * t4;
t1 = m(2) * t75 - t98 * t6 + (t10 * t109 + t105 * t4) * t100;
t3 = [-m(1) * g(1) - t97 * t1 + t99 * t2, t2, t4, t7, t73 * mrSges(5,2) + t118, t12, t17; -m(1) * g(2) + t99 * t1 + t97 * t2, t1, t10, t8, -t72 * mrSges(5,3) + (-t104 * t78 - t108 * t80) * qJD(2) + t115, t13, t16; -m(1) * g(3) + t126, t126, t6, t146, -qJDD(3) * mrSges(5,1) + t72 * mrSges(5,2) - qJD(3) * t80 + t131 * t70 - t116, t120, t114;];
f_new  = t3;
