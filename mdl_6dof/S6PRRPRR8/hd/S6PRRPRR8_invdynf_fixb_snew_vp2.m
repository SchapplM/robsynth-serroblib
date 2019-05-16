% Calculate vector of cutting forces with Newton-Euler
% S6PRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 06:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRPRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:15:28
% EndTime: 2019-05-05 06:15:35
% DurationCPUTime: 3.59s
% Computational Cost: add. (45446->187), mult. (97511->251), div. (0->0), fcn. (72132->14), ass. (0->107)
t149 = -2 * qJD(4);
t109 = qJD(2) ^ 2;
t104 = sin(qJ(2));
t108 = cos(qJ(2));
t97 = sin(pkin(6));
t134 = t108 * t97;
t100 = cos(pkin(6));
t95 = sin(pkin(12));
t98 = cos(pkin(12));
t81 = t95 * g(1) - t98 * g(2);
t136 = t100 * t81;
t82 = -t98 * g(1) - t95 * g(2);
t94 = -g(3) + qJDD(1);
t121 = -t104 * t82 + t108 * t136 + t94 * t134;
t96 = sin(pkin(7));
t44 = t109 * t96 * pkin(9) + qJDD(2) * pkin(2) + t121;
t99 = cos(pkin(7));
t143 = t44 * t99;
t61 = t100 * t94 - t97 * t81;
t148 = t61 * t96 + t143;
t103 = sin(qJ(3));
t132 = qJD(2) * t96;
t124 = t103 * t132;
t90 = t99 * qJD(2) + qJD(3);
t147 = (pkin(3) * t90 + t149) * t124;
t107 = cos(qJ(3));
t146 = t107 * t148;
t130 = qJD(2) * t107;
t125 = t96 * t130;
t135 = t104 * t97;
t127 = t104 * t136 + t108 * t82 + t94 * t135;
t129 = qJDD(2) * t96;
t45 = -t109 * pkin(2) + pkin(9) * t129 + t127;
t128 = t103 * t148 + t107 * t45;
t68 = (-pkin(3) * t107 - qJ(4) * t103) * t132;
t88 = t90 ^ 2;
t89 = t99 * qJDD(2) + qJDD(3);
t145 = t88 * pkin(3) - t89 * qJ(4) - t68 * t125 + t149 * t90 - t128;
t144 = -pkin(3) - pkin(10);
t141 = mrSges(4,1) - mrSges(5,2);
t140 = mrSges(4,3) + mrSges(5,1);
t102 = sin(qJ(5));
t106 = cos(qJ(5));
t42 = t103 * t45;
t117 = -t88 * qJ(4) + t68 * t124 + qJDD(4) + t42;
t133 = t109 * t96 ^ 2;
t72 = (qJD(3) * t130 + qJDD(2) * t103) * t96;
t25 = t72 * pkin(4) + t144 * t89 + (-pkin(10) * t103 * t133 - t143 + (-pkin(4) * qJD(2) * t90 - t61) * t96) * t107 + t117;
t126 = t107 ^ 2 * t133;
t57 = t99 * t61;
t71 = pkin(4) * t124 - t90 * pkin(10);
t73 = -qJD(3) * t124 + t107 * t129;
t27 = -pkin(4) * t126 - t72 * qJ(4) + t57 + t144 * t73 + (-t44 + (-qJ(4) * t107 * t90 - t103 * t71) * qJD(2)) * t96 + t147;
t139 = t102 * t25 + t106 * t27;
t67 = mrSges(5,1) * t124 + t90 * mrSges(5,2);
t138 = t90 * mrSges(4,1) - mrSges(4,3) * t124 - t67;
t69 = (mrSges(5,2) * t107 - mrSges(5,3) * t103) * t132;
t137 = t69 + (-mrSges(4,1) * t107 + mrSges(4,2) * t103) * t132;
t122 = -t96 * t44 + t57;
t101 = sin(qJ(6));
t105 = cos(qJ(6));
t59 = -t102 * t90 - t106 * t125;
t60 = -t102 * t125 + t106 * t90;
t47 = -t59 * pkin(5) - t60 * pkin(11);
t63 = qJDD(5) + t72;
t80 = qJD(5) + t124;
t78 = t80 ^ 2;
t20 = -t78 * pkin(5) + t63 * pkin(11) + t59 * t47 + t139;
t110 = t73 * pkin(4) - pkin(10) * t126 + t90 * t71 - t145;
t39 = -t60 * qJD(5) - t102 * t89 - t106 * t73;
t40 = t59 * qJD(5) - t102 * t73 + t106 * t89;
t21 = (-t59 * t80 - t40) * pkin(11) + (t60 * t80 - t39) * pkin(5) + t110;
t48 = -t101 * t60 + t105 * t80;
t32 = t48 * qJD(6) + t101 * t63 + t105 * t40;
t49 = t101 * t80 + t105 * t60;
t33 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t58 = qJD(6) - t59;
t34 = -t58 * mrSges(7,2) + t48 * mrSges(7,3);
t37 = qJDD(6) - t39;
t17 = m(7) * (-t101 * t20 + t105 * t21) - t32 * mrSges(7,3) + t37 * mrSges(7,1) - t49 * t33 + t58 * t34;
t31 = -t49 * qJD(6) - t101 * t40 + t105 * t63;
t35 = t58 * mrSges(7,1) - t49 * mrSges(7,3);
t18 = m(7) * (t101 * t21 + t105 * t20) + t31 * mrSges(7,3) - t37 * mrSges(7,2) + t48 * t33 - t58 * t35;
t46 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t51 = t80 * mrSges(6,1) - t60 * mrSges(6,3);
t13 = m(6) * t139 - t63 * mrSges(6,2) + t39 * mrSges(6,3) - t101 * t17 + t105 * t18 + t59 * t46 - t80 * t51;
t118 = -t102 * t27 + t106 * t25;
t111 = m(7) * (-t63 * pkin(5) - t78 * pkin(11) + t60 * t47 - t118) - t31 * mrSges(7,1) + t32 * mrSges(7,2) - t48 * t34 + t49 * t35;
t50 = -t80 * mrSges(6,2) + t59 * mrSges(6,3);
t14 = m(6) * t118 + t63 * mrSges(6,1) - t40 * mrSges(6,3) - t60 * t46 + t80 * t50 - t111;
t66 = -mrSges(5,1) * t125 - t90 * mrSges(5,3);
t116 = -t102 * t14 + t106 * t13 + m(5) * (-t73 * pkin(3) + (-t125 * t90 - t72) * qJ(4) + t122 + t147) + t66 * t125 - t72 * mrSges(5,3);
t65 = -t90 * mrSges(4,2) + mrSges(4,3) * t125;
t10 = m(4) * t122 + t72 * mrSges(4,2) - t141 * t73 + (t103 * t138 - t107 * t65) * t132 + t116;
t113 = m(6) * t110 - t39 * mrSges(6,1) + t40 * mrSges(6,2) + t101 * t18 + t105 * t17 - t59 * t50 + t60 * t51;
t112 = -m(5) * t145 + t113;
t11 = m(4) * t128 - t138 * t90 + (-mrSges(4,2) + mrSges(5,3)) * t89 + t140 * t73 + t137 * t125 + t112;
t114 = -m(5) * (-t89 * pkin(3) + t117 - t146) - t102 * t13 - t106 * t14;
t9 = m(4) * (-t42 + t146) + (t65 - t66) * t90 + t141 * t89 - t140 * t72 - t137 * t124 + t114;
t119 = t103 * t11 + t107 * t9;
t4 = m(3) * t121 + qJDD(2) * mrSges(3,1) - t109 * mrSges(3,2) - t96 * t10 + t119 * t99;
t6 = m(3) * t61 + t99 * t10 + t119 * t96;
t8 = m(3) * t127 - t109 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t103 * t9 + t107 * t11;
t123 = m(2) * t94 + t100 * t6 + t4 * t134 + t8 * t135;
t2 = m(2) * t82 - t104 * t4 + t108 * t8;
t1 = m(2) * t81 - t97 * t6 + (t104 * t8 + t108 * t4) * t100;
t3 = [-m(1) * g(1) - t95 * t1 + t98 * t2, t2, t8, t11, t73 * mrSges(5,2) - t124 * t67 + t116, t13, t18; -m(1) * g(2) + t98 * t1 + t95 * t2, t1, t4, t9, -t73 * mrSges(5,1) - t89 * mrSges(5,3) - t125 * t69 - t90 * t67 - t112, t14, t17; -m(1) * g(3) + t123, t123, t6, t10, t72 * mrSges(5,1) + t89 * mrSges(5,2) + t124 * t69 + t90 * t66 - t114, t113, t111;];
f_new  = t3;
