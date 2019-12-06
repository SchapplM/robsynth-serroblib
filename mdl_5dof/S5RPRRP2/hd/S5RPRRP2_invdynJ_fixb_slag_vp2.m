% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:23
% EndTime: 2019-12-05 18:01:30
% DurationCPUTime: 2.39s
% Computational Cost: add. (1896->306), mult. (3271->377), div. (0->0), fcn. (1667->12), ass. (0->136)
t214 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t213 = mrSges(5,2) + mrSges(6,2);
t212 = -mrSges(5,3) - mrSges(6,3);
t211 = Ifges(5,1) + Ifges(6,1);
t209 = Ifges(6,5) + Ifges(5,5);
t208 = Ifges(5,2) + Ifges(6,2);
t207 = Ifges(6,6) + Ifges(5,6);
t126 = -qJ(5) - pkin(7);
t206 = m(6) * t126 + mrSges(4,2);
t121 = qJDD(1) + qJDD(3);
t128 = sin(qJ(3));
t131 = cos(qJ(3));
t124 = sin(pkin(8));
t191 = pkin(1) * t124;
t125 = cos(pkin(8));
t109 = pkin(1) * t125 + pkin(2);
t90 = t109 * qJD(1);
t203 = qJD(3) * t90 + qJDD(1) * t191;
t162 = qJD(1) * t191;
t204 = -qJD(3) * t162 + t109 * qJDD(1);
t17 = t128 * t204 + t131 * t203;
t13 = pkin(7) * t121 + t17;
t205 = qJD(2) * qJD(4) + t13;
t130 = cos(qJ(4));
t202 = (Ifges(5,4) + Ifges(6,4)) * t130;
t127 = sin(qJ(4));
t201 = t127 * t213;
t186 = mrSges(5,1) + mrSges(6,1);
t65 = t109 * t131 - t128 * t191;
t119 = t130 * qJD(2);
t122 = qJD(1) + qJD(3);
t47 = t128 * t90 + t131 * t162;
t37 = pkin(7) * t122 + t47;
t154 = qJ(5) * t122 + t37;
t24 = -t127 * t154 + t119;
t19 = qJD(4) * pkin(4) + t24;
t30 = -t127 * t37 + t119;
t199 = -t30 * mrSges(5,3) - t19 * mrSges(6,3);
t197 = m(5) * pkin(3);
t196 = m(5) * pkin(7);
t190 = pkin(4) * t130;
t167 = qJD(4) * t127;
t5 = t127 * qJDD(2) + t130 * t205 - t167 * t37;
t189 = t130 * t5;
t170 = t122 * t127;
t80 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t170;
t81 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t170;
t184 = t80 + t81;
t169 = t122 * t130;
t82 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t169;
t83 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t169;
t183 = t82 + t83;
t182 = Ifges(5,4) * t127;
t180 = Ifges(6,4) * t127;
t66 = t128 * t109 + t131 * t191;
t63 = pkin(7) + t66;
t177 = -qJ(5) - t63;
t168 = qJD(2) * t127;
t31 = t130 * t37 + t168;
t175 = qJD(4) * t31;
t174 = qJD(4) * t63;
t166 = qJD(4) * t130;
t118 = t130 * qJD(5);
t165 = m(4) + m(5) + m(6);
t123 = qJ(1) + pkin(8);
t161 = pkin(4) * t167;
t160 = m(3) + t165;
t111 = pkin(3) + t190;
t158 = t122 * t167;
t71 = t121 * t130 - t158;
t72 = t121 * t127 + t122 * t166;
t28 = -t71 * mrSges(6,1) + t72 * mrSges(6,2);
t46 = -t128 * t162 + t131 * t90;
t155 = qJD(4) * t126;
t153 = qJD(4) * t177;
t117 = qJ(3) + t123;
t107 = sin(t117);
t108 = cos(t117);
t152 = -t107 * t201 + t108 * t212;
t62 = -pkin(3) - t65;
t25 = t130 * t154 + t168;
t146 = t31 * mrSges(5,3) + t25 * mrSges(6,3);
t145 = g(2) * t108 + g(3) * t107;
t144 = pkin(2) * t165 + mrSges(3,1);
t93 = -mrSges(5,1) * t130 + mrSges(5,2) * t127;
t143 = -mrSges(6,1) * t130 + mrSges(6,2) * t127;
t142 = Ifges(5,2) * t130 + t182;
t141 = Ifges(6,2) * t130 + t180;
t140 = -t127 * t19 + t130 * t25;
t139 = -t127 * t30 + t130 * t31;
t18 = -t128 * t203 + t131 * t204;
t138 = m(6) * t111 + mrSges(4,1) + t197;
t137 = pkin(1) * t160 + mrSges(2,1);
t14 = -pkin(3) * t121 - t18;
t57 = t66 * qJD(3);
t136 = t186 * t130 + t138;
t135 = -t108 * t201 + (t196 - t206 - t212) * t107;
t27 = -t111 * t122 + qJD(5) - t46;
t3 = qJ(5) * t71 + t118 * t122 + t5;
t36 = -pkin(3) * t122 - t46;
t58 = Ifges(6,6) * qJD(4) + t122 * t141;
t59 = Ifges(5,6) * qJD(4) + t122 * t142;
t99 = Ifges(6,4) * t169;
t60 = Ifges(6,1) * t170 + Ifges(6,5) * qJD(4) + t99;
t100 = Ifges(5,4) * t169;
t61 = Ifges(5,1) * t170 + Ifges(5,5) * qJD(4) + t100;
t8 = -pkin(4) * t71 + qJDD(5) + t14;
t134 = t3 * t130 * mrSges(6,3) + t18 * mrSges(4,1) + mrSges(5,3) * t189 + Ifges(4,3) * t121 + t14 * t93 + t8 * t143 + (-t207 * t127 + t209 * t130) * qJD(4) ^ 2 / 0.2e1 - (t59 + t58) * t167 / 0.2e1 + (t130 * t211 - t180 - t182) * t158 / 0.2e1 + (t27 * (mrSges(6,1) * t127 + mrSges(6,2) * t130) + t36 * (mrSges(5,1) * t127 + mrSges(5,2) * t130)) * qJD(4) + (t142 / 0.2e1 + t141 / 0.2e1 + t127 * t214 + t208 * t130 / 0.2e1) * t71 + (t211 * t127 + t202 / 0.2e1 + t130 * t214) * t72 + (t61 + t60 + (-t127 * t208 + t202) * t122) * t166 / 0.2e1 + (t127 * t209 + t130 * t207) * qJDD(4);
t132 = cos(qJ(1));
t129 = sin(qJ(1));
t120 = t130 * qJ(5);
t116 = t130 * qJDD(2);
t114 = cos(t123);
t113 = sin(t123);
t104 = t108 * pkin(7);
t94 = pkin(7) * t130 + t120;
t92 = t126 * t127;
t70 = t93 * t122;
t69 = t143 * t122;
t68 = -qJD(5) * t127 + t130 * t155;
t67 = t127 * t155 + t118;
t56 = t65 * qJD(3);
t51 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t72;
t50 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t72;
t49 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t71;
t48 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t71;
t45 = t62 - t190;
t40 = t57 + t161;
t39 = t130 * t63 + t120;
t38 = t177 * t127;
t29 = -mrSges(5,1) * t71 + mrSges(5,2) * t72;
t11 = (-qJD(5) - t56) * t127 + t130 * t153;
t10 = t127 * t153 + t130 * t56 + t118;
t6 = -t127 * t13 + t116 - t175;
t2 = -t37 * t166 + qJDD(4) * pkin(4) - qJ(5) * t72 + t116 + (-qJD(5) * t122 - t205) * t127;
t1 = [(-t121 * t66 - t122 * t56 - t17) * mrSges(4,2) + (t121 * t65 - t122 * t57) * mrSges(4,1) + (-m(5) * t104 + mrSges(2,2) * t132 + mrSges(3,2) * t114 + t136 * t107 + t108 * t206 + t144 * t113 + t137 * t129 + t152) * g(3) + m(4) * (t17 * t66 + t18 * t65 - t46 * t57 + t47 * t56) + m(6) * (t10 * t25 + t11 * t19 + t2 * t38 + t27 * t40 + t3 * t39 + t45 * t8) + m(5) * (t14 * t62 + t36 * t57) + (-mrSges(2,2) * t129 - mrSges(3,2) * t113 + t136 * t108 + t114 * t144 + t132 * t137 + t135) * g(2) + t11 * t80 + t10 * t82 + t40 * t69 + t57 * t70 + t39 * t48 + t38 * t50 + t62 * t29 + t45 * t28 + t134 + (m(5) * (-t174 * t31 - t30 * t56 - t6 * t63) - t56 * t81 - t63 * t51 - t83 * t174 + (-qJD(4) * t25 - t2) * mrSges(6,3) + (-t6 - t175) * mrSges(5,3)) * t127 + (m(5) * (-t174 * t30 + t31 * t56 + t5 * t63) + t56 * t83 + t63 * t49 - t81 * t174 + t199 * qJD(4)) * t130 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t125 - 0.2e1 * mrSges(3,2) * t124 + m(3) * (t124 ^ 2 + t125 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1); (t50 + t51) * t130 + (t48 + t49) * t127 + (m(3) + m(4)) * qJDD(2) + (-t127 * t184 + t183 * t130) * qJD(4) + m(5) * (qJD(4) * t139 + t127 * t5 + t130 * t6) + m(6) * (qJD(4) * t140 + t127 * t3 + t130 * t2) - t160 * g(1); (-m(5) * (-pkin(3) * t107 + t104) + t107 * mrSges(4,1) + t108 * mrSges(4,2) + t152) * g(3) + (-t6 * t127 + t189 + (-t127 * t31 - t130 * t30) * qJD(4)) * t196 + (t122 * mrSges(4,1) - t69 - t70) * t47 - t14 * t197 + (-t6 * mrSges(5,3) - t2 * mrSges(6,3) - pkin(7) * t51 + t184 * t46 + (pkin(4) * t69 - pkin(7) * t83 - t146) * qJD(4)) * t127 + (pkin(7) * t49 - t183 * t46 + (-pkin(7) * t81 + t199) * qJD(4) + t186 * t145) * t130 + (t108 * t138 + t135) * g(2) + (t122 * t46 - t17) * mrSges(4,2) - t111 * t28 + t68 * t80 + t67 * t82 + t92 * t50 + t94 * t48 - pkin(3) * t29 + t134 - m(5) * (t139 * t46 + t36 * t47) + ((t107 * t111 + t108 * t126) * g(3) - t140 * t46 - t111 * t8 + t19 * t68 + t2 * t92 + t25 * t67 + t3 * t94 + (-t47 + t161) * t27) * m(6); t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t24 * t82 - t30 * t83 + t31 * t81 + t209 * t72 + t207 * t71 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t143 + t93) * g(1) + (t50 + (-g(1) * t130 + t2) * m(6)) * pkin(4) + (-m(6) * (-t19 + t24) + t80) * t25 + ((-t100 / 0.2e1 - t99 / 0.2e1 - t60 / 0.2e1 - t61 / 0.2e1 - t36 * mrSges(5,2) - t27 * mrSges(6,2) + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4) - t199) * t130 + (t58 / 0.2e1 + t59 / 0.2e1 - t36 * mrSges(5,1) - t27 * mrSges(6,1) + t214 * t170 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t27 - t69) * pkin(4) + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t169 + t146) * t127) * t122 + (-g(2) * t107 + g(3) * t108) * ((m(6) * pkin(4) + t186) * t127 + t213 * t130); (t127 * t80 - t130 * t82) * t122 + (-t122 * t140 - t145 + t8) * m(6) + t28;];
tau = t1;
