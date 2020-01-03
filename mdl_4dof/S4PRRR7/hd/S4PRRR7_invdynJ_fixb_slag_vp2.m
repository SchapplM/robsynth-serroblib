% Calculate vector of inverse dynamics joint torques for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:04
% EndTime: 2019-12-31 16:36:13
% DurationCPUTime: 4.36s
% Computational Cost: add. (1435->333), mult. (3573->513), div. (0->0), fcn. (2519->10), ass. (0->172)
t213 = m(5) + m(4);
t100 = cos(qJ(4));
t155 = qJD(3) * t100;
t98 = sin(qJ(3));
t163 = qJD(2) * t98;
t97 = sin(qJ(4));
t69 = -t163 * t97 + t155;
t220 = -t69 / 0.2e1;
t70 = qJD(3) * t97 + t100 * t163;
t219 = -t70 / 0.2e1;
t101 = cos(qJ(3));
t156 = qJD(2) * t101;
t91 = qJD(4) - t156;
t218 = -t91 / 0.2e1;
t126 = mrSges(5,1) * t97 + mrSges(5,2) * t100;
t217 = -t126 - mrSges(4,3) + mrSges(3,2);
t125 = -mrSges(5,1) * t100 + mrSges(5,2) * t97;
t111 = m(5) * pkin(3) - t125;
t127 = mrSges(4,1) * t101 - mrSges(4,2) * t98;
t143 = m(5) * pkin(7) + mrSges(5,3);
t200 = t101 * t111 + t143 * t98 + mrSges(3,1) + t127;
t216 = pkin(2) * t213 + t200;
t149 = qJD(2) * qJD(3);
t74 = qJDD(2) * t101 - t149 * t98;
t215 = t74 / 0.2e1;
t75 = qJDD(2) * t98 + t101 * t149;
t214 = t75 / 0.2e1;
t26 = qJD(4) * t69 + qJDD(3) * t97 + t100 * t75;
t27 = -qJD(4) * t70 + qJDD(3) * t100 - t75 * t97;
t5 = -mrSges(5,1) * t27 + mrSges(5,2) * t26;
t212 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t75 + t5;
t147 = mrSges(4,3) * t163;
t210 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t69 + mrSges(5,2) * t70 + t147;
t152 = qJD(4) * t101;
t102 = cos(qJ(2));
t157 = t101 * t102;
t160 = qJD(4) * t97;
t161 = qJD(3) * t98;
t95 = sin(pkin(4));
t166 = qJD(1) * t95;
t128 = pkin(3) * t98 - pkin(7) * t101;
t73 = t128 * qJD(3);
t78 = -pkin(3) * t101 - pkin(7) * t98 - pkin(2);
t99 = sin(qJ(2));
t209 = -(t100 * t99 - t157 * t97) * t166 - t78 * t160 + t100 * t73 + (-t100 * t152 + t161 * t97) * pkin(6);
t153 = qJD(4) * t100;
t208 = -(t100 * t157 + t97 * t99) * t166 + t78 * t153 + t73 * t97 + (-t152 * t97 - t155 * t98) * pkin(6);
t207 = mrSges(4,1) + t111;
t206 = mrSges(4,2) - t143;
t67 = Ifges(5,4) * t69;
t23 = Ifges(5,1) * t70 + Ifges(5,5) * t91 + t67;
t93 = Ifges(4,4) * t156;
t205 = Ifges(4,1) * t163 + Ifges(4,5) * qJD(3) + t100 * t23 + t93;
t96 = cos(pkin(4));
t165 = qJD(1) * t96;
t141 = t101 * t165;
t150 = qJDD(1) * t96;
t151 = qJDD(1) * t95;
t164 = qJD(2) * t95;
t137 = qJD(1) * t164;
t85 = t102 * t137;
t55 = t99 * t151 + t85;
t49 = qJDD(2) * pkin(6) + t55;
t145 = t99 * t166;
t76 = qJD(2) * pkin(6) + t145;
t10 = qJD(3) * t141 + t101 * t49 + t98 * t150 - t161 * t76;
t46 = t101 * t76 + t165 * t98;
t11 = -t46 * qJD(3) + t101 * t150 - t49 * t98;
t204 = t10 * t101 - t11 * t98;
t39 = qJD(3) * pkin(7) + t46;
t140 = t102 * t166;
t47 = qJD(2) * t78 - t140;
t12 = t100 * t47 - t39 * t97;
t84 = t99 * t137;
t54 = t102 * t151 - t84;
t48 = -qJDD(2) * pkin(2) - t54;
t16 = -pkin(3) * t74 - pkin(7) * t75 + t48;
t8 = qJDD(3) * pkin(7) + t10;
t1 = qJD(4) * t12 + t100 * t8 + t16 * t97;
t13 = t100 * t39 + t47 * t97;
t2 = -qJD(4) * t13 + t100 * t16 - t8 * t97;
t203 = t1 * t100 - t2 * t97;
t187 = Ifges(4,4) * t98;
t122 = Ifges(4,2) * t101 + t187;
t201 = Ifges(4,6) * qJD(3) / 0.2e1 + qJD(2) * t122 / 0.2e1 + Ifges(5,5) * t219 + Ifges(5,6) * t220 + Ifges(5,3) * t218;
t199 = -t2 * mrSges(5,1) + t1 * mrSges(5,2);
t198 = -t213 * pkin(6) + t217;
t103 = qJD(2) ^ 2;
t186 = Ifges(5,4) * t70;
t22 = Ifges(5,2) * t69 + Ifges(5,6) * t91 + t186;
t197 = -t22 / 0.2e1;
t196 = t26 / 0.2e1;
t195 = t27 / 0.2e1;
t68 = qJDD(4) - t74;
t194 = t68 / 0.2e1;
t192 = t70 / 0.2e1;
t9 = -qJDD(3) * pkin(3) - t11;
t188 = t9 * t98;
t185 = Ifges(5,4) * t97;
t94 = sin(pkin(8));
t182 = t94 * t99;
t181 = t95 * t99;
t180 = t97 * t98;
t176 = Ifges(5,4) * t100;
t173 = t100 * t98;
t172 = t101 * Ifges(4,4);
t171 = t101 * t95;
t170 = t101 * t97;
t169 = t102 * t95;
t168 = t94 * t102;
t167 = cos(pkin(8));
t162 = qJD(2) * t99;
t159 = qJD(4) * t98;
t158 = t100 * t101;
t154 = qJD(3) * t101;
t148 = Ifges(5,5) * t26 + Ifges(5,6) * t27 + Ifges(5,3) * t68;
t144 = t95 * t162;
t142 = mrSges(4,3) * t156;
t139 = t102 * t164;
t138 = t97 * t154;
t135 = t95 * t167;
t134 = t167 * t99;
t131 = t167 * t102;
t130 = t149 / 0.2e1;
t129 = t98 * t140;
t124 = Ifges(5,1) * t97 + t176;
t123 = Ifges(5,1) * t100 - t185;
t121 = Ifges(5,2) * t100 + t185;
t120 = -Ifges(5,2) * t97 + t176;
t119 = Ifges(4,5) * t101 - Ifges(4,6) * t98;
t118 = Ifges(5,5) * t97 + Ifges(5,6) * t100;
t117 = Ifges(5,5) * t100 - Ifges(5,6) * t97;
t62 = -t96 * t101 + t181 * t98;
t63 = t171 * t99 + t96 * t98;
t116 = -t100 * t63 + t169 * t97;
t35 = -t100 * t169 - t63 * t97;
t45 = -t76 * t98 + t141;
t77 = -qJD(2) * pkin(2) - t140;
t113 = t77 * (mrSges(4,1) * t98 + mrSges(4,2) * t101);
t112 = t98 * (Ifges(4,1) * t101 - t187);
t109 = t153 * t98 + t138;
t108 = t100 * t154 - t159 * t97;
t107 = Ifges(5,5) * t98 + t101 * t123;
t106 = Ifges(5,6) * t98 + t101 * t120;
t105 = Ifges(5,3) * t98 + t101 * t117;
t80 = -qJD(3) * mrSges(4,2) + t142;
t72 = t128 * qJD(2);
t71 = t127 * qJD(2);
t61 = -t182 * t96 + t131;
t60 = t168 * t96 + t134;
t59 = t134 * t96 + t168;
t58 = -t131 * t96 + t182;
t56 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t74;
t51 = pkin(6) * t158 + t78 * t97;
t50 = -pkin(6) * t170 + t100 * t78;
t44 = mrSges(5,1) * t91 - mrSges(5,3) * t70;
t43 = -mrSges(5,2) * t91 + mrSges(5,3) * t69;
t38 = -qJD(3) * pkin(3) - t45;
t37 = -mrSges(4,1) * t74 + mrSges(4,2) * t75;
t34 = -qJD(3) * t62 + t101 * t139;
t33 = qJD(3) * t63 + t139 * t98;
t32 = t94 * t95 * t98 + t101 * t61;
t30 = t101 * t59 - t135 * t98;
t20 = t100 * t45 + t72 * t97;
t19 = t100 * t72 - t45 * t97;
t15 = -mrSges(5,2) * t68 + mrSges(5,3) * t27;
t14 = mrSges(5,1) * t68 - mrSges(5,3) * t26;
t7 = qJD(4) * t35 + t100 * t34 + t144 * t97;
t6 = qJD(4) * t116 + t100 * t144 - t34 * t97;
t4 = t26 * Ifges(5,1) + t27 * Ifges(5,4) + t68 * Ifges(5,5);
t3 = Ifges(5,4) * t26 + Ifges(5,2) * t27 + Ifges(5,6) * t68;
t17 = [m(2) * qJDD(1) + t35 * t14 - t116 * t15 + t34 * t80 + t7 * t43 + t6 * t44 + t63 * t56 + t212 * t62 + t210 * t33 + ((-mrSges(3,1) * t103 - mrSges(3,2) * qJDD(2) - qJD(2) * t71) * t99 + (mrSges(3,1) * qJDD(2) - mrSges(3,2) * t103 - t37) * t102) * t95 + m(3) * (qJDD(1) * t96 ^ 2 + (t102 * t54 + t55 * t99) * t95) + m(4) * (t10 * t63 - t11 * t62 - t33 * t45 + t34 * t46 + (-t102 * t48 + t162 * t77) * t95) + m(5) * (-t1 * t116 + t12 * t6 + t13 * t7 + t2 * t35 + t33 * t38 + t62 * t9) + (-m(2) - m(3) - t213) * g(3); t208 * t43 + (t1 * t51 + t2 * t50 + (t154 * t38 + t188) * pkin(6) - t129 * t38 + t208 * t13 + t209 * t12) * m(5) + t209 * t44 + t212 * pkin(6) * t98 + t205 * t154 / 0.2e1 - t3 * t180 / 0.2e1 + t4 * t173 / 0.2e1 + t69 * (qJD(3) * t106 - t121 * t159) / 0.2e1 + t91 * (qJD(3) * t105 - t118 * t159) / 0.2e1 + t126 * t188 + (qJD(3) * t107 - t124 * t159) * t192 + t138 * t197 + t71 * t145 + (t198 * t59 + t216 * t58) * g(2) + (t198 * t61 + t216 * t60) * g(1) + (-t213 * (pkin(2) * t169 + pkin(6) * t181) + (-t200 * t102 + t217 * t99) * t95) * g(3) + (-t1 * t180 - t12 * t108 - t13 * t109 - t2 * t173) * mrSges(5,3) + (t12 * mrSges(5,1) - t13 * mrSges(5,2) - t46 * mrSges(4,3) - pkin(6) * t80 - t201) * t161 + (-pkin(2) * t48 + ((-t101 * t45 - t46 * t98) * qJD(3) + t204) * pkin(6) - (t77 * t99 + (t101 * t46 - t45 * t98) * t102) * t166) * m(4) + (-t154 * t45 + t204) * mrSges(4,3) + qJD(3) ^ 2 * t119 / 0.2e1 + qJD(3) * t113 + t38 * (mrSges(5,1) * t109 + mrSges(5,2) * t108) + qJDD(3) * (Ifges(4,5) * t98 + Ifges(4,6) * t101) + (pkin(6) * t154 - t129) * t210 + (-Ifges(5,3) * t194 - Ifges(5,6) * t195 - Ifges(5,5) * t196 - t148 / 0.2e1 + (-Ifges(4,2) * t98 + t172) * t130 + Ifges(4,4) * t214 + Ifges(4,2) * t215 - t80 * t140 + pkin(6) * t56 + t199) * t101 + (Ifges(4,1) * t75 + Ifges(4,4) * t215 + t117 * t194 + t120 * t195 + t123 * t196) * t98 + t50 * t14 + t51 * t15 - pkin(2) * t37 + (t54 + t84) * mrSges(3,1) - t48 * t127 + (-t55 + t85) * mrSges(3,2) + t112 * t130 - (t100 * t22 + t97 * t23) * t159 / 0.2e1 + t172 * t214 + t122 * t215 + Ifges(3,3) * qJDD(2); (-t13 * (-mrSges(5,2) * t98 - mrSges(5,3) * t170) - t12 * (mrSges(5,1) * t98 - mrSges(5,3) * t158) - t113) * qJD(2) + (-m(5) * t38 + t147 - t210) * t46 + (t206 * t32 - t207 * (t171 * t94 - t61 * t98)) * g(1) - t103 * t112 / 0.2e1 - (-Ifges(4,2) * t163 + t205 + t93) * t156 / 0.2e1 + (t206 * t63 + t207 * t62) * g(3) + (t206 * t30 - t207 * (-t101 * t135 - t59 * t98)) * g(2) + t91 * t38 * t126 + t118 * t194 + t121 * t195 + t124 * t196 + t160 * t197 + t23 * t153 / 0.2e1 - t119 * t149 / 0.2e1 + t201 * t163 + (m(5) * ((-t100 * t12 - t13 * t97) * qJD(4) + t203) - t43 * t160 - t44 * t153 + t100 * t15 - t97 * t14) * pkin(7) + (-t12 * t153 - t13 * t160 + t203) * mrSges(5,3) + t100 * t3 / 0.2e1 + Ifges(4,5) * t75 + Ifges(4,6) * t74 - t20 * t43 - t19 * t44 - pkin(3) * t5 - t10 * mrSges(4,2) + t11 * mrSges(4,1) + t9 * t125 + (-t80 + t142) * t45 + (t117 * t91 + t120 * t69 + t123 * t70) * qJD(4) / 0.2e1 - (t91 * t105 + t69 * t106 + t107 * t70) * qJD(2) / 0.2e1 + (t156 * t22 + t4) * t97 / 0.2e1 + (-pkin(3) * t9 - t12 * t19 - t13 * t20) * m(5) + Ifges(4,3) * qJDD(3); -t38 * (mrSges(5,1) * t70 + mrSges(5,2) * t69) + (Ifges(5,1) * t69 - t186) * t219 + t22 * t192 + (Ifges(5,5) * t69 - Ifges(5,6) * t70) * t218 - t12 * t43 + t13 * t44 - g(1) * ((t100 * t60 - t32 * t97) * mrSges(5,1) + (-t100 * t32 - t60 * t97) * mrSges(5,2)) - g(2) * ((t100 * t58 - t30 * t97) * mrSges(5,1) + (-t100 * t30 - t58 * t97) * mrSges(5,2)) - g(3) * (mrSges(5,1) * t35 + mrSges(5,2) * t116) + (t12 * t69 + t13 * t70) * mrSges(5,3) + t148 + (-Ifges(5,2) * t70 + t23 + t67) * t220 - t199;];
tau = t17;
