% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:29
% EndTime: 2019-12-05 16:21:40
% DurationCPUTime: 3.43s
% Computational Cost: add. (2507->311), mult. (6566->451), div. (0->0), fcn. (4617->8), ass. (0->157)
t140 = sin(qJ(3));
t180 = -qJ(4) - pkin(6);
t157 = qJD(3) * t180;
t143 = cos(qJ(3));
t165 = qJD(4) * t143;
t103 = t140 * t157 + t165;
t104 = -qJD(4) * t140 + t143 * t157;
t137 = sin(pkin(9));
t138 = cos(pkin(9));
t115 = t137 * t143 + t138 * t140;
t144 = cos(qJ(2));
t146 = t115 * t144;
t201 = qJD(1) * t146 - t103 * t137 + t104 * t138;
t147 = t137 * t140 - t138 * t143;
t171 = qJD(1) * t144;
t200 = t103 * t138 + t104 * t137 + t147 * t171;
t108 = t147 * qJD(3);
t213 = pkin(7) * t108 + t201;
t107 = t115 * qJD(3);
t212 = -pkin(7) * t107 + t200;
t134 = qJD(3) + qJD(5);
t105 = t147 * qJD(2);
t168 = qJD(2) * t143;
t170 = qJD(2) * t140;
t106 = -t137 * t168 - t138 * t170;
t139 = sin(qJ(5));
t142 = cos(qJ(5));
t156 = -t105 * t142 + t106 * t139;
t52 = Ifges(6,4) * t156;
t61 = -t105 * t139 - t106 * t142;
t18 = Ifges(6,1) * t61 + Ifges(6,5) * t134 + t52;
t184 = Ifges(6,4) * t61;
t167 = qJD(3) * t140;
t141 = sin(qJ(2));
t172 = qJD(1) * t141;
t125 = qJD(2) * pkin(6) + t172;
t164 = qJD(1) * qJD(2);
t158 = t144 * t164;
t75 = -t125 * t167 + t143 * t158;
t63 = (-qJ(4) * t167 + t165) * qJD(2) + t75;
t166 = qJD(3) * t143;
t161 = t125 * t166;
t64 = -t161 + (-qJ(4) * t166 + (-qJD(4) - t171) * t140) * qJD(2);
t25 = -t137 * t63 + t138 * t64;
t95 = qJD(2) * t108;
t13 = pkin(7) * t95 + t25;
t26 = t137 * t64 + t138 * t63;
t94 = qJD(2) * t107;
t14 = -pkin(7) * t94 + t26;
t181 = pkin(7) * t106;
t154 = qJ(4) * qJD(2) + t125;
t99 = t154 * t143;
t77 = t137 * t99;
t98 = t154 * t140;
t83 = qJD(3) * pkin(3) - t98;
t37 = t138 * t83 - t77;
t29 = qJD(3) * pkin(4) + t181 + t37;
t182 = pkin(7) * t105;
t177 = t138 * t99;
t38 = t137 * t83 + t177;
t30 = t38 - t182;
t7 = -t139 * t30 + t142 * t29;
t2 = qJD(5) * t7 + t13 * t139 + t14 * t142;
t22 = qJD(5) * t156 - t139 * t94 - t142 * t95;
t23 = -qJD(5) * t61 + t139 * t95 - t142 * t94;
t8 = t139 * t29 + t142 * t30;
t3 = -qJD(5) * t8 + t13 * t142 - t139 * t14;
t132 = -pkin(3) * t143 - pkin(2);
t109 = qJD(2) * t132 + qJD(4) - t171;
t68 = pkin(4) * t105 + t109;
t211 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t22 + Ifges(6,6) * t23 - (Ifges(6,5) * t156 - Ifges(6,6) * t61) * t134 / 0.2e1 - (-Ifges(6,2) * t61 + t18 + t52) * t156 / 0.2e1 - t68 * (mrSges(6,1) * t61 + mrSges(6,2) * t156) - (Ifges(6,1) * t156 - t184) * t61 / 0.2e1;
t210 = t156 * t7 + t61 * t8;
t208 = -Ifges(4,1) / 0.2e1;
t17 = Ifges(6,2) * t156 + Ifges(6,6) * t134 + t184;
t207 = t17 / 0.2e1;
t206 = -Ifges(4,4) * t168 / 0.2e1;
t122 = t180 * t140;
t123 = t180 * t143;
t69 = t122 * t138 + t123 * t137;
t50 = -pkin(7) * t115 + t69;
t70 = t122 * t137 - t123 * t138;
t51 = -pkin(7) * t147 + t70;
t16 = t139 * t50 + t142 * t51;
t205 = -qJD(5) * t16 - t139 * t212 + t142 * t213;
t15 = -t139 * t51 + t142 * t50;
t204 = qJD(5) * t15 + t139 * t213 + t142 * t212;
t49 = t94 * mrSges(5,1) - mrSges(5,2) * t95;
t6 = -t23 * mrSges(6,1) + mrSges(6,2) * t22;
t203 = t49 + t6;
t24 = -mrSges(6,1) * t156 + mrSges(6,2) * t61;
t62 = mrSges(5,1) * t105 - mrSges(5,2) * t106;
t199 = -t62 - t24;
t130 = pkin(3) * t138 + pkin(4);
t183 = pkin(3) * t137;
t101 = t130 * t142 - t139 * t183;
t41 = t137 * t98 - t177;
t31 = t41 + t182;
t42 = -t138 * t98 - t77;
t32 = t42 + t181;
t198 = qJD(5) * t101 - t139 * t31 - t142 * t32;
t102 = t130 * t139 + t142 * t183;
t197 = -qJD(5) * t102 + t139 * t32 - t142 * t31;
t135 = t140 ^ 2;
t136 = t143 ^ 2;
t192 = t156 / 0.2e1;
t190 = t61 / 0.2e1;
t188 = -t106 / 0.2e1;
t187 = -t107 / 0.2e1;
t186 = -t108 / 0.2e1;
t179 = Ifges(4,4) * t140;
t178 = Ifges(5,4) * t106;
t176 = Ifges(4,5) * qJD(3);
t175 = Ifges(4,6) * qJD(3);
t126 = -qJD(2) * pkin(2) - t171;
t173 = t126 * t141;
t131 = t141 * t164;
t162 = pkin(3) * t167;
t113 = qJD(2) * t162 + t131;
t169 = qJD(2) * t141;
t163 = qJD(2) * qJD(3);
t160 = t176 / 0.2e1;
t159 = -t175 / 0.2e1;
t155 = (t135 + t136) * t125;
t153 = t175 / 0.2e1 + (t143 * Ifges(4,2) + t179) * qJD(2) / 0.2e1 - t126 * mrSges(4,1);
t152 = t170 * t208 + t206 - t176 / 0.2e1 - t126 * mrSges(4,2);
t96 = t115 * t141;
t97 = t147 * t141;
t43 = t139 * t97 - t142 * t96;
t44 = -t139 * t96 - t142 * t97;
t76 = -t140 * t158 - t161;
t150 = -t140 * t76 + t143 * t75;
t66 = -t115 * t139 - t142 * t147;
t67 = t115 * t142 - t139 * t147;
t120 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t170;
t121 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t168;
t149 = t120 * t143 + t121 * t140;
t148 = t120 * t140 - t121 * t143;
t145 = qJD(2) ^ 2;
t112 = (mrSges(4,1) * t140 + mrSges(4,2) * t143) * t163;
t100 = Ifges(5,4) * t105;
t86 = pkin(4) * t147 + t132;
t82 = qJD(3) * mrSges(5,1) + mrSges(5,3) * t106;
t81 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t105;
t74 = pkin(4) * t107 + t162;
t73 = pkin(3) * t170 - pkin(4) * t106;
t65 = pkin(4) * t94 + t113;
t57 = -t106 * Ifges(5,1) + Ifges(5,5) * qJD(3) - t100;
t56 = -t105 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t178;
t46 = -t105 * t144 - t107 * t141;
t45 = -qJD(2) * t146 + t108 * t141;
t40 = mrSges(6,1) * t134 - mrSges(6,3) * t61;
t39 = -mrSges(6,2) * t134 + mrSges(6,3) * t156;
t28 = -qJD(5) * t67 - t107 * t142 + t108 * t139;
t27 = qJD(5) * t66 - t107 * t139 - t108 * t142;
t10 = -qJD(5) * t44 - t139 * t46 + t142 * t45;
t9 = qJD(5) * t43 + t139 * t45 + t142 * t46;
t1 = [t10 * t40 + t9 * t39 + t45 * t82 + t46 * t81 + (-t22 * t43 + t23 * t44) * mrSges(6,3) + (t94 * t97 - t95 * t96) * mrSges(5,3) + (-t145 * mrSges(3,2) - t148 * qJD(2) - t112 - t203) * t144 + (-t145 * mrSges(3,1) - t149 * qJD(3) + (qJD(2) * (-mrSges(4,1) * t143 + mrSges(4,2) * t140) - t199) * qJD(2)) * t141 + m(6) * (t10 * t7 - t144 * t65 + t169 * t68 + t2 * t44 + t3 * t43 + t8 * t9) + m(5) * (t109 * t169 - t113 * t144 - t25 * t96 - t26 * t97 + t37 * t45 + t38 * t46) + m(4) * (t150 * t141 + (t173 + (t155 - t172) * t144) * qJD(2)); t109 * (mrSges(5,1) * t107 - mrSges(5,2) * t108) + (-t108 * t188 - t115 * t95) * Ifges(5,1) + (-t105 * t186 - t107 * t188 - t115 * t94 + t147 * t95) * Ifges(5,4) + (-t107 * t38 + t108 * t37 - t115 * t25 - t147 * t26 + t69 * t95 - t70 * t94) * mrSges(5,3) + t113 * (mrSges(5,1) * t147 + mrSges(5,2) * t115) + (-t105 * t187 + t147 * t94) * Ifges(5,2) + t134 * (Ifges(6,5) * t27 + Ifges(6,6) * t28) / 0.2e1 + t132 * t49 - pkin(2) * t112 + t74 * t24 + t86 * t6 + t150 * mrSges(4,3) + t65 * (-mrSges(6,1) * t66 + mrSges(6,2) * t67) + t68 * (-mrSges(6,1) * t28 + mrSges(6,2) * t27) + t27 * t18 / 0.2e1 + t204 * t39 + (t15 * t3 + t16 * t2 + t65 * t86 + t204 * t8 + t205 * t7 + (-t172 + t74) * t68) * m(6) + t205 * t40 + (t199 * t141 + t148 * t144) * qJD(1) + t200 * t81 + (t113 * t132 + t25 * t69 + t26 * t70 + t200 * t38 + t201 * t37 + (t162 - t172) * t109) * m(5) + t201 * t82 + (-t15 * t22 + t16 * t23 + t2 * t66 - t27 * t7 + t28 * t8 - t3 * t67) * mrSges(6,3) + t28 * t207 + (-0.3e1 / 0.2e1 * t135 + 0.3e1 / 0.2e1 * t136) * Ifges(4,4) * t163 + t57 * t186 + t56 * t187 + (Ifges(5,5) * t186 + Ifges(5,6) * t187 + (-pkin(6) * t120 - t152 + t160) * t143 + (-pkin(6) * t121 + pkin(3) * t62 + t159 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t168 - t153) * t140) * qJD(3) + (t190 * t27 + t22 * t67) * Ifges(6,1) + (t192 * t28 + t23 * t66) * Ifges(6,2) + (t190 * t28 + t192 * t27 + t22 * t66 + t23 * t67) * Ifges(6,4) + (-(t144 * t155 + t173) * qJD(1) - pkin(2) * t131 + pkin(6) * t150) * m(4); -t109 * (-mrSges(5,1) * t106 - mrSges(5,2) * t105) - qJD(3) * (-Ifges(5,5) * t105 + Ifges(5,6) * t106) / 0.2e1 + t106 * (-Ifges(5,1) * t105 + t178) / 0.2e1 + m(5) * (t137 * t26 + t138 * t25) * pkin(3) + (-t37 * t105 - t38 * t106 + (-t137 * t94 + t138 * t95) * pkin(3)) * mrSges(5,3) + (Ifges(5,2) * t106 - t100 + t57) * t105 / 0.2e1 + t61 * t207 - Ifges(5,6) * t94 - Ifges(5,5) * t95 - t73 * t24 - t75 * mrSges(4,2) + t76 * mrSges(4,1) - t42 * t81 - t41 * t82 + t149 * t125 + (-t101 * t22 + t102 * t23 + t210) * mrSges(6,3) + t25 * mrSges(5,1) - t26 * mrSges(5,2) + ((t160 + t206 + t152) * t143 + (t159 + (t179 / 0.2e1 + (Ifges(4,2) / 0.2e1 + t208) * t143) * qJD(2) + (-m(5) * t109 - t62) * pkin(3) + t153) * t140) * qJD(2) + t197 * t40 + (t101 * t3 + t102 * t2 + t197 * t7 + t198 * t8 - t68 * t73) * m(6) + t198 * t39 - m(5) * (t37 * t41 + t38 * t42) + t211 + t56 * t188; t105 * t81 - t106 * t82 - t156 * t39 + t61 * t40 + (-t156 * t8 + t61 * t7 + t65) * m(6) + (t105 * t38 - t106 * t37 + t113) * m(5) + t203; t210 * mrSges(6,3) + t17 * t190 - t7 * t39 + t8 * t40 + t211;];
tauc = t1(:);
