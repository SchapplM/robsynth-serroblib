% Calculate joint inertia matrix for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:44:08
% EndTime: 2019-03-09 06:44:13
% DurationCPUTime: 2.11s
% Computational Cost: add. (4541->402), mult. (12017->562), div. (0->0), fcn. (13519->12), ass. (0->153)
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t209 = t147 ^ 2 + t150 ^ 2;
t142 = sin(pkin(7));
t152 = cos(qJ(3));
t182 = t142 * t152;
t145 = cos(pkin(7));
t148 = sin(qJ(4));
t151 = cos(qJ(4));
t149 = sin(qJ(3));
t183 = t142 * t149;
t93 = t145 * t148 + t151 * t183;
t69 = t147 * t93 + t150 * t182;
t71 = -t147 * t182 + t150 * t93;
t208 = t147 * t69 + t150 * t71;
t207 = 2 * pkin(10);
t206 = mrSges(6,3) + mrSges(7,2);
t141 = sin(pkin(12));
t143 = sin(pkin(6));
t144 = cos(pkin(12));
t146 = cos(pkin(6));
t179 = t145 * t149;
t65 = t146 * t183 + (t141 * t152 + t144 * t179) * t143;
t180 = t143 * t144;
t87 = -t142 * t180 + t145 * t146;
t52 = t148 * t87 + t151 * t65;
t167 = t145 * t180;
t181 = t143 * t141;
t64 = -t146 * t182 + t149 * t181 - t152 * t167;
t41 = t147 * t52 - t64 * t150;
t42 = t147 * t64 + t150 * t52;
t51 = t148 * t65 - t87 * t151;
t12 = Ifges(6,5) * t42 - Ifges(6,6) * t41 + Ifges(6,3) * t51;
t13 = Ifges(7,4) * t42 + Ifges(7,2) * t51 + Ifges(7,6) * t41;
t205 = t12 + t13;
t204 = -m(7) * pkin(5) - mrSges(7,1);
t201 = pkin(1) * t146;
t90 = qJ(2) * t180 + t141 * t201;
t59 = (t142 * t146 + t167) * pkin(9) + t90;
t119 = t144 * t201;
t66 = pkin(2) * t146 + t119 + (-pkin(9) * t145 - qJ(2)) * t181;
t75 = (-pkin(9) * t141 * t142 - pkin(2) * t144 - pkin(1)) * t143;
t34 = -t149 * t59 + (t142 * t75 + t145 * t66) * t152;
t91 = -t151 * t145 + t148 * t183;
t88 = t91 ^ 2;
t203 = 2 * mrSges(3,1);
t202 = 0.2e1 * t146;
t31 = -pkin(3) * t87 - t34;
t18 = pkin(4) * t51 - pkin(11) * t52 + t31;
t46 = -t142 * t66 + t145 * t75;
t28 = pkin(3) * t64 - pkin(10) * t65 + t46;
t35 = t152 * t59 + t66 * t179 + t75 * t183;
t32 = pkin(10) * t87 + t35;
t10 = t148 * t28 + t151 * t32;
t8 = pkin(11) * t64 + t10;
t4 = t147 * t18 + t150 * t8;
t200 = pkin(10) * t151;
t199 = -Ifges(6,3) - Ifges(7,2);
t21 = -mrSges(7,2) * t41 + mrSges(7,3) * t51;
t22 = -mrSges(6,2) * t51 - mrSges(6,3) * t41;
t198 = t21 + t22;
t23 = mrSges(6,1) * t51 - mrSges(6,3) * t42;
t24 = -t51 * mrSges(7,1) + t42 * mrSges(7,2);
t197 = -t23 + t24;
t20 = mrSges(6,1) * t41 + mrSges(6,2) * t42;
t44 = mrSges(5,1) * t64 - mrSges(5,3) * t52;
t196 = -t44 + t20;
t195 = Ifges(6,4) * t147;
t194 = Ifges(6,4) * t150;
t193 = Ifges(7,5) * t147;
t192 = Ifges(7,5) * t150;
t190 = t148 * t91;
t188 = t151 * t93;
t177 = t148 * t150;
t100 = t151 * mrSges(7,1) + mrSges(7,2) * t177;
t99 = -mrSges(6,1) * t151 - mrSges(6,3) * t177;
t187 = t100 - t99;
t178 = t147 * t148;
t101 = -mrSges(7,2) * t178 - mrSges(7,3) * t151;
t98 = mrSges(6,2) * t151 - mrSges(6,3) * t178;
t186 = t101 + t98;
t105 = -mrSges(6,1) * t150 + mrSges(6,2) * t147;
t185 = t105 - mrSges(5,1);
t103 = -pkin(4) * t151 - pkin(11) * t148 - pkin(3);
t77 = t147 * t103 + t150 * t200;
t184 = t103 * t150;
t176 = Ifges(7,4) * t177 + Ifges(7,6) * t178;
t108 = Ifges(6,5) * t147 + Ifges(6,6) * t150;
t175 = Ifges(5,5) * t148 + Ifges(5,6) * t151;
t174 = t209 * pkin(11) ^ 2;
t173 = Ifges(5,5) * t52 - Ifges(5,6) * t51 + Ifges(5,3) * t64;
t172 = Ifges(4,5) * t65 - Ifges(4,6) * t64 + Ifges(4,3) * t87;
t11 = Ifges(7,5) * t42 + Ifges(7,6) * t51 + Ifges(7,3) * t41;
t14 = Ifges(6,4) * t42 - Ifges(6,2) * t41 + Ifges(6,6) * t51;
t171 = t11 / 0.2e1 - t14 / 0.2e1;
t15 = Ifges(7,1) * t42 + Ifges(7,4) * t51 + Ifges(7,5) * t41;
t16 = Ifges(6,1) * t42 - Ifges(6,4) * t41 + Ifges(6,5) * t51;
t170 = t15 / 0.2e1 + t16 / 0.2e1;
t79 = -Ifges(7,6) * t151 + (Ifges(7,3) * t147 + t192) * t148;
t82 = -Ifges(6,6) * t151 + (-Ifges(6,2) * t147 + t194) * t148;
t169 = t79 / 0.2e1 - t82 / 0.2e1;
t83 = -Ifges(7,4) * t151 + (Ifges(7,1) * t150 + t193) * t148;
t84 = -Ifges(6,5) * t151 + (Ifges(6,1) * t150 - t195) * t148;
t168 = t83 / 0.2e1 + t84 / 0.2e1;
t107 = -Ifges(7,3) * t150 + t193;
t110 = Ifges(6,2) * t150 + t195;
t166 = t107 / 0.2e1 - t110 / 0.2e1;
t109 = Ifges(7,4) * t147 - Ifges(7,6) * t150;
t165 = t108 / 0.2e1 + t109 / 0.2e1;
t112 = Ifges(7,1) * t147 - t192;
t113 = Ifges(6,1) * t147 + t194;
t164 = t112 / 0.2e1 + t113 / 0.2e1;
t9 = -t148 * t32 + t151 * t28;
t162 = t208 * pkin(11);
t160 = Ifges(6,5) * t177 - Ifges(6,6) * t178;
t3 = -t147 * t8 + t150 * t18;
t159 = t147 * mrSges(6,1) + t150 * mrSges(6,2);
t158 = t147 * mrSges(7,1) - t150 * mrSges(7,3);
t157 = -pkin(5) * t147 + qJ(6) * t150;
t7 = -pkin(4) * t64 - t9;
t154 = pkin(10) ^ 2;
t140 = t151 ^ 2;
t138 = t148 ^ 2;
t136 = t142 ^ 2;
t134 = t138 * t154;
t125 = t136 * t152 ^ 2;
t117 = mrSges(3,2) * t181;
t114 = Ifges(5,1) * t148 + Ifges(5,4) * t151;
t111 = Ifges(5,4) * t148 + Ifges(5,2) * t151;
t106 = -mrSges(5,1) * t151 + mrSges(5,2) * t148;
t104 = -mrSges(7,1) * t150 - mrSges(7,3) * t147;
t102 = -pkin(5) * t150 - qJ(6) * t147 - pkin(4);
t95 = t159 * t148;
t94 = t158 * t148;
t89 = -qJ(2) * t181 + t119;
t85 = (pkin(10) - t157) * t148;
t81 = -Ifges(7,2) * t151 + t176;
t80 = -Ifges(6,3) * t151 + t160;
t76 = -t147 * t200 + t184;
t74 = -t184 + (pkin(10) * t147 + pkin(5)) * t151;
t73 = -qJ(6) * t151 + t77;
t54 = mrSges(4,1) * t87 - mrSges(4,3) * t65;
t53 = -mrSges(4,2) * t87 - mrSges(4,3) * t64;
t45 = mrSges(4,1) * t64 + mrSges(4,2) * t65;
t43 = -mrSges(5,2) * t64 - mrSges(5,3) * t51;
t33 = mrSges(5,1) * t51 + mrSges(5,2) * t52;
t26 = Ifges(5,1) * t52 - Ifges(5,4) * t51 + Ifges(5,5) * t64;
t25 = Ifges(5,4) * t52 - Ifges(5,2) * t51 + Ifges(5,6) * t64;
t19 = mrSges(7,1) * t41 - mrSges(7,3) * t42;
t5 = pkin(5) * t41 - qJ(6) * t42 + t7;
t2 = -pkin(5) * t51 - t3;
t1 = qJ(6) * t51 + t4;
t6 = [t65 * (Ifges(4,1) * t65 + Ifges(4,5) * t87) + t87 * t172 + (-0.2e1 * Ifges(4,4) * t65 + Ifges(4,2) * t64 - Ifges(4,6) * t87 + t173) * t64 + m(4) * (t34 ^ 2 + t35 ^ 2 + t46 ^ 2) + m(5) * (t10 ^ 2 + t31 ^ 2 + t9 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) + m(3) * (pkin(1) ^ 2 * t143 ^ 2 + t89 ^ 2 + t90 ^ 2) + (-0.2e1 * t90 * mrSges(3,2) + Ifges(3,3) * t146 + t203 * t89) * t146 + (-0.2e1 * pkin(1) * t117 + (-0.2e1 * t89 * mrSges(3,3) + Ifges(3,1) * t181 + Ifges(3,5) * t202) * t141 + (0.2e1 * t90 * mrSges(3,3) + Ifges(3,6) * t202 + (0.2e1 * Ifges(3,4) * t141 + Ifges(3,2) * t144 + pkin(1) * t203) * t143) * t144) * t143 + (-t25 + t205) * t51 + Ifges(2,3) + 0.2e1 * t5 * t19 + 0.2e1 * t7 * t20 + 0.2e1 * t1 * t21 + 0.2e1 * t4 * t22 + 0.2e1 * t3 * t23 + 0.2e1 * t2 * t24 + 0.2e1 * t31 * t33 + 0.2e1 * t10 * t43 + 0.2e1 * t9 * t44 + 0.2e1 * t46 * t45 + t52 * t26 + 0.2e1 * t35 * t53 + 0.2e1 * t34 * t54 + (t11 - t14) * t41 + (t15 + t16) * t42; t145 * t45 + t93 * t43 + t117 + t198 * t71 + t197 * t69 + (-m(3) * pkin(1) - mrSges(3,1) * t144) * t143 + (t19 + t196) * t91 + (t149 * t53 + (-t33 + t54) * t152) * t142 + m(7) * (t1 * t71 + t2 * t69 + t5 * t91) + m(6) * (-t3 * t69 + t4 * t71 + t7 * t91) + m(5) * (t10 * t93 - t182 * t31 - t9 * t91) + m(4) * (t145 * t46 + (t149 * t35 + t152 * t34) * t142); m(3) + m(5) * (t93 ^ 2 + t125 + t88) + m(4) * (t136 * t149 ^ 2 + t145 ^ 2 + t125) + (m(6) + m(7)) * (t69 ^ 2 + t71 ^ 2 + t88); t64 * t175 / 0.2e1 + (pkin(10) * t43 + t10 * mrSges(5,3) - t12 / 0.2e1 - t13 / 0.2e1 + t25 / 0.2e1) * t151 + t168 * t42 + t169 * t41 + m(7) * (t1 * t73 + t2 * t74 + t5 * t85) + m(6) * (pkin(10) * t148 * t7 + t3 * t76 + t4 * t77) + (t80 / 0.2e1 + t81 / 0.2e1 - t111 / 0.2e1) * t51 + t172 + m(5) * (-pkin(3) * t31 + (t10 * t151 - t9 * t148) * pkin(10)) + (-t9 * mrSges(5,3) + t26 / 0.2e1 + t170 * t150 + t171 * t147 + t196 * pkin(10)) * t148 - pkin(3) * t33 + t34 * mrSges(4,1) - t35 * mrSges(4,2) + t73 * t21 + t74 * t24 + t76 * t23 + t77 * t22 + t85 * t19 + t5 * t94 + t7 * t95 + t4 * t98 + t3 * t99 + t2 * t100 + t1 * t101 + t31 * t106 + t52 * t114 / 0.2e1; mrSges(5,3) * t188 + t186 * t71 + t187 * t69 + (t148 * mrSges(5,3) + t94 + t95) * t91 + (-t149 * mrSges(4,2) + (mrSges(4,1) - t106) * t152) * t142 + m(6) * (pkin(10) * t190 - t69 * t76 + t71 * t77) + m(7) * (t69 * t74 + t71 * t73 + t85 * t91) + m(5) * (pkin(3) * t182 + (t188 + t190) * pkin(10)); -0.2e1 * pkin(3) * t106 + 0.2e1 * t74 * t100 + 0.2e1 * t73 * t101 + 0.2e1 * t76 * t99 + 0.2e1 * t77 * t98 + 0.2e1 * t85 * t94 + Ifges(4,3) + (t138 + t140) * mrSges(5,3) * t207 + (-t80 - t81 + t111) * t151 + m(7) * (t73 ^ 2 + t74 ^ 2 + t85 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2 + t134) + m(5) * (pkin(3) ^ 2 + t140 * t154 + t134) + (t95 * t207 + t114 + (t83 + t84) * t150 + (t79 - t82) * t147) * t148; t9 * mrSges(5,1) - t10 * mrSges(5,2) - pkin(4) * t20 + t102 * t19 + t5 * t104 + t7 * t105 + t165 * t51 + t164 * t42 + t166 * t41 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) + pkin(11) * t198 - t171) * t150 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + pkin(11) * t197 + t170) * t147 + m(6) * (-pkin(4) * t7 + (-t147 * t3 + t150 * t4) * pkin(11)) + m(7) * (t102 * t5 + (t1 * t150 + t147 * t2) * pkin(11)) + t173; -t93 * mrSges(5,2) + (t104 + t185) * t91 + m(6) * (-pkin(4) * t91 + t162) + m(7) * (t102 * t91 + t162) + t206 * t208; -pkin(4) * t95 + t85 * t104 + (m(7) * t85 + t94) * t102 + (-pkin(10) * mrSges(5,2) - t165) * t151 + (t73 * mrSges(7,2) + t77 * mrSges(6,3) - t169) * t150 + (t74 * mrSges(7,2) - t76 * mrSges(6,3) + t168) * t147 + (t186 * t150 + t187 * t147 + m(7) * (t147 * t74 + t150 * t73) + m(6) * (-t147 * t76 + t150 * t77)) * pkin(11) + (t164 * t150 + t166 * t147 + (-m(6) * pkin(4) + t185) * pkin(10)) * t148 + t175; -0.2e1 * pkin(4) * t105 + 0.2e1 * t102 * t104 + Ifges(5,3) + (-t107 + t110) * t150 + (t112 + t113) * t147 + m(7) * (t102 ^ 2 + t174) + m(6) * (pkin(4) ^ 2 + t174) + 0.2e1 * t206 * pkin(11) * t209; -pkin(5) * t24 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t21 + t1 * mrSges(7,3) + t3 * mrSges(6,1) - t2 * mrSges(7,1) - t4 * mrSges(6,2) + t205; (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t71 + (-mrSges(6,1) + t204) * t69; -pkin(5) * t100 + m(7) * (-pkin(5) * t74 + qJ(6) * t73) + qJ(6) * t101 + t73 * mrSges(7,3) - t74 * mrSges(7,1) - t77 * mrSges(6,2) + t76 * mrSges(6,1) + t199 * t151 + t160 + t176; t157 * mrSges(7,2) + (m(7) * t157 - t158 - t159) * pkin(11) + t109 + t108; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t199; m(7) * t2 + t24; m(7) * t69; m(7) * t74 + t100; (m(7) * pkin(11) + mrSges(7,2)) * t147; t204; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
