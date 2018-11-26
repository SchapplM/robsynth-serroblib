% Calculate joint inertia matrix for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:41
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_inertiaJ_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:40:50
% EndTime: 2018-11-23 16:40:52
% DurationCPUTime: 2.78s
% Computational Cost: add. (10931->493), mult. (30212->742), div. (0->0), fcn. (35330->16), ass. (0->190)
t248 = 2 * pkin(12);
t177 = sin(pkin(8));
t181 = cos(pkin(8));
t186 = sin(qJ(4));
t190 = cos(qJ(4));
t178 = sin(pkin(7));
t182 = cos(pkin(7));
t183 = cos(pkin(6));
t179 = sin(pkin(6));
t180 = cos(pkin(14));
t214 = t179 * t180;
t129 = -t178 * t214 + t182 * t183;
t176 = sin(pkin(14));
t187 = sin(qJ(3));
t191 = cos(qJ(3));
t209 = t182 * t191;
t216 = t178 * t191;
t94 = t183 * t216 + (-t176 * t187 + t180 * t209) * t179;
t195 = t129 * t177 + t181 * t94;
t109 = (-pkin(10) * t176 * t178 - pkin(2) * t180 - pkin(1)) * t179;
t210 = t182 * t187;
t217 = t178 * t187;
t230 = pkin(1) * t183;
t132 = qJ(2) * t214 + t176 * t230;
t91 = (t178 * t183 + t182 * t214) * pkin(10) + t132;
t157 = t180 * t230;
t215 = t179 * t176;
t96 = pkin(2) * t183 + t157 + (-pkin(10) * t182 - qJ(2)) * t215;
t54 = t109 * t217 + t191 * t91 + t210 * t96;
t39 = pkin(11) * t195 + t54;
t95 = t183 * t217 + (t176 * t191 + t180 * t210) * t179;
t240 = pkin(11) * t95;
t53 = t109 * t216 - t187 * t91 + t209 * t96;
t45 = pkin(3) * t129 - t181 * t240 + t53;
t66 = t109 * t182 - t178 * t96;
t50 = -pkin(3) * t94 - t177 * t240 + t66;
t17 = -t186 * t39 + (t177 * t50 + t181 * t45) * t190;
t184 = sin(qJ(6));
t188 = cos(qJ(6));
t145 = -mrSges(7,1) * t188 + mrSges(7,2) * t184;
t247 = -m(7) * pkin(5) - mrSges(6,1) + t145;
t130 = -t177 * t216 + t181 * t182;
t185 = sin(qJ(5));
t189 = cos(qJ(5));
t211 = t181 * t191;
t219 = t177 * t186;
t99 = t182 * t219 + (t186 * t211 + t187 * t190) * t178;
t71 = -t130 * t189 + t185 * t99;
t246 = t71 ^ 2;
t218 = t177 * t190;
t97 = -t178 * t190 * t211 - t182 * t218 + t186 * t217;
t245 = t97 ^ 2;
t244 = 2 * mrSges(3,1);
t243 = 0.2e1 * t183;
t59 = t186 * t195 + t190 * t95;
t70 = t129 * t181 - t177 * t94;
t43 = t185 * t59 - t189 * t70;
t44 = t185 * t70 + t189 * t59;
t212 = t181 * t190;
t58 = -t129 * t218 + t186 * t95 - t212 * t94;
t23 = Ifges(6,1) * t44 - Ifges(6,4) * t43 + Ifges(6,5) * t58;
t242 = t23 / 0.2e1;
t133 = -t181 * t189 + t185 * t219;
t134 = t181 * t185 + t189 * t219;
t85 = Ifges(6,1) * t134 - Ifges(6,4) * t133 - Ifges(6,5) * t218;
t241 = t85 / 0.2e1;
t222 = Ifges(7,4) * t188;
t123 = -Ifges(7,6) * t189 + (-Ifges(7,2) * t184 + t222) * t185;
t239 = t123 / 0.2e1;
t223 = Ifges(7,4) * t184;
t124 = -Ifges(7,5) * t189 + (Ifges(7,1) * t188 - t223) * t185;
t238 = t124 / 0.2e1;
t147 = Ifges(7,5) * t184 + Ifges(7,6) * t188;
t237 = t147 / 0.2e1;
t149 = Ifges(7,2) * t188 + t223;
t236 = t149 / 0.2e1;
t151 = Ifges(7,1) * t184 + t222;
t235 = t151 / 0.2e1;
t152 = Ifges(6,1) * t185 + Ifges(6,4) * t189;
t234 = t152 / 0.2e1;
t233 = -t184 / 0.2e1;
t232 = t184 / 0.2e1;
t231 = t188 / 0.2e1;
t229 = pkin(3) * t177;
t228 = pkin(12) * t185;
t227 = pkin(12) * t189;
t226 = pkin(13) * t184;
t225 = pkin(13) * t188;
t28 = -t184 * t44 + t188 * t58;
t29 = t184 * t58 + t188 * t44;
t11 = -mrSges(7,1) * t28 + mrSges(7,2) * t29;
t31 = mrSges(6,1) * t58 - mrSges(6,3) * t44;
t224 = t11 - t31;
t213 = t181 * t186;
t18 = t190 * t39 + t213 * t45 + t219 * t50;
t14 = pkin(12) * t70 + t18;
t25 = -t177 * t45 + t181 * t50;
t16 = pkin(4) * t58 - pkin(12) * t59 + t25;
t6 = t14 * t189 + t16 * t185;
t221 = t185 * t71;
t108 = -mrSges(6,1) * t218 - mrSges(6,3) * t134;
t104 = -t134 * t184 - t188 * t218;
t105 = t134 * t188 - t184 * t218;
t67 = -mrSges(7,1) * t104 + mrSges(7,2) * t105;
t220 = -t108 + t67;
t208 = t184 * t185;
t207 = t185 * t188;
t137 = pkin(3) * t213 + pkin(11) * t218;
t120 = pkin(12) * t181 + t137;
t121 = (-pkin(4) * t190 - pkin(12) * t186 - pkin(3)) * t177;
t82 = t120 * t189 + t121 * t185;
t148 = Ifges(6,5) * t185 + Ifges(6,6) * t189;
t206 = t184 ^ 2 + t188 ^ 2;
t8 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t43;
t21 = Ifges(6,5) * t44 - Ifges(6,6) * t43 + Ifges(6,3) * t58;
t32 = Ifges(5,5) * t59 - Ifges(5,6) * t58 + Ifges(5,3) * t70;
t22 = Ifges(6,4) * t44 - Ifges(6,2) * t43 + Ifges(6,6) * t58;
t205 = t8 / 0.2e1 - t22 / 0.2e1;
t204 = Ifges(4,5) * t95 + Ifges(4,6) * t94 + Ifges(4,3) * t129;
t62 = Ifges(7,5) * t105 + Ifges(7,6) * t104 + Ifges(7,3) * t133;
t84 = Ifges(6,4) * t134 - Ifges(6,2) * t133 - Ifges(6,6) * t218;
t203 = t62 / 0.2e1 - t84 / 0.2e1;
t115 = Ifges(5,5) * t219 + Ifges(5,6) * t218 + Ifges(5,3) * t181;
t122 = Ifges(7,5) * t207 - Ifges(7,6) * t208 - Ifges(7,3) * t189;
t150 = Ifges(6,4) * t185 + Ifges(6,2) * t189;
t202 = t122 / 0.2e1 - t150 / 0.2e1;
t4 = pkin(13) * t58 + t6;
t13 = -pkin(4) * t70 - t17;
t7 = pkin(5) * t43 - pkin(13) * t44 + t13;
t1 = -t184 * t4 + t188 * t7;
t2 = t184 * t7 + t188 * t4;
t201 = -t1 * t184 + t188 * t2;
t200 = mrSges(7,1) * t184 + mrSges(7,2) * t188;
t5 = -t14 * t185 + t16 * t189;
t160 = pkin(11) * t219;
t119 = t160 + (-pkin(3) * t190 - pkin(4)) * t181;
t74 = pkin(5) * t133 - pkin(13) * t134 + t119;
t76 = -pkin(13) * t218 + t82;
t51 = -t184 * t76 + t188 * t74;
t52 = t184 * t74 + t188 * t76;
t198 = -t184 * t51 + t188 * t52;
t73 = t130 * t185 + t189 * t99;
t196 = t189 * t73 + t221;
t144 = -pkin(5) * t189 - pkin(13) * t185 - pkin(4);
t113 = t144 * t188 - t184 * t227;
t114 = t144 * t184 + t188 * t227;
t194 = -t113 * t184 + t114 * t188;
t81 = -t120 * t185 + t121 * t189;
t83 = Ifges(6,5) * t134 - Ifges(6,6) * t133 - Ifges(6,3) * t218;
t193 = pkin(12) ^ 2;
t175 = t189 ^ 2;
t173 = t185 ^ 2;
t170 = t173 * t193;
t155 = mrSges(3,2) * t215;
t146 = -mrSges(6,1) * t189 + mrSges(6,2) * t185;
t142 = -mrSges(7,1) * t189 - mrSges(7,3) * t207;
t141 = mrSges(7,2) * t189 - mrSges(7,3) * t208;
t140 = -mrSges(5,2) * t181 + mrSges(5,3) * t218;
t139 = mrSges(5,1) * t181 - mrSges(5,3) * t219;
t138 = t200 * t185;
t136 = pkin(3) * t212 - t160;
t135 = (-mrSges(5,1) * t190 + mrSges(5,2) * t186) * t177;
t131 = -qJ(2) * t215 + t157;
t117 = Ifges(5,5) * t181 + (Ifges(5,1) * t186 + Ifges(5,4) * t190) * t177;
t116 = Ifges(5,6) * t181 + (Ifges(5,4) * t186 + Ifges(5,2) * t190) * t177;
t107 = mrSges(6,2) * t218 - mrSges(6,3) * t133;
t90 = mrSges(6,1) * t133 + mrSges(6,2) * t134;
t80 = mrSges(7,1) * t133 - mrSges(7,3) * t105;
t79 = -mrSges(7,2) * t133 + mrSges(7,3) * t104;
t78 = mrSges(4,1) * t129 - mrSges(4,3) * t95;
t77 = -mrSges(4,2) * t129 + mrSges(4,3) * t94;
t75 = pkin(5) * t218 - t81;
t65 = -mrSges(4,1) * t94 + mrSges(4,2) * t95;
t64 = Ifges(7,1) * t105 + Ifges(7,4) * t104 + Ifges(7,5) * t133;
t63 = Ifges(7,4) * t105 + Ifges(7,2) * t104 + Ifges(7,6) * t133;
t61 = t184 * t97 + t188 * t73;
t60 = -t184 * t73 + t188 * t97;
t47 = mrSges(5,1) * t70 - mrSges(5,3) * t59;
t46 = -mrSges(5,2) * t70 - mrSges(5,3) * t58;
t35 = mrSges(5,1) * t58 + mrSges(5,2) * t59;
t34 = Ifges(5,1) * t59 - Ifges(5,4) * t58 + Ifges(5,5) * t70;
t33 = Ifges(5,4) * t59 - Ifges(5,2) * t58 + Ifges(5,6) * t70;
t30 = -mrSges(6,2) * t58 - mrSges(6,3) * t43;
t24 = mrSges(6,1) * t43 + mrSges(6,2) * t44;
t20 = mrSges(7,1) * t43 - mrSges(7,3) * t29;
t19 = -mrSges(7,2) * t43 + mrSges(7,3) * t28;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t43;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t43;
t3 = -pkin(5) * t58 - t5;
t12 = [(0.2e1 * Ifges(4,4) * t95 + Ifges(4,2) * t94 + Ifges(4,6) * t129) * t94 + (t21 - t33) * t58 + t95 * (Ifges(4,1) * t95 + Ifges(4,5) * t129) + t129 * t204 + (-0.2e1 * mrSges(3,2) * t132 + Ifges(3,3) * t183 + t131 * t244) * t183 + (-0.2e1 * pkin(1) * t155 + (-0.2e1 * mrSges(3,3) * t131 + Ifges(3,1) * t215 + Ifges(3,5) * t243) * t176 + (0.2e1 * t132 * mrSges(3,3) + Ifges(3,6) * t243 + (0.2e1 * Ifges(3,4) * t176 + Ifges(3,2) * t180 + pkin(1) * t244) * t179) * t180) * t179 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t13 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2 + t25 ^ 2) + m(4) * (t53 ^ 2 + t54 ^ 2 + t66 ^ 2) + m(3) * (pkin(1) ^ 2 * t179 ^ 2 + t131 ^ 2 + t132 ^ 2) + (t8 - t22) * t43 + Ifges(2,3) + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t19 + 0.2e1 * t1 * t20 + 0.2e1 * t13 * t24 + t28 * t9 + t29 * t10 + 0.2e1 * t6 * t30 + 0.2e1 * t5 * t31 + 0.2e1 * t25 * t35 + t44 * t23 + 0.2e1 * t18 * t46 + 0.2e1 * t17 * t47 + t59 * t34 + 0.2e1 * t66 * t65 + t70 * t32 + 0.2e1 * t54 * t77 + 0.2e1 * t53 * t78; t130 * t35 + t182 * t65 + t61 * t19 + t60 * t20 + t73 * t30 + t99 * t46 + t155 + (t24 - t47) * t97 + t224 * t71 + (-m(3) * pkin(1) - mrSges(3,1) * t180) * t179 + (t187 * t77 + t191 * t78) * t178 + m(7) * (t1 * t60 + t2 * t61 + t3 * t71) + m(6) * (t13 * t97 - t5 * t71 + t6 * t73) + m(5) * (t130 * t25 - t17 * t97 + t18 * t99) + m(4) * (t182 * t66 + (t187 * t54 + t191 * t53) * t178); m(3) + m(7) * (t60 ^ 2 + t61 ^ 2 + t246) + m(6) * (t73 ^ 2 + t245 + t246) + m(5) * (t130 ^ 2 + t99 ^ 2 + t245) + m(4) * (t182 ^ 2 + (t187 ^ 2 + t191 ^ 2) * t178 ^ 2); t203 * t43 + (t83 / 0.2e1 - t116 / 0.2e1) * t58 + (-pkin(3) * t35 + t186 * t34 / 0.2e1 + (-t21 / 0.2e1 + t33 / 0.2e1) * t190) * t177 + t204 + t205 * t133 + m(7) * (t1 * t51 + t2 * t52 + t3 * t75) + m(6) * (t119 * t13 + t5 * t81 + t6 * t82) + m(5) * (t136 * t17 + t137 * t18 - t229 * t25) + t44 * t241 + t134 * t242 + t51 * t20 + t52 * t19 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + t28 * t63 / 0.2e1 + t29 * t64 / 0.2e1 + t3 * t67 + t75 * t11 + t2 * t79 + t1 * t80 + t81 * t31 + t82 * t30 + t13 * t90 + t104 * t9 / 0.2e1 + t105 * t10 / 0.2e1 + t6 * t107 + t5 * t108 + t70 * t115 / 0.2e1 + t59 * t117 / 0.2e1 + t119 * t24 + t25 * t135 + t136 * t47 + t137 * t46 + t17 * t139 + t18 * t140 + t181 * t32 / 0.2e1; t73 * t107 + t130 * t135 + t99 * t140 + t60 * t80 + t61 * t79 + (-t139 + t90) * t97 + t220 * t71 + (mrSges(4,1) * t191 - mrSges(4,2) * t187) * t178 + m(7) * (t51 * t60 + t52 * t61 + t71 * t75) + m(6) * (t119 * t97 - t71 * t81 + t73 * t82) + m(5) * (-t130 * t229 - t136 * t97 + t137 * t99); t104 * t63 + t105 * t64 + 0.2e1 * t82 * t107 + 0.2e1 * t81 * t108 + t181 * t115 + 0.2e1 * t119 * t90 + t134 * t85 + 0.2e1 * t136 * t139 + 0.2e1 * t137 * t140 + 0.2e1 * t51 * t80 + 0.2e1 * t52 * t79 + 0.2e1 * t75 * t67 + Ifges(4,3) + (t62 - t84) * t133 + (-0.2e1 * pkin(3) * t135 + t186 * t117 + (t116 - t83) * t190) * t177 + m(7) * (t51 ^ 2 + t52 ^ 2 + t75 ^ 2) + m(6) * (t119 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (pkin(3) ^ 2 * t177 ^ 2 + t136 ^ 2 + t137 ^ 2); t32 + m(6) * (-pkin(4) * t13 + (-t185 * t5 + t189 * t6) * pkin(12)) + t202 * t43 + (t6 * mrSges(6,3) + pkin(12) * t30 - t205) * t189 + (-t5 * mrSges(6,3) + pkin(12) * t224 + t10 * t231 + t233 * t9 + t242) * t185 + m(7) * (t1 * t113 + t114 * t2 + t228 * t3) + t17 * mrSges(5,1) - t18 * mrSges(5,2) - pkin(4) * t24 + t113 * t20 + t114 * t19 + t28 * t239 + t29 * t238 + t3 * t138 + t2 * t141 + t1 * t142 + t13 * t146 + t58 * t148 / 0.2e1 + t44 * t234; -t99 * mrSges(5,2) + t71 * t138 + t61 * t141 + t60 * t142 + (-mrSges(5,1) + t146) * t97 + t196 * mrSges(6,3) + m(7) * (pkin(12) * t221 + t113 * t60 + t114 * t61) + m(6) * (-pkin(4) * t97 + pkin(12) * t196); m(7) * (t113 * t51 + t114 * t52 + t228 * t75) + m(6) * (-pkin(4) * t119 + (-t185 * t81 + t189 * t82) * pkin(12)) + t202 * t133 + (-t81 * mrSges(6,3) + pkin(12) * t220 + t231 * t64 + t233 * t63 + t241) * t185 + t115 + (t82 * mrSges(6,3) + pkin(12) * t107 - t203) * t189 - t148 * t218 / 0.2e1 - pkin(4) * t90 + t113 * t80 + t114 * t79 + t104 * t239 + t105 * t238 + t136 * mrSges(5,1) - t137 * mrSges(5,2) + t75 * t138 + t52 * t141 + t51 * t142 + t119 * t146 + t134 * t234; -0.2e1 * pkin(4) * t146 + 0.2e1 * t113 * t142 + 0.2e1 * t114 * t141 + Ifges(5,3) + (-t122 + t150) * t189 + (t173 + t175) * mrSges(6,3) * t248 + m(7) * (t113 ^ 2 + t114 ^ 2 + t170) + m(6) * (pkin(4) ^ 2 + t175 * t193 + t170) + (-t123 * t184 + t124 * t188 + t138 * t248 + t152) * t185; t9 * t231 + t10 * t232 - pkin(5) * t11 + m(7) * (-pkin(5) * t3 + pkin(13) * t201) - t20 * t226 + t19 * t225 + t3 * t145 + t28 * t236 + t43 * t237 + t29 * t235 + t5 * mrSges(6,1) - t6 * mrSges(6,2) + t201 * mrSges(7,3) + t21; -t73 * mrSges(6,2) + (m(7) * pkin(13) + mrSges(7,3)) * (-t184 * t60 + t188 * t61) + t247 * t71; -pkin(5) * t67 + t63 * t231 + t75 * t145 + t64 * t232 - t80 * t226 + t105 * t235 + t133 * t237 + t104 * t236 + m(7) * (-pkin(5) * t75 + pkin(13) * t198) + t79 * t225 - t82 * mrSges(6,2) + t81 * mrSges(6,1) + t198 * mrSges(7,3) + t83; t124 * t232 + t123 * t231 - pkin(5) * t138 + (m(7) * t194 + t188 * t141 - t184 * t142) * pkin(13) + (pkin(12) * t247 + t149 * t233 + t151 * t231) * t185 + (-pkin(12) * mrSges(6,2) - t147 / 0.2e1) * t189 + t194 * mrSges(7,3) + t148; Ifges(6,3) - 0.2e1 * pkin(5) * t145 + t184 * t151 + t188 * t149 + m(7) * (pkin(13) ^ 2 * t206 + pkin(5) ^ 2) + 0.2e1 * t206 * pkin(13) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t60 - mrSges(7,2) * t61; mrSges(7,1) * t51 - mrSges(7,2) * t52 + t62; mrSges(7,1) * t113 - mrSges(7,2) * t114 + t122; -pkin(13) * t200 + t147; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
