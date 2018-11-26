% Calculate joint inertia matrix for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2018-11-23 18:18
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:18:08
% EndTime: 2018-11-23 18:18:10
% DurationCPUTime: 2.24s
% Computational Cost: add. (5222->440), mult. (11560->629), div. (0->0), fcn. (13164->12), ass. (0->168)
t195 = sin(pkin(12));
t197 = cos(pkin(12));
t200 = sin(qJ(4));
t201 = sin(qJ(3));
t204 = cos(qJ(4));
t205 = cos(qJ(3));
t160 = t200 * t201 - t204 * t205;
t161 = t200 * t205 + t201 * t204;
t185 = -pkin(3) * t205 - pkin(2);
t111 = pkin(4) * t160 - qJ(5) * t161 + t185;
t256 = -pkin(10) - pkin(9);
t173 = t256 * t205;
t222 = t256 * t201;
t131 = -t204 * t173 + t200 * t222;
t65 = t111 * t197 - t131 * t195;
t66 = t111 * t195 + t131 * t197;
t216 = -t195 * t65 + t197 * t66;
t196 = sin(pkin(6));
t206 = cos(qJ(2));
t228 = t196 * t206;
t198 = cos(pkin(6));
t202 = sin(qJ(2));
t229 = t196 * t202;
t145 = t198 * t201 + t205 * t229;
t147 = pkin(1) * t198 * t202 + pkin(8) * t228;
t137 = pkin(9) * t198 + t147;
t138 = (-pkin(2) * t206 - pkin(9) * t202 - pkin(1)) * t196;
t85 = -t137 * t201 + t138 * t205;
t64 = -pkin(3) * t228 - pkin(10) * t145 + t85;
t144 = t198 * t205 - t201 * t229;
t86 = t137 * t205 + t138 * t201;
t71 = pkin(10) * t144 + t86;
t35 = t200 * t64 + t204 * t71;
t29 = -qJ(5) * t228 + t35;
t100 = t144 * t200 + t145 * t204;
t176 = pkin(8) * t229;
t245 = pkin(1) * t206;
t136 = t176 + (-pkin(2) - t245) * t198;
t107 = -pkin(3) * t144 + t136;
t99 = -t144 * t204 + t145 * t200;
t38 = pkin(4) * t99 - qJ(5) * t100 + t107;
t10 = t195 * t38 + t197 * t29;
t9 = -t195 * t29 + t197 * t38;
t218 = t10 * t197 - t195 * t9;
t199 = sin(qJ(6));
t203 = cos(qJ(6));
t158 = -t195 * t199 + t197 * t203;
t159 = t195 * t203 + t197 * t199;
t115 = Ifges(7,5) * t159 + Ifges(7,6) * t158;
t247 = t197 / 0.2e1;
t248 = t195 / 0.2e1;
t262 = Ifges(6,5) * t248 + Ifges(6,6) * t247 + t115 / 0.2e1;
t261 = -Ifges(4,5) * t145 - Ifges(4,6) * t144;
t129 = -t173 * t200 - t204 * t222;
t260 = t129 ^ 2;
t114 = -t158 * mrSges(7,1) + mrSges(7,2) * t159;
t259 = 0.2e1 * t114;
t258 = 0.2e1 * t129;
t77 = -t100 * t195 - t197 * t228;
t257 = t77 / 0.2e1;
t116 = Ifges(7,4) * t159 + Ifges(7,2) * t158;
t254 = t116 / 0.2e1;
t117 = Ifges(7,1) * t159 + Ifges(7,4) * t158;
t253 = t117 / 0.2e1;
t252 = t158 / 0.2e1;
t251 = t159 / 0.2e1;
t239 = Ifges(6,4) * t197;
t169 = Ifges(6,1) * t195 + t239;
t249 = t169 / 0.2e1;
t244 = pkin(3) * t204;
t242 = Ifges(4,3) + Ifges(5,3);
t241 = -Ifges(5,5) * t100 + Ifges(5,6) * t99;
t240 = Ifges(6,4) * t195;
t146 = t198 * t245 - t176;
t237 = t146 * mrSges(3,1);
t236 = t147 * mrSges(3,2);
t235 = t158 * mrSges(7,3);
t234 = t159 * mrSges(7,3);
t231 = t161 * t195;
t230 = t161 * t197;
t110 = mrSges(6,1) * t231 + mrSges(6,2) * t230;
t227 = Ifges(5,5) * t161 - Ifges(5,6) * t160;
t226 = Ifges(4,5) * t201 + Ifges(4,6) * t205;
t225 = t195 ^ 2 + t197 ^ 2;
t224 = t201 ^ 2 + t205 ^ 2;
t223 = 0.2e1 * mrSges(7,3);
t78 = t100 * t197 - t195 * t228;
t42 = -t199 * t78 + t203 * t77;
t43 = t199 * t77 + t203 * t78;
t11 = Ifges(7,5) * t43 + Ifges(7,6) * t42 + Ifges(7,3) * t99;
t101 = t159 * t161;
t102 = t158 * t161;
t50 = Ifges(7,5) * t102 - Ifges(7,6) * t101 + Ifges(7,3) * t160;
t221 = Ifges(3,5) * t229 + Ifges(3,6) * t228 + Ifges(3,3) * t198;
t182 = -pkin(5) * t197 - pkin(4);
t46 = -t77 * mrSges(6,1) + mrSges(6,2) * t78;
t15 = -t42 * mrSges(7,1) + mrSges(7,2) * t43;
t57 = t101 * mrSges(7,1) + mrSges(7,2) * t102;
t34 = -t200 * t71 + t204 * t64;
t165 = -t197 * mrSges(6,1) + mrSges(6,2) * t195;
t181 = pkin(3) * t200 + qJ(5);
t220 = t225 * t181;
t30 = pkin(4) * t228 - t34;
t168 = Ifges(6,2) * t197 + t240;
t219 = t116 * t158 + t117 * t159 + t168 * t197 + t169 * t195 + Ifges(5,3);
t48 = -mrSges(6,2) * t99 + mrSges(6,3) * t77;
t49 = mrSges(6,1) * t99 - mrSges(6,3) * t78;
t217 = -t195 * t49 + t197 * t48;
t215 = 0.2e1 * mrSges(6,3) * t225;
t112 = -mrSges(6,2) * t160 - mrSges(6,3) * t231;
t113 = mrSges(6,1) * t160 - mrSges(6,3) * t230;
t214 = t112 * t197 - t113 * t195;
t213 = (mrSges(5,1) * t204 - mrSges(5,2) * t200) * pkin(3);
t212 = t114 + t165;
t12 = Ifges(7,4) * t43 + Ifges(7,2) * t42 + Ifges(7,6) * t99;
t13 = Ifges(7,1) * t43 + Ifges(7,4) * t42 + Ifges(7,5) * t99;
t17 = -pkin(5) * t77 + t30;
t4 = pkin(5) * t99 - pkin(11) * t78 + t9;
t5 = pkin(11) * t77 + t10;
t2 = -t199 * t5 + t203 * t4;
t27 = Ifges(6,4) * t78 + Ifges(6,2) * t77 + Ifges(6,6) * t99;
t28 = Ifges(6,1) * t78 + Ifges(6,4) * t77 + Ifges(6,5) * t99;
t3 = t199 * t4 + t203 * t5;
t211 = t34 * mrSges(5,1) - t35 * mrSges(5,2) + t17 * t114 + t12 * t252 + t13 * t251 + t30 * t165 + t168 * t257 - t2 * t234 + t3 * t235 + t27 * t247 + t28 * t248 + t78 * t249 + t43 * t253 + t42 * t254 - t241 + t262 * t99 + t218 * mrSges(6,3);
t47 = pkin(5) * t160 - pkin(11) * t230 + t65;
t55 = -pkin(11) * t231 + t66;
t18 = -t199 * t55 + t203 * t47;
t19 = t199 * t47 + t203 * t55;
t51 = Ifges(7,4) * t102 - Ifges(7,2) * t101 + Ifges(7,6) * t160;
t52 = Ifges(7,1) * t102 - Ifges(7,4) * t101 + Ifges(7,5) * t160;
t81 = Ifges(6,6) * t160 + (-Ifges(6,2) * t195 + t239) * t161;
t82 = Ifges(6,5) * t160 + (Ifges(6,1) * t197 - t240) * t161;
t92 = pkin(5) * t231 + t129;
t210 = -t131 * mrSges(5,2) - t101 * t254 + t102 * t253 + t92 * t114 - t18 * t234 + t19 * t235 + t81 * t247 + t82 * t248 + t52 * t251 + t51 * t252 - t168 * t231 / 0.2e1 + t230 * t249 + t227 + t262 * t160 + (t165 - mrSges(5,1)) * t129 + t216 * mrSges(6,3);
t188 = t197 * pkin(11);
t184 = -pkin(4) - t244;
t172 = Ifges(4,1) * t201 + Ifges(4,4) * t205;
t171 = Ifges(4,4) * t201 + Ifges(4,2) * t205;
t170 = -mrSges(4,1) * t205 + mrSges(4,2) * t201;
t166 = qJ(5) * t197 + t188;
t164 = (-pkin(11) - qJ(5)) * t195;
t163 = t182 - t244;
t149 = t181 * t197 + t188;
t148 = (-pkin(11) - t181) * t195;
t128 = -mrSges(4,1) * t228 - mrSges(4,3) * t145;
t127 = mrSges(4,2) * t228 + mrSges(4,3) * t144;
t125 = t164 * t199 + t166 * t203;
t124 = t164 * t203 - t166 * t199;
t120 = Ifges(5,1) * t161 - Ifges(5,4) * t160;
t119 = Ifges(5,4) * t161 - Ifges(5,2) * t160;
t118 = mrSges(5,1) * t160 + mrSges(5,2) * t161;
t109 = -mrSges(4,1) * t144 + mrSges(4,2) * t145;
t106 = t148 * t199 + t149 * t203;
t105 = t148 * t203 - t149 * t199;
t91 = Ifges(4,1) * t145 + Ifges(4,4) * t144 - Ifges(4,5) * t228;
t90 = Ifges(4,4) * t145 + Ifges(4,2) * t144 - Ifges(4,6) * t228;
t84 = -mrSges(5,1) * t228 - mrSges(5,3) * t100;
t83 = mrSges(5,2) * t228 - mrSges(5,3) * t99;
t80 = Ifges(6,3) * t160 + (Ifges(6,5) * t197 - Ifges(6,6) * t195) * t161;
t73 = mrSges(7,1) * t160 - mrSges(7,3) * t102;
t72 = -mrSges(7,2) * t160 - mrSges(7,3) * t101;
t56 = mrSges(5,1) * t99 + mrSges(5,2) * t100;
t54 = Ifges(5,1) * t100 - Ifges(5,4) * t99 - Ifges(5,5) * t228;
t53 = Ifges(5,4) * t100 - Ifges(5,2) * t99 - Ifges(5,6) * t228;
t26 = Ifges(6,5) * t78 + Ifges(6,6) * t77 + Ifges(6,3) * t99;
t22 = mrSges(7,1) * t99 - mrSges(7,3) * t43;
t21 = -mrSges(7,2) * t99 + mrSges(7,3) * t42;
t1 = [m(7) * (t17 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t10 ^ 2 + t30 ^ 2 + t9 ^ 2) + m(5) * (t107 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(4) * (t136 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(3) * (pkin(1) ^ 2 * t196 ^ 2 + t146 ^ 2 + t147 ^ 2) + (t11 + t26 - t53) * t99 + (t221 - 0.2e1 * t236 + 0.2e1 * t237) * t198 + ((-0.2e1 * t146 * mrSges(3,3) + Ifges(3,5) * t198 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t202) * t196) * t202 + (0.2e1 * t147 * mrSges(3,3) + Ifges(3,6) * t198 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t202 + (Ifges(3,2) + t242) * t206) * t196 + t241 + t261) * t206) * t196 + Ifges(2,3) + 0.2e1 * t17 * t15 + 0.2e1 * t3 * t21 + 0.2e1 * t2 * t22 + t42 * t12 + t43 * t13 + 0.2e1 * t30 * t46 + 0.2e1 * t10 * t48 + 0.2e1 * t9 * t49 + t77 * t27 + t78 * t28 + 0.2e1 * t35 * t83 + 0.2e1 * t34 * t84 + t100 * t54 + 0.2e1 * t107 * t56 + 0.2e1 * t86 * t127 + 0.2e1 * t85 * t128 + 0.2e1 * t136 * t109 + t144 * t90 + t145 * t91; t81 * t257 + (-t35 * mrSges(5,3) + t11 / 0.2e1 + t26 / 0.2e1 - t53 / 0.2e1) * t160 + (t86 * mrSges(4,3) + pkin(9) * t127 + t90 / 0.2e1) * t205 + (-t85 * mrSges(4,3) - pkin(9) * t128 + t91 / 0.2e1) * t201 + (t46 - t84) * t129 - t236 + t237 + (-t34 * mrSges(5,3) - t195 * t27 / 0.2e1 + t28 * t247 + t54 / 0.2e1) * t161 + m(7) * (t17 * t92 + t18 * t2 + t19 * t3) + m(6) * (t10 * t66 + t129 * t30 + t65 * t9) + m(5) * (t107 * t185 - t129 * t34 + t131 * t35) + (t50 / 0.2e1 + t80 / 0.2e1 - t119 / 0.2e1) * t99 + t221 + m(4) * (-pkin(2) * t136 + (-t201 * t85 + t205 * t86) * pkin(9)) - (t226 + t227) * t228 / 0.2e1 + t19 * t21 + t18 * t22 + t42 * t51 / 0.2e1 + t43 * t52 / 0.2e1 + t17 * t57 + t65 * t49 + t66 * t48 + t3 * t72 + t2 * t73 + t78 * t82 / 0.2e1 + t92 * t15 - t101 * t12 / 0.2e1 + t102 * t13 / 0.2e1 - pkin(2) * t109 + t30 * t110 + t10 * t112 + t9 * t113 + t107 * t118 + t100 * t120 / 0.2e1 + t131 * t83 + t136 * t170 + t144 * t171 / 0.2e1 + t145 * t172 / 0.2e1 + t185 * t56; -0.2e1 * pkin(2) * t170 - t101 * t51 + t102 * t52 + t110 * t258 + 0.2e1 * t66 * t112 + 0.2e1 * t65 * t113 + 0.2e1 * t185 * t118 + t205 * t171 + t201 * t172 + 0.2e1 * t18 * t73 + 0.2e1 * t19 * t72 + 0.2e1 * t92 * t57 + Ifges(3,3) + 0.2e1 * t224 * pkin(9) * mrSges(4,3) + (mrSges(5,3) * t258 - t195 * t81 + t197 * t82 + t120) * t161 + (-0.2e1 * mrSges(5,3) * t131 - t119 + t50 + t80) * t160 + m(7) * (t18 ^ 2 + t19 ^ 2 + t92 ^ 2) + m(6) * (t65 ^ 2 + t66 ^ 2 + t260) + m(5) * (t131 ^ 2 + t185 ^ 2 + t260) + m(4) * (pkin(9) ^ 2 * t224 + pkin(2) ^ 2); m(7) * (t105 * t2 + t106 * t3 + t163 * t17) + t217 * t181 - t242 * t228 + (t200 * t83 + t204 * t84 + m(5) * (t200 * t35 + t204 * t34)) * pkin(3) + m(6) * (t181 * t218 + t184 * t30) + t211 + t85 * mrSges(4,1) - t86 * mrSges(4,2) + t105 * t22 + t106 * t21 + t163 * t15 + t184 * t46 - t261; m(7) * (t105 * t18 + t106 * t19 + t163 * t92) + (m(5) * (-t129 * t204 + t131 * t200) + (-t160 * t200 - t161 * t204) * mrSges(5,3)) * pkin(3) + m(6) * (t129 * t184 + t181 * t216) + t210 + t214 * t181 + (-mrSges(4,1) * t201 - mrSges(4,2) * t205) * pkin(9) + t105 * t73 + t106 * t72 + t163 * t57 + t184 * t110 + t226; t163 * t259 + 0.2e1 * t184 * t165 + Ifges(4,3) + 0.2e1 * t213 + (-t105 * t159 + t106 * t158) * t223 + t181 * t215 + m(7) * (t105 ^ 2 + t106 ^ 2 + t163 ^ 2) + m(6) * (t181 ^ 2 * t225 + t184 ^ 2) + m(5) * (t200 ^ 2 + t204 ^ 2) * pkin(3) ^ 2 + t219; m(6) * (-pkin(4) * t30 + qJ(5) * t218) + t211 - Ifges(5,3) * t228 + m(7) * (t124 * t2 + t125 * t3 + t17 * t182) + t217 * qJ(5) - pkin(4) * t46 + t124 * t22 + t125 * t21 + t182 * t15; t214 * qJ(5) + m(7) * (t124 * t18 + t125 * t19 + t182 * t92) + m(6) * (-pkin(4) * t129 + qJ(5) * t216) + t210 - pkin(4) * t110 + t124 * t73 + t125 * t72 + t182 * t57; (t184 - pkin(4)) * t165 + (t163 + t182) * t114 + t213 + m(7) * (t105 * t124 + t106 * t125 + t163 * t182) + m(6) * (-pkin(4) * t184 + qJ(5) * t220) + ((-t105 - t124) * t159 + (t106 + t125) * t158) * mrSges(7,3) + (qJ(5) * t225 + t220) * mrSges(6,3) + t219; -0.2e1 * pkin(4) * t165 + t182 * t259 + (-t124 * t159 + t125 * t158) * t223 + qJ(5) * t215 + m(7) * (t124 ^ 2 + t125 ^ 2 + t182 ^ 2) + m(6) * (qJ(5) ^ 2 * t225 + pkin(4) ^ 2) + t219; m(6) * t30 + m(7) * t17 + t15 + t46; m(6) * t129 + m(7) * t92 + t110 + t57; m(6) * t184 + m(7) * t163 + t212; -m(6) * pkin(4) + m(7) * t182 + t212; m(6) + m(7); mrSges(7,1) * t2 - mrSges(7,2) * t3 + t11; mrSges(7,1) * t18 - mrSges(7,2) * t19 + t50; mrSges(7,1) * t105 - mrSges(7,2) * t106 + t115; mrSges(7,1) * t124 - mrSges(7,2) * t125 + t115; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
