% Calculate joint inertia matrix for
% S6RRRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:44
% EndTime: 2019-03-10 02:32:52
% DurationCPUTime: 2.88s
% Computational Cost: add. (6336->559), mult. (16401->770), div. (0->0), fcn. (18398->12), ass. (0->194)
t263 = 2 * pkin(11);
t262 = -2 * mrSges(7,3);
t208 = sin(qJ(5));
t212 = cos(qJ(5));
t204 = sin(pkin(7));
t206 = cos(pkin(7));
t207 = cos(pkin(6));
t205 = sin(pkin(6));
t215 = cos(qJ(2));
t239 = t205 * t215;
t136 = -t204 * t239 + t206 * t207;
t209 = sin(qJ(4));
t213 = cos(qJ(4));
t210 = sin(qJ(3));
t211 = sin(qJ(2));
t214 = cos(qJ(3));
t237 = t206 * t215;
t242 = t204 * t210;
t98 = t207 * t242 + (t210 * t237 + t211 * t214) * t205;
t73 = t136 * t209 + t213 * t98;
t229 = t205 * t237;
t240 = t205 * t211;
t241 = t204 * t214;
t97 = -t207 * t241 + t210 * t240 - t214 * t229;
t47 = -t208 * t73 + t212 * t97;
t48 = t208 * t97 + t212 * t73;
t72 = -t136 * t213 + t209 * t98;
t11 = Ifges(7,5) * t48 + Ifges(7,6) * t47 + Ifges(7,3) * t72;
t12 = Ifges(6,5) * t48 + Ifges(6,6) * t47 + Ifges(6,3) * t72;
t261 = t11 + t12;
t138 = t206 * t209 + t213 * t242;
t105 = -t138 * t208 - t212 * t241;
t106 = t138 * t212 - t208 * t241;
t137 = -t213 * t206 + t209 * t242;
t55 = Ifges(7,5) * t106 + Ifges(7,6) * t105 + Ifges(7,3) * t137;
t56 = Ifges(6,5) * t106 + Ifges(6,6) * t105 + Ifges(6,3) * t137;
t260 = t55 + t56;
t159 = -mrSges(6,1) * t212 + mrSges(6,2) * t208;
t259 = -m(6) * pkin(4) + t159;
t254 = pkin(1) * t207;
t182 = t215 * t254;
t104 = pkin(2) * t207 + t182 + (-pkin(10) * t206 - pkin(9)) * t240;
t115 = (-pkin(10) * t204 * t211 - pkin(2) * t215 - pkin(1)) * t205;
t143 = pkin(9) * t239 + t211 * t254;
t92 = (t204 * t207 + t229) * pkin(10) + t143;
t40 = -t210 * t92 + (t104 * t206 + t115 * t204) * t214;
t30 = Ifges(5,1) * t73 - Ifges(5,4) * t72 + Ifges(5,5) * t97;
t257 = t30 / 0.2e1;
t87 = Ifges(5,1) * t138 - Ifges(5,4) * t137 - Ifges(5,5) * t241;
t256 = t87 / 0.2e1;
t170 = Ifges(5,1) * t209 + Ifges(5,4) * t213;
t255 = t170 / 0.2e1;
t35 = -pkin(3) * t136 - t40;
t19 = pkin(4) * t72 - pkin(12) * t73 + t35;
t64 = -t104 * t204 + t206 * t115;
t32 = pkin(3) * t97 - pkin(11) * t98 + t64;
t238 = t206 * t210;
t41 = t104 * t238 + t115 * t242 + t214 * t92;
t36 = pkin(11) * t136 + t41;
t10 = t209 * t32 + t213 * t36;
t8 = pkin(12) * t97 + t10;
t4 = t208 * t19 + t212 * t8;
t253 = pkin(2) * t214;
t252 = pkin(11) * t209;
t251 = pkin(11) * t213;
t250 = -Ifges(6,3) - Ifges(7,3);
t249 = -qJ(6) - pkin(12);
t176 = pkin(10) * t242;
t122 = t176 + (-pkin(3) - t253) * t206;
t74 = pkin(4) * t137 - pkin(12) * t138 + t122;
t142 = pkin(2) * t238 + pkin(10) * t241;
t123 = pkin(11) * t206 + t142;
t124 = (-pkin(3) * t214 - pkin(11) * t210 - pkin(2)) * t204;
t84 = t213 * t123 + t209 * t124;
t76 = -pkin(12) * t241 + t84;
t38 = t208 * t74 + t212 * t76;
t248 = Ifges(6,4) * t208;
t247 = Ifges(6,4) * t212;
t246 = Ifges(7,4) * t208;
t245 = Ifges(7,4) * t212;
t141 = -pkin(9) * t240 + t182;
t244 = t141 * mrSges(3,1);
t243 = t143 * mrSges(3,2);
t236 = t208 * t209;
t235 = t209 * t212;
t144 = mrSges(7,1) * t236 + mrSges(7,2) * t235;
t156 = -pkin(4) * t213 - pkin(12) * t209 - pkin(3);
t117 = t208 * t156 + t212 * t251;
t162 = Ifges(7,5) * t208 + Ifges(7,6) * t212;
t163 = Ifges(6,5) * t208 + Ifges(6,6) * t212;
t164 = Ifges(5,5) * t209 + Ifges(5,6) * t213;
t234 = t208 ^ 2 + t212 ^ 2;
t28 = Ifges(5,5) * t73 - Ifges(5,6) * t72 + Ifges(5,3) * t97;
t52 = Ifges(4,5) * t98 - Ifges(4,6) * t97 + Ifges(4,3) * t136;
t13 = Ifges(7,4) * t48 + Ifges(7,2) * t47 + Ifges(7,6) * t72;
t14 = Ifges(6,4) * t48 + Ifges(6,2) * t47 + Ifges(6,6) * t72;
t233 = t13 / 0.2e1 + t14 / 0.2e1;
t15 = Ifges(7,1) * t48 + Ifges(7,4) * t47 + Ifges(7,5) * t72;
t16 = Ifges(6,1) * t48 + Ifges(6,4) * t47 + Ifges(6,5) * t72;
t232 = t15 / 0.2e1 + t16 / 0.2e1;
t57 = Ifges(7,4) * t106 + Ifges(7,2) * t105 + Ifges(7,6) * t137;
t58 = Ifges(6,4) * t106 + Ifges(6,2) * t105 + Ifges(6,6) * t137;
t231 = t57 / 0.2e1 + t58 / 0.2e1;
t59 = Ifges(7,1) * t106 + Ifges(7,4) * t105 + Ifges(7,5) * t137;
t60 = Ifges(6,1) * t106 + Ifges(6,4) * t105 + Ifges(6,5) * t137;
t230 = t59 / 0.2e1 + t60 / 0.2e1;
t118 = Ifges(4,5) * t242 + Ifges(4,6) * t241 + Ifges(4,3) * t206;
t228 = Ifges(3,5) * t240 + Ifges(3,6) * t239 + Ifges(3,3) * t207;
t127 = -Ifges(7,6) * t213 + (-Ifges(7,2) * t208 + t245) * t209;
t128 = -Ifges(6,6) * t213 + (-Ifges(6,2) * t208 + t247) * t209;
t227 = t127 / 0.2e1 + t128 / 0.2e1;
t129 = -Ifges(7,5) * t213 + (Ifges(7,1) * t212 - t246) * t209;
t130 = -Ifges(6,5) * t213 + (Ifges(6,1) * t212 - t248) * t209;
t226 = t129 / 0.2e1 + t130 / 0.2e1;
t225 = t162 / 0.2e1 + t163 / 0.2e1;
t165 = Ifges(7,2) * t212 + t246;
t166 = Ifges(6,2) * t212 + t248;
t224 = t165 / 0.2e1 + t166 / 0.2e1;
t168 = Ifges(7,1) * t208 + t245;
t169 = Ifges(6,1) * t208 + t247;
t223 = t168 / 0.2e1 + t169 / 0.2e1;
t20 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t3 = t212 * t19 - t208 * t8;
t62 = -t105 * mrSges(7,1) + t106 * mrSges(7,2);
t37 = -t208 * t76 + t212 * t74;
t9 = -t209 * t36 + t213 * t32;
t158 = -t212 * mrSges(7,1) + t208 * mrSges(7,2);
t83 = -t209 * t123 + t124 * t213;
t29 = Ifges(5,4) * t73 - Ifges(5,2) * t72 + Ifges(5,6) * t97;
t222 = t11 / 0.2e1 + t12 / 0.2e1 - t29 / 0.2e1;
t86 = Ifges(5,4) * t138 - Ifges(5,2) * t137 - Ifges(5,6) * t241;
t221 = -t86 / 0.2e1 + t55 / 0.2e1 + t56 / 0.2e1;
t75 = pkin(4) * t241 - t83;
t185 = Ifges(7,5) * t235;
t125 = -Ifges(7,6) * t236 - Ifges(7,3) * t213 + t185;
t186 = Ifges(6,5) * t235;
t126 = -Ifges(6,6) * t236 - Ifges(6,3) * t213 + t186;
t167 = Ifges(5,4) * t209 + Ifges(5,2) * t213;
t220 = -t167 / 0.2e1 + t125 / 0.2e1 + t126 / 0.2e1;
t219 = mrSges(6,1) * t208 + mrSges(6,2) * t212;
t85 = Ifges(5,5) * t138 - Ifges(5,6) * t137 - Ifges(5,3) * t241;
t7 = -pkin(4) * t97 - t9;
t217 = pkin(11) ^ 2;
t203 = t213 ^ 2;
t201 = t209 ^ 2;
t199 = t201 * t217;
t189 = -pkin(5) * t212 - pkin(4);
t161 = t249 * t212;
t160 = -mrSges(5,1) * t213 + mrSges(5,2) * t209;
t157 = t249 * t208;
t155 = (pkin(5) * t208 + pkin(11)) * t209;
t153 = -mrSges(6,1) * t213 - mrSges(6,3) * t235;
t152 = -mrSges(7,1) * t213 - mrSges(7,3) * t235;
t151 = mrSges(6,2) * t213 - mrSges(6,3) * t236;
t150 = mrSges(7,2) * t213 - mrSges(7,3) * t236;
t149 = -mrSges(4,2) * t206 + mrSges(4,3) * t241;
t148 = mrSges(4,1) * t206 - mrSges(4,3) * t242;
t147 = t212 * t156;
t145 = t219 * t209;
t140 = t206 * t253 - t176;
t139 = (-mrSges(4,1) * t214 + mrSges(4,2) * t210) * t204;
t120 = Ifges(4,5) * t206 + (Ifges(4,1) * t210 + Ifges(4,4) * t214) * t204;
t119 = Ifges(4,6) * t206 + (Ifges(4,4) * t210 + Ifges(4,2) * t214) * t204;
t116 = -t208 * t251 + t147;
t111 = -mrSges(5,1) * t241 - mrSges(5,3) * t138;
t110 = mrSges(5,2) * t241 - mrSges(5,3) * t137;
t108 = -qJ(6) * t236 + t117;
t96 = -qJ(6) * t235 + t147 + (-pkin(11) * t208 - pkin(5)) * t213;
t88 = mrSges(5,1) * t137 + mrSges(5,2) * t138;
t82 = mrSges(6,1) * t137 - mrSges(6,3) * t106;
t81 = mrSges(7,1) * t137 - mrSges(7,3) * t106;
t80 = -mrSges(6,2) * t137 + mrSges(6,3) * t105;
t79 = -mrSges(7,2) * t137 + mrSges(7,3) * t105;
t78 = mrSges(4,1) * t136 - mrSges(4,3) * t98;
t77 = -mrSges(4,2) * t136 - mrSges(4,3) * t97;
t63 = -mrSges(6,1) * t105 + mrSges(6,2) * t106;
t61 = mrSges(4,1) * t97 + mrSges(4,2) * t98;
t54 = Ifges(4,1) * t98 - Ifges(4,4) * t97 + Ifges(4,5) * t136;
t53 = Ifges(4,4) * t98 - Ifges(4,2) * t97 + Ifges(4,6) * t136;
t51 = -pkin(5) * t105 + t75;
t50 = mrSges(5,1) * t97 - mrSges(5,3) * t73;
t49 = -mrSges(5,2) * t97 - mrSges(5,3) * t72;
t39 = mrSges(5,1) * t72 + mrSges(5,2) * t73;
t27 = qJ(6) * t105 + t38;
t26 = mrSges(6,1) * t72 - mrSges(6,3) * t48;
t25 = mrSges(7,1) * t72 - mrSges(7,3) * t48;
t24 = -mrSges(6,2) * t72 + mrSges(6,3) * t47;
t23 = -mrSges(7,2) * t72 + mrSges(7,3) * t47;
t22 = pkin(5) * t137 - qJ(6) * t106 + t37;
t21 = -mrSges(6,1) * t47 + mrSges(6,2) * t48;
t5 = -pkin(5) * t47 + t7;
t2 = qJ(6) * t47 + t4;
t1 = pkin(5) * t72 - qJ(6) * t48 + t3;
t6 = [(t15 + t16) * t48 + (t13 + t14) * t47 + ((Ifges(3,5) * t211 + Ifges(3,6) * t215) * t207 + 0.2e1 * (-t141 * t211 + t143 * t215) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t215 + mrSges(3,2) * t211) + t215 * (Ifges(3,4) * t211 + Ifges(3,2) * t215) + t211 * (Ifges(3,1) * t211 + Ifges(3,4) * t215) + m(3) * pkin(1) ^ 2) * t205) * t205 + t136 * t52 + t98 * t54 + 0.2e1 * t41 * t77 + 0.2e1 * t40 * t78 + t73 * t30 + 0.2e1 * t64 * t61 + 0.2e1 * t10 * t49 + 0.2e1 * t9 * t50 + 0.2e1 * t35 * t39 + 0.2e1 * t2 * t23 + 0.2e1 * t4 * t24 + 0.2e1 * t1 * t25 + 0.2e1 * t3 * t26 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) + m(5) * (t10 ^ 2 + t35 ^ 2 + t9 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2 + t64 ^ 2) + (t28 - t53) * t97 + Ifges(2,3) + (-t29 + t261) * t72 + 0.2e1 * t5 * t20 + 0.2e1 * t7 * t21 + m(3) * (t141 ^ 2 + t143 ^ 2) + (t228 - 0.2e1 * t243 + 0.2e1 * t244) * t207; t206 * t52 / 0.2e1 + (-pkin(2) * t61 + t210 * t54 / 0.2e1 + (-t28 / 0.2e1 + t53 / 0.2e1) * t214) * t204 - t243 + t244 + m(7) * (t1 * t22 + t2 * t27 + t5 * t51) + m(6) * (t3 * t37 + t38 * t4 + t7 * t75) + m(5) * (t10 * t84 + t122 * t35 + t83 * t9) + m(4) * (-pkin(2) * t204 * t64 + t140 * t40 + t142 * t41) + t228 + t73 * t256 + t138 * t257 + (-t119 / 0.2e1 + t85 / 0.2e1) * t97 + t64 * t139 + t140 * t78 + t142 * t77 + t40 * t148 + t41 * t149 + t136 * t118 / 0.2e1 + t98 * t120 / 0.2e1 + t122 * t39 + t10 * t110 + t9 * t111 + t35 * t88 + t2 * t79 + t4 * t80 + t1 * t81 + t3 * t82 + t83 * t50 + t84 * t49 + t75 * t21 + t5 * t62 + t7 * t63 + t51 * t20 + t37 * t26 + t38 * t24 + t22 * t25 + t27 * t23 + t221 * t72 + t222 * t137 + t230 * t48 + t231 * t47 + t232 * t106 + t233 * t105; 0.2e1 * t84 * t110 + 0.2e1 * t83 * t111 + t206 * t118 + 0.2e1 * t122 * t88 + t138 * t87 + 0.2e1 * t140 * t148 + 0.2e1 * t142 * t149 + 0.2e1 * t22 * t81 + 0.2e1 * t27 * t79 + 0.2e1 * t37 * t82 + 0.2e1 * t38 * t80 + 0.2e1 * t51 * t62 + 0.2e1 * t75 * t63 + Ifges(3,3) + (t60 + t59) * t106 + (t57 + t58) * t105 + (-t86 + t260) * t137 + (-0.2e1 * pkin(2) * t139 + t210 * t120 + (t119 - t85) * t214) * t204 + m(7) * (t22 ^ 2 + t27 ^ 2 + t51 ^ 2) + m(6) * (t37 ^ 2 + t38 ^ 2 + t75 ^ 2) + m(5) * (t122 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (pkin(2) ^ 2 * t204 ^ 2 + t140 ^ 2 + t142 ^ 2); m(5) * (-pkin(3) * t35 + (t10 * t213 - t209 * t9) * pkin(11)) + t52 + m(7) * (t1 * t96 + t108 * t2 + t155 * t5) + t73 * t255 + t3 * t153 + t155 * t20 + t35 * t160 + t97 * t164 / 0.2e1 + t5 * t144 + t7 * t145 + t2 * t150 + t4 * t151 + t1 * t152 + t116 * t26 + t117 * t24 + t108 * t23 + t96 * t25 - pkin(3) * t39 + t40 * mrSges(4,1) - t41 * mrSges(4,2) + t220 * t72 + (t10 * mrSges(5,3) + pkin(11) * t49 - t222) * t213 + t226 * t48 + t227 * t47 + m(6) * (t116 * t3 + t117 * t4 + t252 * t7) + (t257 - t9 * mrSges(5,3) + t232 * t212 - t233 * t208 + (t21 - t50) * pkin(11)) * t209; m(5) * (-pkin(3) * t122 + (-t83 * t209 + t84 * t213) * pkin(11)) + m(7) * (t108 * t27 + t155 * t51 + t22 * t96) + t118 + t138 * t255 + t37 * t153 + t155 * t62 + t122 * t160 + t140 * mrSges(4,1) - t142 * mrSges(4,2) + t51 * t144 + t75 * t145 + t27 * t150 + t38 * t151 + t22 * t152 + t116 * t82 + t117 * t80 + t108 * t79 - pkin(3) * t88 + t96 * t81 + t220 * t137 + (t84 * mrSges(5,3) + pkin(11) * t110 - t221) * t213 + t226 * t106 + t227 * t105 - t164 * t241 / 0.2e1 + m(6) * (t116 * t37 + t117 * t38 + t252 * t75) + (t256 - t83 * mrSges(5,3) + t230 * t212 - t231 * t208 + (t63 - t111) * pkin(11)) * t209; -0.2e1 * pkin(3) * t160 + 0.2e1 * t108 * t150 + 0.2e1 * t116 * t153 + 0.2e1 * t117 * t151 + 0.2e1 * t155 * t144 + 0.2e1 * t96 * t152 + Ifges(4,3) + (t201 + t203) * mrSges(5,3) * t263 + (-t125 - t126 + t167) * t213 + m(6) * (t116 ^ 2 + t117 ^ 2 + t199) + m(7) * (t108 ^ 2 + t155 ^ 2 + t96 ^ 2) + m(5) * (pkin(3) ^ 2 + t203 * t217 + t199) + (t145 * t263 + t170 + (t129 + t130) * t212 + (-t127 - t128) * t208) * t209; t9 * mrSges(5,1) - t10 * mrSges(5,2) - pkin(4) * t21 + t157 * t25 + t5 * t158 - t161 * t23 + t189 * t20 + t225 * t72 + t223 * t48 + t224 * t47 + m(7) * (t1 * t157 - t161 * t2 + t189 * t5) + (t2 * mrSges(7,3) + t4 * mrSges(6,3) + (m(6) * t4 + t24) * pkin(12) + t233) * t212 + (-t1 * mrSges(7,3) - t3 * mrSges(6,3) + (-m(6) * t3 - t26) * pkin(12) + t232) * t208 + t28 + t259 * t7; t83 * mrSges(5,1) - t84 * mrSges(5,2) - pkin(4) * t63 + t157 * t81 + t51 * t158 - t161 * t79 + t189 * t62 + t225 * t137 + t223 * t106 + t224 * t105 + m(7) * (t157 * t22 - t161 * t27 + t189 * t51) + (t27 * mrSges(7,3) + t38 * mrSges(6,3) + (m(6) * t38 + t80) * pkin(12) + t231) * t212 + (-t22 * mrSges(7,3) - t37 * mrSges(6,3) + (-m(6) * t37 - t82) * pkin(12) + t230) * t208 + t85 + t259 * t75; t189 * t144 + t157 * t152 + t155 * t158 - t161 * t150 - pkin(4) * t145 + m(7) * (-t108 * t161 + t155 * t189 + t157 * t96) - t225 * t213 + (-t213 * mrSges(5,2) + (-mrSges(5,1) + t259) * t209) * pkin(11) + (t108 * mrSges(7,3) + t117 * mrSges(6,3) + t223 * t209 + (m(6) * t117 + t151) * pkin(12) + t227) * t212 + (-t96 * mrSges(7,3) - t116 * mrSges(6,3) - t224 * t209 + (-m(6) * t116 - t153) * pkin(12) + t226) * t208 + t164; -0.2e1 * pkin(4) * t159 + 0.2e1 * t189 * t158 + Ifges(5,3) + 0.2e1 * t234 * pkin(12) * mrSges(6,3) + m(7) * (t157 ^ 2 + t161 ^ 2 + t189 ^ 2) + m(6) * (pkin(12) ^ 2 * t234 + pkin(4) ^ 2) + (t161 * t262 + t165 + t166) * t212 + (t157 * t262 + t168 + t169) * t208; mrSges(6,1) * t3 + mrSges(7,1) * t1 - mrSges(6,2) * t4 - mrSges(7,2) * t2 + (m(7) * t1 + t25) * pkin(5) + t261; mrSges(6,1) * t37 + mrSges(7,1) * t22 - mrSges(6,2) * t38 - mrSges(7,2) * t27 + (m(7) * t22 + t81) * pkin(5) + t260; mrSges(6,1) * t116 + mrSges(7,1) * t96 - mrSges(6,2) * t117 - mrSges(7,2) * t108 + t185 + t186 + t250 * t213 + (-Ifges(6,6) - Ifges(7,6)) * t236 + (m(7) * t96 + t152) * pkin(5); mrSges(7,1) * t157 + mrSges(7,2) * t161 - t219 * pkin(12) + (m(7) * t157 - t208 * mrSges(7,3)) * pkin(5) + t163 + t162; (m(7) * pkin(5) + 0.2e1 * mrSges(7,1)) * pkin(5) - t250; m(7) * t5 + t20; m(7) * t51 + t62; m(7) * t155 + t144; m(7) * t189 + t158; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
