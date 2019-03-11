% Calculate joint inertia matrix for
% S6RRRRRP12
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:45
% EndTime: 2019-03-10 02:59:53
% DurationCPUTime: 3.38s
% Computational Cost: add. (6391->565), mult. (16525->779), div. (0->0), fcn. (18481->12), ass. (0->192)
t205 = sin(qJ(5));
t209 = cos(qJ(5));
t263 = t205 ^ 2 + t209 ^ 2;
t262 = 2 * pkin(11);
t201 = sin(pkin(7));
t203 = cos(pkin(7));
t204 = cos(pkin(6));
t202 = sin(pkin(6));
t212 = cos(qJ(2));
t241 = t202 * t212;
t136 = -t201 * t241 + t203 * t204;
t206 = sin(qJ(4));
t210 = cos(qJ(4));
t207 = sin(qJ(3));
t208 = sin(qJ(2));
t211 = cos(qJ(3));
t239 = t203 * t212;
t244 = t201 * t207;
t96 = t204 * t244 + (t207 * t239 + t208 * t211) * t202;
t71 = t136 * t206 + t210 * t96;
t230 = t202 * t239;
t242 = t202 * t208;
t243 = t201 * t211;
t95 = -t204 * t243 + t207 * t242 - t211 * t230;
t47 = t205 * t71 - t95 * t209;
t48 = t205 * t95 + t209 * t71;
t70 = -t136 * t210 + t206 * t96;
t12 = Ifges(6,5) * t48 - Ifges(6,6) * t47 + Ifges(6,3) * t70;
t13 = Ifges(7,4) * t48 + Ifges(7,2) * t70 + Ifges(7,6) * t47;
t261 = t12 + t13;
t138 = t203 * t206 + t210 * t244;
t103 = t138 * t205 + t209 * t243;
t104 = t138 * t209 - t205 * t243;
t137 = -t210 * t203 + t206 * t244;
t55 = Ifges(6,5) * t104 - Ifges(6,6) * t103 + Ifges(6,3) * t137;
t56 = Ifges(7,4) * t104 + Ifges(7,2) * t137 + Ifges(7,6) * t103;
t260 = t55 + t56;
t256 = pkin(1) * t204;
t179 = t212 * t256;
t102 = pkin(2) * t204 + t179 + (-pkin(10) * t203 - pkin(9)) * t242;
t114 = (-pkin(10) * t201 * t208 - pkin(2) * t212 - pkin(1)) * t202;
t143 = pkin(9) * t241 + t208 * t256;
t91 = (t201 * t204 + t230) * pkin(10) + t143;
t40 = -t207 * t91 + (t102 * t203 + t114 * t201) * t211;
t27 = Ifges(5,1) * t71 - Ifges(5,4) * t70 + Ifges(5,5) * t95;
t259 = t27 / 0.2e1;
t85 = Ifges(5,1) * t138 - Ifges(5,4) * t137 - Ifges(5,5) * t243;
t258 = t85 / 0.2e1;
t167 = Ifges(5,1) * t206 + Ifges(5,4) * t210;
t257 = t167 / 0.2e1;
t34 = -pkin(3) * t136 - t40;
t18 = pkin(4) * t70 - pkin(12) * t71 + t34;
t63 = -t102 * t201 + t203 * t114;
t31 = pkin(3) * t95 - pkin(11) * t96 + t63;
t240 = t203 * t207;
t41 = t102 * t240 + t114 * t244 + t211 * t91;
t35 = pkin(11) * t136 + t41;
t10 = t206 * t31 + t210 * t35;
t8 = pkin(12) * t95 + t10;
t4 = t205 * t18 + t209 * t8;
t255 = pkin(2) * t211;
t254 = pkin(11) * t206;
t253 = pkin(11) * t210;
t252 = -Ifges(7,2) - Ifges(6,3);
t173 = pkin(10) * t244;
t121 = t173 + (-pkin(3) - t255) * t203;
t72 = pkin(4) * t137 - pkin(12) * t138 + t121;
t142 = pkin(2) * t240 + pkin(10) * t243;
t122 = pkin(11) * t203 + t142;
t123 = (-pkin(3) * t211 - pkin(11) * t207 - pkin(2)) * t201;
t82 = t210 * t122 + t206 * t123;
t74 = -pkin(12) * t243 + t82;
t37 = t205 * t72 + t209 * t74;
t251 = Ifges(6,4) * t205;
t250 = Ifges(6,4) * t209;
t249 = Ifges(7,5) * t205;
t248 = Ifges(7,5) * t209;
t141 = -pkin(9) * t242 + t179;
t247 = t141 * mrSges(3,1);
t246 = t143 * mrSges(3,2);
t155 = -pkin(4) * t210 - pkin(12) * t206 - pkin(3);
t245 = t155 * t209;
t238 = t205 * t206;
t237 = t206 * t209;
t151 = t210 * mrSges(7,1) + mrSges(7,2) * t237;
t236 = Ifges(7,4) * t237 + Ifges(7,6) * t238;
t116 = t205 * t155 + t209 * t253;
t160 = Ifges(6,5) * t205 + Ifges(6,6) * t209;
t161 = Ifges(5,5) * t206 + Ifges(5,6) * t210;
t235 = t263 * pkin(12) ^ 2;
t25 = Ifges(5,5) * t71 - Ifges(5,6) * t70 + Ifges(5,3) * t95;
t51 = Ifges(4,5) * t96 - Ifges(4,6) * t95 + Ifges(4,3) * t136;
t11 = Ifges(7,5) * t48 + Ifges(7,6) * t70 + Ifges(7,3) * t47;
t14 = Ifges(6,4) * t48 - Ifges(6,2) * t47 + Ifges(6,6) * t70;
t234 = t11 / 0.2e1 - t14 / 0.2e1;
t15 = Ifges(7,1) * t48 + Ifges(7,4) * t70 + Ifges(7,5) * t47;
t16 = Ifges(6,1) * t48 - Ifges(6,4) * t47 + Ifges(6,5) * t70;
t233 = t15 / 0.2e1 + t16 / 0.2e1;
t54 = Ifges(7,5) * t104 + Ifges(7,6) * t137 + Ifges(7,3) * t103;
t57 = Ifges(6,4) * t104 - Ifges(6,2) * t103 + Ifges(6,6) * t137;
t232 = t54 / 0.2e1 - t57 / 0.2e1;
t58 = Ifges(7,1) * t104 + Ifges(7,4) * t137 + Ifges(7,5) * t103;
t59 = Ifges(6,1) * t104 - Ifges(6,4) * t103 + Ifges(6,5) * t137;
t231 = t58 / 0.2e1 + t59 / 0.2e1;
t117 = Ifges(4,5) * t244 + Ifges(4,6) * t243 + Ifges(4,3) * t203;
t229 = Ifges(3,5) * t242 + Ifges(3,6) * t241 + Ifges(3,3) * t204;
t124 = -Ifges(7,6) * t210 + (Ifges(7,3) * t205 + t248) * t206;
t127 = -Ifges(6,6) * t210 + (-Ifges(6,2) * t205 + t250) * t206;
t228 = t124 / 0.2e1 - t127 / 0.2e1;
t128 = -Ifges(7,4) * t210 + (Ifges(7,1) * t209 + t249) * t206;
t129 = -Ifges(6,5) * t210 + (Ifges(6,1) * t209 - t251) * t206;
t227 = t128 / 0.2e1 + t129 / 0.2e1;
t159 = -Ifges(7,3) * t209 + t249;
t163 = Ifges(6,2) * t209 + t251;
t226 = t159 / 0.2e1 - t163 / 0.2e1;
t162 = Ifges(7,4) * t205 - Ifges(7,6) * t209;
t225 = t160 / 0.2e1 + t162 / 0.2e1;
t165 = Ifges(7,1) * t205 - t248;
t166 = Ifges(6,1) * t205 + t250;
t224 = t165 / 0.2e1 + t166 / 0.2e1;
t24 = -t70 * mrSges(7,1) + t48 * mrSges(7,2);
t80 = -t137 * mrSges(7,1) + t104 * mrSges(7,2);
t9 = -t206 * t35 + t210 * t31;
t81 = -t206 * t122 + t123 * t210;
t26 = Ifges(5,4) * t71 - Ifges(5,2) * t70 + Ifges(5,6) * t95;
t222 = t12 / 0.2e1 + t13 / 0.2e1 - t26 / 0.2e1;
t84 = Ifges(5,4) * t138 - Ifges(5,2) * t137 - Ifges(5,6) * t243;
t221 = -t84 / 0.2e1 + t55 / 0.2e1 + t56 / 0.2e1;
t220 = Ifges(6,5) * t237 - Ifges(6,6) * t238;
t73 = pkin(4) * t243 - t81;
t125 = -Ifges(6,3) * t210 + t220;
t126 = -Ifges(7,2) * t210 + t236;
t164 = Ifges(5,4) * t206 + Ifges(5,2) * t210;
t219 = -t164 / 0.2e1 + t125 / 0.2e1 + t126 / 0.2e1;
t3 = t18 * t209 - t205 * t8;
t218 = t205 * mrSges(6,1) + t209 * mrSges(6,2);
t217 = t205 * mrSges(7,1) - t209 * mrSges(7,3);
t216 = -pkin(5) * t205 + qJ(6) * t209;
t36 = -t205 * t74 + t209 * t72;
t83 = Ifges(5,5) * t138 - Ifges(5,6) * t137 - Ifges(5,3) * t243;
t7 = -pkin(4) * t95 - t9;
t214 = pkin(11) ^ 2;
t200 = t210 ^ 2;
t198 = t206 ^ 2;
t195 = t198 * t214;
t158 = -mrSges(5,1) * t210 + mrSges(5,2) * t206;
t157 = -mrSges(6,1) * t209 + mrSges(6,2) * t205;
t156 = -mrSges(7,1) * t209 - mrSges(7,3) * t205;
t154 = -pkin(5) * t209 - qJ(6) * t205 - pkin(4);
t152 = -mrSges(7,2) * t238 - mrSges(7,3) * t210;
t150 = -mrSges(6,1) * t210 - mrSges(6,3) * t237;
t149 = mrSges(6,2) * t210 - mrSges(6,3) * t238;
t148 = -mrSges(4,2) * t203 + mrSges(4,3) * t243;
t147 = mrSges(4,1) * t203 - mrSges(4,3) * t244;
t145 = t218 * t206;
t144 = t217 * t206;
t140 = t203 * t255 - t173;
t139 = (-mrSges(4,1) * t211 + mrSges(4,2) * t207) * t201;
t130 = (pkin(11) - t216) * t206;
t119 = Ifges(4,5) * t203 + (Ifges(4,1) * t207 + Ifges(4,4) * t211) * t201;
t118 = Ifges(4,6) * t203 + (Ifges(4,4) * t207 + Ifges(4,2) * t211) * t201;
t115 = -t205 * t253 + t245;
t110 = -t245 + (pkin(11) * t205 + pkin(5)) * t210;
t109 = -qJ(6) * t210 + t116;
t108 = -mrSges(5,1) * t243 - mrSges(5,3) * t138;
t107 = mrSges(5,2) * t243 - mrSges(5,3) * t137;
t86 = mrSges(5,1) * t137 + mrSges(5,2) * t138;
t79 = mrSges(6,1) * t137 - mrSges(6,3) * t104;
t78 = -mrSges(6,2) * t137 - mrSges(6,3) * t103;
t77 = -mrSges(7,2) * t103 + mrSges(7,3) * t137;
t76 = mrSges(4,1) * t136 - mrSges(4,3) * t96;
t75 = -mrSges(4,2) * t136 - mrSges(4,3) * t95;
t62 = mrSges(6,1) * t103 + mrSges(6,2) * t104;
t61 = mrSges(7,1) * t103 - mrSges(7,3) * t104;
t60 = mrSges(4,1) * t95 + mrSges(4,2) * t96;
t53 = Ifges(4,1) * t96 - Ifges(4,4) * t95 + Ifges(4,5) * t136;
t52 = Ifges(4,4) * t96 - Ifges(4,2) * t95 + Ifges(4,6) * t136;
t50 = mrSges(5,1) * t95 - mrSges(5,3) * t71;
t49 = -mrSges(5,2) * t95 - mrSges(5,3) * t70;
t39 = mrSges(5,1) * t70 + mrSges(5,2) * t71;
t38 = pkin(5) * t103 - qJ(6) * t104 + t73;
t29 = -pkin(5) * t137 - t36;
t28 = qJ(6) * t137 + t37;
t23 = mrSges(6,1) * t70 - mrSges(6,3) * t48;
t22 = -mrSges(6,2) * t70 - mrSges(6,3) * t47;
t21 = -mrSges(7,2) * t47 + mrSges(7,3) * t70;
t20 = mrSges(6,1) * t47 + mrSges(6,2) * t48;
t19 = mrSges(7,1) * t47 - mrSges(7,3) * t48;
t5 = pkin(5) * t47 - qJ(6) * t48 + t7;
t2 = -pkin(5) * t70 - t3;
t1 = qJ(6) * t70 + t4;
t6 = [(t25 - t52) * t95 + (t16 + t15) * t48 + (t11 - t14) * t47 + t136 * t51 + ((Ifges(3,5) * t208 + Ifges(3,6) * t212) * t204 + 0.2e1 * (-t141 * t208 + t143 * t212) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t212 + mrSges(3,2) * t208) + t208 * (Ifges(3,1) * t208 + Ifges(3,4) * t212) + t212 * (Ifges(3,4) * t208 + Ifges(3,2) * t212) + m(3) * pkin(1) ^ 2) * t202) * t202 + m(3) * (t141 ^ 2 + t143 ^ 2) + t96 * t53 + 0.2e1 * t41 * t75 + 0.2e1 * t40 * t76 + t71 * t27 + 0.2e1 * t63 * t60 + 0.2e1 * t10 * t49 + 0.2e1 * t9 * t50 + 0.2e1 * t34 * t39 + 0.2e1 * t1 * t21 + 0.2e1 * t4 * t22 + 0.2e1 * t3 * t23 + 0.2e1 * t2 * t24 + 0.2e1 * t5 * t19 + 0.2e1 * t7 * t20 + Ifges(2,3) + (t229 - 0.2e1 * t246 + 0.2e1 * t247) * t204 + (-t26 + t261) * t70 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t7 ^ 2) + m(5) * (t10 ^ 2 + t34 ^ 2 + t9 ^ 2) + m(4) * (t40 ^ 2 + t41 ^ 2 + t63 ^ 2); t203 * t51 / 0.2e1 + t63 * t139 + t140 * t76 + t142 * t75 + t40 * t147 + t41 * t148 + t136 * t117 / 0.2e1 + (-pkin(2) * t60 + t207 * t53 / 0.2e1 + (-t25 / 0.2e1 + t52 / 0.2e1) * t211) * t201 + t247 - t246 + m(7) * (t1 * t28 + t2 * t29 + t38 * t5) + m(6) * (t3 * t36 + t37 * t4 + t7 * t73) + m(5) * (t10 * t82 + t121 * t34 + t81 * t9) + m(4) * (-pkin(2) * t201 * t63 + t140 * t40 + t142 * t41) + (-t118 / 0.2e1 + t83 / 0.2e1) * t95 + t229 + t71 * t258 + t138 * t259 + t221 * t70 + t222 * t137 + t96 * t119 / 0.2e1 + t121 * t39 + t10 * t107 + t9 * t108 + t34 * t86 + t1 * t77 + t4 * t78 + t3 * t79 + t2 * t80 + t81 * t50 + t82 * t49 + t73 * t20 + t5 * t61 + t7 * t62 + t36 * t23 + t37 * t22 + t38 * t19 + t28 * t21 + t29 * t24 + t231 * t48 + t232 * t47 + t233 * t104 + t234 * t103; 0.2e1 * t82 * t107 + 0.2e1 * t81 * t108 + t203 * t117 + 0.2e1 * t121 * t86 + t138 * t85 + 0.2e1 * t140 * t147 + 0.2e1 * t142 * t148 + 0.2e1 * t28 * t77 + 0.2e1 * t29 * t80 + 0.2e1 * t36 * t79 + 0.2e1 * t37 * t78 + 0.2e1 * t38 * t61 + 0.2e1 * t73 * t62 + Ifges(3,3) + (t58 + t59) * t104 + (t54 - t57) * t103 + (-t84 + t260) * t137 + (-0.2e1 * pkin(2) * t139 + t207 * t119 + (t118 - t83) * t211) * t201 + m(7) * (t28 ^ 2 + t29 ^ 2 + t38 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2 + t73 ^ 2) + m(5) * (t121 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (pkin(2) ^ 2 * t201 ^ 2 + t140 ^ 2 + t142 ^ 2); m(7) * (t1 * t109 + t110 * t2 + t130 * t5) + t34 * t158 + t95 * t161 / 0.2e1 + t5 * t144 + t7 * t145 + t4 * t149 + t3 * t150 + t2 * t151 + t1 * t152 + m(5) * (-pkin(3) * t34 + (t10 * t210 - t206 * t9) * pkin(11)) + t51 + t71 * t257 + t227 * t48 + t228 * t47 + (t10 * mrSges(5,3) + pkin(11) * t49 - t222) * t210 + t219 * t70 + t115 * t23 + t116 * t22 + t130 * t19 + t109 * t21 + t110 * t24 - pkin(3) * t39 + t40 * mrSges(4,1) - t41 * mrSges(4,2) + m(6) * (t115 * t3 + t116 * t4 + t254 * t7) + (t259 - t9 * mrSges(5,3) + t233 * t209 + t234 * t205 + (t20 - t50) * pkin(11)) * t206; t121 * t158 + t140 * mrSges(4,1) - t142 * mrSges(4,2) + t38 * t144 + t73 * t145 + t37 * t149 + t36 * t150 + t29 * t151 + t28 * t152 + t117 + m(5) * (-pkin(3) * t121 + (-t206 * t81 + t210 * t82) * pkin(11)) + m(7) * (t109 * t28 + t110 * t29 + t130 * t38) + t138 * t257 + t227 * t104 + t228 * t103 + (t82 * mrSges(5,3) + pkin(11) * t107 - t221) * t210 + t219 * t137 + t115 * t79 + t116 * t78 + t130 * t61 + t109 * t77 + t110 * t80 - pkin(3) * t86 - t161 * t243 / 0.2e1 + m(6) * (t115 * t36 + t116 * t37 + t254 * t73) + (t258 - t81 * mrSges(5,3) + t231 * t209 + t232 * t205 + (t62 - t108) * pkin(11)) * t206; -0.2e1 * pkin(3) * t158 + 0.2e1 * t109 * t152 + 0.2e1 * t110 * t151 + 0.2e1 * t115 * t150 + 0.2e1 * t116 * t149 + 0.2e1 * t130 * t144 + Ifges(4,3) + (t198 + t200) * mrSges(5,3) * t262 + (-t126 - t125 + t164) * t210 + m(6) * (t115 ^ 2 + t116 ^ 2 + t195) + m(7) * (t109 ^ 2 + t110 ^ 2 + t130 ^ 2) + m(5) * (pkin(3) ^ 2 + t200 * t214 + t195) + (t145 * t262 + t167 + (t128 + t129) * t209 + (t124 - t127) * t205) * t206; t9 * mrSges(5,1) - t10 * mrSges(5,2) - pkin(4) * t20 + t154 * t19 + t5 * t156 + t7 * t157 + t225 * t70 + t224 * t48 + t226 * t47 + (t1 * mrSges(7,2) + t4 * mrSges(6,3) + (t21 + t22) * pkin(12) - t234) * t209 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + (-t23 + t24) * pkin(12) + t233) * t205 + m(6) * (-pkin(4) * t7 + (-t205 * t3 + t209 * t4) * pkin(12)) + m(7) * (t154 * t5 + (t1 * t209 + t2 * t205) * pkin(12)) + t25; t81 * mrSges(5,1) - t82 * mrSges(5,2) - pkin(4) * t62 + t154 * t61 + t38 * t156 + t73 * t157 + t225 * t137 + t224 * t104 + t226 * t103 + (t28 * mrSges(7,2) + t37 * mrSges(6,3) + (t77 + t78) * pkin(12) - t232) * t209 + (t29 * mrSges(7,2) - t36 * mrSges(6,3) + (-t79 + t80) * pkin(12) + t231) * t205 + m(6) * (-pkin(4) * t73 + (-t205 * t36 + t209 * t37) * pkin(12)) + m(7) * (t154 * t38 + (t205 * t29 + t209 * t28) * pkin(12)) + t83; -pkin(4) * t145 + t154 * t144 + (m(7) * t154 + t156) * t130 + (-pkin(11) * mrSges(5,2) - t225) * t210 + (t109 * mrSges(7,2) + t116 * mrSges(6,3) - t228) * t209 + (t110 * mrSges(7,2) - t115 * mrSges(6,3) + t227) * t205 + ((t149 + t152) * t209 + (-t150 + t151) * t205 + m(7) * (t109 * t209 + t110 * t205) + m(6) * (-t115 * t205 + t116 * t209)) * pkin(12) + (t224 * t209 + t226 * t205 + (-m(6) * pkin(4) - mrSges(5,1) + t157) * pkin(11)) * t206 + t161; -0.2e1 * pkin(4) * t157 + 0.2e1 * t154 * t156 + Ifges(5,3) + (-t159 + t163) * t209 + (t165 + t166) * t205 + m(7) * (t154 ^ 2 + t235) + m(6) * (pkin(4) ^ 2 + t235) + 0.2e1 * (mrSges(7,2) + mrSges(6,3)) * pkin(12) * t263; -pkin(5) * t24 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t21 + t1 * mrSges(7,3) - t4 * mrSges(6,2) + t3 * mrSges(6,1) - t2 * mrSges(7,1) + t261; t36 * mrSges(6,1) - pkin(5) * t80 + m(7) * (-pkin(5) * t29 + qJ(6) * t28) + qJ(6) * t77 + t28 * mrSges(7,3) - t37 * mrSges(6,2) - t29 * mrSges(7,1) + t260; -pkin(5) * t151 + m(7) * (-pkin(5) * t110 + qJ(6) * t109) + qJ(6) * t152 + t109 * mrSges(7,3) - t110 * mrSges(7,1) - t116 * mrSges(6,2) + t115 * mrSges(6,1) + t252 * t210 + t220 + t236; t216 * mrSges(7,2) + (m(7) * t216 - t217 - t218) * pkin(12) + t162 + t160; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t252; m(7) * t2 + t24; m(7) * t29 + t80; m(7) * t110 + t151; (m(7) * pkin(12) + mrSges(7,2)) * t205; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
