% Calculate joint inertia matrix for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR7_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR7_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR7_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:19:54
% EndTime: 2019-03-09 05:19:59
% DurationCPUTime: 2.38s
% Computational Cost: add. (5360->343), mult. (5233->488), div. (0->0), fcn. (5323->10), ass. (0->168)
t256 = Icges(5,3) + Icges(6,3);
t169 = qJ(3) + qJ(4);
t157 = pkin(10) + t169;
t154 = sin(t157);
t155 = cos(t157);
t158 = sin(t169);
t159 = cos(t169);
t255 = Icges(5,5) * t158 + Icges(6,5) * t154 + Icges(5,6) * t159 + Icges(6,6) * t155;
t172 = sin(qJ(1));
t175 = cos(qJ(1));
t254 = t172 * t255 + t175 * t256;
t253 = t172 * t256 - t175 * t255;
t170 = sin(qJ(6));
t173 = cos(qJ(6));
t78 = rSges(7,3) * t154 + (rSges(7,1) * t173 - rSges(7,2) * t170) * t155;
t252 = pkin(5) * t155 + pkin(9) * t154 + t78;
t251 = Icges(5,5) * t159 + Icges(6,5) * t155 - Icges(5,6) * t158 - Icges(6,6) * t154;
t226 = Icges(6,4) * t155;
t120 = -Icges(6,2) * t154 + t226;
t227 = Icges(6,4) * t154;
t121 = Icges(6,1) * t155 - t227;
t228 = Icges(5,4) * t159;
t128 = -Icges(5,2) * t158 + t228;
t229 = Icges(5,4) * t158;
t129 = Icges(5,1) * t159 - t229;
t250 = t120 * t155 + t121 * t154 + t128 * t159 + t129 * t158;
t248 = t172 * t175;
t247 = (rSges(6,1) * t154 + rSges(6,2) * t155) * t175;
t235 = rSges(5,1) * t158;
t246 = (rSges(5,2) * t159 + t235) * t175;
t171 = sin(qJ(3));
t174 = cos(qJ(3));
t245 = (rSges(4,1) * t171 + rSges(4,2) * t174) * t175;
t167 = t172 ^ 2;
t168 = t175 ^ 2;
t176 = -pkin(8) - pkin(7);
t244 = t172 / 0.2e1;
t243 = t175 / 0.2e1;
t242 = pkin(3) * t171;
t241 = pkin(3) * t174;
t240 = pkin(4) * t159;
t239 = pkin(5) * t154;
t75 = Icges(7,3) * t154 + (Icges(7,5) * t173 - Icges(7,6) * t170) * t155;
t77 = Icges(7,5) * t154 + (Icges(7,1) * t173 - Icges(7,4) * t170) * t155;
t238 = t155 * t173 * t77 + t154 * t75;
t212 = t172 * t176 + t175 * t242;
t133 = pkin(4) * t158 + t242;
t166 = -qJ(5) + t176;
t214 = -t133 * t175 - t166 * t172;
t67 = t175 * (t212 + t214);
t237 = t67 + t175 * (rSges(6,3) * t172 - t247);
t124 = t172 * t133;
t219 = t171 * t172;
t151 = pkin(3) * t219;
t71 = t124 - t151 + (-t166 + t176) * t175;
t224 = t155 * t172;
t225 = t154 * t172;
t89 = rSges(6,1) * t225 + rSges(6,2) * t224 + rSges(6,3) * t175;
t236 = -t71 - t89;
t76 = Icges(7,6) * t154 + (Icges(7,4) * t173 - Icges(7,2) * t170) * t155;
t234 = t170 * t76;
t216 = t173 * t175;
t221 = t170 * t172;
t112 = -t154 * t221 + t216;
t218 = t172 * t173;
t220 = t170 * t175;
t113 = t154 * t218 + t220;
t55 = Icges(7,5) * t113 + Icges(7,6) * t112 - Icges(7,3) * t224;
t57 = Icges(7,4) * t113 + Icges(7,2) * t112 - Icges(7,6) * t224;
t59 = Icges(7,1) * t113 + Icges(7,4) * t112 - Icges(7,5) * t224;
t24 = t154 * t55 + (-t170 * t57 + t173 * t59) * t155;
t233 = t24 * t175;
t114 = t154 * t220 + t218;
t115 = -t154 * t216 + t221;
t223 = t155 * t175;
t56 = Icges(7,5) * t115 + Icges(7,6) * t114 + Icges(7,3) * t223;
t58 = Icges(7,4) * t115 + Icges(7,2) * t114 + Icges(7,6) * t223;
t60 = Icges(7,1) * t115 + Icges(7,4) * t114 + Icges(7,5) * t223;
t25 = t154 * t56 + (-t170 * t58 + t173 * t60) * t155;
t232 = t25 * t172;
t231 = Icges(4,4) * t171;
t230 = Icges(4,4) * t174;
t222 = t159 * t172;
t217 = t172 * t174;
t215 = rSges(7,1) * t113 + rSges(7,2) * t112;
t147 = t167 + t168;
t213 = (m(6) + m(7)) * t147;
t122 = rSges(6,1) * t155 - rSges(6,2) * t154;
t146 = pkin(4) * t222;
t80 = t122 * t172 + t146;
t211 = pkin(1) * t175 + qJ(2) * t172;
t199 = -rSges(7,1) * t115 - rSges(7,2) * t114;
t62 = rSges(7,3) * t223 - t199;
t210 = t175 * t62 + t67 + t168 * (pkin(9) * t155 - t239);
t136 = pkin(5) * t225;
t61 = -rSges(7,3) * t224 + t215;
t209 = pkin(9) * t224 - t136 - t61 - t71;
t49 = t172 * t252 + t146;
t161 = t175 * qJ(2);
t208 = t161 - t214;
t99 = rSges(5,2) * t222 + rSges(5,3) * t175 + t172 * t235;
t207 = rSges(4,1) * t219 + rSges(4,2) * t217 + rSges(4,3) * t175;
t206 = (-rSges(7,3) - pkin(9)) * t155;
t18 = t114 * t57 + t115 * t59 + t223 * t55;
t19 = t114 * t58 + t115 * t60 + t223 * t56;
t10 = t172 * t19 + t175 * t18;
t16 = t112 * t57 + t113 * t59 - t224 * t55;
t17 = t112 * t58 + t113 * t60 - t224 * t56;
t29 = t112 * t76 + t113 * t77 - t224 * t75;
t3 = t154 * t29 + (-t16 * t172 + t17 * t175) * t155;
t30 = t114 * t76 + t115 * t77 + t223 * t75;
t4 = t154 * t30 + (-t172 * t18 + t175 * t19) * t155;
t9 = t16 * t175 + t17 * t172;
t205 = t3 * t243 + t4 * t244 + t154 * (t232 + t233) / 0.2e1 - t9 * t224 / 0.2e1 + t10 * t223 / 0.2e1;
t204 = -t240 - t241;
t203 = (t167 * t253 + t10) * t172 + (t9 + t254 * t168 + (t172 * t254 + t175 * t253) * t172) * t175;
t130 = rSges(5,1) * t159 - rSges(5,2) * t158;
t152 = pkin(3) * t217;
t91 = t130 * t172 + t152;
t92 = (-t130 - t241) * t175;
t194 = t172 * t91 - t175 * t92;
t193 = Icges(4,1) * t171 + t230;
t192 = Icges(5,1) * t158 + t228;
t191 = Icges(6,1) * t154 + t226;
t190 = Icges(4,2) * t174 + t231;
t189 = Icges(5,2) * t159 + t229;
t188 = Icges(6,2) * t155 + t227;
t187 = Icges(4,5) * t171 + Icges(4,6) * t174;
t180 = -t175 * t166 + t124 + t211;
t73 = t161 + t245 + (-rSges(4,3) - pkin(1) - pkin(7)) * t172;
t74 = pkin(7) * t175 + t207 + t211;
t179 = m(4) * (t172 * t73 - t175 * t74);
t64 = t161 + t246 + (-rSges(5,3) - pkin(1)) * t172 + t212;
t65 = -t175 * t176 + t151 + t211 + t99;
t178 = m(5) * (t172 * t64 - t175 * t65);
t177 = t233 / 0.2e1 + t232 / 0.2e1 + (-t154 * (Icges(6,6) * t172 - t175 * t188) + t155 * (Icges(6,5) * t172 - t175 * t191) - t158 * (Icges(5,6) * t172 - t175 * t189) + t159 * (Icges(5,5) * t172 - t175 * t192) + t251 * t172 - t250 * t175 + t30) * t244 + (-t154 * (Icges(6,6) * t175 + t172 * t188) + t155 * (Icges(6,5) * t175 + t172 * t191) - t158 * (Icges(5,6) * t175 + t172 * t189) + t159 * (Icges(5,5) * t175 + t172 * t192) + t250 * t172 + t251 * t175 + t29) * t243;
t143 = rSges(2,1) * t175 - rSges(2,2) * t172;
t142 = rSges(4,1) * t174 - rSges(4,2) * t171;
t141 = -rSges(2,1) * t172 - rSges(2,2) * t175;
t118 = t151 + (-pkin(7) - t176) * t175;
t117 = -rSges(3,2) * t175 + rSges(3,3) * t172 + t211;
t116 = rSges(3,3) * t175 + t161 + (rSges(3,2) - pkin(1)) * t172;
t106 = Icges(4,3) * t172 - t175 * t187;
t105 = Icges(4,3) * t175 + t172 * t187;
t101 = t175 * (-pkin(7) * t172 - t212);
t88 = t175 * (rSges(5,3) * t172 - t246);
t81 = (-t122 - t240) * t175;
t70 = (-t122 + t204) * t175;
t69 = t152 + t80;
t63 = -t172 * t207 + (t172 * rSges(4,3) - t245) * t175;
t54 = -t172 * t99 + t88;
t52 = t180 + t89;
t51 = t247 + (-rSges(6,3) - pkin(1)) * t172 + t208;
t50 = (-t240 - t252) * t175;
t44 = (t204 - t252) * t175;
t43 = t152 + t49;
t38 = t101 + t88 + (-t118 - t99) * t172;
t37 = -t154 * t62 + t223 * t78;
t36 = t154 * t61 + t224 * t78;
t35 = t172 * t206 + t136 + t180 + t215;
t34 = -pkin(1) * t172 + (t206 + t239) * t175 + t199 + t208;
t33 = (-t155 * t234 + t238) * t154;
t32 = (-t172 * t62 - t175 * t61) * t155;
t31 = t172 * t236 + t237;
t26 = t101 + (-t118 + t236) * t172 + t237;
t13 = t172 * t209 + t210;
t12 = t101 + (-t118 + t209) * t172 + t210;
t1 = [-t154 * t120 - t158 * t128 + t159 * t129 - t171 * (-Icges(4,2) * t171 + t230) + t174 * (Icges(4,1) * t174 - t231) + Icges(3,1) + Icges(2,3) + (t121 - t234) * t155 + m(7) * (t34 ^ 2 + t35 ^ 2) + m(6) * (t51 ^ 2 + t52 ^ 2) + m(5) * (t64 ^ 2 + t65 ^ 2) + m(4) * (t73 ^ 2 + t74 ^ 2) + m(3) * (t116 ^ 2 + t117 ^ 2) + m(2) * (t141 ^ 2 + t143 ^ 2) + t238; m(7) * (t172 * t34 - t175 * t35) + m(6) * (t172 * t51 - t175 * t52) + t178 + t179 + m(3) * (t116 * t172 - t117 * t175); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1) * t147 + t213; m(7) * (t34 * t43 + t35 * t44) + m(6) * (t51 * t69 + t52 * t70) + m(5) * (t64 * t91 + t65 * t92) + t177 + (-(Icges(4,6) * t175 + t172 * t190) * t171 + (Icges(4,5) * t175 + t172 * t193) * t174) * t243 + (-(Icges(4,6) * t172 - t175 * t190) * t171 + (Icges(4,5) * t172 - t175 * t193) * t174) * t244 + (t168 / 0.2e1 + t167 / 0.2e1) * (Icges(4,5) * t174 - Icges(4,6) * t171) + t142 * t179; m(5) * t194 + m(6) * (t172 * t69 - t175 * t70) + m(7) * (t172 * t43 - t175 * t44) + m(4) * t147 * t142; m(7) * (t12 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t26 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t38 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(4) * (t142 ^ 2 * t147 + t63 ^ 2) + t172 * (t105 * t248 + t167 * t106) + t175 * (t168 * t105 + t106 * t248) + t203; m(7) * (t34 * t49 + t35 * t50) + m(6) * (t51 * t80 + t52 * t81) + t177 + t130 * t178; m(6) * (t172 * t80 - t175 * t81) + m(7) * (t172 * t49 - t175 * t50) + m(5) * t147 * t130; m(7) * (t12 * t13 + t43 * t49 + t44 * t50) + m(6) * (t26 * t31 + t69 * t80 + t70 * t81) + m(5) * (t130 * t194 + t38 * t54) + t203; m(7) * (t13 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t31 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t130 ^ 2 * t147 + t54 ^ 2) + t203; m(7) * (t172 * t35 + t175 * t34) + m(6) * (t172 * t52 + t175 * t51); 0; m(7) * (t172 * t44 + t175 * t43) + m(6) * (t172 * t70 + t175 * t69); m(7) * (t172 * t50 + t175 * t49) + m(6) * (t172 * t81 + t175 * t80); t213; m(7) * (t34 * t37 + t35 * t36) + t33 + ((t25 / 0.2e1 + t30 / 0.2e1) * t175 + (-t24 / 0.2e1 - t29 / 0.2e1) * t172) * t155; m(7) * (t172 * t37 - t175 * t36); m(7) * (t12 * t32 + t36 * t44 + t37 * t43) + t205; m(7) * (t13 * t32 + t36 * t50 + t37 * t49) + t205; m(7) * (t172 * t36 + t175 * t37); t154 * t33 + m(7) * (t32 ^ 2 + t36 ^ 2 + t37 ^ 2) + (-t172 * t3 + t175 * t4 + t154 * (-t24 * t172 + t25 * t175)) * t155;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
