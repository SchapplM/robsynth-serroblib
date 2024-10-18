% Calculate joint inertia matrix for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:07:16
% EndTime: 2024-09-28 18:07:19
% DurationCPUTime: 0.48s
% Computational Cost: add. (24911->442), mult. (30322->649), div. (0->0), fcn. (34780->26), ass. (0->239)
t199 = sin(pkin(11));
t198 = qJ(2) + qJ(3);
t194 = cos(t198);
t240 = pkin(3) * t194;
t275 = t240 * t199;
t202 = cos(pkin(11));
t274 = t240 * t202;
t204 = cos(pkin(5));
t200 = sin(pkin(6));
t266 = pkin(10) * t200;
t236 = t202 * t266;
t164 = -t199 * pkin(4) + t204 * t236;
t237 = t199 * t266;
t250 = t202 * t204;
t166 = pkin(4) * t250 + t237;
t206 = sin(qJ(4));
t210 = cos(qJ(4));
t115 = pkin(3) * t250 + t164 * t206 + t166 * t210;
t207 = sin(qJ(3));
t208 = sin(qJ(2));
t211 = cos(qJ(3));
t212 = cos(qJ(2));
t216 = pkin(3) * t199 - t164 * t210 + t166 * t206;
t189 = t211 * pkin(3) + pkin(2);
t213 = pkin(7) + pkin(8);
t197 = pkin(9) + t213;
t201 = sin(pkin(5));
t234 = t212 * t207 * pkin(3);
t144 = -t201 * t197 + (t208 * t189 + t234) * t204;
t203 = cos(pkin(6));
t177 = t203 * pkin(10) + t197;
t246 = t204 * t208;
t235 = pkin(2) * t246;
t218 = -t201 * t177 - t144 + t235;
t195 = qJ(4) + t198;
t187 = sin(t195);
t188 = cos(t195);
t209 = cos(qJ(5));
t205 = sin(qJ(5));
t247 = t204 * t205;
t230 = t203 * t247;
t254 = t200 * t201;
t232 = t205 * t254;
t245 = t204 * t209;
t249 = t203 * t205;
t110 = (t199 * t209 + t202 * t230) * t188 + (-t199 * t249 + t202 * t245) * t187 - t202 * t232;
t229 = t203 * t245;
t231 = t209 * t254;
t248 = t203 * t209;
t111 = (-t199 * t205 + t202 * t229) * t188 + (-t199 * t248 - t202 * t247) * t187 - t202 * t231;
t252 = t201 * t203;
t257 = t199 * t187;
t133 = -t202 * t252 + (-t188 * t250 + t257) * t200;
t66 = t110 * rSges(6,1) + t111 * rSges(6,2) + t133 * rSges(6,3);
t273 = (t115 * t207 + t216 * t211) * t212 + (t115 * t211 - t216 * t207) * t208 + t218 * t202 - t275 + t66;
t163 = t202 * pkin(4) + t204 * t237;
t255 = t199 * t204;
t165 = pkin(4) * t255 - t236;
t114 = -pkin(3) * t255 - t163 * t206 - t165 * t210;
t217 = t202 * pkin(3) + t163 * t210 - t165 * t206;
t108 = (-t199 * t230 + t202 * t209) * t188 + (-t199 * t245 - t202 * t249) * t187 + t199 * t232;
t109 = (-t199 * t229 - t202 * t205) * t188 + (t199 * t247 - t202 * t248) * t187 + t199 * t231;
t251 = t202 * t187;
t132 = t200 * t251 + (t188 * t200 * t204 + t252) * t199;
t65 = t108 * rSges(6,1) + t109 * rSges(6,2) + t132 * rSges(6,3);
t272 = (t114 * t207 + t217 * t211) * t212 + (t114 * t211 - t217 * t207) * t208 - t218 * t199 - t274 + t65;
t253 = t201 * t202;
t192 = pkin(5) - t198;
t186 = -qJ(4) + t192;
t176 = cos(t186) / 0.2e1;
t191 = pkin(5) + t198;
t185 = qJ(4) + t191;
t182 = cos(t185);
t157 = t176 + t182 / 0.2e1;
t140 = t202 * t157 - t257;
t175 = sin(t185) / 0.2e1;
t181 = sin(t186);
t156 = t175 - t181 / 0.2e1;
t141 = t202 * t156 + t199 * t188;
t89 = Icges(5,5) * t141 + Icges(5,6) * t140 - Icges(5,3) * t253;
t271 = t89 * t253;
t179 = cos(t192) / 0.2e1;
t184 = cos(t191);
t169 = t179 + t184 / 0.2e1;
t193 = sin(t198);
t145 = t202 * t169 - t199 * t193;
t178 = sin(t191) / 0.2e1;
t183 = sin(t192);
t168 = t178 - t183 / 0.2e1;
t146 = t202 * t168 + t199 * t194;
t100 = Icges(4,5) * t146 + Icges(4,6) * t145 - Icges(4,3) * t253;
t147 = -t199 * t169 - t202 * t193;
t148 = -t199 * t168 + t202 * t194;
t256 = t199 * t201;
t101 = Icges(4,5) * t148 + Icges(4,6) * t147 + Icges(4,3) * t256;
t102 = Icges(4,4) * t146 + Icges(4,2) * t145 - Icges(4,6) * t253;
t103 = Icges(4,4) * t148 + Icges(4,2) * t147 + Icges(4,6) * t256;
t104 = Icges(4,1) * t146 + Icges(4,4) * t145 - Icges(4,5) * t253;
t105 = Icges(4,1) * t148 + Icges(4,4) * t147 + Icges(4,5) * t256;
t167 = t178 + t183 / 0.2e1;
t170 = t179 - t184 / 0.2e1;
t134 = Icges(4,5) * t170 + Icges(4,6) * t167 + Icges(4,3) * t204;
t135 = Icges(4,4) * t170 + Icges(4,2) * t167 + Icges(4,6) * t204;
t136 = Icges(4,1) * t170 + Icges(4,4) * t167 + Icges(4,5) * t204;
t155 = t175 + t181 / 0.2e1;
t158 = t176 - t182 / 0.2e1;
t125 = Icges(5,5) * t158 + Icges(5,6) * t155 + Icges(5,3) * t204;
t126 = Icges(5,4) * t158 + Icges(5,2) * t155 + Icges(5,6) * t204;
t127 = Icges(5,1) * t158 + Icges(5,4) * t155 + Icges(5,5) * t204;
t59 = Icges(6,5) * t108 + Icges(6,6) * t109 + Icges(6,3) * t132;
t61 = Icges(6,4) * t108 + Icges(6,2) * t109 + Icges(6,6) * t132;
t63 = Icges(6,1) * t108 + Icges(6,4) * t109 + Icges(6,5) * t132;
t31 = t110 * t63 + t111 * t61 + t133 * t59;
t60 = Icges(6,5) * t110 + Icges(6,6) * t111 + Icges(6,3) * t133;
t62 = Icges(6,4) * t110 + Icges(6,2) * t111 + Icges(6,6) * t133;
t64 = Icges(6,1) * t110 + Icges(6,4) * t111 + Icges(6,5) * t133;
t32 = t110 * t64 + t111 * t62 + t133 * t60;
t138 = t200 * t245 + (-t187 * t205 + t188 * t248) * t201;
t139 = t200 * t247 + (t187 * t209 + t188 * t249) * t201;
t153 = -t188 * t254 + t204 * t203;
t75 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t153;
t76 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t153;
t77 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t153;
t37 = t110 * t77 + t111 * t76 + t133 * t75;
t9 = t37 * t204 + (t199 * t31 - t202 * t32) * t201;
t142 = -t199 * t157 - t251;
t143 = -t199 * t156 + t202 * t188;
t90 = Icges(5,5) * t143 + Icges(5,6) * t142 + Icges(5,3) * t256;
t91 = Icges(5,4) * t141 + Icges(5,2) * t140 - Icges(5,6) * t253;
t92 = Icges(5,4) * t143 + Icges(5,2) * t142 + Icges(5,6) * t256;
t93 = Icges(5,1) * t141 + Icges(5,4) * t140 - Icges(5,5) * t253;
t94 = Icges(5,1) * t143 + Icges(5,4) * t142 + Icges(5,5) * t256;
t267 = -(t140 * t92 + t141 * t94 - t90 * t253) * t256 + (t140 * t91 + t141 * t93 - t271) * t253 - (-t125 * t253 + t140 * t126 + t141 * t127) * t204 - t9;
t270 = -(-t101 * t253 + t145 * t103 + t146 * t105) * t256 + (-t100 * t253 + t145 * t102 + t146 * t104) * t253 - (-t134 * t253 + t145 * t135 + t146 * t136) * t204 + t267;
t269 = t204 * t125;
t265 = t212 * pkin(2);
t264 = -pkin(2) + t189;
t263 = t272 * t204;
t172 = pkin(4) * t210 + t206 * t266 + pkin(3);
t173 = -t206 * pkin(4) + t210 * t266;
t78 = t139 * rSges(6,1) + t138 * rSges(6,2) + t153 * rSges(6,3);
t261 = -(-t197 + t177) * t204 - ((-t173 * t211 + (-pkin(3) + t172) * t207) * t212 + (t172 * t211 + t173 * t207 - t264) * t208) * t201 - t78;
t171 = -t201 * t213 + t235;
t241 = t144 - t171;
t85 = t241 * t202 + t275;
t86 = -t241 * t199 + t274;
t260 = t86 * t253 + t85 * t256;
t81 = t204 * t86;
t96 = t143 * rSges(5,1) + t142 * rSges(5,2) + rSges(5,3) * t256;
t84 = t204 * t96;
t259 = t81 + t84;
t95 = t141 * rSges(5,1) + t140 * rSges(5,2) - rSges(5,3) * t253;
t67 = t96 * t253 + t95 * t256;
t258 = -t85 - t95;
t106 = t146 * rSges(4,1) + t145 * rSges(4,2) - rSges(4,3) * t253;
t107 = t148 * rSges(4,1) + t147 * rSges(4,2) + rSges(4,3) * t256;
t68 = t106 * t256 + t107 * t253;
t244 = t204 * t212;
t226 = pkin(7) * t201 + t171;
t130 = t265 * t199 + t226 * t202;
t131 = -t226 * t199 + t265 * t202;
t243 = t130 * t256 + t131 * t253;
t128 = (t197 - t213) * t204 + (t264 * t208 + t234) * t201;
t129 = t158 * rSges(5,1) + t155 * rSges(5,2) + t204 * rSges(5,3);
t242 = -t128 - t129;
t239 = t81 + t263;
t238 = -t85 - t273;
t233 = -t128 + t261;
t228 = t261 * t201;
t34 = t138 * t61 + t139 * t63 + t153 * t59;
t35 = t138 * t62 + t139 * t64 + t153 * t60;
t40 = t138 * t76 + t139 * t77 + t153 * t75;
t14 = t40 * t204 + (t199 * t34 - t202 * t35) * t201;
t29 = t108 * t63 + t109 * t61 + t132 * t59;
t30 = t108 * t64 + t109 * t62 + t132 * t60;
t36 = t108 * t77 + t109 * t76 + t132 * t75;
t8 = t36 * t204 + (t199 * t29 - t202 * t30) * t201;
t227 = (t14 + ((t155 * t92 + t158 * t94) * t199 - (t155 * t91 + t158 * t93) * t202) * t201 + (t155 * t126 + t158 * t127 + (t90 * t199 - t89 * t202) * t201 + t269) * t204) * t204 + (t8 - (t142 * t91 + t143 * t93) * t253 + (t142 * t126 + t143 * t127) * t204 + (t142 * t92 + t143 * t94 + t90 * t256 + t269 - t271) * t256) * t256;
t17 = t272 * t253 + t273 * t256;
t43 = t67 + t260;
t225 = t242 * t201;
t137 = t170 * rSges(4,1) + t167 * rSges(4,2) + t204 * rSges(4,3);
t154 = t201 * t208 * pkin(2) + (-pkin(7) + t213) * t204;
t224 = (-t137 - t154) * t201;
t223 = t233 * t201;
t222 = (-t154 + t242) * t201;
t11 = t34 * t132 + t35 * t133 + t40 * t153;
t3 = t29 * t132 + t30 * t133 + t36 * t153;
t4 = t31 * t132 + t32 * t133 + t37 * t153;
t221 = t3 * t256 / 0.2e1 + t204 * t11 / 0.2e1 + t153 * t14 / 0.2e1 - t4 * t253 / 0.2e1 + t132 * t8 / 0.2e1 + t133 * t9 / 0.2e1;
t220 = ((t101 * t256 + t147 * t103 + t148 * t105) * t256 - (t100 * t256 + t147 * t102 + t148 * t104) * t253 + (t134 * t256 + t147 * t135 + t148 * t136) * t204) * t256 + t204 * ((t204 * t134 + t167 * t135 + t170 * t136) * t204 + ((t204 * t101 + t167 * t103 + t170 * t105) * t199 - (t204 * t100 + t167 * t102 + t170 * t104) * t202) * t201) + t227;
t16 = t17 + t260;
t219 = (-t154 + t233) * t201;
t215 = t267 * t253 + t227;
t214 = t270 * t253 + t220;
t162 = -t199 * t246 + t202 * t212;
t161 = -t199 * t244 - t202 * t208;
t160 = t199 * t212 + t202 * t246;
t159 = -t199 * t208 + t202 * t244;
t152 = t204 * rSges(3,3) + (rSges(3,1) * t208 + rSges(3,2) * t212) * t201;
t151 = Icges(3,5) * t204 + (Icges(3,1) * t208 + Icges(3,4) * t212) * t201;
t150 = Icges(3,6) * t204 + (Icges(3,4) * t208 + Icges(3,2) * t212) * t201;
t149 = Icges(3,3) * t204 + (Icges(3,5) * t208 + Icges(3,6) * t212) * t201;
t124 = t204 * t131;
t123 = t162 * rSges(3,1) + t161 * rSges(3,2) + rSges(3,3) * t256;
t122 = t160 * rSges(3,1) + t159 * rSges(3,2) - rSges(3,3) * t253;
t121 = Icges(3,1) * t162 + Icges(3,4) * t161 + Icges(3,5) * t256;
t120 = Icges(3,1) * t160 + Icges(3,4) * t159 - Icges(3,5) * t253;
t119 = Icges(3,4) * t162 + Icges(3,2) * t161 + Icges(3,6) * t256;
t118 = Icges(3,4) * t160 + Icges(3,2) * t159 - Icges(3,6) * t253;
t117 = Icges(3,5) * t162 + Icges(3,6) * t161 + Icges(3,3) * t256;
t116 = Icges(3,5) * t160 + Icges(3,6) * t159 - Icges(3,3) * t253;
t99 = t204 * t107;
t88 = -t204 * t122 - t152 * t253;
t87 = t204 * t123 - t152 * t256;
t74 = (t122 * t199 + t123 * t202) * t201;
t73 = -t204 * t106 - t137 * t253;
t72 = -t137 * t256 + t99;
t70 = -t129 * t253 - t204 * t95;
t69 = -t129 * t256 + t84;
t55 = (-t106 - t130) * t204 + t202 * t224;
t54 = t199 * t224 + t124 + t99;
t53 = t243 + t68;
t52 = t202 * t225 + t258 * t204;
t51 = t199 * t225 + t259;
t45 = (-t130 + t258) * t204 + t202 * t222;
t44 = t199 * t222 + t124 + t259;
t42 = t133 * t78 - t153 * t66;
t41 = -t132 * t78 + t153 * t65;
t39 = t43 + t243;
t38 = t132 * t66 - t133 * t65;
t28 = t202 * t228 - t204 * t273;
t27 = t199 * t228 + t263;
t21 = t202 * t223 + t238 * t204;
t20 = t199 * t223 + t239;
t19 = (-t130 + t238) * t204 + t202 * t219;
t18 = t199 * t219 + t124 + t239;
t15 = t16 + t243;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t74 + m(4) * t53 + m(5) * t39 + m(6) * t15; ((t117 * t256 + t161 * t119 + t162 * t121) * t256 + (t149 * t256 + t161 * t150 + t162 * t151) * t204) * t256 + m(6) * (t15 ^ 2 + t18 ^ 2 + t19 ^ 2) + m(5) * (t39 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(3) * (t74 ^ 2 + t87 ^ 2 + t88 ^ 2) + t204 * (t204 ^ 2 * t149 + (((t119 * t212 + t121 * t208) * t199 - (t118 * t212 + t120 * t208) * t202) * t201 + (-t116 * t202 + t117 * t199 + t150 * t212 + t151 * t208) * t204) * t201) + t220 + (-(t159 * t150 + t160 * t151) * t204 + (-t116 * t253 + t159 * t118 + t160 * t120 + t204 * t149) * t253 + (-t116 * t256 + t117 * t253 - t161 * t118 - t159 * t119 - t162 * t120 - t160 * t121) * t256 + t270) * t253; m(4) * t68 + m(5) * t43 + m(6) * t16; m(6) * (t16 * t15 + t20 * t18 + t21 * t19) + m(5) * (t43 * t39 + t51 * t44 + t52 * t45) + m(4) * (t68 * t53 + t72 * t54 + t73 * t55) + t214; m(6) * (t16 ^ 2 + t20 ^ 2 + t21 ^ 2) + m(5) * (t43 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(4) * (t68 ^ 2 + t72 ^ 2 + t73 ^ 2) + t214; m(5) * t67 + m(6) * t17; m(6) * (t17 * t15 + t27 * t18 + t28 * t19) + m(5) * (t67 * t39 + t69 * t44 + t70 * t45) + t215; m(6) * (t17 * t16 + t27 * t20 + t28 * t21) + m(5) * (t67 * t43 + t69 * t51 + t70 * t52) + t215; m(5) * (t67 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t17 ^ 2 + t27 ^ 2 + t28 ^ 2) + t215; m(6) * t38; m(6) * (t38 * t15 + t41 * t18 + t42 * t19) + t221; m(6) * (t38 * t16 + t41 * t20 + t42 * t21) + t221; m(6) * (t38 * t17 + t41 * t27 + t42 * t28) + t221; m(6) * (t38 ^ 2 + t41 ^ 2 + t42 ^ 2) + t132 * t3 + t133 * t4 + t153 * t11;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
