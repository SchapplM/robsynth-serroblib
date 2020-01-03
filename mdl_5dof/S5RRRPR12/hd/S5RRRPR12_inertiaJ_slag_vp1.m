% Calculate joint inertia matrix for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR12_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR12_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR12_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:06
% EndTime: 2019-12-31 21:37:15
% DurationCPUTime: 3.20s
% Computational Cost: add. (13319->485), mult. (28865->701), div. (0->0), fcn. (36662->12), ass. (0->232)
t224 = cos(pkin(5));
t229 = cos(qJ(2));
t230 = cos(qJ(1));
t259 = t229 * t230;
t227 = sin(qJ(2));
t228 = sin(qJ(1));
t264 = t227 * t228;
t206 = -t224 * t264 + t259;
t226 = sin(qJ(3));
t222 = sin(pkin(5));
t273 = cos(qJ(3));
t241 = t222 * t273;
t186 = t206 * t226 - t228 * t241;
t262 = t228 * t222;
t187 = t206 * t273 + t226 * t262;
t223 = cos(pkin(10));
t216 = pkin(4) * t223 + pkin(3);
t225 = -pkin(9) - qJ(4);
t261 = t228 * t229;
t263 = t227 * t230;
t205 = t224 * t261 + t263;
t221 = sin(pkin(10));
t267 = t205 * t221;
t220 = pkin(10) + qJ(5);
t217 = sin(t220);
t218 = cos(t220);
t143 = -t187 * t217 + t205 * t218;
t144 = t187 * t218 + t205 * t217;
t87 = t144 * rSges(6,1) + t143 * rSges(6,2) + t186 * rSges(6,3);
t278 = pkin(4) * t267 - t186 * t225 + t187 * t216 + t87;
t204 = t224 * t263 + t261;
t184 = t204 * t226 + t230 * t241;
t266 = t222 * t227;
t201 = -t224 * t273 + t226 * t266;
t258 = t230 * t222;
t185 = t204 * t273 - t226 * t258;
t203 = -t224 * t259 + t264;
t141 = -t185 * t217 + t203 * t218;
t142 = t185 * t218 + t203 * t217;
t78 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t184;
t80 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t184;
t82 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t184;
t30 = t143 * t80 + t144 * t82 + t186 * t78;
t79 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t186;
t81 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t186;
t83 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t186;
t31 = t143 * t81 + t144 * t83 + t186 * t79;
t202 = t224 * t226 + t227 * t241;
t265 = t222 * t229;
t173 = -t202 * t217 - t218 * t265;
t174 = t202 * t218 - t217 * t265;
t109 = Icges(6,5) * t174 + Icges(6,6) * t173 + Icges(6,3) * t201;
t110 = Icges(6,4) * t174 + Icges(6,2) * t173 + Icges(6,6) * t201;
t111 = Icges(6,1) * t174 + Icges(6,4) * t173 + Icges(6,5) * t201;
t43 = t109 * t186 + t110 * t143 + t111 * t144;
t2 = t184 * t30 + t186 * t31 + t201 * t43;
t277 = t2 / 0.2e1;
t276 = t184 / 0.2e1;
t275 = t186 / 0.2e1;
t274 = t201 / 0.2e1;
t272 = -pkin(3) + t216;
t176 = t184 * qJ(4);
t268 = t203 * t221;
t249 = pkin(4) * t268;
t236 = -rSges(6,1) * t142 - rSges(6,2) * t141;
t86 = rSges(6,3) * t184 - t236;
t271 = -t184 * t225 + t185 * t272 - t176 + t249 + t86;
t136 = t187 * pkin(3) + qJ(4) * t186;
t270 = -t136 + t278;
t152 = -t187 * t221 + t205 * t223;
t153 = t187 * t223 + t267;
t97 = t153 * rSges(5,1) + t152 * rSges(5,2) + t186 * rSges(5,3);
t269 = -t136 - t97;
t260 = t228 * t230;
t112 = rSges(6,1) * t174 + rSges(6,2) * t173 + rSges(6,3) * t201;
t245 = t221 * t265;
t257 = t112 - pkin(4) * t245 + t272 * t202 + (-qJ(4) - t225) * t201;
t182 = -t202 * t221 - t223 * t265;
t183 = t202 * t223 - t245;
t118 = rSges(5,1) * t183 + rSges(5,2) * t182 + rSges(5,3) * t201;
t169 = pkin(3) * t202 + qJ(4) * t201;
t256 = -t118 - t169;
t135 = pkin(3) * t185 + t176;
t255 = t135 * t265 + t203 * t169;
t172 = t206 * pkin(2) + pkin(8) * t205;
t168 = t224 * t172;
t254 = t224 * t136 + t168;
t171 = pkin(2) * t204 + t203 * pkin(8);
t253 = -t135 - t171;
t252 = t171 * t262 + t172 * t258;
t251 = t230 * pkin(1) + pkin(7) * t262;
t52 = t109 * t201 + t110 * t173 + t111 * t174;
t115 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t201;
t116 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t201;
t117 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t201;
t55 = t201 * t115 + t116 * t182 + t117 * t183;
t154 = Icges(4,5) * t202 - Icges(4,6) * t201 - Icges(4,3) * t265;
t155 = Icges(4,4) * t202 - Icges(4,2) * t201 - Icges(4,6) * t265;
t156 = Icges(4,1) * t202 - Icges(4,4) * t201 - Icges(4,5) * t265;
t75 = -t154 * t265 - t201 * t155 + t202 * t156;
t250 = -t52 - t55 - t75;
t248 = -t136 - t270;
t36 = t173 * t80 + t174 * t82 + t201 * t78;
t42 = t109 * t184 + t110 * t141 + t111 * t142;
t247 = t42 / 0.2e1 + t36 / 0.2e1;
t37 = t173 * t81 + t174 * t83 + t201 * t79;
t246 = t43 / 0.2e1 + t37 / 0.2e1;
t244 = -t169 - t257;
t126 = t187 * rSges(4,1) - t186 * rSges(4,2) + t205 * rSges(4,3);
t191 = Icges(3,3) * t224 + (Icges(3,5) * t227 + Icges(3,6) * t229) * t222;
t192 = Icges(3,6) * t224 + (Icges(3,4) * t227 + Icges(3,2) * t229) * t222;
t193 = Icges(3,5) * t224 + (Icges(3,1) * t227 + Icges(3,4) * t229) * t222;
t242 = t224 * t191 + t192 * t265 + t193 * t266;
t165 = t206 * rSges(3,1) - t205 * rSges(3,2) + rSges(3,3) * t262;
t240 = -t228 * pkin(1) + pkin(7) * t258;
t157 = rSges(4,1) * t202 - rSges(4,2) * t201 - rSges(4,3) * t265;
t207 = (pkin(2) * t227 - pkin(8) * t229) * t222;
t239 = t222 * (-t157 - t207);
t238 = t135 * t262 + t136 * t258 + t252;
t237 = t222 * (-t207 + t256);
t235 = t172 + t251;
t234 = t222 * (-t207 + t244);
t233 = -t171 + t240;
t125 = rSges(4,1) * t185 - rSges(4,2) * t184 + rSges(4,3) * t203;
t150 = -t185 * t221 + t203 * t223;
t151 = t185 * t223 + t268;
t96 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t184;
t164 = rSges(3,1) * t204 - rSges(3,2) * t203 - rSges(3,3) * t258;
t90 = Icges(5,5) * t151 + Icges(5,6) * t150 + Icges(5,3) * t184;
t92 = Icges(5,4) * t151 + Icges(5,2) * t150 + Icges(5,6) * t184;
t94 = Icges(5,1) * t151 + Icges(5,4) * t150 + Icges(5,5) * t184;
t39 = t182 * t92 + t183 * t94 + t201 * t90;
t46 = t115 * t184 + t116 * t150 + t117 * t151;
t119 = Icges(4,5) * t185 - Icges(4,6) * t184 + Icges(4,3) * t203;
t121 = Icges(4,4) * t185 - Icges(4,2) * t184 + Icges(4,6) * t203;
t123 = Icges(4,1) * t185 - Icges(4,4) * t184 + Icges(4,5) * t203;
t64 = -t119 * t265 - t121 * t201 + t123 * t202;
t69 = t154 * t203 - t155 * t184 + t156 * t185;
t232 = t64 / 0.2e1 + t39 / 0.2e1 + t69 / 0.2e1 + t46 / 0.2e1 + t247;
t91 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t186;
t93 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t186;
t95 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t186;
t40 = t182 * t93 + t183 * t95 + t201 * t91;
t47 = t115 * t186 + t116 * t152 + t117 * t153;
t120 = Icges(4,5) * t187 - Icges(4,6) * t186 + Icges(4,3) * t205;
t122 = Icges(4,4) * t187 - Icges(4,2) * t186 + Icges(4,6) * t205;
t124 = Icges(4,1) * t187 - Icges(4,4) * t186 + Icges(4,5) * t205;
t65 = -t120 * t265 - t122 * t201 + t124 * t202;
t70 = t154 * t205 - t155 * t186 + t156 * t187;
t231 = t70 / 0.2e1 + t47 / 0.2e1 + t65 / 0.2e1 + t40 / 0.2e1 + t246;
t209 = rSges(2,1) * t230 - rSges(2,2) * t228;
t208 = -rSges(2,1) * t228 - rSges(2,2) * t230;
t194 = rSges(3,3) * t224 + (rSges(3,1) * t227 + rSges(3,2) * t229) * t222;
t163 = Icges(3,1) * t206 - Icges(3,4) * t205 + Icges(3,5) * t262;
t162 = Icges(3,1) * t204 - Icges(3,4) * t203 - Icges(3,5) * t258;
t161 = Icges(3,4) * t206 - Icges(3,2) * t205 + Icges(3,6) * t262;
t160 = Icges(3,4) * t204 - Icges(3,2) * t203 - Icges(3,6) * t258;
t159 = Icges(3,5) * t206 - Icges(3,6) * t205 + Icges(3,3) * t262;
t158 = Icges(3,5) * t204 - Icges(3,6) * t203 - Icges(3,3) * t258;
t146 = t165 + t251;
t145 = -t164 + t240;
t130 = -t164 * t224 - t194 * t258;
t129 = t165 * t224 - t194 * t262;
t127 = t205 * t135;
t114 = t242 * t224;
t107 = (t164 * t228 + t165 * t230) * t222;
t106 = t191 * t262 - t192 * t205 + t193 * t206;
t105 = -t191 * t258 - t192 * t203 + t193 * t204;
t99 = t235 + t126;
t98 = -t125 + t233;
t89 = -t126 * t265 - t157 * t205;
t88 = t125 * t265 + t157 * t203;
t85 = t159 * t224 + (t161 * t229 + t163 * t227) * t222;
t84 = t158 * t224 + (t160 * t229 + t162 * t227) * t222;
t74 = t75 * t224;
t73 = t125 * t205 - t126 * t203;
t72 = (-t125 - t171) * t224 + t230 * t239;
t71 = t126 * t224 + t228 * t239 + t168;
t68 = t235 - t269;
t67 = -t135 + t233 - t96;
t66 = (t125 * t228 + t126 * t230) * t222 + t252;
t63 = t235 + t278;
t62 = -t249 - t185 * t216 + (-rSges(6,3) + t225) * t184 + t233 + t236;
t61 = -t112 * t186 + t201 * t87;
t60 = t112 * t184 - t201 * t86;
t59 = t120 * t205 - t122 * t186 + t124 * t187;
t58 = t119 * t205 - t121 * t186 + t123 * t187;
t57 = t120 * t203 - t122 * t184 + t124 * t185;
t56 = t119 * t203 - t121 * t184 + t123 * t185;
t54 = t55 * t224;
t53 = -t184 * t87 + t186 * t86;
t51 = t52 * t224;
t50 = t205 * t256 + t265 * t269;
t49 = t118 * t203 + t265 * t96 + t255;
t48 = t52 * t201;
t45 = (-t96 + t253) * t224 + t230 * t237;
t44 = t224 * t97 + t228 * t237 + t254;
t41 = t203 * t269 + t205 * t96 + t127;
t38 = (t228 * t96 + t230 * t97) * t222 + t238;
t35 = t152 * t93 + t153 * t95 + t186 * t91;
t34 = t152 * t92 + t153 * t94 + t186 * t90;
t33 = t150 * t93 + t151 * t95 + t184 * t91;
t32 = t150 * t92 + t151 * t94 + t184 * t90;
t29 = t141 * t81 + t142 * t83 + t184 * t79;
t28 = t141 * t80 + t142 * t82 + t184 * t78;
t27 = t205 * t244 + t248 * t265;
t26 = t203 * t257 + t265 * t271 + t255;
t25 = (t253 - t271) * t224 + t230 * t234;
t24 = t224 * t270 + t228 * t234 + t254;
t23 = t74 + (t65 * t228 - t64 * t230) * t222;
t22 = t64 * t203 + t65 * t205 - t265 * t75;
t21 = t203 * t248 + t205 * t271 + t127;
t20 = (t228 * t271 + t230 * t270) * t222 + t238;
t19 = t70 * t224 + (t228 * t59 - t230 * t58) * t222;
t18 = t69 * t224 + (t228 * t57 - t230 * t56) * t222;
t17 = t203 * t58 + t205 * t59 - t265 * t70;
t16 = t203 * t56 + t205 * t57 - t265 * t69;
t15 = t54 + (t40 * t228 - t39 * t230) * t222;
t14 = t39 * t203 + t40 * t205 - t265 * t55;
t13 = t51 + (t37 * t228 - t36 * t230) * t222;
t12 = t36 * t203 + t37 * t205 - t265 * t52;
t11 = t36 * t184 + t37 * t186 + t48;
t10 = t47 * t224 + (t228 * t35 - t230 * t34) * t222;
t9 = t46 * t224 + (t228 * t33 - t230 * t32) * t222;
t8 = t203 * t34 + t205 * t35 - t265 * t47;
t7 = t203 * t32 + t205 * t33 - t265 * t46;
t6 = t43 * t224 + (t228 * t31 - t230 * t30) * t222;
t5 = t42 * t224 + (t228 * t29 - t230 * t28) * t222;
t4 = t203 * t30 + t205 * t31 - t265 * t43;
t3 = t203 * t28 + t205 * t29 - t265 * t42;
t1 = t184 * t28 + t186 * t29 + t201 * t42;
t76 = [Icges(2,3) + m(6) * (t62 ^ 2 + t63 ^ 2) + m(5) * (t67 ^ 2 + t68 ^ 2) + m(4) * (t98 ^ 2 + t99 ^ 2) + m(3) * (t145 ^ 2 + t146 ^ 2) + m(2) * (t208 ^ 2 + t209 ^ 2) + t242 - t250; t51 + t74 + t54 + t114 + m(6) * (t24 * t63 + t25 * t62) + m(5) * (t44 * t68 + t45 * t67) + m(4) * (t71 * t99 + t72 * t98) + m(3) * (t129 * t146 + t130 * t145) + ((-t105 / 0.2e1 - t84 / 0.2e1 - t232) * t230 + (t85 / 0.2e1 + t106 / 0.2e1 + t231) * t228) * t222; (t13 + t23 + t15 + t114) * t224 + m(6) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(5) * (t38 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t66 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(3) * (t107 ^ 2 + t129 ^ 2 + t130 ^ 2) + (-t230 * t5 + t228 * t6 + t228 * t10 + t228 * t19 - t230 * t18 - t230 * t9 + (-t230 * ((-t161 * t203 + t163 * t204) * t228 - (-t160 * t203 + t162 * t204) * t230) + t228 * ((-t161 * t205 + t163 * t206) * t228 - (-t160 * t205 + t162 * t206) * t230) + (-t230 * (t158 * t230 ^ 2 - t159 * t260) + t228 * (t159 * t228 ^ 2 - t158 * t260)) * t222) * t222 + ((-t105 - t84) * t230 + (t106 + t85) * t228) * t224) * t222; t250 * t265 + m(6) * (t26 * t62 + t27 * t63) + m(5) * (t49 * t67 + t50 * t68) + m(4) * (t88 * t98 + t89 * t99) + t231 * t205 + t232 * t203; (t12 / 0.2e1 + t22 / 0.2e1 + t14 / 0.2e1) * t224 + (t6 / 0.2e1 + t10 / 0.2e1 + t19 / 0.2e1) * t205 + (t5 / 0.2e1 + t18 / 0.2e1 + t9 / 0.2e1) * t203 + m(6) * (t20 * t21 + t24 * t27 + t25 * t26) + m(5) * (t38 * t41 + t44 * t50 + t45 * t49) + m(4) * (t66 * t73 + t71 * t89 + t72 * t88) + ((-t3 / 0.2e1 - t16 / 0.2e1 - t7 / 0.2e1) * t230 + (-t13 / 0.2e1 - t23 / 0.2e1 - t15 / 0.2e1) * t229 + (t4 / 0.2e1 + t17 / 0.2e1 + t8 / 0.2e1) * t228) * t222; (-t12 - t14 - t22) * t265 + (t4 + t17 + t8) * t205 + (t3 + t7 + t16) * t203 + m(6) * (t21 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t41 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(4) * (t73 ^ 2 + t88 ^ 2 + t89 ^ 2); m(6) * (t184 * t63 + t186 * t62) + m(5) * (t184 * t68 + t186 * t67); m(6) * (t184 * t24 + t186 * t25 + t20 * t201) + m(5) * (t184 * t44 + t186 * t45 + t201 * t38); m(6) * (t184 * t27 + t186 * t26 + t201 * t21) + m(5) * (t184 * t50 + t186 * t49 + t201 * t41); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (t184 ^ 2 + t186 ^ 2 + t201 ^ 2); t48 + m(6) * (t60 * t62 + t61 * t63) + t246 * t186 + t247 * t184; t13 * t274 + t224 * t11 / 0.2e1 + t5 * t276 + t6 * t275 + m(6) * (t20 * t53 + t24 * t61 + t25 * t60) + (t228 * t277 - t230 * t1 / 0.2e1) * t222; t12 * t274 + t205 * t277 + t203 * t1 / 0.2e1 + t3 * t276 + t4 * t275 + m(6) * (t21 * t53 + t26 * t60 + t27 * t61) - t11 * t265 / 0.2e1; m(6) * (t184 * t61 + t186 * t60 + t201 * t53); m(6) * (t53 ^ 2 + t60 ^ 2 + t61 ^ 2) + t186 * t2 + t184 * t1 + t201 * t11;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t76(1), t76(2), t76(4), t76(7), t76(11); t76(2), t76(3), t76(5), t76(8), t76(12); t76(4), t76(5), t76(6), t76(9), t76(13); t76(7), t76(8), t76(9), t76(10), t76(14); t76(11), t76(12), t76(13), t76(14), t76(15);];
Mq = res;
