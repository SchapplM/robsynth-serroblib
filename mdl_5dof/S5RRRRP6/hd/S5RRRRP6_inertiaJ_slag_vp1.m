% Calculate joint inertia matrix for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:50
% EndTime: 2019-12-31 21:52:59
% DurationCPUTime: 3.26s
% Computational Cost: add. (6040->349), mult. (7700->501), div. (0->0), fcn. (8262->8), ass. (0->180)
t177 = qJ(2) + qJ(3);
t171 = cos(t177);
t182 = cos(qJ(4));
t184 = cos(qJ(1));
t223 = t184 * t182;
t179 = sin(qJ(4));
t181 = sin(qJ(1));
t226 = t181 * t179;
t141 = -t171 * t226 - t223;
t224 = t184 * t179;
t225 = t181 * t182;
t142 = t171 * t225 - t224;
t170 = sin(t177);
t230 = t170 * t181;
t83 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t230;
t85 = Icges(5,5) * t142 + Icges(5,6) * t141 + Icges(5,3) * t230;
t276 = t83 + t85;
t143 = -t171 * t224 + t225;
t144 = t171 * t223 + t226;
t228 = t170 * t184;
t84 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t228;
t86 = Icges(5,5) * t144 + Icges(5,6) * t143 + Icges(5,3) * t228;
t275 = t84 + t86;
t87 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t230;
t89 = Icges(5,4) * t142 + Icges(5,2) * t141 + Icges(5,6) * t230;
t274 = t87 + t89;
t88 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t228;
t90 = Icges(5,4) * t144 + Icges(5,2) * t143 + Icges(5,6) * t228;
t273 = t88 + t90;
t91 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t230;
t93 = Icges(5,1) * t142 + Icges(5,4) * t141 + Icges(5,5) * t230;
t272 = t91 + t93;
t92 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t228;
t94 = Icges(5,1) * t144 + Icges(5,4) * t143 + Icges(5,5) * t228;
t271 = t92 + t94;
t178 = -qJ(5) - pkin(8);
t270 = rSges(6,3) - t178;
t269 = t274 * t141 + t272 * t142 + t276 * t230;
t268 = t273 * t141 + t271 * t142 + t275 * t230;
t267 = t274 * t143 + t272 * t144 + t276 * t228;
t266 = t273 * t143 + t271 * t144 + t275 * t228;
t113 = -Icges(6,3) * t171 + (Icges(6,5) * t182 - Icges(6,6) * t179) * t170;
t115 = -Icges(6,6) * t171 + (Icges(6,4) * t182 - Icges(6,2) * t179) * t170;
t117 = -Icges(6,5) * t171 + (Icges(6,1) * t182 - Icges(6,4) * t179) * t170;
t50 = t113 * t230 + t141 * t115 + t142 * t117;
t114 = -Icges(5,3) * t171 + (Icges(5,5) * t182 - Icges(5,6) * t179) * t170;
t116 = -Icges(5,6) * t171 + (Icges(5,4) * t182 - Icges(5,2) * t179) * t170;
t118 = -Icges(5,5) * t171 + (Icges(5,1) * t182 - Icges(5,4) * t179) * t170;
t51 = t114 * t230 + t141 * t116 + t142 * t118;
t265 = -t50 - t51;
t52 = t113 * t228 + t143 * t115 + t144 * t117;
t53 = t114 * t228 + t143 * t116 + t144 * t118;
t264 = -t52 - t53;
t167 = t182 * pkin(4) + pkin(3);
t227 = t171 * t184;
t263 = t144 * rSges(6,1) + t143 * rSges(6,2) + pkin(4) * t226 + t167 * t227 + t270 * t228;
t262 = t265 * t171 + (t269 * t181 + t268 * t184) * t170;
t261 = t264 * t171 + (t267 * t181 + t266 * t184) * t170;
t260 = t268 * t181 - t269 * t184;
t259 = t266 * t181 - t267 * t184;
t203 = -t142 * rSges(6,1) - t141 * rSges(6,2);
t246 = pkin(8) + t178;
t247 = -pkin(3) + t167;
t244 = -pkin(4) * t224 + (-t170 * t246 + t171 * t247) * t181 + rSges(6,3) * t230 - t203;
t217 = pkin(3) * t227 + pkin(8) * t228;
t258 = -t217 + t263;
t257 = -t113 - t114;
t256 = (-t115 - t116) * t179;
t255 = (t117 + t118) * t170 * t182;
t197 = Icges(4,5) * t171 - Icges(4,6) * t170;
t123 = -Icges(4,3) * t184 + t181 * t197;
t124 = Icges(4,3) * t181 + t184 * t197;
t176 = t184 ^ 2;
t232 = Icges(4,4) * t171;
t199 = -Icges(4,2) * t170 + t232;
t126 = Icges(4,6) * t181 + t184 * t199;
t233 = Icges(4,4) * t170;
t201 = Icges(4,1) * t171 - t233;
t128 = Icges(4,5) * t181 + t184 * t201;
t195 = -t126 * t170 + t128 * t171;
t125 = -Icges(4,6) * t184 + t181 * t199;
t127 = -Icges(4,5) * t184 + t181 * t201;
t196 = t125 * t170 - t127 * t171;
t254 = -t176 * t123 - (t195 * t181 + (-t124 + t196) * t184) * t181 - t260;
t253 = t171 ^ 2;
t175 = t181 ^ 2;
t251 = t181 / 0.2e1;
t250 = -t184 / 0.2e1;
t180 = sin(qJ(2));
t249 = pkin(2) * t180;
t248 = pkin(3) * t171;
t245 = t170 * t256 + t257 * t171 + t255;
t183 = cos(qJ(2));
t242 = rSges(3,1) * t183;
t241 = rSges(3,2) * t180;
t240 = t184 * rSges(3,3);
t41 = -t171 * t83 + (-t179 * t87 + t182 * t91) * t170;
t239 = t41 * t184;
t42 = -t171 * t84 + (-t179 * t88 + t182 * t92) * t170;
t238 = t42 * t181;
t43 = -t171 * t85 + (-t179 * t89 + t182 * t93) * t170;
t237 = t43 * t184;
t44 = -t171 * t86 + (-t179 * t90 + t182 * t94) * t170;
t236 = t44 * t181;
t235 = Icges(3,4) * t180;
t234 = Icges(3,4) * t183;
t185 = -pkin(7) - pkin(6);
t222 = t184 * t185;
t221 = (t246 - rSges(6,3)) * t171 + (rSges(6,1) * t182 - rSges(6,2) * t179 + t247) * t170;
t168 = t183 * pkin(2) + pkin(1);
t159 = t184 * t168;
t174 = t184 * pkin(6);
t220 = t181 * (t222 + t174 + (-pkin(1) + t168) * t181) + t184 * (-t184 * pkin(1) + t159 + (-pkin(6) - t185) * t181);
t190 = rSges(4,1) * t227 - rSges(4,2) * t228 + t181 * rSges(4,3);
t205 = rSges(4,1) * t171 - rSges(4,2) * t170;
t72 = t181 * (-t184 * rSges(4,3) + t181 * t205) + t184 * t190;
t120 = -t171 * rSges(5,3) + (rSges(5,1) * t182 - rSges(5,2) * t179) * t170;
t151 = t170 * pkin(3) - t171 * pkin(8);
t219 = -t120 - t151;
t218 = t175 * (pkin(8) * t170 + t248) + t184 * t217;
t216 = t181 * rSges(3,3) + t184 * t242;
t215 = t175 + t176;
t214 = (t175 * t124 + (t196 * t184 + (-t123 + t195) * t181) * t184 + t259) * t181;
t213 = -t151 - t221;
t98 = t144 * rSges(5,1) + t143 * rSges(5,2) + rSges(5,3) * t228;
t150 = t170 * rSges(4,1) + t171 * rSges(4,2);
t210 = -t150 - t249;
t209 = -t151 - t249;
t208 = -t181 * t185 + t159;
t204 = -t142 * rSges(5,1) - t141 * rSges(5,2);
t96 = rSges(5,3) * t230 - t204;
t45 = t181 * t96 + t184 * t98 + t218;
t207 = -t120 + t209;
t206 = -t241 + t242;
t202 = Icges(3,1) * t183 - t235;
t200 = -Icges(3,2) * t180 + t234;
t198 = Icges(3,5) * t183 - Icges(3,6) * t180;
t147 = Icges(4,2) * t171 + t233;
t148 = Icges(4,1) * t170 + t232;
t192 = -t147 * t170 + t148 * t171;
t191 = t209 - t221;
t22 = t244 * t181 + t258 * t184 + t218;
t188 = t254 * t184 + t214;
t187 = -(t238 - t239 + t236 - t237) * t171 / 0.2e1 + t261 * t251 + t262 * t250 + t260 * t230 / 0.2e1 + t259 * t228 / 0.2e1;
t146 = Icges(4,5) * t170 + Icges(4,6) * t171;
t186 = -t239 / 0.2e1 + t236 / 0.2e1 + t238 / 0.2e1 - t237 / 0.2e1 + (t171 * t126 + t170 * t128 + t181 * t146 + t184 * t192 - t264) * t251 + (t171 * t125 + t170 * t127 - t184 * t146 + t181 * t192 - t265) * t250;
t158 = t184 * rSges(2,1) - t181 * rSges(2,2);
t157 = -t181 * rSges(2,1) - t184 * rSges(2,2);
t156 = t180 * rSges(3,1) + t183 * rSges(3,2);
t132 = Icges(3,3) * t181 + t184 * t198;
t131 = -Icges(3,3) * t184 + t181 * t198;
t122 = t210 * t184;
t121 = t210 * t181;
t108 = t181 * pkin(6) + (pkin(1) - t241) * t184 + t216;
t107 = t240 + t174 + (-pkin(1) - t206) * t181;
t103 = t190 + t208;
t102 = (rSges(4,3) - t185) * t184 + (-t168 - t205) * t181;
t101 = t219 * t184;
t100 = t219 * t181;
t99 = t184 * (-t184 * t241 + t216) + (t181 * t206 - t240) * t181;
t76 = t207 * t184;
t75 = t207 * t181;
t67 = t208 + t98 + t217;
t66 = -t222 + (-t248 - t168 + (-rSges(5,3) - pkin(8)) * t170) * t181 + t204;
t65 = t213 * t184;
t64 = t213 * t181;
t63 = -t120 * t228 - t171 * t98;
t62 = t120 * t230 + t171 * t96;
t61 = t191 * t184;
t60 = t191 * t181;
t59 = t208 + t263;
t58 = (pkin(4) * t179 - t185) * t184 + (-t167 * t171 - t270 * t170 - t168) * t181 + t203;
t55 = t72 + t220;
t54 = (-t181 * t98 + t184 * t96) * t170;
t36 = -t171 * t258 - t221 * t228;
t35 = t171 * t244 + t221 * t230;
t24 = t45 + t220;
t23 = (-t181 * t258 + t184 * t244) * t170;
t21 = t22 + t220;
t1 = [t183 * (Icges(3,2) * t183 + t235) + t180 * (Icges(3,1) * t180 + t234) + Icges(2,3) + (t147 + t257) * t171 + (t148 + t256) * t170 + m(6) * (t58 ^ 2 + t59 ^ 2) + m(5) * (t66 ^ 2 + t67 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2) + m(3) * (t107 ^ 2 + t108 ^ 2) + m(2) * (t157 ^ 2 + t158 ^ 2) + t255; (t175 / 0.2e1 + t176 / 0.2e1) * (Icges(3,5) * t180 + Icges(3,6) * t183) + t186 + m(6) * (t61 * t58 + t60 * t59) + m(5) * (t76 * t66 + t75 * t67) + m(4) * (t122 * t102 + t121 * t103) + m(3) * (-t107 * t184 - t108 * t181) * t156 + (t183 * (Icges(3,6) * t181 + t184 * t200) + t180 * (Icges(3,5) * t181 + t184 * t202)) * t251 + (t183 * (-Icges(3,6) * t184 + t181 * t200) + t180 * (-Icges(3,5) * t184 + t181 * t202)) * t250; m(6) * (t21 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(5) * (t24 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(4) * (t121 ^ 2 + t122 ^ 2 + t55 ^ 2) + t181 * t175 * t132 + m(3) * (t156 ^ 2 * t215 + t99 ^ 2) + t214 + (-t176 * t131 + (-t181 * t131 + t184 * t132) * t181 + t254) * t184; t186 + m(6) * (t65 * t58 + t64 * t59) + m(5) * (t100 * t67 + t101 * t66) + m(4) * (-t102 * t184 - t103 * t181) * t150; m(6) * (t22 * t21 + t64 * t60 + t65 * t61) + m(5) * (t100 * t75 + t101 * t76 + t45 * t24) + m(4) * (t72 * t55 + (-t121 * t181 - t122 * t184) * t150) + t188; m(6) * (t22 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2 + t45 ^ 2) + m(4) * (t150 ^ 2 * t215 + t72 ^ 2) + t188; -t245 * t171 + m(6) * (t35 * t58 + t36 * t59) + m(5) * (t62 * t66 + t63 * t67) + ((t53 / 0.2e1 + t52 / 0.2e1 + t44 / 0.2e1 + t42 / 0.2e1) * t184 + (t41 / 0.2e1 + t51 / 0.2e1 + t50 / 0.2e1 + t43 / 0.2e1) * t181) * t170; m(6) * (t23 * t21 + t35 * t61 + t36 * t60) + m(5) * (t54 * t24 + t62 * t76 + t63 * t75) + t187; m(6) * (t23 * t22 + t35 * t65 + t36 * t64) + m(5) * (t63 * t100 + t62 * t101 + t54 * t45) + t187; m(6) * (t23 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t54 ^ 2 + t62 ^ 2 + t63 ^ 2) + t245 * t253 + (t261 * t184 + t262 * t181 + ((-t42 - t44) * t184 + (-t41 - t43) * t181) * t171) * t170; m(6) * (t181 * t59 + t184 * t58) * t170; m(6) * (-t171 * t21 + (t181 * t60 + t184 * t61) * t170); m(6) * (-t171 * t22 + (t181 * t64 + t184 * t65) * t170); m(6) * (-t171 * t23 + (t181 * t36 + t184 * t35) * t170); m(6) * (t170 ^ 2 * t215 + t253);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
