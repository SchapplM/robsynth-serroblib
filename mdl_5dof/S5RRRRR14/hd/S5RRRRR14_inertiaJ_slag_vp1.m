% Calculate joint inertia matrix for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR14_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR14_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR14_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR14_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:42:21
% EndTime: 2024-09-27 18:42:22
% DurationCPUTime: 0.78s
% Computational Cost: add. (17196->352), mult. (10228->485), div. (0->0), fcn. (9474->22), ass. (0->186)
t199 = qJ(1) + qJ(2);
t191 = sin(t199);
t193 = cos(t199);
t203 = sin(qJ(3));
t201 = cos(pkin(5));
t205 = cos(qJ(3));
t228 = t201 * t205;
t141 = -t191 * t203 + t193 * t228;
t229 = t201 * t203;
t142 = t191 * t205 + t193 * t229;
t200 = sin(pkin(5));
t231 = t193 * t200;
t101 = Icges(4,1) * t142 + Icges(4,4) * t141 - Icges(4,5) * t231;
t137 = Icges(4,3) * t201 + (Icges(4,5) * t203 + Icges(4,6) * t205) * t200;
t138 = Icges(4,6) * t201 + (Icges(4,4) * t203 + Icges(4,2) * t205) * t200;
t139 = Icges(4,5) * t201 + (Icges(4,1) * t203 + Icges(4,4) * t205) * t200;
t97 = Icges(4,5) * t142 + Icges(4,6) * t141 - Icges(4,3) * t231;
t99 = Icges(4,4) * t142 + Icges(4,2) * t141 - Icges(4,6) * t231;
t241 = t201 * t97 + (t101 * t203 + t205 * t99) * t200 - t137 * t231 + t141 * t138 + t142 * t139;
t143 = -t191 * t228 - t193 * t203;
t144 = -t191 * t229 + t193 * t205;
t233 = t191 * t200;
t100 = Icges(4,4) * t144 + Icges(4,2) * t143 + Icges(4,6) * t233;
t102 = Icges(4,1) * t144 + Icges(4,4) * t143 + Icges(4,5) * t233;
t98 = Icges(4,5) * t144 + Icges(4,6) * t143 + Icges(4,3) * t233;
t240 = t201 * t98 + (t100 * t205 + t102 * t203) * t200 + t137 * t233 + t143 * t138 + t144 * t139;
t207 = pkin(8) + pkin(9);
t198 = qJ(3) + qJ(4);
t189 = pkin(5) - t198;
t183 = -qJ(5) + t189;
t173 = cos(t183) / 0.2e1;
t188 = pkin(5) + t198;
t182 = qJ(5) + t188;
t177 = cos(t182);
t151 = t173 + t177 / 0.2e1;
t194 = qJ(5) + t198;
t184 = sin(t194);
t123 = -t191 * t151 - t193 * t184;
t172 = sin(t182) / 0.2e1;
t176 = sin(t183);
t150 = t172 - t176 / 0.2e1;
t185 = cos(t194);
t124 = -t191 * t150 + t193 * t185;
t149 = t172 + t176 / 0.2e1;
t152 = t173 - t177 / 0.2e1;
t121 = t193 * t151 - t191 * t184;
t122 = t193 * t150 + t191 * t185;
t70 = Icges(6,5) * t122 + Icges(6,6) * t121 - Icges(6,3) * t231;
t72 = Icges(6,4) * t122 + Icges(6,2) * t121 - Icges(6,6) * t231;
t74 = Icges(6,1) * t122 + Icges(6,4) * t121 - Icges(6,5) * t231;
t13 = t149 * t72 + t152 * t74 + t201 * t70;
t71 = Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t233;
t73 = Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t233;
t75 = Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t233;
t14 = t149 * t73 + t152 * t75 + t201 * t71;
t109 = Icges(6,5) * t152 + Icges(6,6) * t149 + Icges(6,3) * t201;
t110 = Icges(6,4) * t152 + Icges(6,2) * t149 + Icges(6,6) * t201;
t111 = Icges(6,1) * t152 + Icges(6,4) * t149 + Icges(6,5) * t201;
t25 = t109 * t233 + t123 * t110 + t124 * t111;
t223 = t201 * t109 + t149 * t110 + t152 * t111;
t37 = t223 * t201;
t239 = ((t123 * t73 + t124 * t75 + t71 * t233) * t233 - (t123 * t72 + t124 * t74 + t70 * t233) * t231 + t25 * t201) * t233 + t201 * (t37 + (-t13 * t193 + t14 * t191) * t200);
t204 = sin(qJ(1));
t238 = t204 * pkin(1);
t187 = t205 * pkin(3) + pkin(2);
t76 = t122 * rSges(6,1) + t121 * rSges(6,2) - rSges(6,3) * t231;
t77 = t124 * rSges(6,1) + t123 * rSges(6,2) + rSges(6,3) * t233;
t34 = t77 * t231 + t76 * t233;
t64 = t201 * t77;
t192 = cos(t198);
t160 = pkin(4) * t192 + t187;
t148 = t193 * t160;
t161 = t193 * t187;
t186 = cos(qJ(4)) * pkin(4) + pkin(3);
t197 = pkin(10) + t207;
t224 = pkin(4) * sin(qJ(4)) * t205;
t131 = -t200 * t197 + (t186 * t203 + t224) * t201;
t157 = pkin(3) * t229 - t200 * t207;
t226 = t131 - t157;
t66 = -t226 * t191 + t148 - t161;
t237 = t201 * t66 + t64;
t65 = t226 * t193 + (t160 - t187) * t191;
t236 = -t65 - t76;
t175 = cos(t189) / 0.2e1;
t180 = cos(t188);
t155 = t175 + t180 / 0.2e1;
t190 = sin(t198);
t127 = t193 * t155 - t191 * t190;
t174 = sin(t188) / 0.2e1;
t179 = sin(t189);
t154 = t174 - t179 / 0.2e1;
t128 = t193 * t154 + t191 * t192;
t85 = t128 * rSges(5,1) + t127 * rSges(5,2) - rSges(5,3) * t231;
t129 = -t191 * t155 - t193 * t190;
t130 = -t191 * t154 + t193 * t192;
t86 = t130 * rSges(5,1) + t129 * rSges(5,2) + rSges(5,3) * t233;
t38 = t86 * t231 + t85 * t233;
t171 = pkin(8) * t231;
t232 = t193 * t157;
t107 = t232 + t171 + (-pkin(2) + t187) * t191;
t215 = -t191 * t157 + t161;
t225 = t193 * pkin(2) + pkin(8) * t233;
t108 = t215 - t225;
t235 = t107 * t233 + t108 * t231;
t234 = t191 * t193;
t230 = t200 * t203;
t113 = t152 * rSges(6,1) + t149 * rSges(6,2) + t201 * rSges(6,3);
t227 = -(t197 - t207) * t201 - (t224 + (-pkin(3) + t186) * t203) * t200 - t113;
t153 = t174 + t179 / 0.2e1;
t156 = t175 - t180 / 0.2e1;
t115 = Icges(5,5) * t156 + Icges(5,6) * t153 + Icges(5,3) * t201;
t116 = Icges(5,4) * t156 + Icges(5,2) * t153 + Icges(5,6) * t201;
t117 = Icges(5,1) * t156 + Icges(5,4) * t153 + Icges(5,5) * t201;
t222 = t201 * t115 + t153 * t116 + t156 * t117;
t79 = Icges(5,5) * t128 + Icges(5,6) * t127 - Icges(5,3) * t231;
t81 = Icges(5,4) * t128 + Icges(5,2) * t127 - Icges(5,6) * t231;
t83 = Icges(5,1) * t128 + Icges(5,4) * t127 - Icges(5,5) * t231;
t20 = t153 * t81 + t156 * t83 + t201 * t79;
t80 = Icges(5,5) * t130 + Icges(5,6) * t129 + Icges(5,3) * t233;
t82 = Icges(5,4) * t130 + Icges(5,2) * t129 + Icges(5,6) * t233;
t84 = Icges(5,1) * t130 + Icges(5,4) * t129 + Icges(5,5) * t233;
t21 = t153 * t82 + t156 * t84 + t201 * t80;
t33 = t115 * t233 + t129 * t116 + t130 * t117;
t39 = t222 * t201;
t221 = ((t129 * t82 + t130 * t84 + t80 * t233) * t233 - (t129 * t81 + t130 * t83 + t79 * t233) * t231 + t33 * t201) * t233 + t201 * (t39 + (t191 * t21 - t193 * t20) * t200) + t239;
t220 = t200 * t205 * t138 + t201 * t137 + t139 * t230;
t104 = t144 * rSges(4,1) + t143 * rSges(4,2) + rSges(4,3) * t233;
t219 = t233 / 0.2e1;
t218 = -t231 / 0.2e1;
t8 = t66 * t231 + t65 * t233 + t34;
t159 = t193 * rSges(3,1) - t191 * rSges(3,2);
t217 = t227 * t200;
t118 = t156 * rSges(5,1) + t153 * rSges(5,2) + t201 * rSges(5,3);
t145 = pkin(3) * t230 + (-pkin(8) + t207) * t201;
t216 = (-t118 - t145) * t200;
t24 = -t109 * t231 + t121 * t110 + t122 * t111;
t214 = t37 + (t14 + t25) * t219 + (t13 + t24) * t218;
t213 = (-t145 + t227) * t200;
t2 = (t121 * t73 + t122 * t75 - t71 * t231) * t233 - (t121 * t72 + t122 * t74 - t70 * t231) * t231 + t24 * t201;
t212 = -t2 * t231 + t239;
t158 = -t191 * rSges(3,1) - t193 * rSges(3,2);
t90 = t104 + t225;
t103 = t142 * rSges(4,1) + t141 * rSges(4,2) - rSges(4,3) * t231;
t45 = -t191 * t131 + t148 + t77;
t56 = t215 + t86;
t32 = -t115 * t231 + t127 * t116 + t128 * t117;
t4 = (t127 * t82 + t128 * t84 - t80 * t231) * t233 - (t127 * t81 + t128 * t83 - t79 * t231) * t231 + t32 * t201;
t211 = (-t2 - t4) * t231 + t221;
t210 = t214 + t39 + (t21 + t33) * t219 + (t20 + t32) * t218;
t89 = -t191 * pkin(2) - t103 + t171;
t209 = Icges(3,3) + t220 + t222 + t223;
t55 = -t191 * t187 - t232 - t85;
t44 = -t193 * t131 - t191 * t160 - t76;
t67 = t220 * t201;
t208 = t241 * t218 + t240 * t219 + t210 + t67;
t206 = cos(qJ(1));
t196 = t206 * pkin(1);
t163 = t206 * rSges(2,1) - t204 * rSges(2,2);
t162 = -t204 * rSges(2,1) - t206 * rSges(2,2);
t147 = t159 + t196;
t146 = t158 - t238;
t140 = t201 * rSges(4,3) + (rSges(4,1) * t203 + rSges(4,2) * t205) * t200;
t105 = t201 * t108;
t88 = t196 + t90;
t87 = t89 - t238;
t78 = t201 * t86;
t63 = -t201 * t103 - t140 * t231;
t62 = t201 * t104 - t140 * t233;
t54 = t196 + t56;
t53 = t55 - t238;
t48 = (t103 * t191 + t104 * t193) * t200;
t47 = -t118 * t231 - t201 * t85;
t46 = -t118 * t233 + t78;
t43 = -t113 * t231 - t201 * t76;
t42 = -t113 * t233 + t64;
t41 = t196 + t45;
t40 = t44 - t238;
t27 = (-t107 - t85) * t201 + t193 * t216;
t26 = t191 * t216 + t105 + t78;
t17 = t235 + t38;
t16 = t193 * t217 + t236 * t201;
t15 = t191 * t217 + t237;
t10 = (-t107 + t236) * t201 + t193 * t213;
t9 = t191 * t213 + t105 + t237;
t7 = t8 + t235;
t1 = [Icges(2,3) + m(6) * (t40 ^ 2 + t41 ^ 2) + m(5) * (t53 ^ 2 + t54 ^ 2) + m(4) * (t87 ^ 2 + t88 ^ 2) + m(3) * (t146 ^ 2 + t147 ^ 2) + m(2) * (t162 ^ 2 + t163 ^ 2) + t209; m(6) * (t40 * t44 + t41 * t45) + m(5) * (t53 * t55 + t54 * t56) + m(4) * (t87 * t89 + t88 * t90) + m(3) * (t146 * t158 + t147 * t159) + t209; m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t55 ^ 2 + t56 ^ 2) + m(4) * (t89 ^ 2 + t90 ^ 2) + m(3) * (t158 ^ 2 + t159 ^ 2) + t209; m(4) * (t62 * t88 + t63 * t87) + m(6) * (t10 * t40 + t41 * t9) + m(5) * (t26 * t54 + t27 * t53) + t208; m(6) * (t10 * t44 + t45 * t9) + m(5) * (t26 * t56 + t27 * t55) + m(4) * (t62 * t90 + t63 * t89) + t208; t201 * t67 + (-t193 * t2 - t193 * t4 + (t191 * ((t143 * t100 + t144 * t102) * t191 - (t144 * t101 + t143 * t99) * t193) - t193 * ((t141 * t100 + t142 * t102) * t191 - (t142 * t101 + t141 * t99) * t193) + (t191 * (t191 ^ 2 * t98 - t97 * t234) - t193 * (t193 ^ 2 * t97 - t98 * t234)) * t200) * t200 + (t240 * t191 - t241 * t193) * t201) * t200 + m(6) * (t10 ^ 2 + t7 ^ 2 + t9 ^ 2) + m(5) * (t17 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(4) * (t48 ^ 2 + t62 ^ 2 + t63 ^ 2) + t221; m(6) * (t15 * t41 + t16 * t40) + m(5) * (t46 * t54 + t47 * t53) + t210; m(6) * (t15 * t45 + t16 * t44) + m(5) * (t46 * t56 + t47 * t55) + t210; m(6) * (t10 * t16 + t15 * t9 + t7 * t8) + m(5) * (t17 * t38 + t26 * t46 + t27 * t47) + t211; m(5) * (t38 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t15 ^ 2 + t16 ^ 2 + t8 ^ 2) + t211; m(6) * (t40 * t43 + t41 * t42) + t214; m(6) * (t42 * t45 + t43 * t44) + t214; m(6) * (t10 * t43 + t34 * t7 + t42 * t9) + t212; m(6) * (t15 * t42 + t16 * t43 + t34 * t8) + t212; m(6) * (t34 ^ 2 + t42 ^ 2 + t43 ^ 2) + t212;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
