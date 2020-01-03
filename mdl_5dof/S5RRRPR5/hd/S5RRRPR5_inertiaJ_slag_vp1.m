% Calculate joint inertia matrix for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:02
% EndTime: 2019-12-31 21:13:08
% DurationCPUTime: 2.34s
% Computational Cost: add. (5092->306), mult. (4837->451), div. (0->0), fcn. (5001->10), ass. (0->163)
t239 = Icges(4,3) + Icges(5,3);
t154 = qJ(2) + qJ(3);
t143 = pkin(9) + t154;
t140 = sin(t143);
t141 = cos(t143);
t144 = sin(t154);
t145 = cos(t154);
t238 = Icges(4,5) * t145 + Icges(5,5) * t141 - Icges(4,6) * t144 - Icges(5,6) * t140;
t157 = sin(qJ(1));
t160 = cos(qJ(1));
t237 = t239 * t157 + t238 * t160;
t236 = -t238 * t157 + t239 * t160;
t210 = Icges(4,4) * t144;
t179 = Icges(4,1) * t145 - t210;
t100 = Icges(4,5) * t157 + t179 * t160;
t207 = Icges(5,4) * t141;
t175 = -Icges(5,2) * t140 + t207;
t86 = Icges(5,6) * t157 + t175 * t160;
t208 = Icges(5,4) * t140;
t178 = Icges(5,1) * t141 - t208;
t88 = Icges(5,5) * t157 + t178 * t160;
t209 = Icges(4,4) * t145;
t176 = -Icges(4,2) * t144 + t209;
t98 = Icges(4,6) * t157 + t176 * t160;
t235 = -t100 * t145 + t140 * t86 - t141 * t88 + t144 * t98;
t85 = -Icges(5,6) * t160 + t175 * t157;
t87 = -Icges(5,5) * t160 + t178 * t157;
t97 = -Icges(4,6) * t160 + t176 * t157;
t99 = -Icges(4,5) * t160 + t179 * t157;
t234 = t140 * t85 - t141 * t87 + t144 * t97 - t145 * t99;
t152 = t157 ^ 2;
t233 = t157 * pkin(6);
t204 = t141 * t160;
t205 = t140 * t160;
t158 = cos(qJ(5));
t201 = t157 * t158;
t155 = sin(qJ(5));
t202 = t155 * t160;
t111 = -t141 * t202 + t201;
t200 = t158 * t160;
t203 = t155 * t157;
t112 = t141 * t200 + t203;
t60 = t112 * rSges(6,1) + t111 * rSges(6,2) + rSges(6,3) * t205;
t232 = pkin(4) * t204 + pkin(8) * t205 + t60;
t231 = Icges(4,5) * t144 + Icges(5,5) * t140 + Icges(4,6) * t145 + Icges(5,6) * t141;
t114 = Icges(5,2) * t141 + t208;
t115 = Icges(5,1) * t140 + t207;
t121 = Icges(4,2) * t145 + t210;
t122 = Icges(4,1) * t144 + t209;
t230 = -t114 * t140 + t115 * t141 - t121 * t144 + t122 * t145;
t187 = rSges(4,1) * t145 - rSges(4,2) * t144;
t74 = -rSges(6,3) * t141 + (rSges(6,1) * t158 - rSges(6,2) * t155) * t140;
t229 = -pkin(4) * t140 + pkin(8) * t141 - t74;
t153 = t160 ^ 2;
t109 = -t141 * t203 - t200;
t110 = t141 * t201 - t202;
t206 = t140 * t157;
t53 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t206;
t55 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t206;
t57 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t206;
t15 = t109 * t55 + t110 * t57 + t53 * t206;
t54 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t205;
t56 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t205;
t58 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t205;
t16 = t109 * t56 + t110 * t58 + t54 * t206;
t8 = -t15 * t160 + t157 * t16;
t228 = -t8 + t236 * t153 + (t235 * t157 + (-t234 + t237) * t160) * t157;
t161 = -pkin(7) - pkin(6);
t227 = t157 / 0.2e1;
t226 = -t160 / 0.2e1;
t156 = sin(qJ(2));
t225 = pkin(2) * t156;
t224 = pkin(3) * t144;
t223 = pkin(4) * t141;
t159 = cos(qJ(2));
t142 = t159 * pkin(2) + pkin(1);
t124 = pkin(3) * t145 + t142;
t118 = t160 * t124;
t136 = t160 * t142;
t222 = t160 * (t118 - t136) + (t124 - t142) * t152;
t150 = t160 * pkin(6);
t221 = t157 * (t150 + (-pkin(1) + t142) * t157) + t160 * (-t160 * pkin(1) + t136 - t233);
t167 = t157 * rSges(4,3) + t187 * t160;
t52 = t157 * (-t160 * rSges(4,3) + t187 * t157) + t160 * t167;
t220 = rSges(3,1) * t159;
t218 = rSges(3,2) * t156;
t72 = -Icges(6,6) * t141 + (Icges(6,4) * t158 - Icges(6,2) * t155) * t140;
t216 = t155 * t72;
t215 = t160 * rSges(3,3);
t24 = -t141 * t53 + (-t155 * t55 + t158 * t57) * t140;
t214 = t24 * t160;
t25 = -t141 * t54 + (-t155 * t56 + t158 * t58) * t140;
t213 = t25 * t157;
t212 = Icges(3,4) * t156;
t211 = Icges(3,4) * t159;
t198 = t157 * rSges(3,3) + t160 * t220;
t196 = t152 + t153;
t17 = t111 * t55 + t112 * t57 + t53 * t205;
t18 = t111 * t56 + t112 * t58 + t54 * t205;
t9 = t157 * t18 - t160 * t17;
t195 = (t9 + t237 * t152 + ((-t235 + t236) * t157 + t234 * t160) * t160) * t157;
t123 = rSges(4,1) * t144 + rSges(4,2) * t145;
t194 = -t123 - t225;
t116 = rSges(5,1) * t140 + rSges(5,2) * t141;
t193 = -t116 - t224;
t166 = rSges(5,1) * t204 - rSges(5,2) * t205 + t157 * rSges(5,3);
t186 = rSges(5,1) * t141 - rSges(5,2) * t140;
t28 = t157 * (-rSges(5,3) * t160 + t186 * t157) + t160 * t166 + t222;
t151 = -qJ(4) + t161;
t192 = -t157 * t151 + t118;
t71 = -Icges(6,3) * t141 + (Icges(6,5) * t158 - Icges(6,6) * t155) * t140;
t73 = -Icges(6,5) * t141 + (Icges(6,1) * t158 - Icges(6,4) * t155) * t140;
t29 = t109 * t72 + t110 * t73 + t71 * t206;
t3 = -t141 * t29 + (t15 * t157 + t16 * t160) * t140;
t30 = t111 * t72 + t112 * t73 + t71 * t205;
t4 = -t141 * t30 + (t157 * t17 + t160 * t18) * t140;
t191 = t8 * t206 / 0.2e1 + t3 * t226 + t4 * t227 - t141 * (t213 - t214) / 0.2e1 + t9 * t205 / 0.2e1;
t190 = -t224 + t229;
t189 = -t224 - t225;
t188 = -t218 + t220;
t185 = -rSges(6,1) * t110 - rSges(6,2) * t109;
t59 = rSges(6,3) * t206 - t185;
t12 = t157 * t59 + t222 + t152 * (pkin(8) * t140 + t223) + t232 * t160;
t180 = Icges(3,1) * t159 - t212;
t177 = -Icges(3,2) * t156 + t211;
t174 = Icges(3,5) * t159 - Icges(3,6) * t156;
t165 = -t116 + t189;
t164 = t189 + t229;
t163 = t228 * t160 + t195;
t162 = -t214 / 0.2e1 + t213 / 0.2e1 + (t100 * t144 + t140 * t88 + t141 * t86 + t145 * t98 + t231 * t157 + t230 * t160 + t30) * t227 + (t140 * t87 + t141 * t85 + t144 * t99 + t145 * t97 + t230 * t157 - t231 * t160 + t29) * t226;
t135 = rSges(2,1) * t160 - rSges(2,2) * t157;
t134 = -rSges(2,1) * t157 - rSges(2,2) * t160;
t133 = rSges(3,1) * t156 + rSges(3,2) * t159;
t104 = Icges(3,3) * t157 + t174 * t160;
t103 = -Icges(3,3) * t160 + t174 * t157;
t94 = t194 * t160;
t93 = t194 * t157;
t80 = t233 + (pkin(1) - t218) * t160 + t198;
t79 = t215 + t150 + (-pkin(1) - t188) * t157;
t78 = t193 * t160;
t77 = t193 * t157;
t70 = t165 * t160;
t69 = t165 * t157;
t68 = -t157 * t161 + t136 + t167;
t67 = (rSges(4,3) - t161) * t160 + (-t142 - t187) * t157;
t66 = t140 * t158 * t73;
t63 = t160 * (-t160 * t218 + t198) + (t188 * t157 - t215) * t157;
t62 = t166 + t192;
t61 = (rSges(5,3) - t151) * t160 + (-t124 - t186) * t157;
t49 = t190 * t160;
t48 = t190 * t157;
t43 = t164 * t160;
t42 = t164 * t157;
t37 = t192 + t232;
t36 = -t160 * t151 + (-t223 - t124 + (-rSges(6,3) - pkin(8)) * t140) * t157 + t185;
t35 = -t141 * t60 - t74 * t205;
t34 = t141 * t59 + t74 * t206;
t33 = t52 + t221;
t32 = -t140 * t216 - t141 * t71 + t66;
t31 = (-t157 * t60 + t160 * t59) * t140;
t19 = t28 + t221;
t11 = t12 + t221;
t1 = [t145 * t121 + t144 * t122 + t159 * (Icges(3,2) * t159 + t212) + t156 * (Icges(3,1) * t156 + t211) + Icges(2,3) + t66 + (-t71 + t114) * t141 + (t115 - t216) * t140 + m(6) * (t36 ^ 2 + t37 ^ 2) + m(5) * (t61 ^ 2 + t62 ^ 2) + m(4) * (t67 ^ 2 + t68 ^ 2) + m(3) * (t79 ^ 2 + t80 ^ 2) + m(2) * (t134 ^ 2 + t135 ^ 2); (t153 / 0.2e1 + t152 / 0.2e1) * (Icges(3,5) * t156 + Icges(3,6) * t159) + t162 + ((-Icges(3,6) * t160 + t177 * t157) * t159 + (-Icges(3,5) * t160 + t180 * t157) * t156) * t226 + ((Icges(3,6) * t157 + t177 * t160) * t159 + (Icges(3,5) * t157 + t180 * t160) * t156) * t227 + m(3) * (-t157 * t80 - t160 * t79) * t133 + m(6) * (t36 * t43 + t37 * t42) + m(5) * (t61 * t70 + t62 * t69) + m(4) * (t67 * t94 + t68 * t93); m(6) * (t11 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(4) * (t33 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(5) * (t19 ^ 2 + t69 ^ 2 + t70 ^ 2) + t157 * t152 * t104 + m(3) * (t196 * t133 ^ 2 + t63 ^ 2) + t195 + (-t153 * t103 + (-t157 * t103 + t160 * t104) * t157 + t228) * t160; t162 + m(6) * (t36 * t49 + t37 * t48) + m(5) * (t61 * t78 + t62 * t77) + m(4) * (-t157 * t68 - t160 * t67) * t123; m(6) * (t11 * t12 + t42 * t48 + t43 * t49) + m(4) * (t33 * t52 + (-t157 * t93 - t160 * t94) * t123) + m(5) * (t19 * t28 + t69 * t77 + t70 * t78) + t163; m(6) * (t12 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t28 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(4) * (t196 * t123 ^ 2 + t52 ^ 2) + t163; m(6) * (t157 * t36 - t160 * t37) + m(5) * (t157 * t61 - t160 * t62); m(6) * (t157 * t43 - t160 * t42) + m(5) * (t157 * t70 - t160 * t69); m(6) * (t157 * t49 - t160 * t48) + m(5) * (t157 * t78 - t160 * t77); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t196; m(6) * (t34 * t36 + t35 * t37) - t32 * t141 + ((t30 / 0.2e1 + t25 / 0.2e1) * t160 + (t24 / 0.2e1 + t29 / 0.2e1) * t157) * t140; m(6) * (t11 * t31 + t34 * t43 + t35 * t42) + t191; m(6) * (t12 * t31 + t34 * t49 + t35 * t48) + t191; m(6) * (t157 * t34 - t160 * t35); m(6) * (t31 ^ 2 + t34 ^ 2 + t35 ^ 2) + t141 ^ 2 * t32 + (t160 * t4 + t157 * t3 - t141 * (t157 * t24 + t160 * t25)) * t140;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
