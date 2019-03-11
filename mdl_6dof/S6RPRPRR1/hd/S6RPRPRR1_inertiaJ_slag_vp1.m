% Calculate joint inertia matrix for
% S6RPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:33:53
% EndTime: 2019-03-09 03:33:56
% DurationCPUTime: 1.86s
% Computational Cost: add. (7049->315), mult. (4699->452), div. (0->0), fcn. (4845->12), ass. (0->157)
t228 = Icges(4,3) + Icges(5,3);
t148 = qJ(3) + pkin(11);
t140 = sin(t148);
t142 = cos(t148);
t152 = sin(qJ(3));
t155 = cos(qJ(3));
t227 = Icges(4,5) * t155 + Icges(5,5) * t142 - Icges(4,6) * t152 - Icges(5,6) * t140;
t149 = qJ(1) + pkin(10);
t141 = sin(t149);
t138 = t141 ^ 2;
t226 = pkin(7) * t141;
t144 = qJ(5) + t148;
t136 = cos(t144);
t143 = cos(t149);
t197 = t136 * t143;
t135 = sin(t144);
t198 = t135 * t143;
t151 = sin(qJ(6));
t194 = t143 * t151;
t154 = cos(qJ(6));
t195 = t141 * t154;
t100 = -t136 * t194 + t195;
t193 = t143 * t154;
t196 = t141 * t151;
t101 = t136 * t193 + t196;
t54 = t101 * rSges(7,1) + t100 * rSges(7,2) + rSges(7,3) * t198;
t225 = pkin(5) * t197 + pkin(9) * t198 + t54;
t180 = rSges(5,1) * t142 - rSges(5,2) * t140;
t139 = t143 ^ 2;
t200 = Icges(6,4) * t136;
t167 = -Icges(6,2) * t135 + t200;
t73 = Icges(6,6) * t141 + t143 * t167;
t201 = Icges(6,4) * t135;
t170 = Icges(6,1) * t136 - t201;
t75 = Icges(6,5) * t141 + t143 * t170;
t177 = -t135 * t73 + t136 * t75;
t72 = -Icges(6,6) * t143 + t141 * t167;
t74 = -Icges(6,5) * t143 + t141 * t170;
t178 = t135 * t72 - t136 * t74;
t164 = Icges(6,5) * t136 - Icges(6,6) * t135;
t70 = -Icges(6,3) * t143 + t141 * t164;
t71 = Icges(6,3) * t141 + t143 * t164;
t199 = t135 * t141;
t98 = -t136 * t196 - t193;
t99 = t136 * t195 - t194;
t46 = Icges(7,5) * t99 + Icges(7,6) * t98 + Icges(7,3) * t199;
t48 = Icges(7,4) * t99 + Icges(7,2) * t98 + Icges(7,6) * t199;
t50 = Icges(7,1) * t99 + Icges(7,4) * t98 + Icges(7,5) * t199;
t15 = t199 * t46 + t48 * t98 + t50 * t99;
t47 = Icges(7,5) * t101 + Icges(7,6) * t100 + Icges(7,3) * t198;
t49 = Icges(7,4) * t101 + Icges(7,2) * t100 + Icges(7,6) * t198;
t51 = Icges(7,1) * t101 + Icges(7,4) * t100 + Icges(7,5) * t198;
t16 = t199 * t47 + t49 * t98 + t51 * t99;
t8 = t141 * t16 - t143 * t15;
t224 = -t139 * t70 - (t177 * t141 + (t178 - t71) * t143) * t141 - t8;
t223 = -t227 * t141 + t228 * t143;
t222 = t228 * t141 + t227 * t143;
t221 = t141 / 0.2e1;
t220 = -t143 / 0.2e1;
t17 = t100 * t48 + t101 * t50 + t198 * t46;
t18 = t100 * t49 + t101 * t51 + t198 * t47;
t9 = t141 * t18 - t143 * t17;
t219 = (t138 * t71 + t9 + (t178 * t143 + (t177 - t70) * t141) * t143) * t141;
t153 = sin(qJ(1));
t218 = pkin(1) * t153;
t217 = pkin(3) * t152;
t216 = pkin(5) * t136;
t150 = -qJ(4) - pkin(7);
t137 = t155 * pkin(3) + pkin(2);
t122 = t143 * t137;
t134 = t143 * pkin(7);
t215 = t141 * (t134 + (-pkin(2) + t137) * t141) + t143 * (-pkin(2) * t143 + t122 - t226);
t160 = rSges(6,1) * t197 - rSges(6,2) * t198 + t141 * rSges(6,3);
t179 = rSges(6,1) * t136 - rSges(6,2) * t135;
t41 = t141 * (-rSges(6,3) * t143 + t141 * t179) + t143 * t160;
t214 = rSges(4,1) * t155;
t212 = rSges(4,2) * t152;
t210 = t143 * rSges(4,3);
t85 = -Icges(7,6) * t136 + (Icges(7,4) * t154 - Icges(7,2) * t151) * t135;
t209 = t151 * t85;
t21 = -t136 * t46 + (-t151 * t48 + t154 * t50) * t135;
t208 = t21 * t143;
t22 = -t136 * t47 + (-t151 * t49 + t154 * t51) * t135;
t207 = t22 * t141;
t87 = -rSges(7,3) * t136 + (rSges(7,1) * t154 - rSges(7,2) * t151) * t135;
t206 = -pkin(5) * t135 + pkin(9) * t136 - t87;
t205 = Icges(4,4) * t152;
t204 = Icges(4,4) * t155;
t203 = Icges(5,4) * t140;
t202 = Icges(5,4) * t142;
t191 = t141 * rSges(4,3) + t143 * t214;
t190 = t138 + t139;
t188 = -rSges(5,1) * t140 - rSges(5,2) * t142 - t217;
t119 = pkin(4) * t142 + t137;
t104 = t143 * t119;
t187 = t143 * (t104 - t122) + t215 + (t119 - t137) * t138;
t184 = -rSges(7,1) * t99 - rSges(7,2) * t98;
t53 = rSges(7,3) * t199 - t184;
t23 = t141 * t53 + t138 * (pkin(9) * t135 + t216) + t225 * t143;
t84 = -Icges(7,3) * t136 + (Icges(7,5) * t154 - Icges(7,6) * t151) * t135;
t86 = -Icges(7,5) * t136 + (Icges(7,1) * t154 - Icges(7,4) * t151) * t135;
t28 = t199 * t84 + t85 * t98 + t86 * t99;
t3 = -t136 * t28 + (t141 * t15 + t143 * t16) * t135;
t29 = t100 * t85 + t101 * t86 + t198 * t84;
t4 = -t136 * t29 + (t141 * t17 + t143 * t18) * t135;
t186 = t3 * t220 + t4 * t221 - t136 * (t207 - t208) / 0.2e1 + t8 * t199 / 0.2e1 + t9 * t198 / 0.2e1;
t183 = -pkin(4) * t140 - t217;
t156 = cos(qJ(1));
t146 = t156 * pkin(1);
t147 = -pkin(8) + t150;
t182 = -t141 * t147 + t104 + t146;
t181 = -t212 + t214;
t172 = Icges(4,1) * t155 - t205;
t171 = Icges(5,1) * t142 - t203;
t169 = -Icges(4,2) * t152 + t204;
t168 = -Icges(5,2) * t140 + t202;
t107 = Icges(6,2) * t136 + t201;
t108 = Icges(6,1) * t135 + t200;
t163 = -t107 * t135 + t108 * t136;
t162 = t143 * t224 + t219;
t161 = t141 * rSges(5,3) + t180 * t143;
t109 = rSges(6,1) * t135 + rSges(6,2) * t136;
t159 = -t109 + t183;
t158 = t183 + t206;
t106 = Icges(6,5) * t135 + Icges(6,6) * t136;
t157 = -t208 / 0.2e1 + t207 / 0.2e1 + (t106 * t141 + t135 * t75 + t136 * t73 + t143 * t163 + t29) * t221 + (-t106 * t143 + t135 * t74 + t136 * t72 + t141 * t163 + t28) * t220;
t129 = rSges(2,1) * t156 - t153 * rSges(2,2);
t128 = -t153 * rSges(2,1) - rSges(2,2) * t156;
t127 = rSges(4,1) * t152 + rSges(4,2) * t155;
t103 = rSges(3,1) * t143 - rSges(3,2) * t141 + t146;
t102 = -rSges(3,1) * t141 - rSges(3,2) * t143 - t218;
t89 = t188 * t143;
t88 = t188 * t141;
t65 = t135 * t154 * t86;
t64 = t226 + t146 + (pkin(2) - t212) * t143 + t191;
t63 = t210 - t218 + t134 + (-pkin(2) - t181) * t141;
t62 = t159 * t143;
t61 = t159 * t141;
t58 = -t141 * t150 + t122 + t146 + t161;
t57 = -t218 + (rSges(5,3) - t150) * t143 + (-t137 - t180) * t141;
t56 = t206 * t143;
t55 = t206 * t141;
t52 = t143 * (-t143 * t212 + t191) + (t141 * t181 - t210) * t141;
t45 = t160 + t182;
t44 = -t218 + (rSges(6,3) - t147) * t143 + (-t119 - t179) * t141;
t40 = t158 * t143;
t39 = t158 * t141;
t34 = -t136 * t54 - t198 * t87;
t33 = t136 * t53 + t199 * t87;
t32 = -t135 * t209 - t136 * t84 + t65;
t31 = t182 + t225;
t30 = -t218 - t143 * t147 + (-t216 - t119 + (-rSges(7,3) - pkin(9)) * t135) * t141 + t184;
t27 = t143 * t161 + (-rSges(5,3) * t143 + t141 * t180) * t141 + t215;
t26 = (-t141 * t54 + t143 * t53) * t135;
t14 = t187 + t41;
t11 = t23 + t187;
t1 = [t142 * (Icges(5,2) * t142 + t203) + t140 * (Icges(5,1) * t140 + t202) + t155 * (Icges(4,2) * t155 + t205) + t152 * (Icges(4,1) * t152 + t204) + Icges(2,3) + Icges(3,3) + t65 + (-t84 + t107) * t136 + (t108 - t209) * t135 + m(7) * (t30 ^ 2 + t31 ^ 2) + m(6) * (t44 ^ 2 + t45 ^ 2) + m(5) * (t57 ^ 2 + t58 ^ 2) + m(4) * (t63 ^ 2 + t64 ^ 2) + m(3) * (t102 ^ 2 + t103 ^ 2) + m(2) * (t128 ^ 2 + t129 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); m(4) * (-t141 * t64 - t143 * t63) * t127 + t157 + m(7) * (t30 * t40 + t31 * t39) + m(6) * (t44 * t62 + t45 * t61) + m(5) * (t57 * t89 + t58 * t88) + (t140 * (Icges(5,5) * t141 + t143 * t171) + t142 * (Icges(5,6) * t141 + t143 * t168) + t152 * (Icges(4,5) * t141 + t143 * t172) + t155 * (Icges(4,6) * t141 + t143 * t169)) * t221 + (t140 * (-Icges(5,5) * t143 + t141 * t171) + t142 * (-Icges(5,6) * t143 + t141 * t168) + t152 * (-Icges(4,5) * t143 + t141 * t172) + t155 * (-Icges(4,6) * t143 + t141 * t169)) * t220 + (Icges(4,5) * t152 + Icges(5,5) * t140 + Icges(4,6) * t155 + Icges(5,6) * t142) * (t138 / 0.2e1 + t139 / 0.2e1); m(4) * t52 + m(5) * t27 + m(6) * t14 + m(7) * t11; m(7) * (t11 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t14 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(5) * (t27 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(4) * (t127 ^ 2 * t190 + t52 ^ 2) + t219 + t222 * t141 * t138 + (t223 * t139 + (t141 * t223 + t143 * t222) * t141 + t224) * t143; m(7) * (t141 * t30 - t143 * t31) + m(6) * (t141 * t44 - t143 * t45) + m(5) * (t141 * t57 - t143 * t58); 0; m(7) * (t141 * t40 - t143 * t39) + m(6) * (t141 * t62 - t143 * t61) + m(5) * (t141 * t89 - t143 * t88); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t190; m(7) * (t30 * t56 + t31 * t55) + m(6) * (-t141 * t45 - t143 * t44) * t109 + t157; m(6) * t41 + m(7) * t23; m(7) * (t11 * t23 + t39 * t55 + t40 * t56) + m(6) * (t14 * t41 + (-t141 * t61 - t143 * t62) * t109) + t162; m(7) * (t141 * t56 - t143 * t55); m(6) * (t109 ^ 2 * t190 + t41 ^ 2) + m(7) * (t23 ^ 2 + t55 ^ 2 + t56 ^ 2) + t162; m(7) * (t30 * t33 + t31 * t34) - t32 * t136 + ((t29 / 0.2e1 + t22 / 0.2e1) * t143 + (t21 / 0.2e1 + t28 / 0.2e1) * t141) * t135; m(7) * t26; m(7) * (t11 * t26 + t33 * t40 + t34 * t39) + t186; m(7) * (t141 * t33 - t143 * t34); m(7) * (t23 * t26 + t33 * t56 + t34 * t55) + t186; t136 ^ 2 * t32 + m(7) * (t26 ^ 2 + t33 ^ 2 + t34 ^ 2) + (t143 * t4 + t141 * t3 - t136 * (t141 * t21 + t143 * t22)) * t135;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
