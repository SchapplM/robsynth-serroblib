% Calculate joint inertia matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:46
% EndTime: 2019-12-05 16:47:55
% DurationCPUTime: 3.22s
% Computational Cost: add. (6233->262), mult. (8859->404), div. (0->0), fcn. (9654->8), ass. (0->143)
t243 = Icges(5,1) + Icges(6,1);
t242 = Icges(5,4) + Icges(6,4);
t241 = Icges(6,5) + Icges(5,5);
t240 = Icges(5,2) + Icges(6,2);
t239 = Icges(6,6) + Icges(5,6);
t238 = Icges(5,3) + Icges(6,3);
t157 = qJ(3) + qJ(4);
t150 = sin(t157);
t151 = cos(t157);
t159 = cos(pkin(8));
t158 = sin(pkin(8));
t163 = cos(qJ(2));
t200 = t158 * t163;
t130 = -t150 * t200 - t151 * t159;
t131 = -t150 * t159 + t151 * t200;
t161 = sin(qJ(2));
t201 = t158 * t161;
t237 = t130 * t239 + t131 * t241 + t201 * t238;
t197 = t159 * t163;
t132 = -t150 * t197 + t151 * t158;
t133 = t150 * t158 + t151 * t197;
t198 = t159 * t161;
t236 = t132 * t239 + t133 * t241 + t198 * t238;
t235 = t130 * t240 + t131 * t242 + t201 * t239;
t234 = t132 * t240 + t133 * t242 + t198 * t239;
t233 = t130 * t242 + t131 * t243 + t201 * t241;
t232 = t132 * t242 + t133 * t243 + t198 * t241;
t230 = t238 * t163 + (t150 * t239 - t151 * t241) * t161;
t229 = -t130 * t235 - t131 * t233 - t201 * t237;
t228 = t130 * t234 + t131 * t232 + t201 * t236;
t227 = t132 * t235 + t133 * t233 + t198 * t237;
t226 = t132 * t234 + t133 * t232 + t198 * t236;
t225 = t239 * t163 + (t150 * t240 - t151 * t242) * t161;
t224 = -t241 * t163 + (-t150 * t242 + t151 * t243) * t161;
t223 = t230 * t163;
t222 = t237 * t163 + (t150 * t235 - t151 * t233) * t161;
t221 = -t236 * t163 + (-t150 * t234 + t151 * t232) * t161;
t220 = (t130 * t225 - t131 * t224) * t163 + (t228 * t159 + (t223 - t229) * t158) * t161;
t219 = (t132 * t225 - t133 * t224) * t163 + ((t223 + t226) * t159 + t227 * t158) * t161;
t218 = t158 * t228 + t159 * t229;
t217 = t158 * t226 - t159 * t227;
t187 = pkin(4) * t151;
t165 = qJ(5) * t161 + t163 * t187;
t178 = pkin(4) * t150;
t216 = rSges(6,1) * t131 + rSges(6,2) * t130 + rSges(6,3) * t201 + t158 * t165 - t159 * t178;
t215 = (-qJ(5) - rSges(6,3)) * t163 + (rSges(6,1) * t151 - rSges(6,2) * t150 + t187) * t161;
t154 = t158 ^ 2;
t155 = t159 ^ 2;
t185 = t154 + t155;
t213 = t163 ^ 2;
t214 = t230 * t213 + ((t150 * t225 + t151 * t224) * t163 - t221 * t159 + t222 * t158) * t161;
t212 = t158 / 0.2e1;
t211 = -t159 / 0.2e1;
t210 = -t163 / 0.2e1;
t162 = cos(qJ(3));
t208 = t162 * pkin(3);
t206 = t216 * t198;
t205 = rSges(6,1) * t133 + rSges(6,2) * t132 + rSges(6,3) * t198 + t158 * t178 + t159 * t165;
t167 = pkin(7) * t161 + t163 * t208;
t160 = sin(qJ(3));
t202 = t158 * t160;
t108 = pkin(3) * t202 + t159 * t167;
t92 = rSges(5,1) * t133 + rSges(5,2) * t132 + rSges(5,3) * t198;
t204 = -t108 - t92;
t199 = t159 * t160;
t107 = -pkin(3) * t199 + t158 * t167;
t115 = -pkin(7) * t163 + t161 * t208;
t203 = t107 * t163 + t115 * t201;
t123 = -t163 * rSges(5,3) + (rSges(5,1) * t151 - rSges(5,2) * t150) * t161;
t90 = rSges(5,1) * t131 + rSges(5,2) * t130 + rSges(5,3) * t201;
t63 = t123 * t201 + t163 * t90;
t196 = t160 * t163;
t195 = t162 * t163;
t134 = -Icges(4,3) * t163 + (Icges(4,5) * t162 - Icges(4,6) * t160) * t161;
t192 = t163 * t134;
t190 = -t115 - t123;
t137 = -t163 * rSges(4,3) + (rSges(4,1) * t162 - rSges(4,2) * t160) * t161;
t148 = t161 * pkin(2) - t163 * pkin(6);
t189 = -t137 - t148;
t188 = t185 * (pkin(2) * t163 + pkin(6) * t161);
t184 = -t108 - t205;
t183 = t198 * t219 + t201 * t220;
t182 = -t115 - t215;
t181 = -t148 + t190;
t36 = t163 * t216 + t201 * t215;
t177 = t107 * t158 + t108 * t159 + t188;
t176 = -t148 + t182;
t171 = Icges(3,5) * t163 - Icges(3,6) * t161;
t168 = t163 * t214 + t183;
t166 = t219 * t212 + t220 * t211 + (t158 * t221 + t159 * t222) * t210 + t218 * t201 / 0.2e1 + t217 * t198 / 0.2e1;
t147 = rSges(3,1) * t161 + rSges(3,2) * t163;
t143 = t159 * t195 + t202;
t142 = t158 * t162 - t159 * t196;
t141 = t158 * t195 - t199;
t140 = -t158 * t196 - t159 * t162;
t136 = -Icges(4,5) * t163 + (Icges(4,1) * t162 - Icges(4,4) * t160) * t161;
t135 = -Icges(4,6) * t163 + (Icges(4,4) * t162 - Icges(4,2) * t160) * t161;
t125 = Icges(3,3) * t158 + t159 * t171;
t124 = -Icges(3,3) * t159 + t158 * t171;
t111 = t189 * t159;
t110 = t189 * t158;
t106 = rSges(4,1) * t143 + rSges(4,2) * t142 + rSges(4,3) * t198;
t105 = rSges(4,1) * t141 + rSges(4,2) * t140 + rSges(4,3) * t201;
t104 = Icges(4,1) * t143 + Icges(4,4) * t142 + Icges(4,5) * t198;
t103 = Icges(4,1) * t141 + Icges(4,4) * t140 + Icges(4,5) * t201;
t102 = Icges(4,4) * t143 + Icges(4,2) * t142 + Icges(4,6) * t198;
t101 = Icges(4,4) * t141 + Icges(4,2) * t140 + Icges(4,6) * t201;
t100 = Icges(4,5) * t143 + Icges(4,6) * t142 + Icges(4,3) * t198;
t99 = Icges(4,5) * t141 + Icges(4,6) * t140 + Icges(4,3) * t201;
t94 = t185 * (rSges(3,1) * t163 - rSges(3,2) * t161);
t93 = t107 * t198;
t74 = t90 * t198;
t72 = t181 * t159;
t71 = t181 * t158;
t68 = -t106 * t163 - t137 * t198;
t67 = t105 * t163 + t137 * t201;
t64 = -t123 * t198 - t163 * t92;
t62 = (t105 * t159 - t106 * t158) * t161;
t61 = t176 * t159;
t60 = t176 * t158;
t59 = -t201 * t92 + t74;
t58 = t105 * t158 + t106 * t159 + t188;
t57 = -t163 * t100 + (-t102 * t160 + t104 * t162) * t161;
t56 = -t163 * t99 + (-t101 * t160 + t103 * t162) * t161;
t55 = t163 * t204 + t190 * t198;
t54 = t63 + t203;
t53 = t100 * t198 + t102 * t142 + t104 * t143;
t52 = t101 * t142 + t103 * t143 + t198 * t99;
t51 = t100 * t201 + t102 * t140 + t104 * t141;
t50 = t101 * t140 + t103 * t141 + t201 * t99;
t37 = -t163 * t205 - t198 * t215;
t35 = t201 * t204 + t74 + t93;
t34 = t158 * t90 + t159 * t92 + t177;
t33 = -t201 * t205 + t206;
t32 = t163 * t184 + t182 * t198;
t31 = t36 + t203;
t30 = t158 * t53 - t159 * t52;
t29 = t158 * t51 - t159 * t50;
t26 = t184 * t201 + t206 + t93;
t25 = t158 * t216 + t159 * t205 + t177;
t14 = -(t135 * t142 + t136 * t143) * t163 + (t52 * t158 + (t53 - t192) * t159) * t161;
t13 = -(t135 * t140 + t136 * t141) * t163 + (t51 * t159 + (t50 - t192) * t158) * t161;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t94 + m(4) * t58 + m(5) * t34 + m(6) * t25; m(5) * (t34 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t25 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(4) * (t110 ^ 2 + t111 ^ 2 + t58 ^ 2) + m(3) * (t147 ^ 2 * t185 + t94 ^ 2) + (-t155 * t124 - t218 - t29) * t159 + (t154 * t125 + t30 + (-t158 * t124 + t159 * t125) * t159 + t217) * t158; m(4) * t62 + m(5) * t35 + m(6) * t26; m(5) * (t34 * t35 + t54 * t72 + t55 * t71) + m(6) * (t25 * t26 + t31 * t61 + t32 * t60) + m(4) * (t110 * t68 + t111 * t67 + t58 * t62) + (t29 * t212 + t159 * t30 / 0.2e1) * t161 + t166 + t14 * t212 + t13 * t211 + (t158 * t57 - t159 * t56) * t210; t14 * t198 + t13 * t201 + m(6) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(5) * (t35 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(4) * (t62 ^ 2 + t67 ^ 2 + t68 ^ 2) + t183 + (-t213 * t134 + (-t56 * t158 - t57 * t159 - (t135 * t160 - t136 * t162) * t163) * t161 + t214) * t163; m(5) * t59 + m(6) * t33; m(5) * (t34 * t59 + t63 * t72 + t64 * t71) + m(6) * (t25 * t33 + t36 * t61 + t37 * t60) + t166; m(6) * (t26 * t33 + t31 * t36 + t32 * t37) + m(5) * (t35 * t59 + t54 * t63 + t55 * t64) + t168; m(6) * (t33 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(5) * (t59 ^ 2 + t63 ^ 2 + t64 ^ 2) + t168; -m(6) * t163; m(6) * (-t163 * t25 + (t158 * t60 + t159 * t61) * t161); m(6) * (-t163 * t26 + (t158 * t32 + t159 * t31) * t161); m(6) * (-t163 * t33 + (t158 * t37 + t159 * t36) * t161); m(6) * (t161 ^ 2 * t185 + t213);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
