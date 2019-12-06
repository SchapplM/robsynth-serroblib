% Calculate joint inertia matrix for
% S5PRRRP6
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP6_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP6_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:48
% EndTime: 2019-12-05 16:50:56
% DurationCPUTime: 2.93s
% Computational Cost: add. (5771->257), mult. (8439->397), div. (0->0), fcn. (9335->8), ass. (0->143)
t236 = Icges(5,1) + Icges(6,1);
t225 = -Icges(5,4) + Icges(6,5);
t235 = Icges(6,4) + Icges(5,5);
t234 = Icges(5,2) + Icges(6,3);
t233 = Icges(6,6) - Icges(5,6);
t232 = Icges(5,3) + Icges(6,2);
t151 = qJ(3) + qJ(4);
t149 = sin(t151);
t150 = cos(t151);
t153 = cos(pkin(8));
t152 = sin(pkin(8));
t157 = cos(qJ(2));
t189 = t152 * t157;
t129 = t149 * t189 + t153 * t150;
t130 = -t153 * t149 + t150 * t189;
t155 = sin(qJ(2));
t190 = t152 * t155;
t231 = t234 * t129 + t225 * t130 + t233 * t190;
t186 = t153 * t157;
t131 = t149 * t186 - t152 * t150;
t132 = t152 * t149 + t150 * t186;
t187 = t153 * t155;
t230 = t234 * t131 + t225 * t132 + t233 * t187;
t229 = t233 * t129 + t235 * t130 + t232 * t190;
t228 = t233 * t131 + t235 * t132 + t232 * t187;
t227 = t225 * t129 + t236 * t130 + t235 * t190;
t226 = t225 * t131 + t236 * t132 + t235 * t187;
t224 = t232 * t157 + (-t233 * t149 - t235 * t150) * t155;
t223 = rSges(6,1) + pkin(4);
t222 = -t231 * t129 - t227 * t130 - t229 * t190;
t221 = t230 * t129 + t226 * t130 + t228 * t190;
t220 = t231 * t131 + t227 * t132 + t229 * t187;
t219 = t230 * t131 + t226 * t132 + t228 * t187;
t218 = rSges(6,3) + qJ(5);
t217 = -t233 * t157 + (t234 * t149 + t225 * t150) * t155;
t216 = -t235 * t157 + (t225 * t149 + t236 * t150) * t155;
t215 = t224 * t157;
t214 = t229 * t157 + (-t231 * t149 - t227 * t150) * t155;
t213 = -t228 * t157 + (t230 * t149 + t226 * t150) * t155;
t212 = (-t217 * t129 - t216 * t130) * t157 + (t221 * t153 + (t215 - t222) * t152) * t155;
t211 = (-t217 * t131 - t216 * t132) * t157 + ((t215 + t219) * t153 + t220 * t152) * t155;
t210 = t221 * t152 + t222 * t153;
t209 = t219 * t152 - t220 * t153;
t208 = rSges(6,2) * t190 + t218 * t129 + t223 * t130;
t207 = -t157 * rSges(6,2) + (t218 * t149 + t223 * t150) * t155;
t203 = t153 ^ 2;
t204 = t152 ^ 2;
t205 = t203 + t204;
t202 = t157 ^ 2;
t206 = t224 * t202 + ((t217 * t149 + t216 * t150) * t157 - t213 * t153 + t214 * t152) * t155;
t201 = t152 / 0.2e1;
t200 = -t153 / 0.2e1;
t199 = -t157 / 0.2e1;
t156 = cos(qJ(3));
t198 = t156 * pkin(3);
t196 = t208 * t187;
t195 = rSges(6,2) * t187 + t218 * t131 + t223 * t132;
t160 = pkin(7) * t155 + t198 * t157;
t154 = sin(qJ(3));
t191 = t152 * t154;
t107 = pkin(3) * t191 + t160 * t153;
t88 = t132 * rSges(5,1) - t131 * rSges(5,2) + rSges(5,3) * t187;
t194 = -t107 - t88;
t188 = t153 * t154;
t106 = -pkin(3) * t188 + t160 * t152;
t114 = -pkin(7) * t157 + t198 * t155;
t193 = t157 * t106 + t114 * t190;
t122 = -t157 * rSges(5,3) + (rSges(5,1) * t150 - rSges(5,2) * t149) * t155;
t86 = t130 * rSges(5,1) - t129 * rSges(5,2) + rSges(5,3) * t190;
t63 = t122 * t190 + t157 * t86;
t192 = t149 * t155;
t185 = t154 * t157;
t184 = t156 * t157;
t133 = -Icges(4,3) * t157 + (Icges(4,5) * t156 - Icges(4,6) * t154) * t155;
t181 = t157 * t133;
t180 = -t114 - t122;
t136 = -t157 * rSges(4,3) + (rSges(4,1) * t156 - rSges(4,2) * t154) * t155;
t146 = t155 * pkin(2) - t157 * pkin(6);
t178 = -t136 - t146;
t177 = t205 * (pkin(2) * t157 + pkin(6) * t155);
t176 = -t107 - t195;
t175 = t211 * t187 + t212 * t190;
t174 = -t114 - t207;
t173 = -t146 + t180;
t54 = t208 * t157 + t207 * t190;
t170 = t152 * t106 + t153 * t107 + t177;
t169 = -t146 + t174;
t164 = Icges(3,5) * t157 - Icges(3,6) * t155;
t161 = t206 * t157 + t175;
t159 = t211 * t201 + t212 * t200 + (t213 * t152 + t214 * t153) * t199 + t210 * t190 / 0.2e1 + t209 * t187 / 0.2e1;
t145 = t155 * rSges(3,1) + t157 * rSges(3,2);
t143 = t153 * t184 + t191;
t142 = t152 * t156 - t153 * t185;
t141 = t152 * t184 - t188;
t140 = -t152 * t185 - t153 * t156;
t135 = -Icges(4,5) * t157 + (Icges(4,1) * t156 - Icges(4,4) * t154) * t155;
t134 = -Icges(4,6) * t157 + (Icges(4,4) * t156 - Icges(4,2) * t154) * t155;
t124 = Icges(3,3) * t152 + t164 * t153;
t123 = -Icges(3,3) * t153 + t164 * t152;
t109 = t178 * t153;
t108 = t178 * t152;
t105 = t143 * rSges(4,1) + t142 * rSges(4,2) + rSges(4,3) * t187;
t104 = t141 * rSges(4,1) + t140 * rSges(4,2) + rSges(4,3) * t190;
t101 = Icges(4,1) * t143 + Icges(4,4) * t142 + Icges(4,5) * t187;
t100 = Icges(4,1) * t141 + Icges(4,4) * t140 + Icges(4,5) * t190;
t99 = Icges(4,4) * t143 + Icges(4,2) * t142 + Icges(4,6) * t187;
t98 = Icges(4,4) * t141 + Icges(4,2) * t140 + Icges(4,6) * t190;
t97 = Icges(4,5) * t143 + Icges(4,6) * t142 + Icges(4,3) * t187;
t96 = Icges(4,5) * t141 + Icges(4,6) * t140 + Icges(4,3) * t190;
t91 = t205 * (rSges(3,1) * t157 - rSges(3,2) * t155);
t90 = t106 * t187;
t70 = t86 * t187;
t68 = t173 * t153;
t67 = t173 * t152;
t66 = -t157 * t105 - t136 * t187;
t65 = t157 * t104 + t136 * t190;
t64 = -t122 * t187 - t157 * t88;
t62 = t169 * t153;
t61 = t169 * t152;
t60 = (t104 * t153 - t105 * t152) * t155;
t59 = -t88 * t190 + t70;
t58 = t152 * t104 + t153 * t105 + t177;
t57 = -t157 * t97 + (t101 * t156 - t154 * t99) * t155;
t56 = -t157 * t96 + (t100 * t156 - t154 * t98) * t155;
t55 = -t195 * t157 - t187 * t207;
t53 = t194 * t157 + t180 * t187;
t52 = t63 + t193;
t51 = t143 * t101 + t142 * t99 + t97 * t187;
t50 = t143 * t100 + t142 * t98 + t96 * t187;
t49 = t141 * t101 + t140 * t99 + t97 * t190;
t48 = t141 * t100 + t140 * t98 + t96 * t190;
t35 = t194 * t190 + t70 + t90;
t34 = -t195 * t190 + t196;
t33 = t176 * t157 + t174 * t187;
t32 = t54 + t193;
t31 = t152 * t86 + t153 * t88 + t170;
t30 = t176 * t190 + t196 + t90;
t29 = t208 * t152 + t195 * t153 + t170;
t28 = t51 * t152 - t50 * t153;
t27 = t49 * t152 - t48 * t153;
t14 = -(t142 * t134 + t143 * t135) * t157 + (t50 * t152 + (t51 - t181) * t153) * t155;
t13 = -(t140 * t134 + t141 * t135) * t157 + (t49 * t153 + (t48 - t181) * t152) * t155;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6); m(3) * t91 + m(4) * t58 + m(5) * t31 + m(6) * t29; m(5) * (t31 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t29 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(4) * (t108 ^ 2 + t109 ^ 2 + t58 ^ 2) + m(3) * (t205 * t145 ^ 2 + t91 ^ 2) + (-t203 * t123 - t210 - t27) * t153 + (t204 * t124 + t28 + (-t152 * t123 + t153 * t124) * t153 + t209) * t152; m(4) * t60 + m(5) * t35 + m(6) * t30; t159 + m(5) * (t35 * t31 + t52 * t68 + t53 * t67) + m(6) * (t30 * t29 + t32 * t62 + t33 * t61) + m(4) * (t66 * t108 + t65 * t109 + t60 * t58) + (t153 * t28 / 0.2e1 + t27 * t201) * t155 + t14 * t201 + t13 * t200 + (t57 * t152 - t56 * t153) * t199; t13 * t190 + t14 * t187 + m(6) * (t30 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(5) * (t35 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(4) * (t60 ^ 2 + t65 ^ 2 + t66 ^ 2) + t175 + (-t202 * t133 + (-t56 * t152 - t57 * t153 - (t134 * t154 - t135 * t156) * t157) * t155 + t206) * t157; m(5) * t59 + m(6) * t34; m(5) * (t59 * t31 + t63 * t68 + t64 * t67) + m(6) * (t34 * t29 + t54 * t62 + t55 * t61) + t159; m(6) * (t34 * t30 + t54 * t32 + t55 * t33) + m(5) * (t59 * t35 + t63 * t52 + t64 * t53) + t161; m(6) * (t34 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t59 ^ 2 + t63 ^ 2 + t64 ^ 2) + t161; m(6) * t192; m(6) * (t129 * t61 + t131 * t62 + t29 * t192); m(6) * (t129 * t33 + t131 * t32 + t30 * t192); m(6) * (t129 * t55 + t131 * t54 + t34 * t192); m(6) * (t155 ^ 2 * t149 ^ 2 + t129 ^ 2 + t131 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
