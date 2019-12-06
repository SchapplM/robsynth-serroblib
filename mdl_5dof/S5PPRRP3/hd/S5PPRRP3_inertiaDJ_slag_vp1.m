% Calculate time derivative of joint inertia matrix for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:34
% EndTime: 2019-12-05 15:10:44
% DurationCPUTime: 5.75s
% Computational Cost: add. (9522->517), mult. (27585->757), div. (0->0), fcn. (29682->8), ass. (0->210)
t190 = cos(pkin(8));
t191 = cos(pkin(7));
t193 = sin(qJ(3));
t215 = t191 * t193;
t189 = sin(pkin(7));
t194 = cos(qJ(3));
t216 = t189 * t194;
t177 = t190 * t216 - t215;
t188 = sin(pkin(8));
t192 = sin(qJ(4));
t219 = t188 * t192;
t227 = cos(qJ(4));
t157 = t177 * t227 + t189 * t219;
t214 = t191 * t194;
t217 = t189 * t193;
t176 = t190 * t217 + t214;
t201 = t188 * t227;
t196 = -t177 * t192 + t189 * t201;
t92 = Icges(6,5) * t157 + Icges(6,6) * t176 - Icges(6,3) * t196;
t98 = Icges(5,4) * t157 + Icges(5,2) * t196 + Icges(5,6) * t176;
t250 = -t92 + t98;
t179 = t190 * t214 + t217;
t159 = t179 * t227 + t191 * t219;
t178 = t190 * t215 - t216;
t195 = -t179 * t192 + t191 * t201;
t93 = Icges(6,5) * t159 + Icges(6,6) * t178 - Icges(6,3) * t195;
t99 = Icges(5,4) * t159 + Icges(5,2) * t195 + Icges(5,6) * t178;
t249 = -t93 + t99;
t94 = Icges(5,5) * t157 + Icges(5,6) * t196 + Icges(5,3) * t176;
t96 = Icges(6,4) * t157 + Icges(6,2) * t176 - Icges(6,6) * t196;
t248 = t94 + t96;
t95 = Icges(5,5) * t159 + Icges(5,6) * t195 + Icges(5,3) * t178;
t97 = Icges(6,4) * t159 + Icges(6,2) * t178 - Icges(6,6) * t195;
t247 = t95 + t97;
t100 = Icges(6,1) * t157 + Icges(6,4) * t176 - Icges(6,5) * t196;
t102 = Icges(5,1) * t157 + Icges(5,4) * t196 + Icges(5,5) * t176;
t246 = t100 + t102;
t101 = Icges(6,1) * t159 + Icges(6,4) * t178 - Icges(6,5) * t195;
t103 = Icges(5,1) * t159 + Icges(5,4) * t195 + Icges(5,5) * t178;
t245 = t101 + t103;
t238 = t246 * t157 + t248 * t176 + t250 * t196;
t237 = t245 * t157 + t247 * t176 + t249 * t196;
t236 = t246 * t159 + t248 * t178 + t250 * t195;
t235 = t245 * t159 + t247 * t178 + t249 * t195;
t180 = t190 * t227 + t194 * t219;
t198 = t194 * t201;
t181 = -t190 * t192 + t198;
t218 = t188 * t193;
t133 = Icges(6,5) * t181 + Icges(6,6) * t218 + Icges(6,3) * t180;
t135 = Icges(6,4) * t181 + Icges(6,2) * t218 + Icges(6,6) * t180;
t137 = Icges(6,1) * t181 + Icges(6,4) * t218 + Icges(6,5) * t180;
t62 = -t133 * t196 + t135 * t176 + t137 * t157;
t134 = Icges(5,5) * t181 - Icges(5,6) * t180 + Icges(5,3) * t218;
t136 = Icges(5,4) * t181 - Icges(5,2) * t180 + Icges(5,6) * t218;
t138 = Icges(5,1) * t181 - Icges(5,4) * t180 + Icges(5,5) * t218;
t63 = t134 * t176 + t136 * t196 + t138 * t157;
t244 = t62 + t63;
t64 = -t133 * t195 + t135 * t178 + t137 * t159;
t65 = t134 * t178 + t136 * t195 + t138 * t159;
t243 = t64 + t65;
t242 = rSges(6,1) + pkin(4);
t241 = rSges(6,3) + qJ(5);
t167 = t176 * qJD(3);
t121 = qJD(4) * t157 - t167 * t192;
t122 = qJD(4) * t196 - t167 * t227;
t168 = t177 * qJD(3);
t75 = Icges(6,5) * t122 + Icges(6,6) * t168 + Icges(6,3) * t121;
t79 = Icges(6,4) * t122 + Icges(6,2) * t168 + Icges(6,6) * t121;
t83 = Icges(6,1) * t122 + Icges(6,4) * t168 + Icges(6,5) * t121;
t14 = t100 * t122 + t121 * t92 + t157 * t83 + t168 * t96 + t176 * t79 - t196 * t75;
t169 = t178 * qJD(3);
t123 = qJD(4) * t159 - t169 * t192;
t124 = qJD(4) * t195 - t169 * t227;
t170 = t179 * qJD(3);
t76 = Icges(6,5) * t124 + Icges(6,6) * t170 + Icges(6,3) * t123;
t80 = Icges(6,4) * t124 + Icges(6,2) * t170 + Icges(6,6) * t123;
t84 = Icges(6,1) * t124 + Icges(6,4) * t170 + Icges(6,5) * t123;
t15 = t101 * t122 + t121 * t93 + t157 * t84 + t168 * t97 + t176 * t80 - t196 * t76;
t77 = Icges(5,5) * t122 - Icges(5,6) * t121 + Icges(5,3) * t168;
t81 = Icges(5,4) * t122 - Icges(5,2) * t121 + Icges(5,6) * t168;
t85 = Icges(5,1) * t122 - Icges(5,4) * t121 + Icges(5,5) * t168;
t16 = t102 * t122 - t121 * t98 + t157 * t85 + t168 * t94 + t176 * t77 + t196 * t81;
t78 = Icges(5,5) * t124 - Icges(5,6) * t123 + Icges(5,3) * t170;
t82 = Icges(5,4) * t124 - Icges(5,2) * t123 + Icges(5,6) * t170;
t86 = Icges(5,1) * t124 - Icges(5,4) * t123 + Icges(5,5) * t170;
t17 = t103 * t122 - t121 * t99 + t157 * t86 + t168 * t95 + t176 * t78 + t196 * t82;
t206 = qJD(3) * t194;
t205 = t188 * qJD(3);
t200 = t193 * t205;
t160 = -qJD(4) * t198 + (qJD(4) * t190 + t200) * t192;
t161 = -qJD(4) * t180 - t200 * t227;
t199 = t194 * t205;
t109 = Icges(6,5) * t161 + Icges(6,6) * t199 - Icges(6,3) * t160;
t111 = Icges(6,4) * t161 + Icges(6,2) * t199 - Icges(6,6) * t160;
t113 = Icges(6,1) * t161 + Icges(6,4) * t199 - Icges(6,5) * t160;
t29 = -t109 * t196 + t111 * t176 + t113 * t157 + t121 * t133 + t122 * t137 + t135 * t168;
t110 = Icges(5,5) * t161 + Icges(5,6) * t160 + Icges(5,3) * t199;
t112 = Icges(5,4) * t161 + Icges(5,2) * t160 + Icges(5,6) * t199;
t114 = Icges(5,1) * t161 + Icges(5,4) * t160 + Icges(5,5) * t199;
t30 = t110 * t176 + t112 * t196 + t114 * t157 - t121 * t136 + t122 * t138 + t134 * t168;
t240 = (t244 * t206 + (t29 + t30) * t193) * t188 + (t15 + t17) * t178 + (t14 + t16) * t176 + t237 * t170 + t238 * t168;
t18 = t100 * t124 + t123 * t92 + t159 * t83 + t170 * t96 + t178 * t79 - t195 * t75;
t19 = t101 * t124 + t123 * t93 + t159 * t84 + t170 * t97 + t178 * t80 - t195 * t76;
t20 = t102 * t124 - t123 * t98 + t159 * t85 + t170 * t94 + t178 * t77 + t195 * t81;
t21 = t103 * t124 - t123 * t99 + t159 * t86 + t170 * t95 + t178 * t78 + t195 * t82;
t31 = -t109 * t195 + t111 * t178 + t113 * t159 + t123 * t133 + t124 * t137 + t135 * t170;
t32 = t110 * t178 + t112 * t195 + t114 * t159 - t123 * t136 + t124 * t138 + t134 * t170;
t239 = (t243 * t206 + (t31 + t32) * t193) * t188 + (t19 + t21) * t178 + (t18 + t20) * t176 + t235 * t170 + t236 * t168;
t171 = (-Icges(4,5) * t193 - Icges(4,6) * t194) * t205;
t234 = t190 * t171;
t233 = 2 * m(5);
t232 = 2 * m(6);
t226 = rSges(6,2) * t168 - qJD(5) * t196 + t241 * t121 + t242 * t122;
t225 = rSges(6,2) * t170 - qJD(5) * t195 + t241 * t123 + t242 * t124;
t151 = -pkin(3) * t169 + pkin(6) * t170;
t90 = rSges(5,1) * t124 - rSges(5,2) * t123 + rSges(5,3) * t170;
t224 = -t151 - t90;
t223 = Icges(4,4) * t193;
t222 = Icges(4,4) * t194;
t221 = t188 * t189;
t220 = t188 * t191;
t213 = rSges(6,2) * t176 + t242 * t157 - t241 * t196;
t212 = rSges(6,2) * t178 + t242 * t159 - t241 * t195;
t107 = rSges(5,1) * t159 + rSges(5,2) * t195 + rSges(5,3) * t178;
t154 = pkin(3) * t179 + pkin(6) * t178;
t211 = -t107 - t154;
t210 = rSges(6,2) * t199 + qJD(5) * t180 - t241 * t160 + t242 * t161;
t150 = -pkin(3) * t167 + pkin(6) * t168;
t175 = (-pkin(3) * t193 + pkin(6) * t194) * t205;
t209 = t190 * t150 + t175 * t221;
t208 = rSges(6,2) * t218 + t241 * t180 + t242 * t181;
t153 = pkin(3) * t177 + pkin(6) * t176;
t182 = (pkin(3) * t194 + pkin(6) * t193) * t188;
t207 = t190 * t153 + t182 * t221;
t204 = -t151 - t225;
t202 = -t154 - t212;
t148 = -rSges(4,1) * t167 - rSges(4,2) * t168;
t174 = (-rSges(4,1) * t193 - rSges(4,2) * t194) * t205;
t115 = t148 * t190 + t174 * t221;
t149 = -rSges(4,1) * t169 - rSges(4,2) * t170;
t116 = -t149 * t190 - t174 * t220;
t197 = t115 * t189 - t116 * t191;
t173 = (-Icges(4,1) * t193 - t222) * t205;
t172 = (-Icges(4,2) * t194 - t223) * t205;
t165 = -Icges(4,5) * t190 + (Icges(4,1) * t194 - t223) * t188;
t164 = -Icges(4,6) * t190 + (-Icges(4,2) * t193 + t222) * t188;
t147 = -Icges(4,1) * t169 - Icges(4,4) * t170;
t146 = -Icges(4,1) * t167 - Icges(4,4) * t168;
t145 = -Icges(4,4) * t169 - Icges(4,2) * t170;
t144 = -Icges(4,4) * t167 - Icges(4,2) * t168;
t143 = -Icges(4,5) * t169 - Icges(4,6) * t170;
t142 = -Icges(4,5) * t167 - Icges(4,6) * t168;
t141 = t153 * t220;
t140 = rSges(5,1) * t181 - rSges(5,2) * t180 + rSges(5,3) * t218;
t132 = rSges(4,1) * t179 - rSges(4,2) * t178 + rSges(4,3) * t220;
t131 = rSges(4,1) * t177 - rSges(4,2) * t176 + rSges(4,3) * t221;
t130 = Icges(4,1) * t179 - Icges(4,4) * t178 + Icges(4,5) * t220;
t129 = Icges(4,1) * t177 - Icges(4,4) * t176 + Icges(4,5) * t221;
t128 = Icges(4,4) * t179 - Icges(4,2) * t178 + Icges(4,6) * t220;
t127 = Icges(4,4) * t177 - Icges(4,2) * t176 + Icges(4,6) * t221;
t125 = t150 * t220;
t118 = t161 * rSges(5,1) + t160 * rSges(5,2) + rSges(5,3) * t199;
t105 = rSges(5,1) * t157 + rSges(5,2) * t196 + rSges(5,3) * t176;
t91 = (t148 * t191 - t149 * t189) * t188;
t88 = rSges(5,1) * t122 - rSges(5,2) * t121 + rSges(5,3) * t168;
t72 = t107 * t218 - t140 * t178;
t71 = -t105 * t218 + t140 * t176;
t70 = t134 * t218 - t136 * t180 + t138 * t181;
t69 = t133 * t180 + t135 * t218 + t137 * t181;
t68 = t105 * t178 - t107 * t176;
t67 = t211 * t190 + (-t140 - t182) * t220;
t66 = t105 * t190 + t140 * t221 + t207;
t61 = t141 + (t105 * t191 + t189 * t211) * t188;
t60 = -t178 * t208 + t212 * t218;
t59 = t176 * t208 - t213 * t218;
t58 = t103 * t181 - t180 * t99 + t218 * t95;
t57 = t102 * t181 - t180 * t98 + t218 * t94;
t56 = t101 * t181 + t180 * t93 + t218 * t97;
t55 = t100 * t181 + t180 * t92 + t218 * t96;
t54 = t224 * t190 + (-t118 - t175) * t220;
t53 = t118 * t221 + t190 * t88 + t209;
t52 = t202 * t190 + (-t182 - t208) * t220;
t51 = t190 * t213 + t208 * t221 + t207;
t42 = -t176 * t212 + t178 * t213;
t41 = t125 + (t189 * t224 + t191 * t88) * t188;
t40 = -t178 * t118 - t170 * t140 + (t107 * t206 + t193 * t90) * t188;
t39 = t176 * t118 + t168 * t140 + (-t105 * t206 - t193 * t88) * t188;
t38 = t141 + (t189 * t202 + t191 * t213) * t188;
t37 = -t180 * t112 + t181 * t114 + t160 * t136 + t161 * t138 + (t110 * t193 + t134 * t206) * t188;
t36 = t180 * t109 + t181 * t113 - t160 * t133 + t161 * t137 + (t111 * t193 + t135 * t206) * t188;
t35 = t105 * t170 - t107 * t168 - t176 * t90 + t178 * t88;
t34 = t204 * t190 + (-t175 - t210) * t220;
t33 = t190 * t226 + t210 * t221 + t209;
t28 = t125 + (t189 * t204 + t191 * t226) * t188;
t27 = -t210 * t178 - t208 * t170 + (t193 * t225 + t206 * t212) * t188;
t26 = t210 * t176 + t208 * t168 + (-t193 * t226 - t206 * t213) * t188;
t25 = t161 * t103 + t160 * t99 - t180 * t82 + t181 * t86 + (t193 * t78 + t206 * t95) * t188;
t24 = t161 * t102 + t160 * t98 - t180 * t81 + t181 * t85 + (t193 * t77 + t206 * t94) * t188;
t23 = t161 * t101 - t160 * t93 + t180 * t76 + t181 * t84 + (t193 * t80 + t206 * t97) * t188;
t22 = t161 * t100 - t160 * t92 + t180 * t75 + t181 * t83 + (t193 * t79 + t206 * t96) * t188;
t13 = -t168 * t212 + t170 * t213 - t176 * t225 + t178 * t226;
t12 = -t190 * t37 + (t189 * t24 + t191 * t25) * t188;
t11 = -t190 * t36 + (t189 * t22 + t191 * t23) * t188;
t10 = -t190 * t32 + (t189 * t20 + t191 * t21) * t188;
t9 = -t190 * t31 + (t18 * t189 + t19 * t191) * t188;
t8 = -t190 * t30 + (t16 * t189 + t17 * t191) * t188;
t7 = -t190 * t29 + (t14 * t189 + t15 * t191) * t188;
t6 = t57 * t168 + t58 * t170 + t24 * t176 + t25 * t178 + (t193 * t37 + t206 * t70) * t188;
t5 = t55 * t168 + t56 * t170 + t22 * t176 + t23 * t178 + (t193 * t36 + t206 * t69) * t188;
t1 = [0; 0; 0; m(4) * t91 + m(5) * t41 + m(6) * t28; m(4) * t197 + m(5) * (t189 * t53 - t191 * t54) + m(6) * (t189 * t33 - t191 * t34); (t28 * t38 + t33 * t51 + t34 * t52) * t232 + (t41 * t61 + t53 * t66 + t54 * t67) * t233 - t190 * t11 - t190 * t12 - t190 * (t190 ^ 2 * t171 + (((-t145 * t193 + t147 * t194) * t191 + (-t144 * t193 + t146 * t194) * t189 + ((-t128 * t194 - t130 * t193) * t191 + (-t127 * t194 - t129 * t193) * t189) * qJD(3)) * t188 + (-t142 * t189 - t143 * t191 + t172 * t193 - t173 * t194 + (t164 * t194 + t165 * t193) * qJD(3)) * t190) * t188) + 0.2e1 * m(4) * ((t115 * t131 - t116 * t132) * t190 + ((t131 * t191 - t132 * t189) * t91 + t197 * (-t190 * rSges(4,3) + (rSges(4,1) * t194 - rSges(4,2) * t193) * t188)) * t188) + (t7 + t8 - (-t168 * t164 - t167 * t165 - t176 * t172 + t177 * t173) * t190 + (-t127 * t168 - t129 * t167 + t142 * t221 - t144 * t176 + t146 * t177 - t234) * t221) * t221 + (t10 + t9 - (-t170 * t164 - t169 * t165 - t178 * t172 + t179 * t173) * t190 + (-t128 * t170 - t130 * t169 + t143 * t220 - t145 * t178 + t147 * t179 - t234) * t220 + (-t127 * t170 - t128 * t168 - t129 * t169 - t130 * t167 + t142 * t220 + t143 * t221 - t144 * t178 - t145 * t176 + t146 * t179 + t147 * t177) * t221) * t220; m(5) * t35 + m(6) * t13; m(5) * (t189 * t39 - t191 * t40) + m(6) * (t189 * t26 - t191 * t27); (t10 / 0.2e1 + t9 / 0.2e1) * t178 + (t8 / 0.2e1 + t7 / 0.2e1) * t176 + m(6) * (t13 * t38 + t26 * t51 + t27 * t52 + t28 * t42 + t33 * t59 + t34 * t60) + m(5) * (t35 * t61 + t39 * t66 + t40 * t67 + t41 * t68 + t53 * t71 + t54 * t72) + (-t5 / 0.2e1 - t6 / 0.2e1 + (-t65 / 0.2e1 - t64 / 0.2e1) * t170 + (-t63 / 0.2e1 - t62 / 0.2e1) * t168) * t190 + ((t12 / 0.2e1 + t11 / 0.2e1) * t193 + (((t58 / 0.2e1 + t56 / 0.2e1) * t191 + (t57 / 0.2e1 + t55 / 0.2e1) * t189) * t188 + (-t70 / 0.2e1 - t69 / 0.2e1) * t190) * t206 + (t238 * t189 + t237 * t191) * t168 / 0.2e1 + (t236 * t189 + t235 * t191) * t170 / 0.2e1 + t240 * t189 / 0.2e1 + t239 * t191 / 0.2e1) * t188; (t5 + t6) * t218 + t239 * t178 + t240 * t176 + (t236 * t176 + t235 * t178 + t243 * t218) * t170 + (t238 * t176 + t237 * t178 + t244 * t218) * t168 + ((t69 + t70) * t218 + (t56 + t58) * t178 + (t55 + t57) * t176) * t199 + (t13 * t42 + t26 * t59 + t27 * t60) * t232 + (t35 * t68 + t39 * t71 + t40 * t72) * t233; -m(6) * t160; m(6) * (-t121 * t191 + t123 * t189); m(6) * (t121 * t52 + t123 * t51 - t160 * t38 + t180 * t28 - t195 * t33 - t196 * t34); m(6) * (t121 * t60 + t123 * t59 + t13 * t180 - t160 * t42 - t195 * t26 - t196 * t27); (-t121 * t196 - t123 * t195 - t160 * t180) * t232;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
