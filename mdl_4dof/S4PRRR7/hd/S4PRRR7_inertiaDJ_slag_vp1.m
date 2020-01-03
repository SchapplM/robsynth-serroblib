% Calculate time derivative of joint inertia matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR7_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR7_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:04
% EndTime: 2019-12-31 16:36:15
% DurationCPUTime: 5.76s
% Computational Cost: add. (23809->609), mult. (67543->916), div. (0->0), fcn. (77136->10), ass. (0->251)
t216 = sin(pkin(8));
t218 = cos(pkin(8));
t224 = cos(qJ(2));
t219 = cos(pkin(4));
t222 = sin(qJ(2));
t248 = t219 * t222;
t205 = t216 * t224 + t218 * t248;
t196 = t205 * qJD(2);
t238 = t216 * t248;
t239 = qJD(2) * t224;
t198 = -qJD(2) * t238 + t218 * t239;
t247 = t219 * t224;
t206 = t216 * t247 + t218 * t222;
t217 = sin(pkin(4));
t227 = -t216 * t222 + t218 * t247;
t240 = qJD(2) * t222;
t221 = sin(qJ(3));
t250 = t217 * t221;
t258 = cos(qJ(3));
t185 = t205 * t258 - t218 * t250;
t220 = sin(qJ(4));
t223 = cos(qJ(4));
t151 = -t185 * t220 - t223 * t227;
t152 = t185 * t223 - t220 * t227;
t235 = t217 * t258;
t225 = -t205 * t221 - t218 * t235;
t103 = Icges(5,5) * t152 + Icges(5,6) * t151 - Icges(5,3) * t225;
t105 = Icges(5,4) * t152 + Icges(5,2) * t151 - Icges(5,6) * t225;
t107 = Icges(5,1) * t152 + Icges(5,4) * t151 - Icges(5,5) * t225;
t195 = t227 * qJD(2);
t148 = qJD(3) * t225 + t195 * t258;
t111 = -qJD(4) * t152 - t148 * t220 + t196 * t223;
t112 = qJD(4) * t151 + t148 * t223 + t196 * t220;
t147 = qJD(3) * t185 + t195 * t221;
t71 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t147;
t73 = Icges(5,4) * t112 + Icges(5,2) * t111 + Icges(5,6) * t147;
t75 = Icges(5,1) * t112 + Icges(5,4) * t111 + Icges(5,5) * t147;
t16 = t103 * t147 + t105 * t111 + t107 * t112 + t151 * t73 + t152 * t75 - t225 * t71;
t207 = t218 * t224 - t238;
t187 = t207 * t258 + t216 * t250;
t153 = -t187 * t220 + t206 * t223;
t154 = t187 * t223 + t206 * t220;
t226 = -t207 * t221 + t216 * t235;
t104 = Icges(5,5) * t154 + Icges(5,6) * t153 - Icges(5,3) * t226;
t106 = Icges(5,4) * t154 + Icges(5,2) * t153 - Icges(5,6) * t226;
t108 = Icges(5,1) * t154 + Icges(5,4) * t153 - Icges(5,5) * t226;
t197 = t206 * qJD(2);
t150 = qJD(3) * t226 - t197 * t258;
t113 = -qJD(4) * t154 - t150 * t220 + t198 * t223;
t114 = qJD(4) * t153 + t150 * t223 + t198 * t220;
t149 = qJD(3) * t187 - t197 * t221;
t72 = Icges(5,5) * t114 + Icges(5,6) * t113 + Icges(5,3) * t149;
t74 = Icges(5,4) * t114 + Icges(5,2) * t113 + Icges(5,6) * t149;
t76 = Icges(5,1) * t114 + Icges(5,4) * t113 + Icges(5,5) * t149;
t17 = t104 * t147 + t106 * t111 + t108 * t112 + t151 * t74 + t152 * t76 - t225 * t72;
t209 = t219 * t221 + t222 * t235;
t249 = t217 * t224;
t190 = -t209 * t220 - t223 * t249;
t208 = -t219 * t258 + t222 * t250;
t228 = -t209 * t223 + t220 * t249;
t134 = -Icges(5,5) * t228 + Icges(5,6) * t190 + Icges(5,3) * t208;
t135 = -Icges(5,4) * t228 + Icges(5,2) * t190 + Icges(5,6) * t208;
t136 = -Icges(5,1) * t228 + Icges(5,4) * t190 + Icges(5,5) * t208;
t233 = t217 * t239;
t189 = -qJD(3) * t208 + t233 * t258;
t234 = t217 * t240;
t138 = qJD(4) * t228 - t189 * t220 + t223 * t234;
t139 = qJD(4) * t190 + t189 * t223 + t220 * t234;
t188 = qJD(3) * t209 + t221 * t233;
t96 = Icges(5,5) * t139 + Icges(5,6) * t138 + Icges(5,3) * t188;
t97 = Icges(5,4) * t139 + Icges(5,2) * t138 + Icges(5,6) * t188;
t98 = Icges(5,1) * t139 + Icges(5,4) * t138 + Icges(5,5) * t188;
t32 = t111 * t135 + t112 * t136 + t134 * t147 + t151 * t97 + t152 * t98 - t225 * t96;
t49 = -t103 * t225 + t105 * t151 + t107 * t152;
t50 = -t104 * t225 + t106 * t151 + t108 * t152;
t66 = -t134 * t225 + t135 * t151 + t136 * t152;
t3 = -t16 * t227 + t17 * t206 + t49 * t196 + t50 * t198 + (-t224 * t32 + t240 * t66) * t217;
t115 = Icges(4,5) * t148 - Icges(4,6) * t147 + Icges(4,3) * t196;
t117 = Icges(4,4) * t148 - Icges(4,2) * t147 + Icges(4,6) * t196;
t119 = Icges(4,1) * t148 - Icges(4,4) * t147 + Icges(4,5) * t196;
t126 = Icges(4,5) * t185 + Icges(4,6) * t225 - Icges(4,3) * t227;
t128 = Icges(4,4) * t185 + Icges(4,2) * t225 - Icges(4,6) * t227;
t130 = Icges(4,1) * t185 + Icges(4,4) * t225 - Icges(4,5) * t227;
t41 = -t115 * t227 + t117 * t225 + t119 * t185 + t126 * t196 - t128 * t147 + t130 * t148;
t116 = Icges(4,5) * t150 - Icges(4,6) * t149 + Icges(4,3) * t198;
t118 = Icges(4,4) * t150 - Icges(4,2) * t149 + Icges(4,6) * t198;
t120 = Icges(4,1) * t150 - Icges(4,4) * t149 + Icges(4,5) * t198;
t127 = Icges(4,5) * t187 + Icges(4,6) * t226 + Icges(4,3) * t206;
t129 = Icges(4,4) * t187 + Icges(4,2) * t226 + Icges(4,6) * t206;
t131 = Icges(4,1) * t187 + Icges(4,4) * t226 + Icges(4,5) * t206;
t42 = -t116 * t227 + t118 * t225 + t120 * t185 + t127 * t196 - t129 * t147 + t131 * t148;
t140 = Icges(4,5) * t189 - Icges(4,6) * t188 + Icges(4,3) * t234;
t141 = Icges(4,4) * t189 - Icges(4,2) * t188 + Icges(4,6) * t234;
t142 = Icges(4,1) * t189 - Icges(4,4) * t188 + Icges(4,5) * t234;
t164 = Icges(4,5) * t209 - Icges(4,6) * t208 - Icges(4,3) * t249;
t165 = Icges(4,4) * t209 - Icges(4,2) * t208 - Icges(4,6) * t249;
t166 = Icges(4,1) * t209 - Icges(4,4) * t208 - Icges(4,5) * t249;
t53 = -t140 * t227 + t141 * t225 + t142 * t185 - t147 * t165 + t148 * t166 + t164 * t196;
t79 = -t126 * t227 + t128 * t225 + t130 * t185;
t80 = -t127 * t227 + t129 * t225 + t131 * t185;
t91 = -t164 * t227 + t165 * t225 + t166 * t185;
t273 = t79 * t196 + t80 * t198 - t41 * t227 + t42 * t206 + (-t224 * t53 + t240 * t91) * t217 + t3;
t18 = t103 * t149 + t105 * t113 + t107 * t114 + t153 * t73 + t154 * t75 - t226 * t71;
t19 = t104 * t149 + t106 * t113 + t108 * t114 + t153 * t74 + t154 * t76 - t226 * t72;
t33 = t113 * t135 + t114 * t136 + t134 * t149 + t153 * t97 + t154 * t98 - t226 * t96;
t51 = -t103 * t226 + t105 * t153 + t107 * t154;
t52 = -t104 * t226 + t106 * t153 + t108 * t154;
t67 = -t134 * t226 + t135 * t153 + t136 * t154;
t4 = -t18 * t227 + t19 * t206 + t51 * t196 + t52 * t198 + (-t224 * t33 + t240 * t67) * t217;
t43 = t115 * t206 + t117 * t226 + t119 * t187 + t126 * t198 - t128 * t149 + t130 * t150;
t44 = t116 * t206 + t118 * t226 + t120 * t187 + t127 * t198 - t129 * t149 + t131 * t150;
t54 = t140 * t206 + t141 * t226 + t142 * t187 - t149 * t165 + t150 * t166 + t164 * t198;
t81 = t126 * t206 + t128 * t226 + t130 * t187;
t82 = t127 * t206 + t129 * t226 + t131 * t187;
t92 = t164 * t206 + t165 * t226 + t166 * t187;
t272 = t81 * t196 + t82 * t198 - t43 * t227 + t44 * t206 + (-t224 * t54 + t240 * t92) * t217 + t4;
t271 = 2 * m(4);
t270 = 2 * m(5);
t269 = t147 / 0.2e1;
t268 = t149 / 0.2e1;
t267 = -t225 / 0.2e1;
t266 = -t226 / 0.2e1;
t265 = t188 / 0.2e1;
t264 = t196 / 0.2e1;
t263 = t198 / 0.2e1;
t262 = t208 / 0.2e1;
t261 = t216 / 0.2e1;
t260 = -t218 / 0.2e1;
t259 = t219 / 0.2e1;
t77 = rSges(5,1) * t112 + rSges(5,2) * t111 + rSges(5,3) * t147;
t257 = pkin(3) * t148 + pkin(7) * t147 + t77;
t78 = rSges(5,1) * t114 + rSges(5,2) * t113 + rSges(5,3) * t149;
t256 = pkin(3) * t150 + pkin(7) * t149 + t78;
t99 = rSges(5,1) * t139 + rSges(5,2) * t138 + rSges(5,3) * t188;
t255 = pkin(3) * t189 + pkin(7) * t188 + t99;
t254 = Icges(3,4) * t222;
t253 = Icges(3,4) * t224;
t252 = t216 * t217;
t251 = t217 * t218;
t109 = rSges(5,1) * t152 + rSges(5,2) * t151 - rSges(5,3) * t225;
t246 = pkin(3) * t185 - pkin(7) * t225 + t109;
t110 = rSges(5,1) * t154 + rSges(5,2) * t153 - rSges(5,3) * t226;
t245 = pkin(3) * t187 - pkin(7) * t226 + t110;
t137 = -rSges(5,1) * t228 + rSges(5,2) * t190 + rSges(5,3) * t208;
t244 = pkin(3) * t209 + pkin(7) * t208 + t137;
t178 = pkin(2) * t195 + pkin(6) * t196;
t179 = -pkin(2) * t197 + pkin(6) * t198;
t243 = t178 * t252 + t179 * t251;
t181 = pkin(2) * t205 - pkin(6) * t227;
t182 = pkin(2) * t207 + pkin(6) * t206;
t242 = t181 * t252 + t182 * t251;
t241 = qJD(2) * t217;
t143 = rSges(4,1) * t189 - rSges(4,2) * t188 + rSges(4,3) * t234;
t203 = (pkin(2) * t224 + pkin(6) * t222) * t241;
t232 = (-t143 - t203) * t217;
t167 = rSges(4,1) * t209 - rSges(4,2) * t208 - rSges(4,3) * t249;
t210 = (pkin(2) * t222 - pkin(6) * t224) * t217;
t231 = (-t167 - t210) * t217;
t230 = (-t203 - t255) * t217;
t229 = (-t210 - t244) * t217;
t202 = (rSges(3,1) * t224 - rSges(3,2) * t222) * t241;
t201 = (Icges(3,1) * t224 - t254) * t241;
t200 = (-Icges(3,2) * t222 + t253) * t241;
t199 = (Icges(3,5) * t224 - Icges(3,6) * t222) * t241;
t194 = t219 * rSges(3,3) + (rSges(3,1) * t222 + rSges(3,2) * t224) * t217;
t193 = Icges(3,5) * t219 + (Icges(3,1) * t222 + t253) * t217;
t192 = Icges(3,6) * t219 + (Icges(3,2) * t224 + t254) * t217;
t180 = t219 * t182;
t177 = -rSges(3,1) * t197 - rSges(3,2) * t198;
t176 = rSges(3,1) * t195 - rSges(3,2) * t196;
t175 = -Icges(3,1) * t197 - Icges(3,4) * t198;
t174 = Icges(3,1) * t195 - Icges(3,4) * t196;
t173 = -Icges(3,4) * t197 - Icges(3,2) * t198;
t172 = Icges(3,4) * t195 - Icges(3,2) * t196;
t171 = -Icges(3,5) * t197 - Icges(3,6) * t198;
t170 = Icges(3,5) * t195 - Icges(3,6) * t196;
t163 = rSges(3,1) * t207 - rSges(3,2) * t206 + rSges(3,3) * t252;
t162 = rSges(3,1) * t205 + rSges(3,2) * t227 - rSges(3,3) * t251;
t161 = Icges(3,1) * t207 - Icges(3,4) * t206 + Icges(3,5) * t252;
t160 = Icges(3,1) * t205 + Icges(3,4) * t227 - Icges(3,5) * t251;
t159 = Icges(3,4) * t207 - Icges(3,2) * t206 + Icges(3,6) * t252;
t158 = Icges(3,4) * t205 + Icges(3,2) * t227 - Icges(3,6) * t251;
t157 = t219 * t179;
t133 = rSges(4,1) * t187 + rSges(4,2) * t226 + rSges(4,3) * t206;
t132 = rSges(4,1) * t185 + rSges(4,2) * t225 - rSges(4,3) * t227;
t125 = (t176 * t216 + t177 * t218) * t217;
t122 = rSges(4,1) * t150 - rSges(4,2) * t149 + rSges(4,3) * t198;
t121 = rSges(4,1) * t148 - rSges(4,2) * t147 + rSges(4,3) * t196;
t102 = -t133 * t249 - t167 * t206;
t101 = t132 * t249 - t167 * t227;
t100 = -t164 * t249 - t165 * t208 + t166 * t209;
t95 = t132 * t206 + t133 * t227;
t94 = (-t132 - t181) * t219 + t218 * t231;
t93 = t133 * t219 + t216 * t231 + t180;
t90 = (t132 * t216 + t133 * t218) * t217 + t242;
t89 = t110 * t208 + t137 * t226;
t88 = -t109 * t208 - t137 * t225;
t87 = -t127 * t249 - t129 * t208 + t131 * t209;
t86 = -t126 * t249 - t128 * t208 + t130 * t209;
t85 = (-t121 - t178) * t219 + t218 * t232;
t84 = t122 * t219 + t216 * t232 + t157;
t83 = t134 * t208 + t135 * t190 - t136 * t228;
t70 = -t109 * t226 + t110 * t225;
t69 = -t206 * t244 - t245 * t249;
t68 = -t227 * t244 + t246 * t249;
t65 = (t121 * t216 + t122 * t218) * t217 + t243;
t64 = -t206 * t143 - t198 * t167 + (-t122 * t224 + t133 * t240) * t217;
t63 = -t227 * t143 + t196 * t167 + (t121 * t224 - t132 * t240) * t217;
t62 = (-t181 - t246) * t219 + t218 * t229;
t61 = t216 * t229 + t219 * t245 + t180;
t60 = -t208 * t141 + t209 * t142 - t188 * t165 + t189 * t166 + (-t140 * t224 + t164 * t240) * t217;
t59 = t206 * t246 + t227 * t245;
t58 = t104 * t208 + t106 * t190 - t108 * t228;
t57 = t103 * t208 + t105 * t190 - t107 * t228;
t56 = t121 * t206 + t122 * t227 + t132 * t198 - t133 * t196;
t55 = (t216 * t246 + t218 * t245) * t217 + t242;
t48 = -t208 * t118 + t209 * t120 - t188 * t129 + t189 * t131 + (-t116 * t224 + t127 * t240) * t217;
t47 = -t208 * t117 + t209 * t119 - t188 * t128 + t189 * t130 + (-t115 * t224 + t126 * t240) * t217;
t46 = (-t178 - t257) * t219 + t218 * t230;
t45 = t216 * t230 + t219 * t256 + t157;
t40 = t110 * t188 - t137 * t149 + t208 * t78 + t226 * t99;
t39 = -t109 * t188 + t137 * t147 - t208 * t77 - t225 * t99;
t38 = (t216 * t257 + t218 * t256) * t217 + t243;
t37 = t134 * t188 + t135 * t138 + t136 * t139 + t190 * t97 + t208 * t96 - t228 * t98;
t36 = -t255 * t206 - t244 * t198 + (-t224 * t256 + t240 * t245) * t217;
t35 = -t255 * t227 + t244 * t196 + (t224 * t257 - t240 * t246) * t217;
t34 = t109 * t149 - t110 * t147 + t225 * t78 - t226 * t77;
t31 = t219 * t83 + (t216 * t58 - t218 * t57) * t217;
t30 = t206 * t58 - t227 * t57 - t249 * t83;
t29 = t208 * t83 - t225 * t57 - t226 * t58;
t28 = -t196 * t245 + t198 * t246 + t206 * t257 + t227 * t256;
t27 = t219 * t67 + (t216 * t52 - t218 * t51) * t217;
t26 = t219 * t66 + (t216 * t50 - t218 * t49) * t217;
t25 = t206 * t52 - t227 * t51 - t249 * t67;
t24 = t206 * t50 - t227 * t49 - t249 * t66;
t23 = t208 * t67 - t225 * t51 - t226 * t52;
t22 = t208 * t66 - t225 * t49 - t226 * t50;
t21 = t104 * t188 + t106 * t138 + t108 * t139 + t190 * t74 + t208 * t72 - t228 * t76;
t20 = t103 * t188 + t105 * t138 + t107 * t139 + t190 * t73 + t208 * t71 - t228 * t75;
t15 = t219 * t60 + (t216 * t48 - t218 * t47) * t217;
t14 = t219 * t54 + (t216 * t44 - t218 * t43) * t217;
t13 = t219 * t53 + (t216 * t42 - t218 * t41) * t217;
t12 = t86 * t196 + t87 * t198 - t47 * t227 + t48 * t206 + (t100 * t240 - t224 * t60) * t217;
t9 = t219 * t37 + (-t20 * t218 + t21 * t216) * t217;
t8 = t219 * t33 + (-t18 * t218 + t19 * t216) * t217;
t7 = t219 * t32 + (-t16 * t218 + t17 * t216) * t217;
t6 = t57 * t196 + t58 * t198 - t20 * t227 + t21 * t206 + (-t224 * t37 + t240 * t83) * t217;
t5 = t147 * t57 + t149 * t58 + t188 * t83 - t20 * t225 + t208 * t37 - t21 * t226;
t2 = t147 * t51 + t149 * t52 - t18 * t225 + t188 * t67 - t19 * t226 + t208 * t33;
t1 = t147 * t49 + t149 * t50 - t16 * t225 - t17 * t226 + t188 * t66 + t208 * t32;
t10 = [0; m(3) * t125 + m(4) * t65 + m(5) * t38; t8 * t252 - t7 * t251 - t13 * t251 + t14 * t252 + ((-t159 * t198 - t161 * t197 + t171 * t252 - t173 * t206 + t175 * t207) * t252 - (-t158 * t198 - t160 * t197 + t170 * t252 - t172 * t206 + t174 * t207) * t251 + (-t198 * t192 - t197 * t193 + t199 * t252 - t206 * t200 + t207 * t201) * t219) * t252 - ((-t159 * t196 + t161 * t195 - t171 * t251 + t173 * t227 + t175 * t205) * t252 - (-t158 * t196 + t160 * t195 - t170 * t251 + t172 * t227 + t174 * t205) * t251 + (-t196 * t192 + t195 * t193 - t199 * t251 + t200 * t227 + t205 * t201) * t219) * t251 + t219 * t9 + (t38 * t55 + t45 * t61 + t46 * t62) * t270 + t219 * t15 + (t65 * t90 + t84 * t93 + t85 * t94) * t271 + 0.2e1 * m(3) * ((-t162 * t219 - t194 * t251) * (-t176 * t219 - t202 * t251) + (t163 * t219 - t194 * t252) * (t177 * t219 - t202 * t252) + (t162 * t216 + t163 * t218) * t217 * t125) + t219 * (t219 ^ 2 * t199 + (((t173 * t224 + t175 * t222) * t216 - (t172 * t224 + t174 * t222) * t218 + ((-t159 * t222 + t161 * t224) * t216 - (-t158 * t222 + t160 * t224) * t218) * qJD(2)) * t217 + (-t170 * t218 + t171 * t216 + t200 * t224 + t201 * t222 + (-t192 * t222 + t193 * t224) * qJD(2)) * t219) * t217); m(4) * t56 + m(5) * t28; t26 * t264 + t27 * t263 + (t8 / 0.2e1 + t14 / 0.2e1) * t206 - (t7 / 0.2e1 + t13 / 0.2e1) * t227 + m(5) * (t28 * t55 + t35 * t62 + t36 * t61 + t38 * t59 + t45 * t69 + t46 * t68) + m(4) * (t101 * t85 + t102 * t84 + t56 * t90 + t63 * t94 + t64 * t93 + t65 * t95) + (t6 / 0.2e1 + t12 / 0.2e1 + t92 * t263 + t91 * t264) * t219 + ((t216 * t82 - t218 * t81) * t263 + (t216 * t80 - t218 * t79) * t264 + (-t9 / 0.2e1 - t15 / 0.2e1) * t224 + (t31 / 0.2e1 + t100 * t259 + (t216 * t87 - t218 * t86) * t217 / 0.2e1) * t240 + t272 * t261 + t273 * t260) * t217; (-t12 - t6) * t249 + t272 * t206 - t273 * t227 + (t206 * t82 - t227 * t81 - t249 * t92 + t25) * t198 + (t206 * t80 - t227 * t79 - t249 * t91 + t24) * t196 + (-t100 * t249 + t206 * t87 - t227 * t86 + t30) * t234 + (t28 * t59 + t35 * t68 + t36 * t69) * t270 + (t101 * t63 + t102 * t64 + t56 * t95) * t271; m(5) * t34; t31 * t265 + t9 * t262 + t5 * t259 + t26 * t269 + t7 * t267 + t27 * t268 + t8 * t266 + m(5) * (t34 * t55 + t38 * t70 + t39 * t62 + t40 * t61 + t45 * t89 + t46 * t88) + (t1 * t260 + t2 * t261) * t217; m(5) * (t28 * t70 + t34 * t59 + t35 * t88 + t36 * t89 + t39 * t68 + t40 * t69) + t30 * t265 + t6 * t262 + t25 * t268 + t4 * t266 + t23 * t263 + t206 * t2 / 0.2e1 + t22 * t264 - t227 * t1 / 0.2e1 + t24 * t269 + t3 * t267 + (t29 * t240 / 0.2e1 - t224 * t5 / 0.2e1) * t217; (t34 * t70 + t39 * t88 + t40 * t89) * t270 + t149 * t23 - t226 * t2 + t147 * t22 - t225 * t1 + t188 * t29 + t208 * t5;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
