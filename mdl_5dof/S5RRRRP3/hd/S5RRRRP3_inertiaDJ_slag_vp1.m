% Calculate time derivative of joint inertia matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:16
% DurationCPUTime: 4.72s
% Computational Cost: add. (8370->333), mult. (6088->456), div. (0->0), fcn. (4408->8), ass. (0->186)
t172 = cos(qJ(4));
t170 = sin(qJ(4));
t245 = Icges(6,5) * t170;
t247 = Icges(5,4) * t170;
t285 = t245 - t247 + (-Icges(5,2) - Icges(6,3)) * t172;
t244 = Icges(6,5) * t172;
t246 = Icges(5,4) * t172;
t284 = -t244 + t246 + (Icges(5,1) + Icges(6,1)) * t170;
t195 = Icges(6,3) * t170 + t244;
t198 = -Icges(5,2) * t170 + t246;
t199 = Icges(6,1) * t172 + t245;
t200 = Icges(5,1) * t172 - t247;
t283 = (t199 + t200) * t170 - (t195 - t198) * t172;
t272 = rSges(6,3) + qJ(5);
t278 = rSges(6,1) + pkin(4);
t274 = t170 * t272 + t172 * t278;
t275 = t285 * t170 + t284 * t172;
t169 = qJ(1) + qJ(2);
t166 = qJ(3) + t169;
t161 = sin(t166);
t229 = qJD(4) * t170;
t220 = t161 * t229;
t227 = qJD(5) * t170;
t228 = qJD(4) * t172;
t280 = -(t228 * t272 + t227) * t161 + t278 * t220;
t162 = cos(t166);
t168 = qJD(1) + qJD(2);
t163 = qJD(3) + t168;
t241 = t162 * t163;
t140 = Icges(5,5) * t170 + Icges(5,6) * t172;
t141 = Icges(6,4) * t170 - Icges(6,6) * t172;
t277 = t140 + t141;
t196 = Icges(5,5) * t172 - Icges(5,6) * t170;
t197 = Icges(6,4) * t172 + Icges(6,6) * t170;
t273 = t275 * t163 + (-t196 - t197) * qJD(4);
t190 = t198 * t162;
t78 = Icges(5,6) * t161 + t190;
t192 = t200 * t162;
t82 = Icges(5,5) * t161 + t192;
t201 = t170 * t78 - t172 * t82;
t271 = t161 * t201;
t187 = t195 * t162;
t72 = Icges(6,6) * t161 + t187;
t191 = t199 * t162;
t80 = Icges(6,4) * t161 + t191;
t205 = t170 * t72 + t172 * t80;
t270 = t161 * t205;
t77 = -Icges(5,6) * t162 + t161 * t198;
t81 = -Icges(5,5) * t162 + t161 * t200;
t203 = t170 * t77 - t172 * t81;
t269 = t162 * t203;
t71 = -Icges(6,6) * t162 + t161 * t195;
t79 = -Icges(6,4) * t162 + t161 * t199;
t207 = t170 * t71 + t172 * t79;
t268 = t162 * t207;
t265 = 2 * m(3);
t264 = 2 * m(4);
t263 = 2 * m(5);
t262 = 2 * m(6);
t251 = rSges(5,2) * t170;
t252 = rSges(5,1) * t172;
t124 = (-t251 + t252) * qJD(4);
t258 = m(5) * t124;
t148 = rSges(5,1) * t170 + rSges(5,2) * t172;
t257 = m(5) * t148;
t171 = sin(qJ(1));
t256 = pkin(1) * t171;
t164 = sin(t169);
t255 = pkin(2) * t164;
t152 = t162 * rSges(6,2);
t254 = t161 * t274 - t152;
t150 = t161 * rSges(6,2);
t239 = t162 * t172;
t240 = t162 * t170;
t253 = t278 * t239 + t272 * t240 + t150;
t250 = pkin(1) * qJD(1);
t149 = t161 * rSges(5,3);
t232 = t278 * t170 - t272 * t172;
t70 = t232 * t162;
t249 = t163 * t70;
t243 = t161 * t163;
t242 = t161 * t172;
t238 = t164 * t168;
t165 = cos(t169);
t237 = t165 * t168;
t236 = -qJD(4) * t274 + qJD(5) * t172;
t134 = t161 * t251;
t235 = rSges(5,3) * t241 + t163 * t134;
t233 = t162 * rSges(5,3) + t134;
t231 = t162 * pkin(3) + t161 * pkin(8);
t230 = t161 ^ 2 + t162 ^ 2;
t226 = pkin(2) * t238;
t225 = pkin(2) * t237;
t224 = rSges(5,2) * t240;
t223 = t171 * t250;
t173 = cos(qJ(1));
t222 = t173 * t250;
t221 = -t161 * rSges(5,2) * t228 - rSges(5,1) * t220 - t163 * t224;
t219 = t162 * t229;
t218 = t162 * t228;
t215 = -pkin(3) - t252;
t113 = t165 * rSges(3,1) - rSges(3,2) * t164;
t106 = t162 * rSges(4,1) - rSges(4,2) * t161;
t98 = -rSges(3,1) * t237 + rSges(3,2) * t238;
t92 = -rSges(4,1) * t241 + rSges(4,2) * t243;
t160 = pkin(2) * t165;
t94 = t106 + t160;
t112 = -rSges(3,1) * t164 - rSges(3,2) * t165;
t105 = -rSges(4,1) * t161 - rSges(4,2) * t162;
t208 = t170 * t79 - t172 * t71;
t206 = t170 * t80 - t172 * t72;
t204 = t170 * t81 + t172 * t77;
t202 = t170 * t82 + t172 * t78;
t86 = rSges(5,1) * t239 + t149 - t224;
t97 = t112 * t168;
t91 = t105 * t163;
t189 = t197 * t162;
t188 = t196 * t162;
t58 = t231 + t253;
t93 = t105 - t255;
t64 = t86 + t231;
t68 = t92 - t225;
t40 = t160 + t58;
t155 = t162 * pkin(8);
t63 = t161 * t215 + t155 + t233;
t186 = -pkin(3) - t274;
t62 = t160 + t64;
t185 = rSges(6,2) * t241 + t162 * t227 + t272 * t218 - t219 * t278;
t184 = t186 * t161;
t181 = Icges(6,2) * t163 - qJD(4) * t141;
t178 = Icges(5,3) * t163 - qJD(4) * t140;
t67 = t91 - t226;
t61 = t63 - t255;
t177 = t283 * qJD(4) + t284 * t228 + t285 * t229;
t57 = t152 + t155 + t184;
t176 = (-t273 * t161 + (-t201 + t205) * qJD(4) - t283 * t243) * t161 / 0.2e1 - (t273 * t162 + (-t203 + t207) * qJD(4)) * t162 / 0.2e1 + (t275 * t161 - t277 * t162 + t204 + t208) * t243 / 0.2e1 + (t277 * t161 + t275 * t162 + t202 + t206) * t241 / 0.2e1 - ((-t187 + t190) * t172 + (t191 + t192) * t170) * t241 / 0.2e1;
t30 = (t215 * t162 + (-rSges(5,3) - pkin(8)) * t161) * t163 - t221;
t39 = t57 - t255;
t28 = t30 - t225;
t125 = pkin(8) * t241;
t29 = -rSges(5,2) * t218 - pkin(3) * t243 + t125 + (-t163 * t242 - t219) * rSges(5,1) + t235;
t14 = t163 * t184 + t125 + t185;
t27 = t29 - t226;
t12 = t14 - t226;
t15 = ((-rSges(6,2) - pkin(8)) * t161 + t186 * t162) * t163 + t280;
t13 = t15 - t225;
t167 = t173 * pkin(1);
t100 = t113 + t167;
t99 = t112 - t256;
t90 = t98 - t222;
t89 = t97 - t223;
t88 = t167 + t94;
t87 = t93 - t256;
t84 = rSges(5,1) * t242 - t233;
t76 = Icges(6,2) * t161 + t189;
t75 = -Icges(6,2) * t162 + t161 * t197;
t74 = Icges(5,3) * t161 + t188;
t73 = -Icges(5,3) * t162 + t161 * t196;
t69 = t232 * t161;
t66 = t68 - t222;
t65 = t67 - t223;
t60 = t167 + t62;
t59 = t61 - t256;
t50 = t161 * t181 + t163 * t189;
t49 = t162 * t181 - t197 * t243;
t48 = t161 * t178 + t163 * t188;
t47 = t162 * t178 - t196 * t243;
t38 = t167 + t40;
t37 = t39 - t256;
t32 = t161 * t236 - t249;
t31 = t162 * t236 + t232 * t243;
t26 = t28 - t222;
t25 = t27 - t223;
t24 = t161 * t74 - t162 * t201;
t23 = t161 * t73 - t269;
t22 = t161 * t76 + t162 * t205;
t21 = t161 * t75 + t268;
t20 = -t162 * t74 - t271;
t19 = -t161 * t203 - t162 * t73;
t18 = -t162 * t76 + t270;
t17 = t161 * t207 - t162 * t75;
t16 = t161 * t254 + t162 * t253;
t11 = t13 - t222;
t10 = t12 - t223;
t1 = (t163 * t254 + t185) * t162 + ((t150 - t253) * t163 - t280) * t161;
t2 = [(t10 * t38 + t11 * t37) * t262 + (t25 * t60 + t26 * t59) * t263 + (t65 * t88 + t66 * t87) * t264 + (t100 * t89 + t90 * t99) * t265 + t177; m(6) * (t10 * t40 + t11 * t39 + t12 * t38 + t13 * t37) + m(5) * (t25 * t62 + t26 * t61 + t27 * t60 + t28 * t59) + m(4) * (t65 * t94 + t66 * t93 + t67 * t88 + t68 * t87) + m(3) * (t100 * t97 + t112 * t90 + t113 * t89 + t98 * t99) + t177; (t12 * t40 + t13 * t39) * t262 + (t27 * t62 + t28 * t61) * t263 + (t67 * t94 + t68 * t93) * t264 + (t112 * t98 + t113 * t97) * t265 + t177; m(6) * (t10 * t58 + t11 * t57 + t14 * t38 + t15 * t37) + m(5) * (t25 * t64 + t26 * t63 + t29 * t60 + t30 * t59) + m(4) * (t105 * t66 + t106 * t65 + t87 * t92 + t88 * t91) + t177; m(6) * (t12 * t58 + t13 * t57 + t14 * t40 + t15 * t39) + m(5) * (t27 * t64 + t28 * t63 + t29 * t62 + t30 * t61) + m(4) * (t105 * t68 + t106 * t67 + t91 * t94 + t92 * t93) + t177; (t14 * t58 + t15 * t57) * t262 + (t29 * t64 + t30 * t63) * t263 + (t105 * t92 + t106 * t91) * t264 + t177; m(6) * (-t10 * t69 - t11 * t70 + t31 * t37 + t32 * t38) + ((-t163 * t60 - t26) * t162 + (t163 * t59 - t25) * t161) * t257 + (-t161 * t60 - t162 * t59) * t258 + t176; m(6) * (-t12 * t69 - t13 * t70 + t31 * t39 + t32 * t40) + ((-t163 * t62 - t28) * t162 + (t163 * t61 - t27) * t161) * t257 + (-t161 * t62 - t162 * t61) * t258 + t176; m(6) * (-t14 * t69 - t15 * t70 + t31 * t57 + t32 * t58) + ((-t163 * t64 - t30) * t162 + (t163 * t63 - t29) * t161) * t257 + (-t161 * t64 - t162 * t63) * t258 + t176; ((t161 * t84 + t162 * t86) * (((-t86 + t149) * t163 + t221) * t161 + (-qJD(4) * t148 * t162 + t163 * t84 + t235) * t162) + t230 * t148 * t124) * t263 + t161 * ((t161 * t47 + (t23 + t271) * t163) * t161 + (t24 * t163 + (t228 * t77 + t229 * t81) * t162 + (-t202 * qJD(4) - t163 * t203 - t48) * t161) * t162) - t162 * ((t162 * t48 + (t20 + t269) * t163) * t162 + (t19 * t163 + (-t228 * t78 - t229 * t82) * t161 + (t204 * qJD(4) - t163 * t201 - t47) * t162) * t161) + (t1 * t16 - t31 * t70 - t32 * t69) * t262 + t161 * ((t161 * t49 + (t21 - t270) * t163) * t161 + (t22 * t163 + (-t228 * t71 + t229 * t79) * t162 + (-t206 * qJD(4) + t163 * t207 - t50) * t161) * t162) - t162 * ((t162 * t50 + (t18 - t268) * t163) * t162 + (t17 * t163 + (t228 * t72 - t229 * t80) * t161 + (t208 * qJD(4) + t163 * t205 - t49) * t162) * t161) + ((-t17 - t19) * t162 + (t18 + t20) * t161) * t243 + ((-t21 - t23) * t162 + (t22 + t24) * t161) * t241; m(6) * ((t161 * t38 + t162 * t37) * t228 + ((t163 * t38 + t11) * t162 + (-t163 * t37 + t10) * t161) * t170); m(6) * ((t161 * t40 + t162 * t39) * t228 + ((t163 * t40 + t13) * t162 + (-t163 * t39 + t12) * t161) * t170); m(6) * ((t161 * t58 + t162 * t57) * t228 + ((t163 * t58 + t15) * t162 + (-t163 * t57 + t14) * t161) * t170); m(6) * ((-t1 + (-t161 * t69 - t162 * t70) * qJD(4)) * t172 + (qJD(4) * t16 + (-t163 * t69 + t31) * t162 + (t32 + t249) * t161) * t170); (-0.1e1 + t230) * t170 * t228 * t262;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
