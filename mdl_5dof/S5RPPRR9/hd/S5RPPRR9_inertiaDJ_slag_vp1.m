% Calculate time derivative of joint inertia matrix for
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:28
% DurationCPUTime: 4.81s
% Computational Cost: add. (9269->436), mult. (22757->648), div. (0->0), fcn. (25735->8), ass. (0->210)
t248 = sin(pkin(8));
t249 = cos(pkin(8));
t261 = sin(qJ(1));
t262 = cos(qJ(1));
t135 = -t248 * t261 - t249 * t262;
t130 = t135 * qJD(1);
t136 = t262 * t248 - t261 * t249;
t131 = t136 * qJD(1);
t163 = cos(qJ(4));
t161 = sin(qJ(4));
t227 = qJD(4) * t161;
t180 = t131 * t163 + t135 * t227;
t178 = rSges(5,1) * t180 + t130 * rSges(5,3);
t226 = qJD(4) * t163;
t240 = t131 * t161;
t181 = t135 * t226 - t240;
t182 = -t130 * t163 + t136 * t227;
t242 = t130 * t161;
t183 = t136 * t226 + t242;
t254 = t131 * rSges(5,3);
t205 = -rSges(5,1) * t163 + rSges(5,2) * t161;
t253 = t135 * rSges(5,3);
t88 = t136 * t205 - t253;
t237 = t135 * t163;
t238 = t135 * t161;
t89 = -rSges(5,1) * t237 + rSges(5,2) * t238 + rSges(5,3) * t136;
t27 = t130 * t88 + t136 * (rSges(5,1) * t182 + t254) - t131 * t89 + t135 * t178 + (t135 * t181 + t136 * t183) * rSges(5,2);
t270 = 2 * m(5);
t280 = t27 * t270;
t147 = -rSges(5,1) * t161 - rSges(5,2) * t163;
t279 = t147 ^ 2 * t270;
t229 = pkin(1) * t262 + qJ(2) * t261;
t219 = t262 * pkin(2) + t229;
t278 = -rSges(3,1) * t262 - rSges(3,3) * t261 - t229;
t246 = Icges(5,4) * t163;
t190 = Icges(5,2) * t161 - t246;
t84 = -Icges(5,6) * t135 + t136 * t190;
t247 = Icges(5,4) * t161;
t192 = -Icges(5,1) * t163 + t247;
t86 = -Icges(5,5) * t135 + t136 * t192;
t196 = t161 * t84 - t163 * t86;
t85 = Icges(5,6) * t136 + t135 * t190;
t87 = Icges(5,5) * t136 + t135 * t192;
t195 = t161 * t85 - t163 * t87;
t273 = t196 * t135;
t274 = t195 * t136;
t188 = -Icges(5,5) * t163 + Icges(5,6) * t161;
t82 = -Icges(5,3) * t135 + t136 * t188;
t83 = Icges(5,3) * t136 + t135 * t188;
t277 = -t135 * t83 + t136 * t82 + t273 + t274;
t160 = sin(qJ(5));
t162 = cos(qJ(5));
t169 = qJD(5) * t136 + t180;
t224 = qJD(5) * t163;
t207 = t135 * t224 + t130;
t59 = -t160 * t169 + t162 * t207;
t60 = t160 * t207 + t162 * t169;
t39 = t60 * rSges(6,1) + t59 * rSges(6,2) - rSges(6,3) * t181;
t276 = pkin(4) * t180 - pkin(7) * t181 + t39;
t258 = pkin(4) * t163;
t209 = pkin(7) * t161 + t258;
t275 = t136 * t209;
t141 = t205 * qJD(4);
t272 = t141 * t147 * t270;
t263 = rSges(6,3) + pkin(7);
t271 = t161 * t263 + pkin(3) + t258;
t269 = 2 * m(6);
t268 = t131 / 0.2e1;
t267 = -t135 / 0.2e1;
t266 = -t136 / 0.2e1;
t265 = t136 / 0.2e1;
t264 = t163 / 0.2e1;
t260 = m(5) * t141;
t259 = m(5) * t147;
t257 = t161 * pkin(4);
t234 = t160 * t163;
t100 = -t135 * t162 + t136 * t234;
t232 = t162 * t163;
t101 = -t135 * t160 - t136 * t232;
t236 = t136 * t161;
t69 = Icges(6,4) * t101 + Icges(6,2) * t100 - Icges(6,6) * t236;
t71 = Icges(6,1) * t101 + Icges(6,4) * t100 - Icges(6,5) * t236;
t198 = t160 * t69 - t162 * t71;
t170 = -qJD(5) * t135 + t182;
t206 = t136 * t224 + t131;
t57 = -t160 * t170 + t162 * t206;
t58 = t160 * t206 + t162 * t170;
t32 = Icges(6,5) * t58 + Icges(6,6) * t57 - Icges(6,3) * t183;
t34 = Icges(6,4) * t58 + Icges(6,2) * t57 - Icges(6,6) * t183;
t36 = Icges(6,1) * t58 + Icges(6,4) * t57 - Icges(6,5) * t183;
t67 = Icges(6,5) * t101 + Icges(6,6) * t100 - Icges(6,3) * t236;
t9 = (qJD(4) * t198 + t32) * t163 + (-qJD(4) * t67 + t160 * t34 - t162 * t36 + (t160 * t71 + t162 * t69) * qJD(5)) * t161;
t256 = t9 * t135;
t102 = t135 * t234 + t136 * t162;
t103 = -t135 * t232 + t136 * t160;
t70 = Icges(6,4) * t103 + Icges(6,2) * t102 - Icges(6,6) * t238;
t72 = Icges(6,1) * t103 + Icges(6,4) * t102 - Icges(6,5) * t238;
t197 = t160 * t70 - t162 * t72;
t33 = Icges(6,5) * t60 + Icges(6,6) * t59 - Icges(6,3) * t181;
t35 = Icges(6,4) * t60 + Icges(6,2) * t59 - Icges(6,6) * t181;
t37 = Icges(6,1) * t60 + Icges(6,4) * t59 - Icges(6,5) * t181;
t68 = Icges(6,5) * t103 + Icges(6,6) * t102 - Icges(6,3) * t238;
t10 = (qJD(4) * t197 + t33) * t163 + (-qJD(4) * t68 + t160 * t35 - t162 * t37 + (t160 * t72 + t162 * t70) * qJD(5)) * t161;
t255 = t10 * t136;
t204 = -t101 * rSges(6,1) - t100 * rSges(6,2);
t73 = -rSges(6,3) * t236 - t204;
t252 = t73 - t275;
t74 = rSges(6,1) * t103 + rSges(6,2) * t102 - rSges(6,3) * t238;
t251 = -pkin(4) * t237 - pkin(7) * t238 + t74;
t203 = -rSges(6,1) * t162 + rSges(6,2) * t160;
t225 = qJD(5) * t161;
t97 = (rSges(6,1) * t160 + rSges(6,2) * t162) * t225 + (-rSges(6,3) * t161 + t163 * t203) * qJD(4);
t250 = qJD(4) * t209 - t97;
t245 = Icges(6,4) * t160;
t244 = Icges(6,4) * t162;
t189 = Icges(6,2) * t160 - t244;
t116 = Icges(6,6) * t163 + t161 * t189;
t243 = t116 * t160;
t187 = -Icges(6,5) * t162 + Icges(6,6) * t160;
t115 = Icges(6,3) * t163 + t161 * t187;
t233 = t161 * t115;
t118 = rSges(6,3) * t163 + t161 * t203;
t208 = -pkin(7) * t163 + t257;
t231 = t118 - t208;
t157 = t262 * qJ(2);
t230 = qJD(1) * t157 + qJD(2) * t261;
t228 = qJD(4) * t136;
t223 = t261 * pkin(1);
t29 = t161 * t198 + t163 * t67;
t191 = -Icges(6,1) * t162 + t245;
t117 = Icges(6,5) * t163 + t161 * t191;
t49 = t100 * t116 + t101 * t117 - t136 * t233;
t222 = t29 / 0.2e1 + t49 / 0.2e1;
t30 = t161 * t197 + t163 * t68;
t50 = t102 * t116 + t103 * t117 - t135 * t233;
t221 = t30 / 0.2e1 + t50 / 0.2e1;
t220 = t117 * t232;
t210 = rSges(6,1) * t58 + rSges(6,2) * t57;
t23 = t100 * t69 + t101 * t71 - t236 * t67;
t24 = t100 * t70 + t101 * t72 - t236 * t68;
t202 = -t135 * t24 - t136 * t23;
t25 = t102 * t69 + t103 * t71 - t238 * t67;
t26 = t102 * t70 + t103 * t72 - t238 * t68;
t201 = -t135 * t26 - t136 * t25;
t200 = -t135 * t30 - t136 * t29;
t199 = -t135 * t73 + t136 * t74;
t186 = pkin(3) - t205;
t185 = -t161 * t32 - t226 * t67;
t184 = -t161 * t33 - t226 * t68;
t179 = -pkin(3) * t135 + pkin(6) * t136 + t219;
t177 = -pkin(2) * t261 - t223;
t173 = t157 + t177;
t94 = (Icges(6,5) * t160 + Icges(6,6) * t162) * t225 + (-Icges(6,3) * t161 + t163 * t187) * qJD(4);
t95 = (Icges(6,2) * t162 + t245) * t225 + (-Icges(6,6) * t161 + t163 * t189) * qJD(4);
t96 = (Icges(6,1) * t160 + t244) * t225 + (-Icges(6,5) * t161 + t163 * t191) * qJD(4);
t172 = t163 * t94 + t226 * t243 + (t116 * t225 - t161 * t96) * t162 + (t117 * t225 + t161 * t95) * t160;
t171 = pkin(6) * t135 + t173;
t168 = -rSges(3,1) * t261 + rSges(3,3) * t262 - t223;
t167 = qJD(1) * t177 + t230;
t155 = qJD(2) * t262;
t166 = -qJD(1) * t219 + t155;
t165 = pkin(3) * t131 + pkin(6) * t130 + t167;
t164 = -pkin(6) * t131 + t166;
t121 = t157 + t168;
t107 = qJD(1) * t278 + t155;
t106 = qJD(1) * t168 + t230;
t93 = -rSges(4,1) * t135 - rSges(4,2) * t136 + t219;
t92 = rSges(4,1) * t136 - rSges(4,2) * t135 + t173;
t81 = t231 * t135;
t80 = t231 * t136;
t79 = rSges(4,1) * t130 + rSges(4,2) * t131 + t166;
t78 = rSges(4,1) * t131 - rSges(4,2) * t130 + t167;
t77 = t115 * t163 + (-t117 * t162 + t243) * t161;
t76 = t179 + t89;
t75 = t136 * t186 + t171 + t253;
t62 = Icges(5,5) * t180 + Icges(5,6) * t181 + Icges(5,3) * t130;
t61 = Icges(5,5) * t182 + Icges(5,6) * t183 + Icges(5,3) * t131;
t54 = t118 * t238 + t163 * t74;
t53 = -t118 * t236 - t163 * t73;
t52 = t131 * t231 + t135 * t250;
t51 = -t130 * t231 + t136 * t250;
t48 = t130 * t186 + t147 * t228 + t164 - t254;
t47 = rSges(5,2) * t181 + t165 + t178;
t46 = t179 + t251;
t45 = t136 * t271 + t171 + t204;
t40 = t199 * t161;
t38 = -rSges(6,3) * t183 + t210;
t31 = ((-t220 - t233) * qJD(4) + t172) * t163;
t28 = t135 * t251 + t136 * t252;
t22 = (t163 * t263 - t257) * t228 + t271 * t130 + t164 - t210;
t21 = t165 + t276;
t20 = (qJD(4) * t118 * t135 + t39) * t163 + (-qJD(4) * t74 - t118 * t131 + t135 * t97) * t161;
t19 = (-t118 * t228 - t38) * t163 + (qJD(4) * t73 - t118 * t130 - t136 * t97) * t161;
t18 = t102 * t95 + t103 * t96 - t115 * t181 + t59 * t116 + t60 * t117 - t238 * t94;
t17 = t100 * t95 + t101 * t96 - t115 * t183 + t57 * t116 + t58 * t117 - t236 * t94;
t16 = -t135 * t25 + t136 * t26;
t15 = -t135 * t23 + t136 * t24;
t14 = t161 * t201 + t163 * t50;
t13 = t161 * t202 + t163 * t49;
t12 = t199 * t226 + (t130 * t74 + t131 * t73 - t135 * t38 + t136 * t39) * t161;
t11 = t276 * t135 - t251 * t131 + (t208 * t228 + t38) * t136 - (-t252 + t275) * t130;
t8 = t102 * t35 + t103 * t37 + t135 * t184 + t240 * t68 + t59 * t70 + t60 * t72;
t7 = t102 * t34 + t103 * t36 + t135 * t185 + t240 * t67 + t59 * t69 + t60 * t71;
t6 = t100 * t35 + t101 * t37 + t136 * t184 - t242 * t68 + t57 * t70 + t58 * t72;
t5 = t100 * t34 + t101 * t36 + t136 * t185 - t242 * t67 + t57 * t69 + t58 * t71;
t4 = t130 * t26 + t131 * t25 - t135 * t7 + t136 * t8;
t3 = t130 * t24 + t131 * t23 - t135 * t5 + t136 * t6;
t2 = (qJD(4) * t201 + t18) * t163 + (-qJD(4) * t50 - t130 * t25 + t131 * t26 - t135 * t8 - t136 * t7) * t161;
t1 = (qJD(4) * t202 + t17) * t163 + (-qJD(4) * t49 - t130 * t23 + t131 * t24 - t135 * t6 - t136 * t5) * t161;
t41 = [(t21 * t46 + t22 * t45) * t269 - qJD(4) * t220 + (t47 * t76 + t48 * t75) * t270 + 0.2e1 * m(4) * (t78 * t93 + t79 * t92) + 0.2e1 * m(3) * (-t106 * t278 + t107 * t121) + t172 + (Icges(5,1) * t161 - t190 + t246) * t226 + (-Icges(5,2) * t163 - t115 - t192 - t247) * t227; m(6) * (t261 * t22 - t262 * t21 + (t261 * t46 + t262 * t45) * qJD(1)) + m(5) * (t261 * t48 - t262 * t47 + (t261 * t76 + t262 * t75) * qJD(1)) + m(4) * (t261 * t79 - t262 * t78 + (t261 * t93 + t262 * t92) * qJD(1)) + m(3) * (t261 * t107 - t262 * t106 + (t121 * t262 - t261 * t278) * qJD(1)); 0; 0; 0; 0; -t256 / 0.2e1 + t255 / 0.2e1 + (t135 ^ 2 / 0.2e1 + t136 ^ 2 / 0.2e1) * t188 * qJD(4) + (-t161 * t86 / 0.2e1 - t163 * t84 / 0.2e1 + t222) * t131 - (t161 * t87 / 0.2e1 + t85 * t264 - t221) * t130 + m(6) * (-t21 * t80 - t22 * t81 + t45 * t52 + t46 * t51) + (-t135 * t75 - t136 * t76) * t260 + (-t130 * t76 + t131 * t75 - t135 * t48 - t136 * t47) * t259 + (t130 * t136 - t131 * t135) * (-Icges(5,5) * t161 - Icges(5,6) * t163) + (qJD(4) * t196 - t161 * (Icges(5,1) * t182 + Icges(5,4) * t183 + Icges(5,5) * t131) - t163 * (Icges(5,4) * t182 + Icges(5,2) * t183 + Icges(5,6) * t131) + t17) * t267 + (qJD(4) * t195 - t161 * (Icges(5,1) * t180 + Icges(5,4) * t181 + Icges(5,5) * t130) - t163 * (Icges(5,4) * t180 + Icges(5,2) * t181 + Icges(5,6) * t130) + t18) * t265; m(6) * (t52 * t261 - t51 * t262 + (-t261 * t80 - t262 * t81) * qJD(1)) + (-t135 * t261 + t136 * t262) * t260 + (t261 * t131 + t262 * t130 + (-t135 * t262 - t136 * t261) * qJD(1)) * t259; -m(5) * t27 - m(6) * t11; (t11 * t28 - t51 * t80 - t52 * t81) * t269 + t130 * t16 + t131 * t15 + (t88 * t280 + t4 + (t136 * t62 + t272) * t136 + (-t274 + t277) * t131 + (0.2e1 * t135 * t195 + 0.3e1 * t136 * t83 + t279) * t130) * t136 + (t89 * t280 - t3 + (-t135 * t61 + t272) * t135 + (t135 * t62 - t136 * t61) * t136 + (t273 - t277 + (-t195 - t82) * t136) * t130 + (0.3e1 * t135 * t82 - t279 + (-t196 - t83) * t136) * t131) * t135; m(6) * (t19 * t45 + t20 * t46 + t21 * t54 + t22 * t53) + t31 + (-t135 * t221 - t136 * t222) * t226 + (-t77 * qJD(4) + (-t9 / 0.2e1 - t17 / 0.2e1) * t136 + (-t10 / 0.2e1 - t18 / 0.2e1) * t135 + t221 * t131 - t222 * t130) * t161; m(6) * (t19 * t261 - t20 * t262 + (t261 * t54 + t262 * t53) * qJD(1)); -m(6) * t12; m(6) * (t11 * t40 + t12 * t28 - t19 * t81 - t20 * t80 + t51 * t54 + t52 * t53) + t130 * t14 / 0.2e1 + t2 * t265 + t13 * t268 + t1 * t267 + (t30 * t130 + t29 * t131 + t255 - t256) * t264 + (t15 * t266 + t16 * t267) * t226 + (t16 * t268 + t4 * t267 - t130 * t15 / 0.2e1 + t3 * t266 - qJD(4) * (-t135 * t29 + t136 * t30) / 0.2e1) * t161; (t12 * t40 + t19 * t53 + t20 * t54) * t269 + (t31 + (-t13 * t136 - t135 * t14 + t163 * t200) * qJD(4)) * t163 + (t131 * t14 - t135 * t2 - t130 * t13 - t136 * t1 + t163 * (-t10 * t135 - t29 * t130 + t30 * t131 - t9 * t136) + (-t161 * t200 - 0.2e1 * t77 * t163) * qJD(4)) * t161;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t41(1), t41(2), t41(4), t41(7), t41(11); t41(2), t41(3), t41(5), t41(8), t41(12); t41(4), t41(5), t41(6), t41(9), t41(13); t41(7), t41(8), t41(9), t41(10), t41(14); t41(11), t41(12), t41(13), t41(14), t41(15);];
Mq = res;
