% Calculate time derivative of joint inertia matrix for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:47
% EndTime: 2019-12-05 15:44:56
% DurationCPUTime: 3.61s
% Computational Cost: add. (10821->279), mult. (10909->466), div. (0->0), fcn. (10123->10), ass. (0->172)
t163 = qJD(2) + qJD(4);
t165 = sin(pkin(8));
t161 = t165 ^ 2;
t166 = cos(pkin(8));
t162 = t166 ^ 2;
t289 = t161 + t162;
t283 = t163 * t289;
t164 = qJ(2) + pkin(9);
t156 = sin(t164);
t157 = cos(t164);
t169 = sin(qJ(2));
t171 = cos(qJ(2));
t293 = (-Icges(3,5) * t169 - Icges(4,5) * t156 - Icges(3,6) * t171 - Icges(4,6) * t157) * qJD(2);
t158 = qJ(4) + t164;
t153 = sin(t158);
t154 = cos(t158);
t292 = Icges(5,5) * t153 + Icges(5,6) * t154;
t182 = t163 * t292;
t102 = t165 * t182;
t103 = t166 * t182;
t273 = t292 * t283;
t170 = cos(qJ(5));
t241 = t166 * t170;
t168 = sin(qJ(5));
t244 = t165 * t168;
t133 = -t154 * t244 - t241;
t242 = t166 * t168;
t243 = t165 * t170;
t134 = t154 * t243 - t242;
t247 = t154 * t163;
t246 = t163 * t165;
t232 = t154 * t246;
t250 = t153 * t163;
t234 = t168 * t250;
t83 = -qJD(5) * t134 + t165 * t234;
t233 = t170 * t250;
t84 = qJD(5) * t133 - t165 * t233;
t47 = Icges(6,5) * t84 + Icges(6,6) * t83 + Icges(6,3) * t232;
t249 = t153 * t165;
t68 = Icges(6,5) * t134 + Icges(6,6) * t133 + Icges(6,3) * t249;
t192 = t153 * t47 + t247 * t68;
t49 = Icges(6,4) * t84 + Icges(6,2) * t83 + Icges(6,6) * t232;
t51 = Icges(6,1) * t84 + Icges(6,4) * t83 + Icges(6,5) * t232;
t70 = Icges(6,4) * t134 + Icges(6,2) * t133 + Icges(6,6) * t249;
t72 = Icges(6,1) * t134 + Icges(6,4) * t133 + Icges(6,5) * t249;
t13 = t133 * t49 + t134 * t51 + t165 * t192 + t70 * t83 + t72 * t84;
t245 = t163 * t166;
t231 = t154 * t245;
t136 = t154 * t241 + t244;
t85 = -qJD(5) * t136 + t166 * t234;
t135 = -t154 * t242 + t243;
t86 = qJD(5) * t135 - t166 * t233;
t48 = Icges(6,5) * t86 + Icges(6,6) * t85 + Icges(6,3) * t231;
t248 = t153 * t166;
t69 = Icges(6,5) * t136 + Icges(6,6) * t135 + Icges(6,3) * t248;
t191 = t153 * t48 + t247 * t69;
t50 = Icges(6,4) * t86 + Icges(6,2) * t85 + Icges(6,6) * t231;
t52 = Icges(6,1) * t86 + Icges(6,4) * t85 + Icges(6,5) * t231;
t71 = Icges(6,4) * t136 + Icges(6,2) * t135 + Icges(6,6) * t248;
t73 = Icges(6,1) * t136 + Icges(6,4) * t135 + Icges(6,5) * t248;
t14 = t133 * t50 + t134 * t52 + t165 * t191 + t71 * t83 + t73 * t84;
t9 = -t13 * t166 + t14 * t165;
t287 = t162 * t102 - (t166 * t103 - t273) * t165 - t9;
t286 = t293 * t165;
t285 = t293 * t166;
t276 = qJD(2) * t289;
t272 = 2 * m(5);
t271 = 2 * m(6);
t15 = t135 * t49 + t136 * t51 + t166 * t192 + t70 * t85 + t72 * t86;
t16 = t135 * t50 + t136 * t52 + t166 * t191 + t71 * t85 + t73 * t86;
t10 = -t15 * t166 + t16 * t165;
t270 = (-t161 * t103 + t10 + (t165 * t102 - t273) * t166) * t165;
t269 = pkin(2) * t169;
t268 = pkin(3) * t157;
t159 = t171 * pkin(2);
t146 = rSges(5,1) * t153 + rSges(5,2) * t154;
t60 = t146 * t283;
t218 = rSges(5,1) * t154 - rSges(5,2) * t153;
t65 = t289 * t218;
t265 = pkin(2) * qJD(2);
t201 = Icges(6,5) * t170 - Icges(6,6) * t168;
t264 = t154 * (t201 * t247 + (Icges(6,3) * t163 + (-Icges(6,5) * t168 - Icges(6,6) * t170) * qJD(5)) * t153);
t30 = t135 * t70 + t136 * t72 + t248 * t68;
t263 = t165 * t30;
t29 = t133 * t71 + t134 * t73 + t249 * t69;
t262 = t166 * t29;
t97 = -Icges(6,3) * t154 + t153 * t201;
t261 = t97 * t154;
t221 = pkin(4) * t154 + pkin(7) * t153;
t217 = rSges(6,1) * t170 - rSges(6,2) * t168;
t64 = t217 * t247 + (rSges(6,3) * t163 + (-rSges(6,1) * t168 - rSges(6,2) * t170) * qJD(5)) * t153;
t260 = -t221 * t163 - t64;
t253 = Icges(6,4) * t168;
t252 = Icges(6,4) * t170;
t108 = -rSges(6,3) * t154 + t153 * t217;
t147 = pkin(4) * t153 - pkin(7) * t154;
t240 = -t108 - t147;
t239 = t289 * t159;
t235 = t169 * t265;
t238 = t289 * t235;
t149 = rSges(4,1) * t156 + rSges(4,2) * t157;
t229 = -t149 - t269;
t55 = rSges(6,1) * t84 + rSges(6,2) * t83 + rSges(6,3) * t232;
t56 = rSges(6,1) * t86 + rSges(6,2) * t85 + rSges(6,3) * t231;
t27 = -t147 * t283 + t165 * t55 + t166 * t56;
t226 = t268 * t289 + t239;
t74 = rSges(6,1) * t134 + rSges(6,2) * t133 + rSges(6,3) * t249;
t75 = rSges(6,1) * t136 + rSges(6,2) * t135 + rSges(6,3) * t248;
t35 = t165 * t74 + t166 * t75 + t221 * t289;
t222 = -pkin(3) * t156 - t269;
t225 = -t238 + t289 * (t222 * qJD(2) + t235);
t219 = rSges(4,1) * t157 - rSges(4,2) * t156;
t224 = -t219 * qJD(2) - t171 * t265;
t152 = t169 * rSges(3,1) + rSges(3,2) * t171;
t214 = -t168 * t70 + t170 * t72;
t33 = t153 * t214 - t154 * t68;
t213 = -t168 * t71 + t170 * t73;
t34 = t153 * t213 - t154 * t69;
t216 = t33 * t165 + t34 * t166;
t215 = -t165 * t75 + t166 * t74;
t202 = -Icges(6,2) * t168 + t252;
t98 = -Icges(6,6) * t154 + t153 * t202;
t206 = Icges(6,1) * t170 - t253;
t99 = -Icges(6,5) * t154 + t153 * t206;
t212 = t168 * t98 - t170 * t99;
t194 = t287 * t166 + t270;
t193 = -t146 + t222;
t11 = (t163 * t214 - t47) * t154 + (t163 * t68 - t168 * t49 + t170 * t51 + (-t168 * t72 - t170 * t70) * qJD(5)) * t153;
t12 = (t163 * t213 - t48) * t154 + (t163 * t69 - t168 * t50 + t170 * t52 + (-t168 * t73 - t170 * t71) * qJD(5)) * t153;
t28 = t133 * t70 + t134 * t72 + t249 * t68;
t37 = t133 * t98 + t134 * t99 + t249 * t97;
t62 = t202 * t247 + (Icges(6,6) * t163 + (-Icges(6,2) * t170 - t253) * qJD(5)) * t153;
t63 = t206 * t247 + (Icges(6,5) * t163 + (-Icges(6,1) * t168 - t252) * qJD(5)) * t153;
t3 = (-t133 * t62 - t134 * t63 - t83 * t98 - t84 * t99 + (t262 + (t28 - t261) * t165) * t163) * t154 + (t14 * t166 + t37 * t163 + (t13 - t264) * t165) * t153;
t31 = t135 * t71 + t136 * t73 + t248 * t69;
t38 = t135 * t98 + t136 * t99 + t248 * t97;
t4 = (-t135 * t62 - t136 * t63 - t85 * t98 - t86 * t99 + (t263 + (t31 - t261) * t166) * t163) * t154 + (t15 * t165 + t38 * t163 + (t16 - t264) * t166) * t153;
t189 = t165 * t4 / 0.2e1 - t166 * t3 / 0.2e1 + (t165 * t34 - t166 * t33) * t250 / 0.2e1 - t154 * (-t11 * t166 + t12 * t165) / 0.2e1 + t9 * t249 / 0.2e1 + t10 * t248 / 0.2e1 + (t165 * (t165 * t29 - t166 * t28) + t166 * (t165 * t31 - t166 * t30)) * t247 / 0.2e1;
t57 = -t149 * t276 - t238;
t188 = t57 * t219;
t187 = (-t159 - t268) * qJD(2);
t185 = t222 + t240;
t121 = t218 * t163;
t173 = -t121 + t187;
t172 = t187 + t260;
t101 = t224 * t166;
t100 = t224 * t165;
t88 = t193 * t166;
t87 = t193 * t165;
t80 = t152 * t276;
t79 = t173 * t166;
t78 = t173 * t165;
t77 = t240 * t166;
t76 = t240 * t165;
t59 = t185 * t166;
t58 = t185 * t165;
t54 = t260 * t166;
t53 = t260 * t165;
t44 = -t108 * t248 - t154 * t75;
t43 = t108 * t249 + t154 * t74;
t42 = t172 * t166;
t41 = t172 * t165;
t40 = -t153 * t212 - t261;
t39 = t215 * t153;
t36 = t225 - t60;
t32 = t226 + t65;
t26 = (-t108 * t245 - t56) * t154 + (t163 * t75 - t166 * t64) * t153;
t25 = (t108 * t246 + t55) * t154 + (-t163 * t74 + t165 * t64) * t153;
t24 = t35 + t226;
t23 = t215 * t247 + (-t165 * t56 + t166 * t55) * t153;
t22 = t225 + t27;
t1 = [0; -m(3) * t80 + m(4) * t57 + m(5) * t36 + m(6) * t22; (t22 * t24 + t41 * t58 + t42 * t59) * t271 + (t32 * t36 + t78 * t87 + t79 * t88) * t272 + 0.2e1 * m(4) * (t239 * t57 + (t100 * t229 + t165 * t188) * t165) + t270 + 0.2e1 * m(3) * (qJD(2) * t152 - t80) * t289 * (rSges(3,1) * t171 - rSges(3,2) * t169) + t285 * t165 * t161 + (-t286 * t162 + 0.2e1 * (t101 * t229 + t166 * t188) * m(4) + (-t286 * t165 + t285 * t166) * t165 + t287) * t166; 0; m(6) * (t165 * t42 - t166 * t41) + m(5) * (t165 * t79 - t166 * t78) + m(4) * (-t100 * t166 + t101 * t165); 0; -m(5) * t60 + m(6) * t27; m(6) * (t35 * t22 + t27 * t24 + t76 * t41 + t42 * t77 + t53 * t58 + t54 * t59) + m(5) * (-t60 * t32 + t65 * t36 + (-t165 * t78 - t166 * t79) * t146 + (-t165 * t87 - t166 * t88) * t121) + t194; m(6) * (t165 * t54 - t166 * t53); (t121 * t146 * t289 - t65 * t60) * t272 + (t35 * t27 + t76 * t53 + t54 * t77) * t271 + t194; m(6) * t23; m(6) * (t22 * t39 + t23 * t24 + t25 * t59 + t26 * t58 + t41 * t44 + t42 * t43) + t189; m(6) * (t165 * t25 - t166 * t26); m(6) * (t23 * t35 + t25 * t77 + t26 * t76 + t39 * t27 + t43 * t54 + t44 * t53) + t189; (t23 * t39 + t25 * t43 + t26 * t44) * t271 + (-t154 * t38 + (t166 * t31 + t263) * t153) * t231 + t4 * t248 + (-t154 * t37 + (t165 * t28 + t262) * t153) * t232 + t3 * t249 + (t153 * t216 - t154 * t40) * t250 - t154 * ((t264 + (t154 * t212 + t216) * t163) * t154 + (t12 * t166 + t11 * t165 - (t163 * t97 - t168 * t62 + t170 * t63 + (-t168 * t99 - t170 * t98) * qJD(5)) * t154 + t40 * t163) * t153);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
