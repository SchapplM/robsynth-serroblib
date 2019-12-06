% Calculate time derivative of joint inertia matrix for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:39
% EndTime: 2019-12-05 15:06:51
% DurationCPUTime: 6.29s
% Computational Cost: add. (8499->388), mult. (12575->594), div. (0->0), fcn. (11851->6), ass. (0->201)
t149 = pkin(8) + qJ(3);
t146 = cos(t149);
t151 = cos(pkin(7));
t153 = sin(qJ(4));
t212 = t151 * t153;
t150 = sin(pkin(7));
t154 = cos(qJ(4));
t213 = t150 * t154;
t138 = -t146 * t212 + t213;
t211 = t151 * t154;
t214 = t150 * t153;
t139 = t146 * t211 + t214;
t145 = sin(t149);
t217 = t145 * t151;
t92 = Icges(6,5) * t139 + Icges(6,6) * t138 + Icges(6,3) * t217;
t94 = Icges(5,5) * t139 + Icges(5,6) * t138 + Icges(5,3) * t217;
t276 = t92 + t94;
t96 = Icges(6,4) * t139 + Icges(6,2) * t138 + Icges(6,6) * t217;
t98 = Icges(5,4) * t139 + Icges(5,2) * t138 + Icges(5,6) * t217;
t286 = t96 + t98;
t100 = Icges(6,1) * t139 + Icges(6,4) * t138 + Icges(6,5) * t217;
t102 = Icges(5,1) * t139 + Icges(5,4) * t138 + Icges(5,5) * t217;
t284 = t100 + t102;
t136 = -t146 * t214 - t211;
t137 = t146 * t213 - t212;
t218 = t145 * t150;
t95 = Icges(6,4) * t137 + Icges(6,2) * t136 + Icges(6,6) * t218;
t97 = Icges(5,4) * t137 + Icges(5,2) * t136 + Icges(5,6) * t218;
t281 = t95 + t97;
t101 = Icges(5,1) * t137 + Icges(5,4) * t136 + Icges(5,5) * t218;
t99 = Icges(6,1) * t137 + Icges(6,4) * t136 + Icges(6,5) * t218;
t280 = t101 + t99;
t283 = t286 * t136 + t284 * t137 + t276 * t218;
t91 = Icges(6,5) * t137 + Icges(6,6) * t136 + Icges(6,3) * t218;
t93 = Icges(5,5) * t137 + Icges(5,6) * t136 + Icges(5,3) * t218;
t282 = t91 + t93;
t173 = Icges(6,5) * t154 - Icges(6,6) * t153;
t114 = -Icges(6,3) * t146 + t145 * t173;
t174 = Icges(5,5) * t154 - Icges(5,6) * t153;
t115 = -Icges(5,3) * t146 + t145 * t174;
t279 = t114 + t115;
t221 = Icges(6,4) * t154;
t175 = -Icges(6,2) * t153 + t221;
t116 = -Icges(6,6) * t146 + t145 * t175;
t223 = Icges(5,4) * t154;
t176 = -Icges(5,2) * t153 + t223;
t117 = -Icges(5,6) * t146 + t145 * t176;
t275 = t116 + t117;
t222 = Icges(6,4) * t153;
t178 = Icges(6,1) * t154 - t222;
t118 = -Icges(6,5) * t146 + t145 * t178;
t224 = Icges(5,4) * t153;
t179 = Icges(5,1) * t154 - t224;
t119 = -Icges(5,5) * t146 + t145 * t179;
t278 = t118 + t119;
t274 = -t153 * t286 + t284 * t154;
t273 = -t281 * t153 + t280 * t154;
t258 = t281 * t136 + t280 * t137 + t282 * t218;
t272 = t275 * t136 + t278 * t137 + t279 * t218;
t204 = qJD(4) * t145;
t85 = (-Icges(6,2) * t154 - t222) * t204 + (Icges(6,6) * t145 + t146 * t175) * qJD(3);
t86 = (-Icges(5,2) * t154 - t224) * t204 + (Icges(5,6) * t145 + t146 * t176) * qJD(3);
t271 = -t85 - t86;
t87 = (-Icges(6,1) * t153 - t221) * t204 + (Icges(6,5) * t145 + t146 * t178) * qJD(3);
t88 = (-Icges(5,1) * t153 - t223) * t204 + (Icges(5,5) * t145 + t146 * t179) * qJD(3);
t270 = t87 + t88;
t215 = t146 * t115;
t216 = t146 * t114;
t267 = -t215 - t216;
t266 = t283 * t151;
t234 = t146 * ((-Icges(5,5) * t153 - Icges(5,6) * t154) * t204 + (Icges(5,3) * t145 + t146 * t174) * qJD(3));
t235 = t146 * ((-Icges(6,5) * t153 - Icges(6,6) * t154) * t204 + (Icges(6,3) * t145 + t146 * t173) * qJD(3));
t265 = -t235 - t234;
t264 = t273 * t145 - t146 * t282;
t263 = t274 * t145 - t276 * t146;
t148 = t151 ^ 2;
t253 = t150 ^ 2 + t148;
t262 = t275 * t153 - t154 * t278;
t206 = qJD(3) * t145;
t196 = t153 * t206;
t109 = -qJD(4) * t137 + t150 * t196;
t195 = t154 * t206;
t110 = qJD(4) * t136 - t150 * t195;
t205 = qJD(3) * t146;
t194 = t150 * t205;
t61 = Icges(6,5) * t110 + Icges(6,6) * t109 + Icges(6,3) * t194;
t166 = t145 * t61 + t205 * t91;
t65 = Icges(6,4) * t110 + Icges(6,2) * t109 + Icges(6,6) * t194;
t69 = Icges(6,1) * t110 + Icges(6,4) * t109 + Icges(6,5) * t194;
t14 = t109 * t95 + t110 * t99 + t136 * t65 + t137 * t69 + t150 * t166;
t111 = -qJD(4) * t139 + t151 * t196;
t112 = qJD(4) * t138 - t151 * t195;
t193 = t151 * t205;
t62 = Icges(6,5) * t112 + Icges(6,6) * t111 + Icges(6,3) * t193;
t165 = t145 * t62 + t205 * t92;
t66 = Icges(6,4) * t112 + Icges(6,2) * t111 + Icges(6,6) * t193;
t70 = Icges(6,1) * t112 + Icges(6,4) * t111 + Icges(6,5) * t193;
t15 = t100 * t110 + t109 * t96 + t136 * t66 + t137 * t70 + t150 * t165;
t63 = Icges(5,5) * t110 + Icges(5,6) * t109 + Icges(5,3) * t194;
t164 = t145 * t63 + t205 * t93;
t67 = Icges(5,4) * t110 + Icges(5,2) * t109 + Icges(5,6) * t194;
t71 = Icges(5,1) * t110 + Icges(5,4) * t109 + Icges(5,5) * t194;
t16 = t101 * t110 + t109 * t97 + t136 * t67 + t137 * t71 + t150 * t164;
t64 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t193;
t163 = t145 * t64 + t205 * t94;
t68 = Icges(5,4) * t112 + Icges(5,2) * t111 + Icges(5,6) * t193;
t72 = Icges(5,1) * t112 + Icges(5,4) * t111 + Icges(5,5) * t193;
t17 = t102 * t110 + t109 * t98 + t136 * t68 + t137 * t72 + t150 * t163;
t261 = (-t109 * t275 - t110 * t278 + t271 * t136 - t270 * t137) * t146 + ((t15 + t17) * t151 + (t14 + t16 + t265) * t150) * t145 + (((t267 + t258) * t150 + t266) * t146 + t272 * t145) * qJD(3);
t260 = (-t273 * qJD(3) + t61 + t63) * t146 + ((-t69 - t71) * t154 + (t65 + t67) * t153 + (t280 * t153 + t281 * t154) * qJD(4) - t282 * qJD(3)) * t145;
t259 = (t274 * qJD(3) - t62 - t64) * t146 + ((t70 + t72) * t154 + (-t66 - t68) * t153 + (-t284 * t153 - t154 * t286) * qJD(4) + t276 * qJD(3)) * t145;
t36 = t100 * t139 + t138 * t96 + t217 * t92;
t38 = t102 * t139 + t138 * t98 + t217 * t94;
t257 = t36 + t38;
t256 = t262 * t145 - t267;
t252 = t264 * t150 + t263 * t151;
t251 = qJD(3) * (rSges(4,1) * t145 + rSges(4,2) * t146);
t240 = pkin(4) * t154;
t249 = qJ(5) * t146 - t145 * t240;
t247 = 2 * m(5);
t246 = 2 * m(6);
t243 = t150 / 0.2e1;
t203 = qJD(4) * t153;
t199 = pkin(4) * t203;
t155 = t249 * qJD(3) + qJD(5) * t145 - t146 * t199;
t202 = qJD(4) * t154;
t198 = pkin(4) * t202;
t238 = rSges(6,1) * t110 + rSges(6,2) * t109 + rSges(6,3) * t194 + t150 * t155 - t151 * t198;
t237 = rSges(6,1) * t112 + rSges(6,2) * t111 + rSges(6,3) * t193 + t150 * t198 + t151 * t155;
t157 = qJ(5) * t145 + t146 * t240;
t187 = rSges(6,1) * t154 - rSges(6,2) * t153;
t236 = -qJD(5) * t146 - t145 * t199 + (-rSges(6,1) * t153 - rSges(6,2) * t154) * t204 + (rSges(6,3) * t145 + t146 * t187 + t157) * qJD(3);
t35 = t138 * t95 + t139 * t99 + t217 * t91;
t233 = t150 * t35;
t37 = t101 * t139 + t138 * t97 + t217 * t93;
t232 = t150 * t37;
t229 = rSges(6,1) * t137 + rSges(6,2) * t136 + rSges(6,3) * t218 - pkin(4) * t212 + t150 * t157;
t228 = rSges(6,1) * t139 + rSges(6,2) * t138 + rSges(6,3) * t217 + pkin(4) * t214 + t151 * t157;
t190 = pkin(3) * t146 + pkin(6) * t145;
t141 = t190 * qJD(3);
t188 = rSges(5,1) * t154 - rSges(5,2) * t153;
t90 = (-rSges(5,1) * t153 - rSges(5,2) * t154) * t204 + (rSges(5,3) * t145 + t146 * t188) * qJD(3);
t227 = -t141 - t90;
t123 = -t146 * rSges(5,3) + t145 * t188;
t220 = t123 * t146;
t210 = -t146 * rSges(6,3) + t145 * t187 - t249;
t143 = pkin(3) * t145 - pkin(6) * t146;
t209 = t253 * qJD(3) * t143;
t208 = -t123 - t143;
t207 = t253 * t190;
t201 = -t141 - t236;
t200 = m(6) * t206;
t197 = -t143 - t210;
t192 = t151 * t210;
t191 = t210 * t150;
t104 = rSges(5,1) * t137 + rSges(5,2) * t136 + rSges(5,3) * t218;
t106 = rSges(5,1) * t139 + rSges(5,2) * t138 + rSges(5,3) * t217;
t172 = t104 * t151 - t106 * t150;
t158 = qJD(3) * (-Icges(4,5) * t145 - Icges(4,6) * t146);
t156 = -t150 * t228 + t151 * t229;
t130 = t150 * t158;
t108 = t208 * t151;
t107 = t208 * t150;
t80 = t253 * t251;
t78 = t227 * t151;
t77 = t227 * t150;
t76 = rSges(5,1) * t112 + rSges(5,2) * t111 + rSges(5,3) * t193;
t74 = rSges(5,1) * t110 + rSges(5,2) * t109 + rSges(5,3) * t194;
t58 = t197 * t151;
t57 = t197 * t150;
t56 = -t106 * t146 - t123 * t217;
t55 = t104 * t146 + t123 * t218;
t52 = t172 * t145;
t51 = t115 * t217 + t117 * t138 + t119 * t139;
t50 = t114 * t217 + t116 * t138 + t118 * t139;
t47 = t201 * t151;
t46 = t201 * t150;
t45 = t104 * t150 + t106 * t151 + t207;
t40 = -t145 * t192 - t146 * t228;
t39 = t145 * t191 + t146 * t229;
t30 = t150 * t74 + t151 * t76 - t209;
t29 = t156 * t145;
t28 = -t90 * t217 - t146 * t76 + (t106 * t145 - t151 * t220) * qJD(3);
t27 = t90 * t218 + t146 * t74 + (-t104 * t145 + t150 * t220) * qJD(3);
t26 = t150 * t229 + t151 * t228 + t207;
t25 = (-t150 * t76 + t151 * t74) * t145 + t172 * t205;
t24 = t150 * t238 + t151 * t237 - t209;
t23 = -t237 * t146 - t236 * t217 + (t145 * t228 - t146 * t192) * qJD(3);
t22 = t238 * t146 + t236 * t218 + (-t145 * t229 + t146 * t191) * qJD(3);
t21 = t102 * t112 + t111 * t98 + t138 * t68 + t139 * t72 + t151 * t163;
t20 = t101 * t112 + t111 * t97 + t138 * t67 + t139 * t71 + t151 * t164;
t19 = t100 * t112 + t111 * t96 + t138 * t66 + t139 * t70 + t151 * t165;
t18 = t111 * t95 + t112 * t99 + t138 * t65 + t139 * t69 + t151 * t166;
t9 = (-t150 * t237 + t151 * t238) * t145 + t156 * t205;
t8 = t150 * t21 - t151 * t20;
t7 = t150 * t19 - t151 * t18;
t6 = t150 * t17 - t151 * t16;
t5 = -t14 * t151 + t15 * t150;
t4 = -(t111 * t117 + t112 * t119 + t138 * t86 + t139 * t88) * t146 + (t20 * t150 + (t21 - t234) * t151) * t145 + (t51 * t145 + (t232 + (t38 - t215) * t151) * t146) * qJD(3);
t3 = -(t111 * t116 + t112 * t118 + t138 * t85 + t139 * t87) * t146 + (t18 * t150 + (t19 - t235) * t151) * t145 + (t50 * t145 + (t233 + (t36 - t216) * t151) * t146) * qJD(3);
t1 = [0; 0; 0; -m(4) * t80 + m(5) * t30 + m(6) * t24; m(5) * (t150 * t78 - t151 * t77) + m(6) * (t150 * t47 - t151 * t46); (t24 * t26 + t46 * t57 + t47 * t58) * t246 + (t107 * t77 + t108 * t78 + t45 * t30) * t247 + 0.2e1 * m(4) * (t251 - t80) * t253 * (rSges(4,1) * t146 - rSges(4,2) * t145) + (-t148 * t130 - t5 - t6) * t151 + (t7 + t8 + (-t150 * t130 + t253 * t158) * t151) * t150; m(5) * t25 + m(6) * t9; m(5) * (t150 * t27 - t151 * t28) + m(6) * (t150 * t22 - t151 * t23); t4 * t243 + t3 * t243 + m(6) * (t22 * t58 + t23 * t57 + t24 * t29 + t26 * t9 + t39 * t47 + t40 * t46) + m(5) * (t107 * t28 + t108 * t27 + t25 * t45 + t52 * t30 + t55 * t78 + t56 * t77) + ((t8 / 0.2e1 + t7 / 0.2e1) * t151 + (t6 / 0.2e1 + t5 / 0.2e1) * t150) * t145 + (((t283 * t150 - t258 * t151) * t243 + ((-t35 - t37) * t151 + t257 * t150) * t151 / 0.2e1) * t146 + (t263 * t150 - t264 * t151) * t145 / 0.2e1) * qJD(3) - (t259 * t150 + t260 * t151) * t146 / 0.2e1 - t261 * t151 / 0.2e1; (t22 * t39 + t23 * t40 + t29 * t9) * t246 + (t25 * t52 + t27 * t55 + t28 * t56) * t247 + t261 * t218 + (t3 + t4) * t217 + (t252 * t206 + (t258 * t150 + t266) * t194 + (t257 * t151 + t232 + t233) * t193) * t145 + (t265 * t146 + t256 * t206 - t272 * t194 + (-t51 - t50) * t193 + ((t271 * t153 + t270 * t154 - t202 * t275 - t203 * t278) * t146 - t259 * t151 + t260 * t150) * t145 + ((-t262 * t146 - t252) * t146 + (t279 * t146 + t256) * t145) * qJD(3)) * t146; t200; 0; m(6) * (-t146 * t24 + (t150 * t46 + t151 * t47) * t145 + (t145 * t26 + (t150 * t57 + t151 * t58) * t146) * qJD(3)); m(6) * (-t146 * t9 + (t150 * t23 + t151 * t22) * t145 + (t145 * t29 + (t150 * t40 + t151 * t39) * t146) * qJD(3)); 0.2e1 * (-0.1e1 + t253) * t146 * t200;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
