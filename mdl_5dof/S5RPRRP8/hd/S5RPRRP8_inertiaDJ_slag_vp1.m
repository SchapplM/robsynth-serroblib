% Calculate time derivative of joint inertia matrix for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:08
% EndTime: 2019-12-31 18:47:18
% DurationCPUTime: 5.78s
% Computational Cost: add. (4677->345), mult. (10896->501), div. (0->0), fcn. (10816->6), ass. (0->171)
t328 = Icges(5,1) + Icges(6,1);
t327 = Icges(5,4) - Icges(6,5);
t326 = Icges(6,4) + Icges(5,5);
t175 = cos(qJ(4));
t301 = t175 / 0.2e1;
t302 = -t175 / 0.2e1;
t174 = sin(qJ(4));
t303 = t174 / 0.2e1;
t304 = -t174 / 0.2e1;
t314 = Icges(5,6) * t302 + Icges(6,6) * t301 + t326 * t304;
t325 = t301 * Icges(5,6) + Icges(6,6) * t302 + t326 * t303 - t314;
t261 = Icges(6,5) * t174;
t210 = -Icges(6,1) * t175 - t261;
t263 = Icges(5,4) * t174;
t211 = -Icges(5,1) * t175 + t263;
t324 = (t210 + t211) * qJD(4);
t260 = Icges(6,5) * t175;
t262 = Icges(5,4) * t175;
t243 = -t328 * t174 + t260 - t262;
t239 = qJD(4) * t174;
t320 = rSges(6,1) + pkin(4);
t323 = -qJD(5) * t174 + t320 * t239;
t319 = rSges(6,3) + qJ(5);
t269 = sin(qJ(3));
t270 = sin(qJ(1));
t271 = cos(qJ(3));
t272 = cos(qJ(1));
t143 = -t270 * t269 - t272 * t271;
t312 = qJD(1) - qJD(3);
t101 = t312 * t143;
t322 = t101 * t174;
t257 = t101 * t175;
t144 = t272 * t269 - t270 * t271;
t102 = t312 * t144;
t321 = t102 * t174;
t254 = t102 * t175;
t296 = -qJD(4) / 0.2e1;
t207 = -Icges(5,5) * t175 + Icges(5,6) * t174;
t208 = -Icges(6,4) * t175 - Icges(6,6) * t174;
t318 = (t207 + t208) * qJD(4);
t240 = t272 * pkin(1) + t270 * qJ(2);
t235 = t272 * pkin(2) + t240;
t206 = -Icges(6,3) * t174 - t260;
t137 = t206 * qJD(4);
t209 = Icges(5,2) * t174 - t262;
t140 = t209 * qJD(4);
t153 = -Icges(5,2) * t175 - t263;
t280 = t153 * t239 - (-t137 + t140) * t175 - t324 * t174;
t238 = qJD(4) * t175;
t199 = t144 * t238 + t322;
t198 = t144 * t239 - t257;
t197 = t143 * t238 - t321;
t196 = t143 * t239 + t254;
t317 = -rSges(3,1) * t272 - rSges(3,3) * t270 - t240;
t150 = Icges(6,3) * t175 - t261;
t248 = t174 * t150;
t316 = -t248 / 0.2e1 + t153 * t303 + t243 * t302;
t315 = t137 * t303 + t140 * t304 + t324 * t301 + ((-t150 + t153) * t175 + t243 * t174) * t296;
t179 = t101 * rSges(6,2) + t323 * t143 - t197 * t319 + t254 * t320;
t180 = t102 * rSges(6,2) + t323 * t144 - t199 * t319 - t257 * t320;
t311 = -(t175 * t243 + t248) * qJD(4) + t280;
t295 = qJD(4) / 0.2e1;
t67 = -Icges(5,6) * t143 + t144 * t209;
t71 = -Icges(5,5) * t143 + t144 * t211;
t217 = t174 * t67 - t175 * t71;
t294 = t143 * t217;
t61 = -Icges(6,6) * t143 + t144 * t206;
t69 = -Icges(6,4) * t143 + t144 * t210;
t219 = t174 * t61 + t175 * t69;
t293 = t143 * t219;
t68 = Icges(5,6) * t144 + t143 * t209;
t72 = Icges(5,5) * t144 + t143 * t211;
t216 = t174 * t68 - t175 * t72;
t292 = t144 * t216;
t62 = Icges(6,6) * t144 + t143 * t206;
t70 = Icges(6,4) * t144 + t143 * t210;
t218 = t174 * t62 + t175 * t70;
t291 = t144 * t218;
t288 = t143 ^ 2 + t144 ^ 2;
t287 = t101 * t144 - t102 * t143;
t194 = rSges(5,1) * t198 + t102 * rSges(5,3);
t286 = -rSges(5,2) * t199 - t194;
t195 = rSges(5,1) * t196 + t101 * rSges(5,3);
t285 = rSges(5,2) * t197 + t195;
t278 = -(t301 * t62 + t302 * t68 + (t70 + t72) * t304 + t314 * t144 + t316 * t143) * t101 - (t301 * t61 + t302 * t67 + (t69 + t71) * t304 + t316 * t144 - t314 * t143) * t102;
t277 = 2 * m(4);
t276 = 2 * m(5);
t275 = 2 * m(6);
t158 = -t174 * rSges(5,1) - rSges(5,2) * t175;
t268 = m(5) * t158;
t250 = t144 * t175;
t251 = t144 * t174;
t267 = t143 * rSges(6,2) + t250 * t320 + t251 * t319;
t252 = t143 * t175;
t253 = t143 * t174;
t266 = -t144 * rSges(6,2) + t252 * t320 + t253 * t319;
t265 = t102 * pkin(3) + t101 * pkin(7);
t264 = -t101 * pkin(3) + t102 * pkin(7);
t247 = -qJD(5) * t175 + (t174 * t319 + t175 * t320) * qJD(4);
t246 = t144 * pkin(3) + t143 * pkin(7);
t245 = t143 * pkin(3) - t144 * pkin(7);
t242 = -t174 * t320 + t175 * t319;
t170 = t272 * qJ(2);
t241 = qJD(1) * t170 + qJD(2) * t270;
t236 = t270 * pkin(1);
t60 = -t102 * rSges(4,1) + t101 * rSges(4,2);
t59 = -t101 * rSges(4,1) - t102 * rSges(4,2);
t103 = -t144 * rSges(4,1) + t143 * rSges(4,2);
t104 = rSges(4,1) * t143 + rSges(4,2) * t144;
t203 = t217 * t296 + t219 * t295 + (Icges(6,5) * t198 - Icges(6,3) * t199) * t302 + (Icges(5,4) * t198 + Icges(5,2) * t199) * t301 + t318 * t143 / 0.2e1 + (t328 * t198 + t327 * t199) * t303 + t315 * t144 - t316 * t101 + t325 * t102;
t202 = t216 * t296 + t218 * t295 + (Icges(6,5) * t196 - Icges(6,3) * t197) * t302 + t301 * (Icges(5,4) * t196 + Icges(5,2) * t197) - t318 * t144 / 0.2e1 + (t328 * t196 + t327 * t197) * t303 + t315 * t143 + t316 * t102 + t325 * t101;
t74 = -rSges(5,1) * t250 + rSges(5,2) * t251 - t143 * rSges(5,3);
t76 = -rSges(5,1) * t252 + rSges(5,2) * t253 + t144 * rSges(5,3);
t193 = -pkin(2) * t270 - t236;
t44 = t245 + t266;
t188 = -t143 * t270 + t144 * t272;
t54 = -t76 + t245;
t187 = t170 + t193;
t186 = t187 + t246;
t185 = -rSges(3,1) * t270 + rSges(3,3) * t272 - t236;
t184 = qJD(1) * t193 + t241;
t181 = t184 + t265;
t168 = qJD(2) * t272;
t178 = -qJD(1) * t235 + t168;
t177 = t178 - t264;
t176 = t270 * t102 + t272 * t101 + (-t143 * t272 - t144 * t270) * qJD(1);
t148 = (-rSges(5,1) * t175 + rSges(5,2) * t174) * qJD(4);
t111 = t170 + t185;
t106 = qJD(1) * t317 + t168;
t105 = qJD(1) * t185 + t241;
t80 = -t104 + t235;
t79 = -t103 + t187;
t78 = t242 * t143;
t77 = t242 * t144;
t66 = Icges(6,2) * t144 + t143 * t208;
t65 = -Icges(6,2) * t143 + t144 * t208;
t64 = Icges(5,3) * t144 + t143 * t207;
t63 = -Icges(5,3) * t143 + t144 * t207;
t53 = t74 - t246;
t52 = t178 - t59;
t51 = t184 - t60;
t50 = -t54 + t235;
t49 = t186 - t74;
t43 = -t246 - t267;
t42 = -t44 + t235;
t41 = t186 + t267;
t40 = t102 * t242 + t143 * t247;
t39 = -t101 * t242 + t144 * t247;
t32 = Icges(6,4) * t196 + Icges(6,2) * t101 - Icges(6,6) * t197;
t31 = Icges(6,4) * t198 + Icges(6,2) * t102 - Icges(6,6) * t199;
t30 = Icges(5,5) * t196 + Icges(5,6) * t197 + Icges(5,3) * t101;
t29 = Icges(5,5) * t198 + Icges(5,6) * t199 + Icges(5,3) * t102;
t26 = t143 * t216 + t144 * t64;
t25 = t144 * t63 + t294;
t24 = -t143 * t218 + t144 * t66;
t23 = t144 * t65 - t293;
t22 = -t143 * t64 + t292;
t21 = -t143 * t63 + t144 * t217;
t20 = -t143 * t66 - t291;
t19 = -t143 * t65 - t144 * t219;
t18 = -t265 - t285;
t17 = t264 - t286;
t16 = -t143 * t266 - t144 * t267;
t15 = t177 + t286;
t14 = t181 + t285;
t9 = -t179 - t265;
t8 = t180 + t264;
t7 = t177 - t180;
t6 = t179 + t181;
t1 = -t101 * t267 + t102 * t266 + t143 * t179 + t144 * t180;
t2 = [(t41 * t7 + t42 * t6) * t275 + (t14 * t50 + t15 * t49) * t276 - t150 * t239 + (t51 * t80 + t52 * t79) * t277 + 0.2e1 * m(3) * (-t105 * t317 + t106 * t111) - t243 * t238 + t280; m(6) * (t270 * t7 - t272 * t6 + (t270 * t42 + t272 * t41) * qJD(1)) + m(5) * (t270 * t15 - t272 * t14 + (t270 * t50 + t272 * t49) * qJD(1)) + m(4) * (t270 * t52 - t272 * t51 + (t270 * t80 + t272 * t79) * qJD(1)) + m(3) * (t270 * t106 - t272 * t105 + (t111 * t272 - t270 * t317) * qJD(1)); 0; m(6) * (t41 * t8 + t42 * t9 + t43 * t7 + t44 * t6) + m(5) * (t14 * t54 + t15 * t53 + t17 * t49 + t18 * t50) + m(4) * (t103 * t52 + t104 * t51 + t59 * t79 + t60 * t80) - t311; m(4) * (t59 * t270 - t60 * t272 + (t103 * t272 + t104 * t270) * qJD(1)) + m(5) * (t17 * t270 - t18 * t272 + (t270 * t54 + t272 * t53) * qJD(1)) + m(6) * (t8 * t270 - t9 * t272 + (t270 * t44 + t272 * t43) * qJD(1)); (t43 * t8 + t44 * t9) * t275 + (t17 * t53 + t18 * t54) * t276 + (t103 * t59 + t104 * t60) * t277 + t311; m(6) * (t39 * t42 + t40 * t41 - t6 * t77 - t7 * t78) + (-t101 * t50 + t102 * t49) * t268 + (m(5) * (-t14 * t158 - t148 * t50) - t202) * t144 + (m(5) * (-t148 * t49 - t15 * t158) + t203) * t143 - t278; m(6) * (t40 * t270 - t39 * t272 + (-t270 * t77 - t272 * t78) * qJD(1)) + m(5) * t188 * t148 + t176 * t268; m(6) * (t39 * t44 + t40 * t43 - t77 * t9 - t78 * t8) + (-t101 * t54 + t102 * t53) * t268 + (m(5) * (-t148 * t54 - t158 * t18) + t202) * t144 + (m(5) * (-t148 * t53 - t158 * t17) - t203) * t143 + t278; ((t143 * t76 + t144 * t74) * (t101 * t74 + t144 * t194 - t102 * t76 + t143 * t195 + (t143 * t197 + t144 * t199) * rSges(5,2)) + (t288 * t148 + t158 * t287) * t158) * t276 + t101 * (-t25 * t143 + t26 * t144) + t144 * ((t101 * t64 + t144 * t30) * t144 + t26 * t101 + (t25 - t292) * t102 + (-t101 * t63 - t144 * t29 - t254 * t71 + t321 * t67) * t143) + t102 * (-t21 * t143 + t22 * t144) - t143 * (-(t102 * t63 - t143 * t29) * t143 + t21 * t102 - (-t22 + t294) * t101 + (t102 * t64 - t143 * t30 - t257 * t72 + t322 * t68) * t144) + (t16 * t1 - t39 * t77 - t40 * t78) * t275 + t101 * (-t23 * t143 + t24 * t144) + t144 * ((t101 * t66 + t144 * t32) * t144 + t24 * t101 + (t23 + t291) * t102 + (-t101 * t65 - t144 * t31 - t254 * t69 - t321 * t61) * t143) + t102 * (-t19 * t143 + t20 * t144) - t143 * (-(t102 * t65 - t143 * t31) * t143 + t19 * t102 - (-t20 - t293) * t101 + (t102 * t66 - t143 * t32 - t257 * t70 - t322 * t62) * t144); m(6) * ((-t143 * t41 - t144 * t42) * t238 + (-t101 * t42 + t102 * t41 - t143 * t7 - t144 * t6) * t174); m(6) * (t174 * t176 + t188 * t238); m(6) * ((-t143 * t43 - t144 * t44) * t238 + (-t101 * t44 + t102 * t43 - t143 * t8 - t144 * t9) * t174); m(6) * ((t1 + (t143 * t78 + t144 * t77) * qJD(4)) * t175 + (-qJD(4) * t16 + t101 * t77 - t102 * t78 - t143 * t40 - t144 * t39) * t174); ((-0.1e1 + t288) * t238 + t287 * t174) * t174 * t275;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
