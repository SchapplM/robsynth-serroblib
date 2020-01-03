% Calculate time derivative of joint inertia matrix for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:20
% DurationCPUTime: 6.68s
% Computational Cost: add. (1800->373), mult. (4794->530), div. (0->0), fcn. (3674->4), ass. (0->192)
t148 = sin(qJ(2));
t150 = cos(qJ(2));
t249 = Icges(4,6) * t150;
t257 = Icges(3,4) * t150;
t297 = t249 + t257 + (Icges(3,1) + Icges(4,2)) * t148;
t230 = qJD(2) * t150;
t248 = Icges(5,6) * t148;
t250 = Icges(4,6) * t148;
t258 = Icges(3,4) * t148;
t296 = -t248 + t250 + t258 + (Icges(3,2) + Icges(5,2) + Icges(4,3)) * t150;
t263 = rSges(5,3) + qJ(4);
t295 = t148 * t263;
t149 = sin(qJ(1));
t274 = rSges(5,1) + pkin(3);
t294 = t149 * t274;
t151 = cos(qJ(1));
t293 = t151 * t274;
t292 = t297 * qJD(2);
t291 = t296 * t230;
t290 = t149 / 0.2e1;
t289 = -t151 / 0.2e1;
t288 = -qJD(1) / 0.2e1;
t287 = qJD(1) / 0.2e1;
t229 = qJD(2) * t151;
t218 = t148 * t229;
t234 = qJD(1) * t149;
t286 = t150 * t234 + t218;
t182 = Icges(4,2) * t150 - t250;
t281 = Icges(4,4) * t151 + t149 * t182;
t178 = -Icges(4,3) * t148 + t249;
t282 = Icges(4,5) * t151 + t149 * t178;
t191 = -t148 * t282 + t150 * t281;
t168 = t191 * t151;
t66 = Icges(4,5) * t149 - t151 * t178;
t70 = Icges(4,4) * t149 - t151 * t182;
t192 = t148 * t66 - t150 * t70;
t169 = t192 * t149;
t175 = Icges(5,3) * t150 + t248;
t65 = -Icges(5,5) * t151 + t149 * t175;
t247 = Icges(5,6) * t150;
t179 = Icges(5,2) * t148 + t247;
t69 = -Icges(5,4) * t151 + t149 * t179;
t193 = t148 * t69 + t150 * t65;
t170 = t193 * t151;
t64 = Icges(5,5) * t149 + t151 * t175;
t68 = Icges(5,4) * t149 + t151 * t179;
t194 = t148 * t68 + t150 * t64;
t171 = t194 * t149;
t187 = -Icges(3,2) * t148 + t257;
t61 = Icges(3,6) * t149 + t151 * t187;
t189 = Icges(3,1) * t150 - t258;
t63 = Icges(3,5) * t149 + t151 * t189;
t195 = t148 * t61 - t150 * t63;
t172 = t195 * t149;
t60 = -Icges(3,6) * t151 + t149 * t187;
t62 = -Icges(3,5) * t151 + t149 * t189;
t196 = t148 * t60 - t150 * t62;
t173 = t196 * t151;
t240 = t150 * t151;
t284 = t149 * rSges(4,1) - rSges(4,2) * t240;
t242 = t148 * t151;
t283 = -rSges(3,2) * t242 + t149 * rSges(3,3);
t183 = Icges(3,5) * t150 - Icges(3,6) * t148;
t58 = -Icges(3,3) * t151 + t149 * t183;
t184 = Icges(5,4) * t148 + Icges(5,5) * t150;
t73 = -Icges(5,1) * t151 + t149 * t184;
t185 = Icges(4,4) * t150 - Icges(4,5) * t148;
t280 = Icges(4,1) * t151 + t149 * t185;
t279 = 2 * m(3);
t278 = 2 * m(4);
t277 = 2 * m(5);
t276 = m(4) / 0.2e1;
t275 = m(5) / 0.2e1;
t114 = rSges(3,1) * t148 + rSges(3,2) * t150;
t273 = m(3) * t114;
t272 = pkin(2) * t148;
t271 = pkin(2) * t150;
t203 = qJ(3) * t148 + t271;
t86 = t203 * t149;
t131 = pkin(2) * t240;
t87 = qJ(3) * t242 + t131;
t270 = t149 * t86 + t151 * t87;
t206 = -rSges(4,2) * t150 + rSges(4,3) * t148;
t85 = qJD(2) * t203 - qJD(3) * t150;
t269 = -t206 * qJD(2) - t85;
t268 = rSges(4,1) * t151;
t267 = rSges(4,2) * t148;
t266 = rSges(3,3) * t151;
t265 = -rSges(5,2) - qJ(3);
t264 = -rSges(4,3) - qJ(3);
t262 = rSges(5,2) * t242 + t263 * t240 + t294;
t204 = rSges(5,2) * t148 + rSges(5,3) * t150;
t241 = t149 * t150;
t261 = qJ(4) * t241 + t149 * t204 - t293;
t59 = Icges(3,3) * t149 + t151 * t183;
t245 = qJD(1) * t59;
t72 = Icges(5,1) * t149 + t151 * t184;
t244 = qJD(1) * t72;
t74 = Icges(4,1) * t149 - t151 * t185;
t243 = qJD(1) * t74;
t112 = -qJ(3) * t150 + t272;
t205 = rSges(4,3) * t150 + t267;
t239 = -t112 + t205;
t217 = t150 * t229;
t228 = qJD(3) * t148;
t238 = qJ(3) * t217 + t151 * t228;
t220 = t148 * t234;
t233 = qJD(1) * t151;
t237 = rSges(3,2) * t220 + rSges(3,3) * t233;
t236 = t151 * pkin(1) + t149 * pkin(5);
t235 = t149 ^ 2 + t151 ^ 2;
t232 = qJD(2) * t148;
t231 = qJD(2) * t149;
t227 = qJD(4) * t150;
t226 = -pkin(2) - t263;
t225 = -0.1e1 + t235;
t122 = t231 * t272;
t224 = t149 * (qJD(1) * t131 + t149 * t228 - t122 + (t148 * t233 + t230 * t149) * qJ(3)) + t151 * (-pkin(2) * t286 - qJ(3) * t220 + t238) + t86 * t233;
t135 = pkin(5) * t233;
t221 = t135 + t238;
t57 = t239 * t151;
t216 = -rSges(5,2) * t150 + t295;
t215 = rSges(4,1) * t233 + t286 * rSges(4,2) + rSges(4,3) * t217;
t214 = rSges(5,2) * t217 + t151 * t227 + t274 * t233;
t213 = t236 + t87;
t212 = (-t185 / 0.2e1 + t184 / 0.2e1 + t183 / 0.2e1) * qJD(2);
t211 = t263 * t232;
t153 = t148 * t265 + t150 * t226 - pkin(1);
t152 = t153 * t149;
t3 = qJD(1) * t152 + t218 * t226 + t214 + t221;
t4 = t122 + (-t228 - t227 + (t150 * t265 + t295) * qJD(2)) * t149 + ((-pkin(5) - t274) * t149 + t153 * t151) * qJD(1);
t210 = t149 * t3 + t151 * t4;
t209 = -t112 - t216;
t208 = t148 * t264 - pkin(1);
t207 = rSges(3,1) * t150 - rSges(3,2) * t148;
t142 = t151 * pkin(5);
t25 = t142 + t152 + t293;
t26 = t213 + t262;
t190 = t149 * t26 + t151 * t25;
t78 = rSges(3,1) * t240 + t283;
t80 = rSges(4,3) * t242 + t284;
t174 = -pkin(1) - t207;
t52 = t209 * t151;
t167 = qJD(2) * t114;
t164 = qJD(2) * (Icges(4,4) * t148 + Icges(4,5) * t150);
t163 = qJD(2) * (Icges(5,4) * t150 - Icges(5,5) * t148);
t162 = qJD(2) * (-Icges(3,5) * t148 - Icges(3,6) * t150);
t5 = t149 * t261 + t151 * t262 + t270;
t155 = -qJ(4) * t230 - t204 * qJD(2) - qJD(4) * t148 - t85;
t88 = t112 * t234;
t8 = t151 * t155 + t216 * t234 + t88;
t9 = qJD(1) * t52 + t149 * t155;
t157 = qJD(2) * t5 + t149 * t9 + t151 * t8;
t154 = (rSges(4,2) - pkin(2)) * t150 + t208;
t100 = t207 * qJD(2);
t82 = t149 * t206 - t268;
t77 = t149 * t207 - t266;
t56 = t239 * t149;
t55 = t78 + t236;
t54 = t149 * t174 + t142 + t266;
t53 = t225 * t148 * t230;
t51 = t209 * t149;
t50 = qJD(1) * t280 + t151 * t164;
t49 = t149 * t164 + t243;
t48 = -qJD(1) * t73 + t151 * t163;
t47 = t149 * t163 + t244;
t34 = t149 * t162 + t245;
t33 = -qJD(1) * t58 + t151 * t162;
t32 = t80 + t213;
t31 = t149 * t154 + t142 + t268;
t28 = t114 * t231 + ((-rSges(3,3) - pkin(5)) * t149 + t174 * t151) * qJD(1);
t27 = -rSges(3,1) * t286 - rSges(3,2) * t217 - pkin(1) * t234 + t135 + t237;
t24 = qJD(1) * t57 + t149 * t269;
t23 = t151 * t269 - t205 * t234 + t88;
t22 = t149 * t191 + t151 * t280;
t21 = -t151 * t74 + t169;
t20 = t149 * t193 - t151 * t73;
t19 = -t151 * t72 + t171;
t18 = -t149 * t280 + t168;
t17 = t149 * t74 + t192 * t151;
t16 = t149 * t73 + t170;
t15 = t149 * t72 + t194 * t151;
t14 = t149 * t59 - t195 * t151;
t13 = t149 * t58 - t173;
t12 = -t151 * t59 - t172;
t11 = -t149 * t196 - t151 * t58;
t10 = t149 * t82 + t151 * t80 + t270;
t7 = t122 + (-t228 + (t150 * t264 - t267) * qJD(2)) * t149 + ((-rSges(4,1) - pkin(5)) * t149 + t154 * t151) * qJD(1);
t6 = -pkin(2) * t218 + (t208 - t271) * t234 + t215 + t221;
t2 = (qJD(1) * t82 + t215) * t151 + (t205 * t231 + (-t80 - t87 + t284) * qJD(1)) * t149 + t224;
t1 = (qJD(1) * t261 - t151 * t211 + t214) * t151 + ((rSges(5,2) * qJD(2) + qJD(4)) * t241 + (-t262 - t87 + t294) * qJD(1) - t149 * t211) * t149 + t224;
t29 = [(t27 * t55 + t28 * t54) * t279 + (t31 * t7 + t32 * t6) * t278 + (t25 * t4 + t26 * t3) * t277 + (t189 + t182 + t175 - t296) * t232 + (Icges(5,3) * t148 + t178 - t179 + t187 - t247 + t297) * t230; m(5) * (t25 * t8 + t26 * t9 + t3 * t51 + t4 * t52) + m(4) * (t23 * t31 + t24 * t32 + t56 * t6 + t57 * t7) + (m(3) * (-t100 * t54 - t114 * t28) + t212 * t151 + (t61 * t288 + (t66 + t68) * t287) * t150 + t290 * t291) * t151 + (m(3) * (-t100 * t55 - t114 * t27) + t212 * t149 + (t69 * t287 + (t282 + t60) * t288) * t150 + t289 * t291) * t149 + ((t70 * t287 + t292 * t290 + (t64 + t63) * t288) * t151 + (t292 * t289 + (t62 + t65 + t281) * t288) * t149) * t148 + (-t170 / 0.2e1 - t168 / 0.2e1 + t173 / 0.2e1 - t172 / 0.2e1 + t171 / 0.2e1 + t169 / 0.2e1) * qJD(2) + ((-t55 * t273 + (t61 / 0.2e1 - t68 / 0.2e1 - t66 / 0.2e1) * t150 + (t63 / 0.2e1 + t64 / 0.2e1 - t70 / 0.2e1) * t148) * t151 + (t54 * t273 + (-t69 / 0.2e1 + t282 / 0.2e1 + t60 / 0.2e1) * t150 + (t65 / 0.2e1 + t281 / 0.2e1 + t62 / 0.2e1) * t148) * t149) * qJD(1); (t10 * t2 + t23 * t57 + t24 * t56) * t278 + (t1 * t5 + t51 * t9 + t52 * t8) * t277 - t151 * ((t151 * t47 + (t19 - t170) * qJD(1)) * t151 + (t20 * qJD(1) + (t230 * t68 - t232 * t64 + t244) * t149 + (-t48 + (t148 * t65 - t150 * t69) * qJD(2) + t194 * qJD(1)) * t151) * t149) + t149 * ((t149 * t48 + (t16 - t171) * qJD(1)) * t149 + (t15 * qJD(1) + (-t230 * t69 + t232 * t65) * t151 + (-t47 + (-t148 * t64 + t150 * t68) * qJD(2) + (t193 + t72) * qJD(1)) * t149) * t151) + t149 * ((t149 * t33 + (t13 + t172) * qJD(1)) * t149 + (t14 * qJD(1) + (t230 * t60 + t232 * t62) * t151 + (-t34 + (-t148 * t63 - t150 * t61) * qJD(2) + (-t196 + t59) * qJD(1)) * t149) * t151) + t149 * ((t149 * t50 + (t18 - t169) * qJD(1)) * t149 + (t17 * qJD(1) + (t230 * t282 + t232 * t281) * t151 + (-t49 + (t148 * t70 + t150 * t66) * qJD(2) + (t191 + t74) * qJD(1)) * t149) * t151) - t151 * ((t151 * t49 + (t21 - t168) * qJD(1)) * t151 + (t22 * qJD(1) + (t230 * t66 + t232 * t70 + t243) * t149 + (-t50 + (t148 * t281 + t150 * t282) * qJD(2) + t192 * qJD(1)) * t151) * t149) - t151 * ((t151 * t34 + (t12 + t173) * qJD(1)) * t151 + (t11 * qJD(1) + (-t230 * t61 - t232 * t63 + t245) * t149 + (-t33 + (t148 * t62 + t150 * t60) * qJD(2) - t195 * qJD(1)) * t151) * t149) + ((t149 * t77 + t151 * t78) * ((qJD(1) * t77 - t151 * t167 + t237) * t151 + (-t149 * t167 + (-t78 + t283) * qJD(1)) * t149) + t235 * t114 * t100) * t279 + ((-t11 - t20 - t22) * t151 + (t12 + t19 + t21) * t149) * t234 + ((-t13 - t16 - t18) * t151 + (t14 + t15 + t17) * t149) * t233; 0.2e1 * ((t149 * t32 + t151 * t31) * t276 + t190 * t275) * t230 + 0.2e1 * ((t149 * t6 + t151 * t7 + t233 * t32 - t234 * t31) * t276 + (t233 * t26 - t234 * t25 + t210) * t275) * t148; 0.2e1 * ((t229 * t57 + t231 * t56 - t2) * t276 + (t229 * t52 + t231 * t51 - t1) * t275) * t150 + 0.2e1 * ((qJD(2) * t10 + t149 * t24 + t151 * t23 + t233 * t56 - t234 * t57) * t276 + (t233 * t51 - t234 * t52 + t157) * t275) * t148; 0.4e1 * (t276 + t275) * t53; m(5) * (-t190 * t232 + ((-t149 * t25 + t151 * t26) * qJD(1) + t210) * t150); m(5) * ((t1 + (-t149 * t51 - t151 * t52) * qJD(2)) * t148 + ((-t149 * t52 + t151 * t51) * qJD(1) + t157) * t150); m(5) * (-t148 ^ 2 + t150 ^ 2) * t225 * qJD(2); -0.2e1 * m(5) * t53;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t29(1), t29(2), t29(4), t29(7); t29(2), t29(3), t29(5), t29(8); t29(4), t29(5), t29(6), t29(9); t29(7), t29(8), t29(9), t29(10);];
Mq = res;
