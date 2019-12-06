% Calculate time derivative of joint inertia matrix for
% S5PRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:38
% EndTime: 2019-12-05 16:41:46
% DurationCPUTime: 3.57s
% Computational Cost: add. (5364->209), mult. (4206->298), div. (0->0), fcn. (3113->6), ass. (0->138)
t291 = Icges(5,4) - Icges(6,5);
t290 = Icges(5,1) + Icges(6,1);
t289 = Icges(5,2) + Icges(6,3);
t148 = cos(qJ(4));
t288 = t291 * t148;
t147 = sin(qJ(4));
t287 = t291 * t147;
t286 = Icges(6,4) + Icges(5,5);
t285 = Icges(5,6) - Icges(6,6);
t284 = t289 * t147 - t288;
t283 = t290 * t148 - t287;
t282 = -t289 * t148 - t287;
t281 = t290 * t147 + t288;
t145 = pkin(8) + qJ(2);
t144 = qJ(3) + t145;
t141 = cos(t144);
t280 = t284 * t141;
t279 = t283 * t141;
t274 = -t285 * t147 + t286 * t148;
t278 = t283 * t147 - t284 * t148;
t140 = sin(t144);
t272 = -t284 * t140 - t285 * t141;
t271 = t285 * t140 - t280;
t270 = t283 * t140 - t286 * t141;
t269 = t286 * t140 + t279;
t260 = t272 * t147 - t270 * t148;
t277 = t260 * t140;
t276 = Icges(6,2) + Icges(5,3);
t257 = t286 * t147 + t285 * t148;
t240 = rSges(6,3) + qJ(5);
t258 = rSges(6,1) + pkin(4);
t252 = t147 * t240 + t148 * t258;
t268 = t274 * t141;
t255 = t147 * t282 + t281 * t148;
t203 = qJD(4) * t147;
t195 = t140 * t203;
t201 = qJD(5) * t147;
t202 = qJD(4) * t148;
t266 = -(t240 * t202 + t201) * t140 + t258 * t195;
t265 = t274 * t140 - t276 * t141;
t264 = t276 * t140 + t268;
t146 = qJD(2) + qJD(3);
t211 = t141 * t146;
t239 = t271 * t147 - t269 * t148;
t263 = t239 * t140;
t262 = t239 * t141;
t261 = -t257 * qJD(4) + t276 * t146;
t254 = -t269 * t147 - t271 * t148;
t253 = -t270 * t147 - t272 * t148;
t251 = -qJD(4) * t274 + t255 * t146;
t250 = t265 * t140 - t263 + (-t260 - t264) * t141;
t249 = t140 ^ 2;
t248 = t141 ^ 2;
t128 = t147 * rSges(5,1) + rSges(5,2) * t148;
t129 = t140 * rSges(5,3);
t210 = t141 * t147;
t199 = rSges(5,2) * t210;
t196 = -t140 * rSges(5,2) * t202 - rSges(5,1) * t195 - t146 * t199;
t223 = rSges(5,2) * t147;
t115 = t140 * t223;
t218 = rSges(5,3) * t211 + t146 * t115;
t207 = t141 * rSges(5,3) + t115;
t212 = t140 * t148;
t75 = rSges(5,1) * t212 - t207;
t209 = t141 * t148;
t77 = rSges(5,1) * t209 + t129 - t199;
t2 = ((-t77 + t129) * t146 + t196) * t140 + (-qJD(4) * t128 * t141 + t146 * t75 + t218) * t141;
t234 = 2 * m(5);
t247 = t2 * t234;
t246 = t265 * t141 + t277;
t243 = t264 * t140 - t262;
t213 = t140 * t146;
t242 = t261 * t141 - t213 * t274;
t241 = -t261 * t140 - t268 * t146;
t235 = 2 * m(4);
t233 = 2 * m(6);
t224 = rSges(5,1) * t148;
t102 = (-t223 + t224) * qJD(4);
t229 = m(5) * t102;
t228 = m(5) * t128;
t142 = sin(t145);
t227 = pkin(2) * t142;
t132 = t141 * rSges(6,2);
t226 = t252 * t140 - t132;
t130 = t140 * rSges(6,2);
t225 = t258 * t209 + t240 * t210 + t130;
t222 = pkin(2) * qJD(2);
t206 = t258 * t147 - t240 * t148;
t59 = t206 * t141;
t221 = t146 * t59;
t219 = -t252 * qJD(4) + qJD(5) * t148;
t205 = t141 * pkin(3) + t140 * pkin(7);
t204 = t248 + t249;
t200 = m(6) * t203;
t198 = t142 * t222;
t143 = cos(t145);
t197 = t143 * t222;
t194 = t141 * t203;
t193 = t141 * t202;
t190 = -pkin(3) - t224;
t86 = t141 * rSges(4,1) - rSges(4,2) * t140;
t79 = -rSges(4,1) * t211 + rSges(4,2) * t213;
t85 = -rSges(4,1) * t140 - rSges(4,2) * t141;
t78 = t85 * t146;
t41 = t205 + t225;
t57 = t77 + t205;
t135 = t141 * pkin(7);
t56 = t140 * t190 + t135 + t207;
t161 = -pkin(3) - t252;
t160 = rSges(6,2) * t211 + t141 * t201 + t240 * t193 - t194 * t258;
t159 = t161 * t140;
t152 = t278 * qJD(4) + t281 * t202 + t203 * t282;
t40 = t132 + t135 + t159;
t151 = (-t239 * qJD(4) - t251 * t140 - t213 * t278) * t140 / 0.2e1 - (-qJD(4) * t260 + t251 * t141) * t141 / 0.2e1 + (t255 * t140 - t257 * t141 - t253) * t213 / 0.2e1 + (t257 * t140 + t255 * t141 - t254) * t211 / 0.2e1 - (t147 * t279 - t148 * t280) * t211 / 0.2e1;
t27 = (t190 * t141 + (-rSges(5,3) - pkin(7)) * t140) * t146 - t196;
t112 = pkin(7) * t211;
t26 = -rSges(5,2) * t193 - pkin(3) * t213 + t112 + (-t146 * t212 - t194) * rSges(5,1) + t218;
t13 = t146 * t159 + t112 + t160;
t14 = ((-rSges(6,2) - pkin(7)) * t140 + t161 * t141) * t146 + t266;
t139 = pkin(2) * t143;
t81 = t139 + t86;
t80 = t85 - t227;
t61 = t79 - t197;
t60 = t78 - t198;
t58 = t206 * t140;
t55 = t139 + t57;
t54 = t56 - t227;
t35 = t139 + t41;
t34 = t40 - t227;
t29 = t140 * t219 - t221;
t28 = t141 * t219 + t206 * t213;
t25 = t27 - t197;
t24 = t26 - t198;
t15 = t140 * t226 + t141 * t225;
t12 = t14 - t197;
t11 = t13 - t198;
t1 = (t146 * t226 + t160) * t141 + ((t130 - t225) * t146 - t266) * t140;
t3 = [0; 0; (t11 * t35 + t12 * t34) * t233 + (t24 * t55 + t25 * t54) * t234 + (t60 * t81 + t61 * t80) * t235 + t152; 0; m(6) * (t11 * t41 + t12 * t40 + t13 * t35 + t14 * t34) + m(5) * (t24 * t57 + t25 * t56 + t26 * t55 + t27 * t54) + m(4) * (t60 * t86 + t61 * t85 + t78 * t81 + t79 * t80) + t152; (t13 * t41 + t14 * t40) * t233 + (t26 * t57 + t27 * t56) * t234 + (t78 * t86 + t79 * t85) * t235 + t152; m(5) * t2 + m(6) * t1; ((-t146 * t55 - t25) * t141 + (t146 * t54 - t24) * t140) * t228 + (-t140 * t55 - t141 * t54) * t229 + m(6) * (-t11 * t58 - t12 * t59 + t28 * t34 + t29 * t35) + t151; ((-t146 * t57 - t27) * t141 + (t146 * t56 - t26) * t140) * t228 + m(6) * (-t13 * t58 - t14 * t59 + t28 * t40 + t29 * t41) + (-t140 * t57 - t141 * t56) * t229 + t151; t204 * t128 * t102 * t234 + (t1 * t15 - t28 * t59 - t29 * t58) * t233 + (t246 * t213 + t241 * t248 + t77 * t247 + (-t260 * t141 - t250) * t211) * t141 + (t75 * t247 + t242 * t249 + t243 * t211 + (t241 * t140 + t242 * t141 + (t269 * t140 + t270 * t141) * t203 + (t271 * t140 + t272 * t141) * t202 + (t254 * t140 + t253 * t141) * qJD(4) + (t243 + t246 + t262 - t277) * t146) * t141 + (t250 + t263) * t213) * t140; t200; m(6) * ((t140 * t35 + t141 * t34) * t202 + ((t146 * t35 + t12) * t141 + (-t146 * t34 + t11) * t140) * t147); m(6) * ((t140 * t41 + t141 * t40) * t202 + ((t146 * t41 + t14) * t141 + (-t146 * t40 + t13) * t140) * t147); m(6) * ((-t1 + (-t140 * t58 - t141 * t59) * qJD(4)) * t148 + (qJD(4) * t15 + (-t146 * t58 + t28) * t141 + (t29 + t221) * t140) * t147); 0.2e1 * (-0.1e1 + t204) * t148 * t200;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
