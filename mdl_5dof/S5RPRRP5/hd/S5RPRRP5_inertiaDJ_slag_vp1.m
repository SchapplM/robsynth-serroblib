% Calculate time derivative of joint inertia matrix for
% S5RPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:48
% DurationCPUTime: 3.51s
% Computational Cost: add. (5414->211), mult. (4288->299), div. (0->0), fcn. (3165->8), ass. (0->135)
t297 = Icges(5,4) - Icges(6,5);
t296 = Icges(5,1) + Icges(6,1);
t295 = Icges(5,2) + Icges(6,3);
t150 = cos(qJ(4));
t294 = t297 * t150;
t148 = sin(qJ(4));
t293 = t297 * t148;
t292 = Icges(6,4) + Icges(5,5);
t291 = Icges(5,6) - Icges(6,6);
t290 = t295 * t148 - t294;
t289 = t296 * t150 - t293;
t288 = -t295 * t150 - t293;
t287 = t296 * t148 + t294;
t147 = qJ(1) + pkin(8);
t144 = qJ(3) + t147;
t141 = cos(t144);
t286 = t290 * t141;
t285 = t289 * t141;
t280 = -t291 * t148 + t292 * t150;
t284 = t289 * t148 - t290 * t150;
t140 = sin(t144);
t278 = -t291 * t140 + t286;
t277 = -t290 * t140 - t291 * t141;
t276 = t289 * t140 - t292 * t141;
t275 = t292 * t140 + t285;
t266 = t277 * t148 - t276 * t150;
t283 = t266 * t140;
t282 = Icges(6,2) + Icges(5,3);
t263 = t292 * t148 + t291 * t150;
t246 = rSges(6,3) + qJ(5);
t264 = rSges(6,1) + pkin(4);
t258 = t148 * t246 + t150 * t264;
t274 = t280 * t141;
t261 = t288 * t148 + t287 * t150;
t209 = qJD(4) * t148;
t203 = t140 * t209;
t207 = qJD(5) * t148;
t208 = qJD(4) * t150;
t272 = -(t246 * t208 + t207) * t140 + t264 * t203;
t271 = t280 * t140 - t282 * t141;
t270 = t282 * t140 + t274;
t146 = qJD(1) + qJD(3);
t218 = t141 * t146;
t245 = -t278 * t148 - t275 * t150;
t269 = t245 * t140;
t268 = t245 * t141;
t267 = -t263 * qJD(4) + t282 * t146;
t260 = -t275 * t148 + t278 * t150;
t259 = -t276 * t148 - t277 * t150;
t257 = -qJD(4) * t280 + t261 * t146;
t256 = t271 * t140 - t269 + (-t266 - t270) * t141;
t255 = t140 ^ 2;
t254 = t141 ^ 2;
t128 = rSges(5,1) * t148 + rSges(5,2) * t150;
t129 = t140 * rSges(5,3);
t217 = t141 * t148;
t205 = rSges(5,2) * t217;
t204 = -t140 * rSges(5,2) * t208 - rSges(5,1) * t203 - t146 * t205;
t229 = rSges(5,2) * t148;
t115 = t140 * t229;
t225 = rSges(5,3) * t218 + t146 * t115;
t214 = t141 * rSges(5,3) + t115;
t219 = t140 * t150;
t75 = rSges(5,1) * t219 - t214;
t216 = t141 * t150;
t77 = rSges(5,1) * t216 + t129 - t205;
t2 = ((-t77 + t129) * t146 + t204) * t140 + (-qJD(4) * t128 * t141 + t146 * t75 + t225) * t141;
t240 = 2 * m(5);
t253 = t2 * t240;
t252 = t271 * t141 + t283;
t249 = t270 * t140 - t268;
t220 = t140 * t146;
t248 = t267 * t141 - t220 * t280;
t247 = -t267 * t140 - t274 * t146;
t210 = pkin(2) * cos(t147) + cos(qJ(1)) * pkin(1);
t241 = 2 * m(4);
t239 = 2 * m(6);
t230 = rSges(5,1) * t150;
t102 = (-t229 + t230) * qJD(4);
t235 = m(5) * t102;
t234 = m(5) * t128;
t132 = t141 * rSges(6,2);
t232 = t258 * t140 - t132;
t130 = t140 * rSges(6,2);
t231 = t264 * t216 + t246 * t217 + t130;
t213 = t264 * t148 - t246 * t150;
t61 = t213 * t141;
t228 = t146 * t61;
t226 = -t258 * qJD(4) + qJD(5) * t150;
t212 = t141 * pkin(3) + t140 * pkin(7);
t211 = t254 + t255;
t206 = m(6) * t209;
t202 = t141 * t209;
t201 = t141 * t208;
t198 = -pkin(3) - t230;
t86 = t141 * rSges(4,1) - rSges(4,2) * t140;
t81 = -rSges(4,1) * t218 + rSges(4,2) * t220;
t197 = -pkin(2) * sin(t147) - sin(qJ(1)) * pkin(1);
t85 = -rSges(4,1) * t140 - rSges(4,2) * t141;
t173 = t197 * qJD(1);
t172 = t210 * qJD(1);
t80 = t85 * t146;
t41 = t212 + t231;
t57 = t77 + t212;
t135 = t141 * pkin(7);
t56 = t140 * t198 + t135 + t214;
t164 = -pkin(3) - t258;
t163 = rSges(6,2) * t218 + t141 * t207 + t246 * t201 - t202 * t264;
t162 = t164 * t140;
t155 = t284 * qJD(4) + t287 * t208 + t288 * t209;
t40 = t132 + t135 + t162;
t154 = (-t245 * qJD(4) - t257 * t140 - t284 * t220) * t140 / 0.2e1 - (-qJD(4) * t266 + t257 * t141) * t141 / 0.2e1 + (t261 * t140 - t263 * t141 - t259) * t220 / 0.2e1 + (t263 * t140 + t261 * t141 - t260) * t218 / 0.2e1 - (t285 * t148 - t286 * t150) * t218 / 0.2e1;
t27 = (t198 * t141 + (-rSges(5,3) - pkin(7)) * t140) * t146 - t204;
t112 = pkin(7) * t218;
t26 = -rSges(5,2) * t201 - pkin(3) * t220 + t112 + (-t146 * t219 - t202) * rSges(5,1) + t225;
t13 = t146 * t162 + t112 + t163;
t14 = ((-rSges(6,2) - pkin(7)) * t140 + t164 * t141) * t146 + t272;
t79 = t86 + t210;
t78 = t197 + t85;
t60 = t213 * t140;
t59 = -t172 + t81;
t58 = t80 + t173;
t55 = t57 + t210;
t54 = t197 + t56;
t35 = t41 + t210;
t34 = t197 + t40;
t29 = t140 * t226 - t228;
t28 = t141 * t226 + t213 * t220;
t25 = -t172 + t27;
t24 = t173 + t26;
t15 = t140 * t232 + t141 * t231;
t12 = -t172 + t14;
t11 = t173 + t13;
t1 = (t146 * t232 + t163) * t141 + ((t130 - t231) * t146 - t272) * t140;
t3 = [(t11 * t35 + t12 * t34) * t239 + (t24 * t55 + t25 * t54) * t240 + (t58 * t79 + t59 * t78) * t241 + t155; 0; 0; m(6) * (t11 * t41 + t12 * t40 + t13 * t35 + t14 * t34) + m(5) * (t24 * t57 + t25 * t56 + t26 * t55 + t27 * t54) + m(4) * (t58 * t86 + t59 * t85 + t78 * t81 + t79 * t80) + t155; 0; (t13 * t41 + t14 * t40) * t239 + (t26 * t57 + t27 * t56) * t240 + (t80 * t86 + t81 * t85) * t241 + t155; m(6) * (-t11 * t60 - t12 * t61 + t28 * t34 + t29 * t35) + ((-t146 * t55 - t25) * t141 + (t146 * t54 - t24) * t140) * t234 + (-t140 * t55 - t141 * t54) * t235 + t154; m(5) * t2 + m(6) * t1; ((-t146 * t57 - t27) * t141 + (t146 * t56 - t26) * t140) * t234 + m(6) * (-t13 * t60 - t14 * t61 + t28 * t40 + t29 * t41) + (-t140 * t57 - t141 * t56) * t235 + t154; t211 * t128 * t102 * t240 + (t1 * t15 - t28 * t61 - t29 * t60) * t239 + (t252 * t220 + t247 * t254 + t77 * t253 + (-t266 * t141 - t256) * t218) * t141 + (t75 * t253 + t248 * t255 + t249 * t218 + (t247 * t140 + t248 * t141 + (t275 * t140 + t276 * t141) * t209 + (-t278 * t140 + t277 * t141) * t208 + (t260 * t140 + t259 * t141) * qJD(4) + (t249 + t252 + t268 - t283) * t146) * t141 + (t256 + t269) * t220) * t140; m(6) * ((t140 * t35 + t141 * t34) * t208 + ((t146 * t35 + t12) * t141 + (-t146 * t34 + t11) * t140) * t148); t206; m(6) * ((t140 * t41 + t141 * t40) * t208 + ((t146 * t41 + t14) * t141 + (-t146 * t40 + t13) * t140) * t148); m(6) * ((-t1 + (-t140 * t60 - t141 * t61) * qJD(4)) * t150 + (qJD(4) * t15 + (-t146 * t60 + t28) * t141 + (t29 + t228) * t140) * t148); 0.2e1 * (-0.1e1 + t211) * t150 * t206;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
