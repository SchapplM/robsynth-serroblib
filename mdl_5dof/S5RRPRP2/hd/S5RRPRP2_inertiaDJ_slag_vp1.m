% Calculate time derivative of joint inertia matrix for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:33
% DurationCPUTime: 3.60s
% Computational Cost: add. (5592->225), mult. (4404->311), div. (0->0), fcn. (3233->8), ass. (0->150)
t307 = Icges(5,4) - Icges(6,5);
t306 = Icges(5,1) + Icges(6,1);
t305 = Icges(5,2) + Icges(6,3);
t160 = cos(qJ(4));
t304 = t307 * t160;
t158 = sin(qJ(4));
t303 = t307 * t158;
t302 = Icges(6,4) + Icges(5,5);
t301 = Icges(5,6) - Icges(6,6);
t300 = t305 * t158 - t304;
t299 = t306 * t160 - t303;
t298 = -t305 * t160 - t303;
t297 = t306 * t158 + t304;
t157 = qJ(1) + qJ(2);
t152 = pkin(8) + t157;
t150 = cos(t152);
t296 = t300 * t150;
t295 = t299 * t150;
t290 = -t301 * t158 + t302 * t160;
t294 = t299 * t158 - t300 * t160;
t149 = sin(t152);
t288 = -t301 * t149 + t296;
t287 = -t300 * t149 - t301 * t150;
t286 = t299 * t149 - t302 * t150;
t285 = t302 * t149 + t295;
t276 = t287 * t158 - t286 * t160;
t293 = t276 * t149;
t292 = Icges(6,2) + Icges(5,3);
t273 = t302 * t158 + t301 * t160;
t256 = rSges(6,3) + qJ(5);
t274 = rSges(6,1) + pkin(4);
t268 = t158 * t256 + t160 * t274;
t284 = t290 * t150;
t271 = t298 * t158 + t297 * t160;
t218 = qJD(4) * t158;
t209 = t149 * t218;
t217 = qJD(4) * t160;
t275 = qJD(5) * t158 + t256 * t217;
t282 = -t275 * t149 + t274 * t209;
t281 = t290 * t149 - t292 * t150;
t280 = t292 * t149 + t284;
t156 = qJD(1) + qJD(2);
t225 = t150 * t156;
t254 = -t288 * t158 - t285 * t160;
t279 = t254 * t149;
t278 = t254 * t150;
t277 = -t273 * qJD(4) + t292 * t156;
t153 = sin(t157);
t154 = cos(t157);
t102 = t154 * rSges(3,1) - rSges(3,2) * t153;
t270 = -t285 * t158 + t288 * t160;
t269 = -t286 * t158 - t287 * t160;
t267 = -qJD(4) * t290 + t271 * t156;
t266 = t281 * t149 - t279 + (-t276 - t280) * t150;
t265 = t149 ^ 2;
t264 = t150 ^ 2;
t138 = t149 * rSges(5,3);
t235 = rSges(5,2) * t158;
t123 = t149 * t235;
t137 = rSges(5,1) * t158 + rSges(5,2) * t160;
t167 = -qJD(4) * t137 * t150 + rSges(5,3) * t225 + t156 * t123;
t224 = t150 * t158;
t214 = rSges(5,2) * t224;
t211 = -t149 * rSges(5,2) * t217 - rSges(5,1) * t209 - t156 * t214;
t221 = t150 * rSges(5,3) + t123;
t237 = rSges(5,1) * t160;
t77 = t149 * t237 - t221;
t223 = t150 * t160;
t79 = rSges(5,1) * t223 + t138 - t214;
t2 = ((-t79 + t138) * t156 + t211) * t149 + (t156 * t77 + t167) * t150;
t248 = 2 * m(5);
t263 = t2 * t248;
t262 = t281 * t150 + t293;
t259 = t280 * t149 - t278;
t226 = t149 * t156;
t258 = t277 * t150 - t226 * t290;
t257 = -t277 * t149 - t284 * t156;
t151 = pkin(2) * t154;
t255 = -t150 * rSges(4,1) - t151;
t250 = 2 * m(3);
t249 = 2 * m(4);
t247 = 2 * m(6);
t110 = (-t235 + t237) * qJD(4);
t243 = m(5) * t110;
t242 = m(5) * t137;
t159 = sin(qJ(1));
t241 = pkin(1) * t159;
t240 = pkin(2) * t153;
t141 = t150 * rSges(6,2);
t239 = t268 * t149 - t141;
t139 = t149 * rSges(6,2);
t238 = t274 * t223 + t256 * t224 + t139;
t234 = pkin(1) * qJD(1);
t220 = t274 * t158 - t256 * t160;
t61 = t220 * t150;
t233 = t156 * t61;
t231 = -t268 * qJD(4) + qJD(5) * t160;
t219 = t264 + t265;
t215 = m(6) * t218;
t213 = t159 * t234;
t161 = cos(qJ(1));
t212 = t161 * t234;
t210 = t150 * pkin(3) + t149 * pkin(7) + t151;
t205 = -pkin(3) - t237;
t89 = t102 * t156;
t85 = -rSges(4,2) * t149 - t255;
t101 = -rSges(3,1) * t153 - rSges(3,2) * t154;
t88 = t101 * t156;
t84 = -rSges(4,1) * t149 - rSges(4,2) * t150 - t240;
t176 = t149 * t205 - t240;
t63 = rSges(4,2) * t226 + t255 * t156;
t37 = t210 + t238;
t175 = -pkin(3) - t268;
t57 = t79 + t210;
t62 = t84 * t156;
t174 = rSges(6,2) * t225 + (-t218 * t274 + t275) * t150;
t144 = t150 * pkin(7);
t56 = t144 + t176 + t221;
t166 = t294 * qJD(4) + t297 * t217 + t298 * t218;
t165 = (-t254 * qJD(4) - t267 * t149 - t294 * t226) * t149 / 0.2e1 - (-qJD(4) * t276 + t267 * t150) * t150 / 0.2e1 + (t271 * t149 - t273 * t150 - t269) * t226 / 0.2e1 + (t273 * t149 + t271 * t150 - t270) * t225 / 0.2e1 - (t295 * t158 - t296 * t160) * t225 / 0.2e1;
t164 = t149 * t175 - t240;
t36 = t141 + t144 + t164;
t27 = (-t151 + t205 * t150 + (-rSges(5,3) - pkin(7)) * t149) * t156 - t211;
t120 = pkin(7) * t225;
t26 = t156 * t176 + t120 + t167;
t13 = t156 * t164 + t120 + t174;
t14 = (-t151 + (-rSges(6,2) - pkin(7)) * t149 + t175 * t150) * t156 + t282;
t155 = t161 * pkin(1);
t91 = t102 + t155;
t90 = t101 - t241;
t83 = -t89 - t212;
t82 = t88 - t213;
t81 = t155 + t85;
t80 = t84 - t241;
t60 = t220 * t149;
t59 = t63 - t212;
t58 = t62 - t213;
t55 = t155 + t57;
t54 = t56 - t241;
t35 = t155 + t37;
t34 = t36 - t241;
t29 = t149 * t231 - t233;
t28 = t150 * t231 + t220 * t226;
t25 = t27 - t212;
t24 = t26 - t213;
t15 = t149 * t239 + t150 * t238;
t12 = t14 - t212;
t11 = t13 - t213;
t1 = (t156 * t239 + t174) * t150 + ((t139 - t238) * t156 - t282) * t149;
t3 = [(t11 * t35 + t12 * t34) * t247 + (t24 * t55 + t25 * t54) * t248 + (t58 * t81 + t59 * t80) * t249 + (t82 * t91 + t83 * t90) * t250 + t166; m(6) * (t11 * t37 + t12 * t36 + t13 * t35 + t14 * t34) + m(5) * (t24 * t57 + t25 * t56 + t26 * t55 + t27 * t54) + m(4) * (t58 * t85 + t59 * t84 + t62 * t81 + t63 * t80) + m(3) * (t101 * t83 + t102 * t82 + t88 * t91 - t89 * t90) + t166; (t13 * t37 + t14 * t36) * t247 + (t26 * t57 + t27 * t56) * t248 + (t62 * t85 + t63 * t84) * t249 + (-t101 * t89 + t102 * t88) * t250 + t166; 0; 0; 0; ((-t156 * t55 - t25) * t150 + (t156 * t54 - t24) * t149) * t242 + m(6) * (-t11 * t60 - t12 * t61 + t28 * t34 + t29 * t35) + (-t149 * t55 - t150 * t54) * t243 + t165; ((-t156 * t57 - t27) * t150 + (t156 * t56 - t26) * t149) * t242 + (-t149 * t57 - t150 * t56) * t243 + m(6) * (-t13 * t60 - t14 * t61 + t28 * t36 + t29 * t37) + t165; m(5) * t2 + m(6) * t1; t219 * t137 * t110 * t248 + (t1 * t15 - t28 * t61 - t29 * t60) * t247 + (t262 * t226 + t257 * t264 + t79 * t263 + (-t276 * t150 - t266) * t225) * t150 + (t77 * t263 + t258 * t265 + t259 * t225 + (t257 * t149 + t258 * t150 + (t285 * t149 + t286 * t150) * t218 + (-t288 * t149 + t287 * t150) * t217 + (t270 * t149 + t269 * t150) * qJD(4) + (t259 + t262 + t278 - t293) * t156) * t150 + (t266 + t279) * t226) * t149; m(6) * ((t149 * t35 + t150 * t34) * t217 + ((t156 * t35 + t12) * t150 + (-t156 * t34 + t11) * t149) * t158); m(6) * ((t149 * t37 + t150 * t36) * t217 + ((t156 * t37 + t14) * t150 + (-t156 * t36 + t13) * t149) * t158); t215; m(6) * ((-t1 + (-t149 * t60 - t150 * t61) * qJD(4)) * t160 + (qJD(4) * t15 + (-t156 * t60 + t28) * t150 + (t29 + t233) * t149) * t158); 0.2e1 * (-0.1e1 + t219) * t160 * t215;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
