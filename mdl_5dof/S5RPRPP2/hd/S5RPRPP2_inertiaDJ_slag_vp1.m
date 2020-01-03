% Calculate time derivative of joint inertia matrix for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:50
% DurationCPUTime: 4.58s
% Computational Cost: add. (3821->315), mult. (4866->420), div. (0->0), fcn. (3686->6), ass. (0->181)
t316 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t145 = sin(qJ(3));
t147 = cos(qJ(3));
t173 = Icges(6,5) * t147 + Icges(6,6) * t145;
t176 = Icges(4,5) * t147 - Icges(4,6) * t145;
t179 = Icges(5,4) * t147 + Icges(5,6) * t145;
t315 = t173 - t176 - t179;
t144 = qJ(1) + pkin(7);
t141 = sin(t144);
t295 = qJ(5) + rSges(6,3);
t210 = t295 * t141;
t142 = cos(t144);
t234 = t142 * t147;
t235 = t142 * t145;
t268 = rSges(6,1) + pkin(4);
t264 = rSges(6,2) * t235 + t268 * t234 - t210;
t246 = Icges(5,5) * t147;
t249 = Icges(6,4) * t147;
t252 = Icges(4,4) * t147;
t313 = -t246 - t249 + t252 + (Icges(4,1) + Icges(5,1) + Icges(6,1)) * t145;
t209 = t295 * t142;
t215 = t268 * t145;
t211 = -rSges(6,2) * t147 + t215;
t307 = t316 * t141 - t315 * t142;
t175 = Icges(5,3) * t145 + t246;
t61 = -Icges(5,6) * t142 + t141 * t175;
t178 = Icges(6,2) * t145 + t249;
t65 = Icges(6,6) * t142 + t141 * t178;
t312 = t61 + t65;
t62 = Icges(5,6) * t141 + t142 * t175;
t66 = -Icges(6,6) * t141 + t142 * t178;
t311 = t62 + t66;
t181 = -Icges(4,2) * t145 + t252;
t70 = Icges(4,6) * t141 + t142 * t181;
t253 = Icges(4,4) * t145;
t187 = Icges(4,1) * t147 - t253;
t76 = Icges(4,5) * t141 + t142 * t187;
t189 = t145 * t70 - t147 * t76;
t250 = Icges(6,4) * t145;
t183 = Icges(6,1) * t147 + t250;
t72 = -Icges(6,5) * t141 + t142 * t183;
t191 = t145 * t66 + t147 * t72;
t247 = Icges(5,5) * t145;
t185 = Icges(5,1) * t147 + t247;
t74 = Icges(5,4) * t141 + t142 * t185;
t193 = t145 * t62 + t147 * t74;
t276 = t189 - t191 - t193;
t310 = t276 * t142;
t308 = t315 * t141 + t316 * t142;
t71 = Icges(6,5) * t142 + t141 * t183;
t73 = -Icges(5,4) * t142 + t141 * t185;
t75 = -Icges(4,5) * t142 + t141 * t187;
t306 = t73 + t71 + t75;
t305 = t74 + t72 + t76;
t304 = ((Icges(4,6) - Icges(5,6) + Icges(6,6)) * t147 + (Icges(5,4) + Icges(4,5) - Icges(6,5)) * t145) * qJD(3);
t69 = -Icges(4,6) * t142 + t141 * t181;
t190 = t145 * t69 - t147 * t75;
t192 = t145 * t65 + t147 * t71;
t194 = t145 * t61 + t147 * t73;
t302 = t190 - t192 - t194;
t300 = t141 / 0.2e1;
t299 = -t142 / 0.2e1;
t297 = -qJD(1) / 0.2e1;
t296 = qJD(1) / 0.2e1;
t165 = t189 * t141;
t166 = t190 * t142;
t167 = t191 * t141;
t168 = t192 * t142;
t169 = t193 * t141;
t170 = t194 * t142;
t294 = -t308 * t141 - t307 * t142 - t165 - t166 + t167 + t168 + t169 + t170;
t293 = t307 * qJD(1);
t292 = t305 * t141 + t306 * t142;
t291 = (-t70 + t311) * t141 + (-t69 + t312) * t142;
t139 = t141 ^ 2;
t140 = t142 ^ 2;
t274 = 2 * m(4);
t122 = rSges(4,1) * t145 + rSges(4,2) * t147;
t164 = t122 * qJD(3);
t226 = qJD(1) * t145;
t214 = t141 * t226;
t227 = qJD(1) * t142;
t151 = rSges(4,2) * t214 + rSges(4,3) * t227 - t142 * t164;
t287 = -rSges(4,2) * t235 + t141 * rSges(4,3);
t262 = rSges(4,1) * t147;
t205 = -rSges(4,2) * t145 + t262;
t258 = rSges(4,3) * t142;
t80 = t141 * t205 - t258;
t83 = rSges(4,1) * t234 + t287;
t3 = (qJD(1) * t80 + t151) * t142 + (-t141 * t164 + (-t83 + t287) * qJD(1)) * t141;
t290 = t274 * t3;
t228 = qJD(1) * t141;
t286 = t302 * t141 - t308 * t142;
t283 = t307 * t141 - t310;
t282 = t308 * qJD(1) - t304 * t142;
t281 = t304 * t141 - t293;
t273 = 2 * m(5);
t272 = 2 * m(6);
t271 = m(5) / 0.2e1;
t270 = m(6) / 0.2e1;
t269 = -rSges(5,1) - pkin(3);
t267 = m(4) * t122;
t266 = sin(qJ(1)) * pkin(1);
t143 = cos(qJ(1)) * pkin(1);
t203 = rSges(6,1) * t147 + rSges(6,2) * t145;
t265 = t209 + (pkin(4) * t147 + t203) * t141;
t202 = pkin(3) * t147 + qJ(4) * t145;
t86 = t202 * t141;
t130 = pkin(3) * t234;
t87 = qJ(4) * t235 + t130;
t263 = t141 * t86 + t142 * t87;
t261 = rSges(5,1) * t145;
t260 = rSges(5,2) * t142;
t135 = t141 * rSges(5,2);
t257 = -rSges(6,2) - qJ(4);
t256 = -rSges(5,3) - qJ(4);
t204 = rSges(5,1) * t147 + rSges(5,3) * t145;
t88 = qJD(3) * t202 - qJD(4) * t147;
t254 = -t204 * qJD(3) - t88;
t222 = qJD(3) * t147;
t212 = t142 * t222;
t221 = qJD(4) * t145;
t233 = qJ(4) * t212 + t142 * t221;
t232 = rSges(5,2) * t227 + rSges(5,3) * t212;
t119 = pkin(3) * t145 - qJ(4) * t147;
t121 = -rSges(5,3) * t147 + t261;
t231 = -t119 - t121;
t229 = t139 + t140;
t225 = qJD(3) * t141;
t224 = qJD(3) * t142;
t223 = qJD(3) * t145;
t220 = -pkin(3) - t268;
t108 = t141 * pkin(3) * t223;
t213 = t142 * t223;
t219 = t141 * (qJD(1) * t130 + t141 * t221 - t108 + (t141 * t222 + t142 * t226) * qJ(4)) + t142 * (-qJ(4) * t214 + (-t147 * t228 - t213) * pkin(3) + t233) + t86 * t227;
t133 = pkin(6) * t227;
t217 = t133 + t233;
t82 = rSges(5,1) * t234 + rSges(5,3) * t235 + t135;
t216 = t142 * pkin(2) + t141 * pkin(6) + t143;
t58 = t231 * t142;
t208 = (-t173 / 0.2e1 + t179 / 0.2e1 + t176 / 0.2e1) * qJD(3);
t207 = -t119 - t211;
t137 = t142 * pkin(6);
t152 = t145 * t257 + t147 * t220 - pkin(2);
t149 = t141 * t152 - t209 - t266;
t24 = t137 + t149;
t188 = t216 + t87;
t25 = t188 + t264;
t201 = t141 * t25 + t142 * t24;
t180 = Icges(4,2) * t147 + t253;
t172 = -pkin(4) * t222 - t203 * qJD(3) - t88;
t171 = -pkin(2) - t205;
t55 = t207 * t142;
t160 = qJD(3) * t180;
t153 = t145 * t256 + t147 * t269 - pkin(2);
t150 = t141 * t153 - t266;
t106 = rSges(6,2) * t212;
t103 = t205 * qJD(3);
t91 = t119 * t228;
t79 = t141 * t204 - t260;
t57 = t231 * t141;
t54 = t207 * t141;
t53 = t83 + t216;
t52 = t141 * t171 + t137 + t258 - t266;
t31 = t188 + t82;
t30 = t137 + t150 + t260;
t29 = t122 * t225 + (-t143 + (-rSges(4,3) - pkin(6)) * t141 + t171 * t142) * qJD(1);
t28 = t133 + (-t266 + (-pkin(2) - t262) * t141) * qJD(1) + t151;
t27 = qJD(1) * t58 + t141 * t254;
t26 = t121 * t228 + t142 * t254 + t91;
t23 = qJD(1) * t55 + t141 * t172;
t22 = t142 * t172 + t211 * t228 + t91;
t9 = t141 * t79 + t142 * t82 + t263;
t8 = t108 + (-t221 + (t147 * t256 + t261) * qJD(3)) * t141 + (-t143 + (-rSges(5,2) - pkin(6)) * t141 + t153 * t142) * qJD(1);
t7 = qJD(1) * t150 + t213 * t269 + t217 + t232;
t6 = t141 * t265 + t142 * t264 + t263;
t5 = -qJD(5) * t142 + t108 + (-t221 + (t147 * t257 + t215) * qJD(3)) * t141 + (-t143 + (-pkin(6) + t295) * t141 + t152 * t142) * qJD(1);
t4 = qJD(1) * t149 - qJD(5) * t141 + t213 * t220 + t106 + t217;
t2 = t142 * t232 + (-t121 * t139 - t140 * t261) * qJD(3) + (t142 * t79 + (-t82 - t87 + t135) * t141) * qJD(1) + t219;
t1 = t142 * t106 + (-t211 * t139 - t140 * t215) * qJD(3) + ((-t209 + t265) * t142 + (-t87 - t210 - t264) * t141) * qJD(1) + t219;
t10 = [(t24 * t5 + t25 * t4) * t272 + (t30 * t8 + t31 * t7) * t273 + (t28 * t53 + t29 * t52) * t274 + (-t180 + t183 + t185 + t187 + t247 + t250 + (-Icges(6,2) - Icges(5,3)) * t147) * t223 + (t181 - t175 - t178 + t313) * t222; 0; 0; m(5) * (t26 * t30 + t27 * t31 + t57 * t7 + t58 * t8) + m(6) * (t22 * t24 + t23 * t25 + t4 * t54 + t5 * t55) + (m(4) * (-t103 * t52 - t122 * t29) + t208 * t142 + (t160 * t300 + t311 * t296 + t70 * t297) * t147) * t142 + (m(4) * (-t103 * t53 - t122 * t28) + t208 * t141 + (t160 * t299 + t312 * t296 + t69 * t297) * t147) * t141 + ((t306 * t141 + t305 * t142) * t297 + (t299 * t141 + t300 * t142) * t313 * qJD(3)) * t145 + (-t170 / 0.2e1 + t166 / 0.2e1 - t168 / 0.2e1 + t167 / 0.2e1 + t169 / 0.2e1 - t165 / 0.2e1) * qJD(3) + ((-t53 * t267 + (t70 / 0.2e1 - t66 / 0.2e1 - t62 / 0.2e1) * t147 + (t76 / 0.2e1 + t72 / 0.2e1 + t74 / 0.2e1) * t145) * t142 + (t52 * t267 + (t69 / 0.2e1 - t65 / 0.2e1 - t61 / 0.2e1) * t147 + (t75 / 0.2e1 + t71 / 0.2e1 + t73 / 0.2e1) * t145) * t141) * qJD(1); m(4) * t3 + m(5) * t2 + m(6) * t1; (t2 * t9 + t26 * t58 + t27 * t57) * t273 + (t1 * t6 + t22 * t55 + t23 * t54) * t272 + t229 * t122 * t103 * t274 + (t281 * t140 + t286 * t228 + t83 * t290 + (-t142 * t302 - t294) * t227) * t142 + (t80 * t290 + t282 * t139 + t283 * t227 + ((t281 - t293) * t141 + t282 * t142 + t292 * t223 - t291 * t222 + (-t292 * t145 + t291 * t147) * qJD(3) + ((-t302 + t307) * t141 + t310 + t283 + t286) * qJD(1)) * t142 + (t141 * t276 + t294) * t228) * t141; 0.2e1 * (t201 * t270 + (t141 * t31 + t142 * t30) * t271) * t222 + 0.2e1 * ((t141 * t4 + t142 * t5 + t227 * t25 - t228 * t24) * t270 + (t141 * t7 + t142 * t8 + t227 * t31 - t228 * t30) * t271) * t145; (m(5) + m(6)) * t223; 0.2e1 * ((t224 * t58 + t225 * t57 - t2) * t271 + (t224 * t55 + t225 * t54 - t1) * t270) * t147 + 0.2e1 * ((qJD(3) * t9 + t141 * t27 + t142 * t26 + t227 * t57 - t228 * t58) * t271 + (qJD(3) * t6 + t141 * t23 + t142 * t22 + t227 * t54 - t228 * t55) * t270) * t145; 0.4e1 * (t271 + t270) * (-0.1e1 + t229) * t145 * t222; m(6) * (-qJD(1) * t201 - t141 * t5 + t142 * t4); 0; m(6) * (-t141 * t22 + t142 * t23 + (-t141 * t54 - t142 * t55) * qJD(1)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
