% Calculate time derivative of joint inertia matrix for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:40:54
% DurationCPUTime: 4.41s
% Computational Cost: add. (3769->308), mult. (4782->419), div. (0->0), fcn. (3632->4), ass. (0->179)
t310 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t144 = sin(qJ(3));
t145 = cos(qJ(3));
t169 = Icges(6,5) * t145 + Icges(6,6) * t144;
t172 = Icges(4,5) * t145 - Icges(4,6) * t144;
t175 = Icges(5,4) * t145 + Icges(5,6) * t144;
t309 = t169 - t172 - t175;
t143 = pkin(7) + qJ(2);
t141 = sin(t143);
t289 = qJ(5) + rSges(6,3);
t205 = t289 * t141;
t142 = cos(t143);
t230 = t142 * t145;
t231 = t142 * t144;
t262 = rSges(6,1) + pkin(4);
t259 = rSges(6,2) * t231 + t262 * t230 - t205;
t242 = Icges(5,5) * t145;
t245 = Icges(6,4) * t145;
t248 = Icges(4,4) * t145;
t307 = -t242 - t245 + t248 + (Icges(4,1) + Icges(5,1) + Icges(6,1)) * t144;
t204 = t289 * t142;
t210 = t262 * t144;
t206 = -rSges(6,2) * t145 + t210;
t301 = -t310 * t141 + t309 * t142;
t171 = Icges(5,3) * t144 + t242;
t61 = -Icges(5,6) * t142 + t141 * t171;
t174 = Icges(6,2) * t144 + t245;
t65 = Icges(6,6) * t142 + t141 * t174;
t306 = t61 + t65;
t62 = Icges(5,6) * t141 + t142 * t171;
t66 = -Icges(6,6) * t141 + t142 * t174;
t305 = t62 + t66;
t177 = -Icges(4,2) * t144 + t248;
t70 = Icges(4,6) * t141 + t142 * t177;
t249 = Icges(4,4) * t144;
t183 = Icges(4,1) * t145 - t249;
t76 = Icges(4,5) * t141 + t142 * t183;
t184 = t144 * t70 - t145 * t76;
t246 = Icges(6,4) * t144;
t179 = Icges(6,1) * t145 + t246;
t72 = -Icges(6,5) * t141 + t142 * t179;
t186 = t144 * t66 + t145 * t72;
t243 = Icges(5,5) * t144;
t181 = Icges(5,1) * t145 + t243;
t74 = Icges(5,4) * t141 + t142 * t181;
t188 = t144 * t62 + t145 * t74;
t270 = t184 - t186 - t188;
t304 = t270 * t142;
t302 = t309 * t141 + t310 * t142;
t71 = Icges(6,5) * t142 + t141 * t179;
t73 = -Icges(5,4) * t142 + t141 * t181;
t75 = -Icges(4,5) * t142 + t141 * t183;
t300 = t71 + t75 + t73;
t299 = t72 + t76 + t74;
t298 = ((Icges(4,6) - Icges(5,6) + Icges(6,6)) * t145 + (Icges(5,4) + Icges(4,5) - Icges(6,5)) * t144) * qJD(3);
t69 = -Icges(4,6) * t142 + t141 * t177;
t185 = t144 * t69 - t145 * t75;
t187 = t144 * t65 + t145 * t71;
t189 = t144 * t61 + t145 * t73;
t296 = t185 - t187 - t189;
t294 = t141 / 0.2e1;
t293 = -t142 / 0.2e1;
t291 = -qJD(2) / 0.2e1;
t290 = qJD(2) / 0.2e1;
t161 = t184 * t141;
t162 = t185 * t142;
t163 = t186 * t141;
t164 = t187 * t142;
t165 = t188 * t141;
t166 = t189 * t142;
t288 = -t302 * t141 + t301 * t142 - t161 - t162 + t163 + t164 + t165 + t166;
t287 = t301 * qJD(2);
t286 = t299 * t141 + t300 * t142;
t285 = (t70 - t305) * t141 + (t69 - t306) * t142;
t139 = t141 ^ 2;
t140 = t142 ^ 2;
t268 = 2 * m(4);
t122 = t144 * rSges(4,1) + rSges(4,2) * t145;
t160 = qJD(3) * t122;
t220 = qJD(2) * t144;
t209 = t141 * t220;
t221 = qJD(2) * t142;
t227 = rSges(4,2) * t209 + rSges(4,3) * t221;
t281 = -rSges(4,2) * t231 + t141 * rSges(4,3);
t200 = rSges(4,1) * t145 - rSges(4,2) * t144;
t255 = t142 * rSges(4,3);
t80 = t141 * t200 - t255;
t83 = rSges(4,1) * t230 + t281;
t3 = (qJD(2) * t80 - t142 * t160 + t227) * t142 + (-t141 * t160 + (-t83 + t281) * qJD(2)) * t141;
t284 = t268 * t3;
t222 = qJD(2) * t141;
t280 = t296 * t141 - t302 * t142;
t277 = -t301 * t141 - t304;
t276 = t302 * qJD(2) - t298 * t142;
t275 = t298 * t141 + t287;
t267 = 2 * m(5);
t266 = 2 * m(6);
t265 = m(5) / 0.2e1;
t264 = m(6) / 0.2e1;
t263 = -rSges(5,1) - pkin(3);
t261 = m(4) * t122;
t198 = rSges(6,1) * t145 + rSges(6,2) * t144;
t260 = t204 + (pkin(4) * t145 + t198) * t141;
t197 = pkin(3) * t145 + qJ(4) * t144;
t86 = t197 * t141;
t130 = pkin(3) * t230;
t87 = qJ(4) * t231 + t130;
t258 = t141 * t86 + t142 * t87;
t135 = t141 * rSges(5,2);
t256 = t142 * rSges(5,2);
t254 = t144 * rSges(5,1);
t253 = -rSges(6,2) - qJ(4);
t252 = -rSges(5,3) - qJ(4);
t199 = rSges(5,1) * t145 + rSges(5,3) * t144;
t88 = qJD(3) * t197 - qJD(4) * t145;
t250 = -t199 * qJD(3) - t88;
t216 = qJD(3) * t145;
t207 = t142 * t216;
t215 = qJD(4) * t144;
t229 = qJ(4) * t207 + t142 * t215;
t228 = rSges(5,2) * t221 + rSges(5,3) * t207;
t119 = t144 * pkin(3) - qJ(4) * t145;
t121 = -rSges(5,3) * t145 + t254;
t226 = -t119 - t121;
t224 = t142 * pkin(2) + t141 * pkin(6);
t223 = t139 + t140;
t219 = qJD(3) * t141;
t218 = qJD(3) * t142;
t217 = qJD(3) * t144;
t214 = -pkin(3) - t262;
t108 = t141 * pkin(3) * t217;
t208 = t142 * t217;
t150 = -t145 * t222 - t208;
t213 = t141 * (qJD(2) * t130 + t141 * t215 - t108 + (t141 * t216 + t142 * t220) * qJ(4)) + t142 * (pkin(3) * t150 - qJ(4) * t209 + t229) + t86 * t221;
t133 = pkin(6) * t221;
t211 = t133 + t229;
t82 = rSges(5,1) * t230 + rSges(5,3) * t231 + t135;
t58 = t226 * t142;
t203 = t224 + t87;
t202 = (-t169 / 0.2e1 + t175 / 0.2e1 + t172 / 0.2e1) * qJD(3);
t201 = -t119 - t206;
t137 = t142 * pkin(6);
t148 = t144 * t253 + t145 * t214 - pkin(2);
t146 = t141 * t148 - t204;
t26 = t137 + t146;
t27 = t203 + t259;
t196 = t141 * t27 + t142 * t26;
t176 = Icges(4,2) * t145 + t249;
t168 = -pkin(4) * t216 - t198 * qJD(3) - t88;
t167 = -pkin(2) - t200;
t53 = t201 * t142;
t156 = qJD(3) * t176;
t149 = t144 * t252 + t145 * t263 - pkin(2);
t147 = t149 * t141;
t106 = rSges(6,2) * t207;
t103 = t200 * qJD(3);
t91 = t119 * t222;
t79 = t141 * t199 - t256;
t57 = t226 * t141;
t55 = t83 + t224;
t54 = t141 * t167 + t137 + t255;
t52 = t201 * t141;
t33 = t203 + t82;
t32 = t137 + t147 + t256;
t29 = t122 * t219 + ((-rSges(4,3) - pkin(6)) * t141 + t167 * t142) * qJD(2);
t28 = rSges(4,1) * t150 - rSges(4,2) * t207 - pkin(2) * t222 + t133 + t227;
t25 = qJD(2) * t58 + t141 * t250;
t24 = t121 * t222 + t142 * t250 + t91;
t23 = qJD(2) * t53 + t141 * t168;
t22 = t142 * t168 + t206 * t222 + t91;
t9 = t141 * t79 + t142 * t82 + t258;
t8 = t108 + (-t215 + (t145 * t252 + t254) * qJD(3)) * t141 + ((-rSges(5,2) - pkin(6)) * t141 + t149 * t142) * qJD(2);
t7 = qJD(2) * t147 + t208 * t263 + t211 + t228;
t6 = -qJD(5) * t142 + t108 + (-t215 + (t145 * t253 + t210) * qJD(3)) * t141 + ((-pkin(6) + t289) * t141 + t148 * t142) * qJD(2);
t5 = qJD(2) * t146 - qJD(5) * t141 + t208 * t214 + t106 + t211;
t4 = t141 * t260 + t142 * t259 + t258;
t2 = t142 * t228 + (-t121 * t139 - t140 * t254) * qJD(3) + (t142 * t79 + (-t82 - t87 + t135) * t141) * qJD(2) + t213;
t1 = t142 * t106 + (-t206 * t139 - t140 * t210) * qJD(3) + ((-t204 + t260) * t142 + (-t87 - t205 - t259) * t141) * qJD(2) + t213;
t10 = [0; 0; (t32 * t8 + t33 * t7) * t267 + (t26 * t6 + t27 * t5) * t266 + (t28 * t55 + t29 * t54) * t268 + (-t176 + t179 + t181 + t183 + t243 + t246 + (-Icges(6,2) - Icges(5,3)) * t145) * t217 + (t177 - t171 - t174 + t307) * t216; m(4) * t3 + m(5) * t2 + m(6) * t1; m(5) * (t24 * t32 + t25 * t33 + t57 * t7 + t58 * t8) + m(6) * (t22 * t26 + t23 * t27 + t5 * t52 + t53 * t6) + (m(4) * (-t103 * t54 - t122 * t29) + t202 * t142 + (t156 * t294 + t305 * t290 + t70 * t291) * t145) * t142 + (m(4) * (-t103 * t55 - t122 * t28) + t202 * t141 + (t156 * t293 + t306 * t290 + t69 * t291) * t145) * t141 + ((t300 * t141 + t299 * t142) * t291 + (t141 * t293 + t142 * t294) * t307 * qJD(3)) * t144 + (-t164 / 0.2e1 - t166 / 0.2e1 + t162 / 0.2e1 + t163 / 0.2e1 + t165 / 0.2e1 - t161 / 0.2e1) * qJD(3) + ((-t55 * t261 + (t70 / 0.2e1 - t66 / 0.2e1 - t62 / 0.2e1) * t145 + (t76 / 0.2e1 + t72 / 0.2e1 + t74 / 0.2e1) * t144) * t142 + (t54 * t261 + (t69 / 0.2e1 - t65 / 0.2e1 - t61 / 0.2e1) * t145 + (t75 / 0.2e1 + t71 / 0.2e1 + t73 / 0.2e1) * t144) * t141) * qJD(2); (t1 * t4 + t22 * t53 + t23 * t52) * t266 + (t2 * t9 + t24 * t58 + t25 * t57) * t267 + t223 * t122 * t103 * t268 + (t275 * t140 + t280 * t222 + t83 * t284 + (-t142 * t296 - t288) * t221) * t142 + (t80 * t284 + t276 * t139 + t277 * t221 + ((t275 + t287) * t141 + t276 * t142 + t286 * t217 + t285 * t216 + (-t286 * t144 - t285 * t145) * qJD(3) + ((-t296 - t301) * t141 + t304 + t277 + t280) * qJD(2)) * t142 + (t270 * t141 + t288) * t222) * t141; (m(5) + m(6)) * t217; 0.2e1 * ((t141 * t33 + t142 * t32) * t265 + t196 * t264) * t216 + 0.2e1 * ((t141 * t7 + t142 * t8 + t221 * t33 - t222 * t32) * t265 + (t141 * t5 + t142 * t6 + t221 * t27 - t222 * t26) * t264) * t144; 0.2e1 * ((t218 * t53 + t219 * t52 - t1) * t264 + (t218 * t58 + t219 * t57 - t2) * t265) * t145 + 0.2e1 * ((qJD(3) * t4 + t141 * t23 + t142 * t22 + t221 * t52 - t222 * t53) * t264 + (qJD(3) * t9 + t141 * t25 + t142 * t24 + t221 * t57 - t222 * t58) * t265) * t144; 0.4e1 * (t265 + t264) * (-0.1e1 + t223) * t144 * t216; 0; m(6) * (-qJD(2) * t196 - t141 * t6 + t142 * t5); m(6) * (-t141 * t22 + t142 * t23 + (-t141 * t52 - t142 * t53) * qJD(2)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
