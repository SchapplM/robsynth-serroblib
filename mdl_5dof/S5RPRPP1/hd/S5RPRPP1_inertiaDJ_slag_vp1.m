% Calculate time derivative of joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:08:47
% DurationCPUTime: 5.09s
% Computational Cost: add. (4752->308), mult. (4568->413), div. (0->0), fcn. (3453->8), ass. (0->190)
t326 = Icges(6,4) + Icges(5,5);
t325 = -Icges(5,6) + Icges(6,6);
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t180 = Icges(4,5) * t151 - Icges(4,6) * t149;
t146 = qJ(3) + pkin(8);
t141 = sin(t146);
t143 = cos(t146);
t322 = t325 * t141 + t326 * t143;
t324 = t180 + t322;
t323 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t252 = Icges(6,5) * t143;
t257 = Icges(5,4) * t143;
t321 = -t252 + t257 + (Icges(5,1) + Icges(6,1)) * t141;
t261 = rSges(6,3) + qJ(5);
t280 = rSges(6,1) + pkin(4);
t307 = t141 * t261 + t143 * t280;
t295 = t141 * t280;
t318 = t143 * t261;
t239 = t295 - t318;
t147 = qJ(1) + pkin(7);
t142 = sin(t147);
t144 = cos(t147);
t314 = t323 * t142 + t324 * t144;
t253 = Icges(6,5) * t141;
t187 = Icges(6,1) * t143 + t253;
t68 = -Icges(6,4) * t144 + t142 * t187;
t258 = Icges(5,4) * t141;
t189 = Icges(5,1) * t143 - t258;
t70 = -Icges(5,5) * t144 + t142 * t189;
t320 = t68 + t70;
t69 = Icges(6,4) * t142 + t144 * t187;
t71 = Icges(5,5) * t142 + t144 * t189;
t319 = t69 + t71;
t259 = Icges(4,4) * t151;
t185 = -Icges(4,2) * t149 + t259;
t81 = Icges(4,6) * t142 + t144 * t185;
t260 = Icges(4,4) * t149;
t191 = Icges(4,1) * t151 - t260;
t83 = Icges(4,5) * t142 + t144 * t191;
t192 = t149 * t81 - t151 * t83;
t183 = -Icges(5,2) * t141 + t257;
t67 = Icges(5,6) * t142 + t144 * t183;
t198 = t141 * t67 - t143 * t71;
t178 = Icges(6,3) * t141 + t252;
t61 = Icges(6,6) * t142 + t144 * t178;
t200 = t141 * t61 + t143 * t69;
t285 = t192 + t198 - t200;
t317 = t285 * t144;
t316 = t321 * qJD(3);
t315 = t324 * t142 - t323 * t144;
t313 = (-Icges(4,5) * t149 - Icges(4,6) * t151 - t141 * t326 + t325 * t143) * qJD(3);
t80 = -Icges(4,6) * t144 + t142 * t185;
t82 = -Icges(4,5) * t144 + t142 * t191;
t193 = t149 * t80 - t151 * t82;
t66 = -Icges(5,6) * t144 + t142 * t183;
t199 = t141 * t66 - t143 * t70;
t60 = -Icges(6,6) * t144 + t142 * t178;
t201 = t141 * t60 + t143 * t68;
t312 = t193 + t199 - t201;
t311 = t142 / 0.2e1;
t310 = -t144 / 0.2e1;
t309 = -qJD(1) / 0.2e1;
t308 = qJD(1) / 0.2e1;
t168 = t192 * t142;
t169 = t193 * t144;
t170 = t198 * t142;
t171 = t199 * t144;
t172 = t200 * t142;
t173 = t201 * t144;
t306 = t315 * t142 - t314 * t144 - t168 - t169 - t170 - t171 + t172 + t173;
t305 = t314 * qJD(1);
t304 = t319 * t142 + t320 * t144;
t303 = (t61 - t67) * t142 + (t60 - t66) * t144;
t139 = t142 ^ 2;
t140 = t144 ^ 2;
t105 = rSges(5,1) * t141 + rSges(5,2) * t143;
t166 = qJD(3) * t105;
t132 = qJD(4) * t142;
t148 = -qJ(4) - pkin(6);
t228 = qJD(3) * t149;
t221 = pkin(3) * t228;
t235 = qJD(1) * t142;
t218 = qJD(4) * t144 + t142 * t221 + t148 * t235;
t234 = qJD(1) * t144;
t274 = -pkin(6) - t148;
t138 = pkin(3) * t151 + pkin(2);
t275 = pkin(2) - t138;
t276 = pkin(6) * t142;
t294 = t142 * t275;
t137 = t144 * pkin(6);
t240 = t144 * t148;
t58 = t137 + t240 - t294;
t225 = t142 * ((-t144 * t275 - t276) * qJD(1) - t218) + t144 * (-t144 * t221 + t132 + (t144 * t274 + t294) * qJD(1)) + t58 * t234;
t267 = rSges(5,2) * t141;
t238 = rSges(5,3) * t234 + t235 * t267;
t134 = t142 * rSges(5,3);
t242 = t141 * t144;
t293 = -rSges(5,2) * t242 + t134;
t119 = t144 * t138;
t59 = -pkin(2) * t144 + t142 * t274 + t119;
t269 = rSges(5,1) * t143;
t208 = -t267 + t269;
t73 = -rSges(5,3) * t144 + t142 * t208;
t241 = t143 * t144;
t75 = rSges(5,1) * t241 + t293;
t2 = (qJD(1) * t73 - t144 * t166 + t238) * t144 + (-t142 * t166 + (-t59 - t75 + t293) * qJD(1)) * t142 + t225;
t282 = 2 * m(5);
t302 = t2 * t282;
t283 = 2 * m(4);
t135 = t142 * rSges(4,3);
t127 = rSges(4,1) * t149 + rSges(4,2) * t151;
t167 = qJD(3) * t127;
t233 = qJD(1) * t149;
t217 = t142 * t233;
t154 = rSges(4,2) * t217 + rSges(4,3) * t234 - t144 * t167;
t268 = rSges(4,2) * t149;
t223 = t144 * t268;
t270 = rSges(4,1) * t151;
t209 = -t268 + t270;
t266 = rSges(4,3) * t144;
t84 = t142 * t209 - t266;
t237 = t144 * t270 + t135;
t85 = -t223 + t237;
t4 = (qJD(1) * t84 + t154) * t144 + (-t142 * t167 + (-t223 - t85 + t135) * qJD(1)) * t142;
t301 = t283 * t4;
t278 = sin(qJ(1)) * pkin(1);
t292 = t137 - t278;
t291 = t312 * t142 + t315 * t144;
t288 = t314 * t142 - t317;
t287 = -t315 * qJD(1) + t313 * t144;
t286 = -t313 * t142 - t305;
t281 = 2 * m(6);
t279 = m(4) * t127;
t277 = pkin(3) * t149;
t145 = cos(qJ(1)) * pkin(1);
t273 = t142 * t58 + t144 * t59;
t272 = -rSges(6,2) * t144 + t307 * t142;
t136 = t142 * rSges(6,2);
t271 = t280 * t241 + t261 * t242 + t136;
t236 = t139 + t140;
t232 = qJD(3) * t141;
t231 = qJD(3) * t142;
t230 = qJD(3) * t143;
t229 = qJD(3) * t144;
t227 = qJD(3) * t151;
t226 = qJD(5) * t141;
t224 = m(6) * t232;
t220 = pkin(3) * t227;
t219 = t322 * qJD(3) / 0.2e1;
t215 = t280 * qJD(3);
t214 = -t105 - t277;
t213 = rSges(6,2) * t234 + t144 * t226 + t229 * t318;
t212 = -t239 - t277;
t211 = -t142 * t148 + t119 + t145;
t77 = t214 * t144;
t155 = -t138 - t307;
t153 = t142 * t155 - t278;
t23 = (rSges(6,2) - t148) * t144 + t153;
t24 = t211 + t271;
t197 = t142 * t24 + t144 * t23;
t51 = t212 * t142;
t52 = t212 * t144;
t196 = t142 * t51 + t144 * t52;
t182 = Icges(5,2) * t143 + t258;
t176 = -t307 * qJD(3) + qJD(5) * t143 - t220;
t175 = -pkin(2) - t209;
t174 = -t138 - t208;
t161 = qJD(3) * t182;
t122 = pkin(3) * t217;
t112 = t209 * qJD(3);
t96 = t208 * qJD(3);
t76 = t214 * t142;
t54 = t276 + t145 + (pkin(2) - t268) * t144 + t237;
t53 = t142 * t175 + t266 + t292;
t44 = t211 + t75;
t43 = -t278 + (rSges(5,3) - t148) * t144 + t174 * t142;
t28 = -t105 * t234 - t142 * t96 + (-t142 * t227 - t144 * t233) * pkin(3);
t27 = t105 * t235 + t122 + (-t96 - t220) * t144;
t26 = t127 * t231 + (-t145 + (-rSges(4,3) - pkin(6)) * t142 + t175 * t144) * qJD(1);
t25 = ((-pkin(2) - t270) * t142 + t292) * qJD(1) + t154;
t18 = t105 * t231 + (t144 * t174 - t134 - t145) * qJD(1) + t218;
t17 = t132 + qJD(3) * t77 + (-t278 - t240 + (-t138 - t269) * t142) * qJD(1) + t238;
t16 = qJD(1) * t52 + t142 * t176;
t15 = t144 * t176 + t235 * t239 + t122;
t6 = (t239 * qJD(3) - t226) * t142 + (t144 * t155 - t136 - t145) * qJD(1) + t218;
t5 = t132 + (-t277 - t295) * t229 + (t153 - t240) * qJD(1) + t213;
t3 = t142 * t272 + t144 * t271 + t273;
t1 = (qJD(1) * t272 - t215 * t242 + t213) * t144 + t225 + ((-t59 + t136 - t271) * qJD(1) + (t261 * t230 + (qJD(5) - t215) * t141) * t142) * t142;
t7 = [(t25 * t54 + t26 * t53) * t283 + (t17 * t44 + t18 * t43) * t282 + (t23 * t6 + t24 * t5) * t281 + (-Icges(4,2) * t151 + t191 - t260) * t228 + (Icges(4,1) * t149 + t185 + t259) * t227 + (-Icges(6,3) * t143 - t182 + t187 + t189 + t253) * t232 + (t183 - t178 + t321) * t230; 0; 0; m(5) * (t17 * t76 + t18 * t77 + t27 * t43 + t28 * t44) + m(6) * (t15 * t23 + t16 * t24 + t5 * t51 + t52 * t6) + (m(4) * (-t112 * t53 - t127 * t26) + (t161 * t311 + t61 * t308 + t67 * t309) * t143 + (t319 * t309 + t316 * t311) * t141 + t219 * t144) * t144 + (m(4) * (-t112 * t54 - t127 * t25) + (t161 * t310 + t60 * t308 + t66 * t309) * t143 + (t320 * t309 + t316 * t310) * t141 + t219 * t142) * t142 + ((t139 / 0.2e1 + t140 / 0.2e1) * t180 + t169 / 0.2e1 - t168 / 0.2e1 - t173 / 0.2e1 + t171 / 0.2e1 + t172 / 0.2e1 - t170 / 0.2e1) * qJD(3) + ((-t54 * t279 + (-t61 / 0.2e1 + t67 / 0.2e1) * t143 + (t69 / 0.2e1 + t71 / 0.2e1) * t141) * t144 + (t53 * t279 + (-t60 / 0.2e1 + t66 / 0.2e1) * t143 + (t68 / 0.2e1 + t70 / 0.2e1) * t141) * t142) * qJD(1); m(4) * t4 + m(5) * t2 + m(6) * t1; (t1 * t3 + t15 * t52 + t16 * t51) * t281 + (t273 * t2 + t77 * t27 + t76 * t28) * t282 + t236 * t127 * t112 * t283 + (t286 * t140 + t291 * t235 + t85 * t301 + t75 * t302 + (-t144 * t312 - t306) * t234) * t144 + (t73 * t302 + t84 * t301 + t287 * t139 + t288 * t234 + ((t227 * t81 + t228 * t83 + t286 - t305) * t142 + (t227 * t80 + t228 * t82 + t287) * t144 + t304 * t232 - t303 * t230 + ((-t149 * t83 - t151 * t81) * t142 + (-t149 * t82 - t151 * t80) * t144 + t303 * t143 - t304 * t141) * qJD(3) + ((-t312 + t314) * t142 + t317 + t288 + t291) * qJD(1)) * t144 + (t285 * t142 + t306) * t235) * t142; m(5) * (t142 * t18 - t144 * t17 + (t142 * t44 + t144 * t43) * qJD(1)) + m(6) * (qJD(1) * t197 + t142 * t6 - t144 * t5); 0; m(6) * (qJD(1) * t196 + t142 * t15 - t144 * t16) + m(5) * (t142 * t27 - t144 * t28 + (t142 * t76 + t144 * t77) * qJD(1)); 0; m(6) * (t197 * t230 + (t142 * t5 + t144 * t6 + (-t142 * t23 + t144 * t24) * qJD(1)) * t141); t224; m(6) * ((qJD(3) * t196 - t1) * t143 + (qJD(3) * t3 + t142 * t16 + t144 * t15 + (-t142 * t52 + t144 * t51) * qJD(1)) * t141); 0; 0.2e1 * (-0.1e1 + t236) * t143 * t224;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
