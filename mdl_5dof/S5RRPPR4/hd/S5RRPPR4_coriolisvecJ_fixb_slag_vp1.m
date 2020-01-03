% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:46
% DurationCPUTime: 4.39s
% Computational Cost: add. (8996->336), mult. (9162->419), div. (0->0), fcn. (9230->8), ass. (0->201)
t283 = qJ(1) + qJ(2);
t258 = sin(t283);
t259 = cos(t283);
t300 = sin(pkin(8));
t301 = cos(pkin(8));
t123 = -t258 * t300 - t259 * t301;
t124 = -t258 * t301 + t259 * t300;
t168 = sin(qJ(5));
t170 = cos(qJ(5));
t298 = Icges(6,4) * t170;
t222 = -Icges(6,2) * t168 + t298;
t69 = Icges(6,6) * t123 + t124 * t222;
t299 = Icges(6,4) * t168;
t224 = Icges(6,1) * t170 - t299;
t72 = Icges(6,5) * t123 + t124 * t224;
t337 = t168 * t69 - t170 * t72;
t220 = Icges(6,5) * t170 - Icges(6,6) * t168;
t66 = Icges(6,3) * t123 + t124 * t220;
t24 = t123 * t66 - t124 * t337;
t68 = Icges(6,3) * t124 - t123 * t220;
t338 = t124 * t68;
t35 = -t168 * t72 - t170 * t69;
t167 = qJD(1) + qJD(2);
t106 = t123 * t167;
t107 = t124 * t167;
t271 = qJD(5) * t168;
t206 = t107 * t170 + t123 * t271;
t270 = qJD(5) * t170;
t207 = -t107 * t168 + t123 * t270;
t44 = t206 * rSges(6,1) + rSges(6,2) * t207 + t106 * rSges(6,3);
t336 = t107 * pkin(4) + t106 * pkin(7) + t44;
t163 = t259 * pkin(3);
t275 = t259 * pkin(2) + t258 * qJ(3);
t324 = t163 + t275;
t225 = -t123 * rSges(5,1) - t124 * rSges(5,2) + t324;
t219 = Icges(6,5) * t168 + Icges(6,6) * t170;
t83 = t219 * t123;
t82 = t219 * t124;
t209 = t106 * t168 + t124 * t270;
t208 = -t106 * t170 + t124 * t271;
t169 = sin(qJ(1));
t304 = pkin(1) * qJD(1);
t267 = t169 * t304;
t127 = rSges(3,1) * t258 + rSges(3,2) * t259;
t285 = t127 * t167;
t108 = -t267 - t285;
t288 = t123 * t170;
t289 = t123 * t168;
t77 = -rSges(6,1) * t288 + rSges(6,2) * t289 + t124 * rSges(6,3);
t185 = -t123 * pkin(4) + pkin(7) * t124 + t324 + t77;
t335 = 0.2e1 * qJD(5);
t221 = Icges(6,2) * t170 + t299;
t223 = Icges(6,1) * t168 + t298;
t218 = -t168 * t221 + t170 * t223;
t52 = t124 * t218 + t83;
t332 = t167 * t52;
t243 = t258 * pkin(3);
t244 = t258 * pkin(2);
t196 = -t244 - t243;
t274 = qJD(5) * t123;
t158 = t259 * qJ(3);
t125 = t244 - t158;
t213 = -t125 - t243;
t155 = qJD(3) * t258;
t232 = rSges(6,1) * t168 + rSges(6,2) * t170;
t272 = qJD(5) * t232;
t264 = t123 * t272;
t237 = t155 + t264;
t119 = t123 * pkin(7);
t321 = t124 * pkin(4) + t119;
t287 = t124 * t168;
t303 = t123 * rSges(6,3);
t249 = -rSges(6,2) * t287 + t303;
t286 = t124 * t170;
t76 = -rSges(6,1) * t286 - t249;
t30 = -t267 + (t213 - t76 + t321) * t167 + t237;
t156 = qJD(3) * t259;
t171 = cos(qJ(1));
t266 = t171 * t304;
t241 = -t156 + t266;
t319 = t124 * t272 + t167 * t185;
t31 = t241 + t319;
t230 = t123 * t77 + t124 * t76;
t32 = qJD(5) * t230 - qJD(4);
t305 = rSges(6,1) * t170;
t233 = -rSges(6,2) * t168 + t305;
t75 = t124 * t233 + t303;
t331 = (t196 * t31 - t30 * t324) * t167 - t32 * (t75 + t76) * t274;
t43 = t208 * rSges(6,1) + t209 * rSges(6,2) + t107 * rSges(6,3);
t330 = t106 * pkin(4) - t107 * pkin(7) - t43;
t234 = t124 * rSges(5,1) - t123 * rSges(5,2);
t282 = t107 * rSges(5,1) - t106 * rSges(5,2);
t329 = -t167 * t234 + t282;
t328 = t167 * t225;
t327 = t167 * t243;
t276 = t259 * rSges(4,1) + t258 * rSges(4,3);
t323 = t275 + t276;
t326 = t167 * t323;
t160 = t259 * rSges(4,3);
t240 = t258 * rSges(4,1);
t126 = t240 - t160;
t239 = t167 * t259;
t147 = rSges(4,3) * t239;
t325 = t167 * t126 + t147;
t121 = t167 * t125;
t139 = t167 * t158;
t280 = t139 + t155;
t322 = t280 - t155 + t121;
t320 = t280 + t336 + (-t321 - t75) * t167;
t307 = t221 * t124 - t72;
t309 = -t223 * t124 - t69;
t318 = t168 * t307 + t170 * t309;
t166 = t167 ^ 2;
t317 = t106 / 0.2e1;
t316 = t107 / 0.2e1;
t312 = t169 * pkin(1);
t165 = t171 * pkin(1);
t74 = Icges(6,5) * t124 - t123 * t224;
t311 = -t123 * t68 - t74 * t286;
t310 = -t74 * t288 + t338;
t71 = Icges(6,6) * t124 - t123 * t222;
t308 = -t223 * t123 + t71;
t306 = t221 * t123 + t74;
t302 = t168 * t71;
t284 = t220 * t167;
t279 = -t221 + t224;
t278 = -t222 - t223;
t273 = qJD(5) * t124;
t172 = qJD(1) ^ 2;
t269 = t172 * t312;
t268 = t172 * t165;
t254 = t274 / 0.2e1;
t253 = -t273 / 0.2e1;
t247 = -t259 / 0.2e1;
t246 = t258 / 0.2e1;
t245 = t167 * t156 - t268;
t242 = t155 - t267;
t238 = t167 * t258;
t235 = t106 * rSges(5,1) + t107 * rSges(5,2);
t231 = -t123 * t30 - t124 * t31;
t226 = -t170 * t74 + t302;
t217 = -t269 + (-pkin(2) * t238 + t155 + t280) * t167;
t26 = -t124 * t66 + t288 * t72 - t289 * t69;
t216 = -t121 + t242;
t215 = t139 + t242;
t214 = t156 + t235;
t130 = t259 * rSges(3,1) - t258 * rSges(3,2);
t25 = t287 * t71 + t311;
t204 = qJD(5) * (-t123 * t24 + t124 * t25);
t27 = t289 * t71 + t310;
t202 = (-t123 * t26 + t124 * t27) * qJD(5);
t115 = rSges(3,1) * t239 - rSges(3,2) * t238;
t200 = t306 * t168 + t308 * t170;
t198 = t106 * t76 - t107 * t77 + t123 * t44 + t124 * t43;
t197 = (t168 * t278 + t170 * t279) * t167;
t195 = -t240 - t244;
t192 = -t163 * t166 + t245;
t191 = t158 + t196;
t188 = -t166 * t243 + t217;
t187 = -pkin(3) * t238 + t216;
t53 = t123 * t218 - t82;
t49 = t53 * t167;
t10 = t49 + t202;
t132 = t222 * qJD(5);
t133 = t224 * qJD(5);
t15 = -qJD(5) * t337 - t168 * (Icges(6,1) * t208 + Icges(6,4) * t209 + Icges(6,5) * t107) - t170 * (Icges(6,4) * t208 + Icges(6,2) * t209 + Icges(6,6) * t107);
t16 = qJD(5) * t226 - t168 * (Icges(6,1) * t206 + Icges(6,4) * t207 + Icges(6,5) * t106) - t170 * (Icges(6,4) * t206 + Icges(6,2) * t207 + Icges(6,6) * t106);
t131 = t220 * qJD(5);
t181 = -t132 * t168 + t133 * t170 + (-t168 * t223 - t170 * t221) * qJD(5);
t22 = t106 * t218 - t107 * t219 + t123 * t131 + t124 * t181;
t23 = -t106 * t219 - t107 * t218 + t123 * t181 - t124 * t131;
t36 = t168 * t74 + t170 * t71;
t9 = t204 + t332;
t186 = (t132 * t170 + t133 * t168) * t167 + t49 * t254 + (t23 + t16) * t273 / 0.2e1 - (t22 + t15 + t10) * t274 / 0.2e1 + (t9 - t332) * t253 + (t218 * t167 + ((t25 - t26 - t311) * t123 + t310 * t124) * t254 + (-t36 + t53) * t317 + (t52 - t35) * t316 + ((t26 + (-t226 + t66) * t124) * t124 + (t27 + (t66 - t302) * t123 + t338 - t310) * t123) * t253) * qJD(5);
t182 = t158 + t160 + t195;
t180 = t156 + t330;
t178 = t191 + t234;
t176 = t119 + (pkin(4) + t305) * t124 + t191 + t249;
t57 = (t213 + t234) * t167 + t242;
t58 = t241 + t328;
t174 = (t196 * t58 - t324 * t57) * t167;
t80 = (-t125 - t126) * t167 + t242;
t81 = t241 + t326;
t173 = (t81 * t195 - t323 * t80) * t167;
t134 = t233 * qJD(5);
t109 = t130 * t167 + t266;
t100 = t167 * t275 - t156;
t96 = -t115 * t167 - t268;
t95 = -t167 * t285 - t269;
t89 = t232 * t123;
t88 = t232 * t124;
t63 = -t167 * t100 - t166 * t276 + t245;
t62 = t167 * (-rSges(4,1) * t238 + t147) + t217;
t48 = (-t100 + t235) * t167 + t192;
t47 = t167 * t282 + t188;
t38 = Icges(6,5) * t206 + Icges(6,6) * t207 + Icges(6,3) * t106;
t37 = Icges(6,5) * t208 + Icges(6,6) * t209 + Icges(6,3) * t107;
t19 = (-t107 * t232 + t123 * t134) * qJD(5) + (-t100 + t330) * t167 + t192;
t18 = t336 * t167 + (t106 * t232 + t124 * t134) * qJD(5) + t188;
t12 = t198 * qJD(5);
t1 = [t186 + m(3) * (t96 * (-t127 - t312) + t95 * (t165 + t130) + (-t115 - t266 + t109) * t108) + (t19 * (t176 - t312) + t30 * (t180 - t266) + t18 * (t165 + t185) + (-t187 + t30 - t264 - t267 + t320) * t31 + t331) * m(6) + (t48 * (t178 - t312) + t57 * (t214 - t266) + t47 * (t165 + t225) + t174 + (t215 + t57 - t187 + t329) * t58) * m(5) + (t63 * (t182 - t312) - t80 * t241 + t62 * (t165 + t323) + t173 + (t215 - t216 + t80 + t325) * t81) * m(4); t186 + (t19 * t176 + t18 * t185 + (t121 - t237 + t320 + t327) * t31 + (-t156 + t180 + t319) * t30 + t331) * m(6) + (t48 * t178 + t47 * t225 + t174 + (t322 + t327 + t329) * t58 + (-t156 + t214 + t328) * t57) * m(5) + (t80 * t326 + t63 * t182 + t62 * t323 + t173 + (t322 + t325) * t81) * m(4) + (-t108 * t115 - t109 * t285 - t127 * t96 + t130 * t95 - (-t108 * t130 - t109 * t127) * t167) * m(3); 0.2e1 * (t18 * t247 + t19 * t246) * m(6) + 0.2e1 * (t246 * t48 + t247 * t47) * m(5) + 0.2e1 * (t246 * t63 + t247 * t62) * m(4); -m(6) * t12; t167 * (-t106 * t36 - t107 * t35 - t123 * t15 + t124 * t16) / 0.2e1 + ((t83 * t273 - t284) * t124 + (t197 + (-t318 * t123 + (-t82 + t200) * t124) * qJD(5)) * t123) * t253 + ((t82 * t274 + t284) * t123 + (t197 + (t200 * t124 + (-t318 - t83) * t123) * qJD(5)) * t124) * t254 - t167 * ((t279 * t168 - t278 * t170) * t167 + ((t123 * t307 - t124 * t306) * t170 + (-t123 * t309 + t124 * t308) * t168) * qJD(5)) / 0.2e1 - (t167 * t22 + (t106 * t25 + t107 * t24 - t123 * (-t106 * t337 - t107 * t66 - t123 * t37) + t124 * (t106 * t226 + t107 * t68 - t123 * t38)) * t335) * t123 / 0.2e1 + (t167 * t23 + (t106 * t27 + t107 * t26 - t123 * (-t106 * t66 + t107 * t337 + t124 * t37) + t124 * (t106 * t68 - t107 * t226 + t124 * t38)) * t335) * t124 / 0.2e1 + (t10 + t202) * t317 + (t9 + t204) * t316 + (t12 * t230 + t32 * t198 - t231 * t134 - (-t106 * t31 + t107 * t30 - t123 * t19 - t124 * t18) * t232 - (-t30 * t88 + t31 * t89) * t167 - (t32 * (t123 * t89 + t124 * t88) - t231 * t233) * qJD(5)) * m(6);];
tauc = t1(:);
