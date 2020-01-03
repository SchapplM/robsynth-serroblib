% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:27
% DurationCPUTime: 4.89s
% Computational Cost: add. (4803->342), mult. (8367->423), div. (0->0), fcn. (8460->8), ass. (0->190)
t261 = sin(pkin(7));
t262 = cos(pkin(7));
t289 = sin(qJ(1));
t290 = cos(qJ(1));
t129 = -t261 * t289 - t262 * t290;
t130 = t290 * t261 - t289 * t262;
t160 = pkin(8) + qJ(5);
t149 = sin(t160);
t150 = cos(t160);
t259 = Icges(6,4) * t150;
t201 = -Icges(6,2) * t149 + t259;
t56 = Icges(6,6) * t129 + t130 * t201;
t260 = Icges(6,4) * t149;
t203 = Icges(6,1) * t150 - t260;
t59 = Icges(6,5) * t129 + t130 * t203;
t306 = t149 * t56 - t150 * t59;
t199 = Icges(6,5) * t150 - Icges(6,6) * t149;
t53 = Icges(6,3) * t129 + t130 * t199;
t18 = t129 * t53 - t306 * t130;
t55 = Icges(6,3) * t130 - t129 * t199;
t309 = t130 * t55;
t163 = -pkin(6) - qJ(4);
t234 = t289 * pkin(2);
t276 = rSges(6,2) * t149;
t212 = rSges(6,1) * t150 - t276;
t62 = t129 * rSges(6,3) + t130 * t212;
t308 = -t129 * t163 - t234 + t62;
t200 = Icges(6,2) * t150 + t260;
t202 = Icges(6,1) * t149 + t259;
t197 = -t149 * t200 + t150 * t202;
t198 = Icges(6,5) * t149 + Icges(6,6) * t150;
t71 = t198 * t129;
t40 = t130 * t197 + t71;
t307 = qJD(1) * t40;
t29 = -t149 * t59 - t150 * t56;
t70 = t198 * t130;
t179 = t290 * pkin(1) + t289 * qJ(2);
t133 = qJD(1) * t179;
t152 = qJD(2) * t290;
t107 = -t152 + t133;
t120 = t129 * pkin(3);
t219 = qJ(4) * t130 - t120;
t228 = qJD(1) * t290;
t243 = -pkin(2) * t228 + t152;
t231 = t133 - t243;
t240 = t129 * qJD(4);
t299 = qJD(1) * t219 + t231 - t240;
t154 = t290 * qJ(2);
t235 = t289 * pkin(1);
t134 = t235 - t154;
t151 = qJD(2) * t289;
t245 = qJ(2) * t228 + t151;
t298 = qJD(1) * t134 - t151 + t245;
t162 = cos(pkin(8));
t305 = rSges(5,2) * sin(pkin(8)) - rSges(5,1) * t162;
t304 = 2 * qJD(5);
t115 = t129 * qJ(4);
t210 = t130 * pkin(3) + t115;
t303 = -t134 + t210;
t174 = t130 * rSges(4,1) - t129 * rSges(4,2) - t234;
t302 = -t134 + t174;
t180 = t290 * rSges(3,1) + t289 * rSges(3,3);
t301 = t179 + t180;
t68 = t130 * rSges(5,3) + t305 * t129;
t28 = qJD(1) * t68 + t299;
t113 = qJD(4) * t130;
t242 = qJD(5) * t130;
t110 = t129 * qJD(1);
t211 = rSges(6,1) * t149 + rSges(6,2) * t150;
t239 = t129 * qJD(5);
t249 = t113 + t151;
t143 = pkin(4) * t162 + pkin(3);
t287 = pkin(3) - t143;
t16 = t211 * t239 + (-t130 * t287 - t115 + t303 + t308) * qJD(1) + t249;
t108 = t212 * qJD(5);
t111 = t130 * qJD(1);
t145 = qJD(1) * t152;
t158 = t290 * pkin(2);
t164 = qJD(1) ^ 2;
t195 = -t158 * t164 + t145;
t190 = qJD(4) * t110 + t195;
t253 = t111 * t163;
t98 = t111 * qJ(4);
t267 = t110 * pkin(3) - t107 + t240 - t98;
t238 = t149 * qJD(5);
t188 = -t110 * t150 + t130 * t238;
t189 = t110 * t149 + t150 * t242;
t273 = t111 * rSges(6,3);
t37 = rSges(6,1) * t188 + rSges(6,2) * t189 + t273;
t9 = (t108 * t129 - t111 * t211) * qJD(5) + (-t110 * t287 + t253 + t267 - t37 + t98) * qJD(1) + t190;
t297 = t110 * t16 + t130 * t9;
t220 = t111 * pkin(3) + t110 * qJ(4);
t296 = -qJD(1) * t210 + t298;
t281 = t200 * t130 - t59;
t283 = -t202 * t130 - t56;
t295 = t149 * t281 + t150 * t283;
t294 = t110 / 0.2e1;
t293 = t111 / 0.2e1;
t61 = Icges(6,5) * t130 - t129 * t203;
t269 = t150 * t61;
t286 = -t129 * t55 - t130 * t269;
t250 = t129 * t150;
t285 = -t61 * t250 + t309;
t266 = -t129 * t143 - t130 * t163;
t251 = t129 * t149;
t265 = -rSges(6,1) * t250 + t130 * rSges(6,3);
t64 = rSges(6,2) * t251 + t265;
t284 = -t219 + t266 + t64;
t58 = Icges(6,6) * t130 - t129 * t201;
t282 = -t202 * t129 + t58;
t280 = t200 * t129 + t61;
t279 = -t110 * t163 + t111 * t143;
t274 = t111 * rSges(5,3);
t272 = t129 * rSges(5,3);
t270 = t149 * t58;
t268 = t111 * rSges(4,1) - t110 * rSges(4,2);
t227 = qJD(1) * t289;
t264 = (-pkin(1) * t227 + t151 + t245) * qJD(1);
t254 = t111 * t150;
t263 = rSges(6,1) * t254 + t110 * rSges(6,3);
t255 = t111 * t149;
t248 = -t129 * rSges(4,1) - t130 * rSges(4,2);
t247 = -t200 + t203;
t246 = -t201 - t202;
t241 = t199 * qJD(1);
t237 = -t290 / 0.2e1;
t236 = t289 / 0.2e1;
t233 = t289 * rSges(3,1);
t232 = t240 + t243;
t230 = t158 + t179;
t229 = t129 * t238;
t224 = -t242 / 0.2e1;
t221 = t239 / 0.2e1;
t214 = t110 * rSges(4,1) + t111 * rSges(4,2);
t17 = t284 * qJD(1) + t211 * t242 + t299;
t209 = t129 * t16 + t130 * t17;
t208 = t129 * t64 - t130 * t62;
t204 = -t269 + t270;
t20 = -t130 * t53 + t250 * t59 - t251 * t56;
t196 = t110 * rSges(5,3) - t305 * t111;
t194 = pkin(3) - t305;
t191 = -t164 * t234 + t264;
t187 = t150 * t239 - t255;
t186 = t229 + t254;
t19 = t130 * t270 + t286;
t184 = (-t129 * t18 + t130 * t19) * qJD(5);
t21 = t251 * t58 + t285;
t183 = (-t129 * t20 + t130 * t21) * qJD(5);
t182 = -t235 - t234;
t181 = -t233 - t235;
t177 = t149 * t280 + t150 * t282;
t176 = t154 + t182;
t175 = qJD(1) * (t220 + t113) + qJD(4) * t111 + t191;
t38 = rSges(6,1) * t229 + rSges(6,2) * t187 + t263;
t172 = -t110 * t62 - t111 * t64 + t129 * t38 + t130 * t37;
t171 = (t149 * t246 + t150 * t247) * qJD(1);
t105 = t201 * qJD(5);
t106 = t203 * qJD(5);
t166 = -t105 * t149 + t106 * t150 + (-t149 * t202 - t150 * t200) * qJD(5);
t165 = -t130 * t305 - t234 + t272;
t156 = t290 * rSges(3,3);
t147 = rSges(3,3) * t228;
t135 = t233 - t156;
t104 = t199 * qJD(5);
t86 = t151 + (-t134 - t135) * qJD(1);
t79 = -qJD(1) * t107 - t164 * t180 + t145;
t78 = qJD(1) * (-rSges(3,1) * t227 + t147) + t264;
t77 = t211 * t129;
t76 = t211 * t130;
t65 = qJD(1) * t302 + t151;
t49 = (-t107 + t214) * qJD(1) + t195;
t48 = qJD(1) * t268 + t191;
t41 = t129 * t197 - t70;
t39 = t41 * qJD(1);
t32 = Icges(6,5) * t186 + Icges(6,6) * t187 + Icges(6,3) * t110;
t31 = Icges(6,5) * t188 + Icges(6,6) * t189 + Icges(6,3) * t111;
t30 = t149 * t61 + t150 * t58;
t27 = (t165 + t303) * qJD(1) + t249;
t24 = qJD(5) * t208 - qJD(3);
t23 = (-t110 * t305 + t267 - t274) * qJD(1) + t190;
t22 = qJD(1) * t196 + t175;
t13 = -t104 * t130 - t110 * t198 - t111 * t197 + t129 * t166;
t12 = t104 * t129 + t110 * t197 - t111 * t198 + t130 * t166;
t11 = qJD(5) * t204 - t149 * (Icges(6,1) * t186 + Icges(6,4) * t187 + Icges(6,5) * t110) - t150 * (Icges(6,4) * t186 + Icges(6,2) * t187 + Icges(6,6) * t110);
t10 = -qJD(5) * t306 - t149 * (Icges(6,1) * t188 + Icges(6,4) * t189 + Icges(6,5) * t111) - t150 * (Icges(6,4) * t188 + Icges(6,2) * t189 + Icges(6,6) * t111);
t8 = (t108 * t130 + t110 * t211) * qJD(5) + (-t220 + t38 + t279) * qJD(1) + t175;
t7 = t172 * qJD(5);
t6 = t39 + t183;
t5 = t184 + t307;
t1 = [(t105 * t150 + t106 * t149) * qJD(1) + t39 * t221 + (t11 + t13) * t242 / 0.2e1 + (t9 * t176 + t8 * (t230 + t265 + t266) + t297 * (t143 + t212) + (-rSges(6,2) * t255 + t263 + t279 + t296) * t17 + (t9 * (rSges(6,3) - t163) + t8 * t276) * t129 + (t232 + t253 - t273 + t299) * t16 + ((-t130 * t143 + t182 + t210 - t308) * t17 + (-t179 + t284) * t16) * qJD(1)) * m(6) + (t23 * (t115 + t176 + t272) + t22 * (-t120 + t68 + t230) + (t22 * qJ(4) + t194 * t23) * t130 + (t196 + (-t165 + t182) * qJD(1) + t296 + t220) * t28 + (t194 * t110 - t133 + t232 - t274 + t28 - t98) * t27) * m(5) + (t49 * t302 + t48 * (t230 + t248) + (t214 - t231) * t65 + (t268 + t65 + (-t174 + t182) * qJD(1) + t298) * (qJD(1) * t248 + t231)) * m(4) + (t79 * (t154 + t156 + t181) + t78 * t301 + (-qJD(1) * t301 + t152) * t86 + (t147 + t86 + (t135 + t181) * qJD(1) + t298) * (qJD(1) * t180 + t107)) * m(3) + (t5 - t307) * t224 - (t10 + t12 + t6) * t239 / 0.2e1 + (t197 * qJD(1) + ((t19 - t20 - t286) * t129 + t285 * t130) * t221 + (-t30 + t41) * t294 + (-t29 + t40) * t293 + ((t20 + (-t204 + t53) * t130) * t130 + (t21 + (t53 - t270) * t129 + t309 - t285) * t129) * t224) * qJD(5); 0.2e1 * (t236 * t9 + t237 * t8) * m(6) + 0.2e1 * (t22 * t237 + t23 * t236) * m(5) + 0.2e1 * (t236 * t49 + t237 * t48) * m(4) + 0.2e1 * (t236 * t79 + t237 * t78) * m(3); -m(6) * t7; (t110 * t27 + t111 * t28 - t129 * t22 + t130 * t23 - (t129 * t27 + t130 * t28) * qJD(1)) * m(5) + (-qJD(1) * t209 + t111 * t17 - t129 * t8 + t297) * m(6); qJD(1) * (-t10 * t129 + t11 * t130 - t30 * t110 - t29 * t111) / 0.2e1 + ((t71 * t242 - t241) * t130 + (t171 + (-t295 * t129 + (-t70 + t177) * t130) * qJD(5)) * t129) * t224 + ((t70 * t239 + t241) * t129 + (t171 + (t177 * t130 + (-t295 - t71) * t129) * qJD(5)) * t130) * t221 - qJD(1) * ((t247 * t149 - t246 * t150) * qJD(1) + ((t129 * t281 - t130 * t280) * t150 + (-t129 * t283 + t130 * t282) * t149) * qJD(5)) / 0.2e1 - (qJD(1) * t12 + (-(-t110 * t306 - t111 * t53 - t129 * t31) * t129 + t110 * t19 + t111 * t18 + t130 * (t110 * t204 + t111 * t55 - t129 * t32)) * t304) * t129 / 0.2e1 + (t13 * qJD(1) + (t110 * t21 + t111 * t20 - t129 * (-t110 * t53 + t111 * t306 + t130 * t31) + t130 * (t110 * t55 - t111 * t204 + t130 * t32)) * t304) * t130 / 0.2e1 + (t6 + t183) * t294 + (t5 + t184) * t293 + (t7 * t208 + t24 * t172 + t209 * t108 - (-t110 * t17 + t111 * t16 - t129 * t9 - t130 * t8) * t211 - (-t16 * t76 + t17 * t77) * qJD(1) - (t24 * (t129 * t77 + t130 * t76) + t209 * t212) * qJD(5)) * m(6);];
tauc = t1(:);
