% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% m [6x1]
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:12:13
% EndTime: 2022-01-20 09:12:27
% DurationCPUTime: 6.88s
% Computational Cost: add. (8366->405), mult. (8328->563), div. (0->0), fcn. (7746->10), ass. (0->209)
t168 = pkin(9) + qJ(5);
t163 = sin(t168);
t165 = cos(t168);
t171 = sin(pkin(8));
t173 = cos(pkin(8));
t169 = qJ(1) + pkin(7);
t166 = cos(t169);
t164 = sin(t169);
t266 = t164 * t173;
t104 = t163 * t266 + t165 * t166;
t105 = -t166 * t163 + t165 * t266;
t267 = t164 * t171;
t59 = Icges(6,5) * t105 - Icges(6,6) * t104 + Icges(6,3) * t267;
t93 = Icges(6,4) * t105;
t63 = Icges(6,2) * t104 - Icges(6,6) * t267 - t93;
t92 = Icges(6,4) * t104;
t65 = Icges(6,1) * t105 + Icges(6,5) * t267 - t92;
t23 = (t163 * t63 + t165 * t65) * t171 - t173 * t59;
t190 = (-rSges(6,1) * t163 - rSges(6,2) * t165) * t171;
t114 = qJD(5) * t190;
t245 = qJD(5) * t171;
t311 = t114 * t245;
t263 = t166 * t173;
t264 = t166 * t171;
t127 = pkin(3) * t263 + qJ(4) * t264;
t132 = t166 * pkin(2) + t164 * qJ(3);
t176 = cos(qJ(1));
t167 = t176 * pkin(1);
t297 = t167 + t132;
t211 = t127 + t297;
t172 = cos(pkin(9));
t261 = t172 * t173;
t170 = sin(pkin(9));
t262 = t170 * t173;
t268 = t164 * t170;
t302 = -(t166 * t261 + t268) * rSges(5,1) - (t164 * t172 - t166 * t262) * rSges(5,2);
t310 = rSges(5,3) * t264 + t211 - t302;
t309 = t104 * t63 + t105 * t65;
t106 = -t163 * t263 + t164 * t165;
t107 = t163 * t164 + t165 * t263;
t308 = -t106 * t63 + t107 * t65;
t177 = qJD(1) ^ 2;
t305 = t59 * t264;
t197 = rSges(4,1) * t263 - rSges(4,2) * t264 + t164 * rSges(4,3);
t304 = t197 + t297;
t265 = t166 * t170;
t303 = rSges(5,2) * (t164 * t262 + t166 * t172) - rSges(5,1) * (t164 * t261 - t265);
t153 = -qJD(5) * t173 + qJD(1);
t103 = -rSges(6,3) * t173 + (rSges(6,1) * t165 - rSges(6,2) * t163) * t171;
t228 = t103 * t245;
t206 = rSges(6,1) * t105 - rSges(6,2) * t104;
t68 = rSges(6,3) * t267 + t206;
t301 = -t153 * t68 + t164 * t228;
t16 = t305 + t308;
t298 = t16 - t305;
t218 = t166 * rSges(3,1) - rSges(3,2) * t164;
t296 = t167 + t218;
t61 = Icges(6,5) * t107 + Icges(6,6) * t106 + Icges(6,3) * t264;
t272 = Icges(6,4) * t107;
t64 = Icges(6,2) * t106 + Icges(6,6) * t264 + t272;
t94 = Icges(6,4) * t106;
t67 = Icges(6,1) * t107 + Icges(6,5) * t264 + t94;
t15 = -t104 * t64 + t105 * t67 + t61 * t267;
t294 = t164 * (-Icges(6,2) * t105 + t65 - t92) + t166 * (-Icges(6,2) * t107 + t67 + t94);
t249 = qJD(1) * t171;
t232 = t164 * t249;
t76 = qJD(1) * t104 - qJD(5) * t107;
t77 = -qJD(1) * t105 + qJD(5) * t106;
t285 = t77 * rSges(6,1) + t76 * rSges(6,2);
t42 = -rSges(6,3) * t232 + t285;
t78 = qJD(1) * t106 - qJD(5) * t105;
t79 = qJD(1) * t107 - qJD(5) * t104;
t210 = -rSges(6,1) * t79 - rSges(6,2) * t78;
t231 = t166 * t249;
t43 = rSges(6,3) * t231 - t210;
t203 = -t164 * t42 + t166 * t43;
t70 = t107 * rSges(6,1) + t106 * rSges(6,2) + rSges(6,3) * t264;
t9 = ((-t164 * t68 - t166 * t70) * qJD(1) + t203) * t245;
t293 = m(6) * t9;
t292 = t164 / 0.2e1;
t291 = -t166 / 0.2e1;
t175 = sin(qJ(1));
t290 = pkin(1) * t175;
t289 = pkin(3) * t173;
t288 = pkin(4) * t170;
t284 = t303 * qJD(1);
t247 = qJD(4) * t171;
t139 = t166 * t247;
t242 = t177 * t290;
t244 = qJD(1) * qJD(3);
t251 = qJD(1) * t164;
t154 = qJD(3) * t164;
t250 = qJD(1) * t166;
t253 = qJ(3) * t250 + t154;
t200 = t164 * t244 + qJD(1) * (-pkin(2) * t251 + t253) - t242;
t205 = qJ(4) * t171 + t289;
t186 = t200 + (-t205 * t251 + 0.2e1 * t139) * qJD(1);
t174 = -pkin(6) - qJ(4);
t238 = qJD(1) * t288;
t257 = t166 * t238 + t174 * t232;
t162 = pkin(4) * t172 + pkin(3);
t269 = t162 * t173;
t10 = -t166 * t311 + t153 * t42 + ((t228 + (t205 - t269) * qJD(1)) * t164 + t257) * qJD(1) + t186;
t280 = t10 * t166;
t241 = t177 * t167;
t212 = t166 * t244 - t241;
t225 = (pkin(3) - t162) * t173;
t246 = qJD(5) * t103;
t260 = qJ(4) + t174;
t155 = qJD(3) * t166;
t108 = qJD(1) * t132 - t155;
t230 = t164 * t247;
t273 = -t205 * t250 - t108 - t230;
t11 = t164 * t311 - t153 * t43 + ((-t238 - t247) * t164 + (qJD(1) * t225 + (qJD(1) * t260 + t246) * t171) * t166 + t273) * qJD(1) + t212;
t279 = t11 * t164;
t14 = t267 * t59 + t309;
t278 = t14 * t164;
t100 = -Icges(6,3) * t173 + (Icges(6,5) * t165 - Icges(6,6) * t163) * t171;
t270 = Icges(6,4) * t165;
t101 = -Icges(6,6) * t173 + (-Icges(6,2) * t163 + t270) * t171;
t271 = Icges(6,4) * t163;
t102 = -Icges(6,5) * t173 + (Icges(6,1) * t165 - t271) * t171;
t31 = t100 * t267 - t101 * t104 + t102 * t105;
t277 = t153 * t31;
t125 = t205 * t164;
t157 = t166 * qJ(3);
t130 = pkin(2) * t164 - t157;
t222 = -t130 - t290;
t255 = t139 + t154;
t243 = pkin(4) * t265;
t75 = t243 + (t171 * t260 + t225) * t164;
t21 = (-t125 + t222 + t75) * qJD(1) + t255 + t301;
t276 = t21 * t166;
t275 = -rSges(5,3) - qJ(4);
t274 = -rSges(6,3) + t174;
t124 = (-Icges(6,1) * t163 - t270) * t171;
t259 = -t101 + t124;
t123 = (-Icges(6,2) * t165 - t271) * t171;
t258 = t102 + t123;
t256 = rSges(4,2) * t232 + rSges(4,3) * t250;
t254 = rSges(4,2) * t267 + t166 * rSges(4,3);
t252 = -qJD(1) * t130 + t154;
t248 = qJD(4) * t164;
t17 = t106 * t64 + t107 * t67 + t61 * t264;
t240 = rSges(4,1) * t266;
t235 = t139 + t253;
t233 = -pkin(2) - t289;
t226 = -rSges(4,1) * t173 - pkin(2);
t224 = -t245 / 0.2e1;
t223 = t245 / 0.2e1;
t220 = t157 - t290;
t219 = -pkin(2) - t269;
t217 = -qJD(1) * t125 + t139 + t252;
t216 = t164 * t224;
t215 = t164 * t223;
t214 = t166 * t224;
t213 = t166 * t223;
t131 = rSges(3,1) * t164 + rSges(3,2) * t166;
t207 = t302 * qJD(1);
t196 = pkin(4) * t268 + t162 * t263 - t174 * t264;
t22 = t153 * t70 - t155 + (-t166 * t246 + t248) * t171 + (t196 + t211 - t127) * qJD(1);
t204 = t21 * t164 - t166 * t22;
t202 = -t164 * t70 + t166 * t68;
t201 = t164 * (-Icges(6,5) * t104 - Icges(6,6) * t105) + t166 * (Icges(6,5) * t106 - Icges(6,6) * t107);
t199 = qJD(1) * t216;
t198 = qJD(1) * t213;
t191 = t171 * t201;
t189 = (t15 * t166 + t278) * t171;
t188 = (t16 * t164 + t166 * t17) * t171;
t122 = (-Icges(6,5) * t163 - Icges(6,6) * t165) * t171;
t184 = -t290 + t303;
t183 = (Icges(6,1) * t106 - t272 - t64) * t166 + (-Icges(6,1) * t104 + t63 - t93) * t164;
t181 = -rSges(5,3) * t267 + t184;
t36 = Icges(6,5) * t77 + Icges(6,6) * t76 - Icges(6,3) * t232;
t37 = Icges(6,5) * t79 + Icges(6,6) * t78 + Icges(6,3) * t231;
t38 = Icges(6,4) * t77 + Icges(6,2) * t76 - Icges(6,6) * t232;
t39 = Icges(6,4) * t79 + Icges(6,2) * t78 + Icges(6,6) * t231;
t40 = Icges(6,1) * t77 + Icges(6,4) * t76 - Icges(6,5) * t232;
t41 = Icges(6,1) * t79 + Icges(6,4) * t78 + Icges(6,5) * t231;
t180 = ((t106 * t39 + t107 * t41 - t63 * t76 + t65 * t77 + (t166 * t37 - t251 * t59) * t171) * t164 + t166 * (t106 * t38 + t107 * t40 + t64 * t76 + t67 * t77 + (t166 * t36 - t251 * t61) * t171) + (t16 * t166 - t17 * t164) * qJD(1)) * t171;
t179 = (t164 * (-t104 * t39 + t105 * t41 - t63 * t78 + t65 * t79 + (t164 * t37 + t250 * t59) * t171) + t166 * (-t104 * t38 + t105 * t40 + t64 * t78 + t67 * t79 + (t164 * t36 + t250 * t61) * t171) + (t14 * t166 - t15 * t164) * qJD(1)) * t171;
t24 = -t173 * t61 + (-t163 * t64 + t165 * t67) * t171;
t7 = -t173 * t37 + (-t163 * t39 + t165 * t41 + (-t163 * t65 + t165 * t63) * qJD(5)) * t171;
t8 = -t173 * t36 + (-t163 * t38 + t165 * t40 + (-t163 * t67 - t165 * t64) * qJD(5)) * t171;
t178 = (t164 * t7 + t166 * t8 + (-t164 * t24 + t166 * t23) * qJD(1)) * t171;
t113 = qJD(5) * t124;
t112 = qJD(5) * t123;
t111 = qJD(5) * t122;
t91 = t240 - t254;
t87 = rSges(6,1) * t106 - rSges(6,2) * t107;
t86 = -rSges(6,1) * t104 - rSges(6,2) * t105;
t72 = qJD(1) * t304 - t155;
t71 = t154 + (t222 - t91) * qJD(1);
t53 = (-qJD(1) * t197 - t108) * qJD(1) + t212;
t52 = qJD(1) * (-qJD(1) * t240 + t256) + t200;
t45 = qJD(1) * t310 - t155 + t230;
t44 = (-t125 - t130 + t181) * qJD(1) + t255;
t32 = t100 * t264 + t101 * t106 + t102 * t107;
t30 = t32 * t153;
t29 = ((-rSges(5,3) * t250 - t248) * t171 + t207 + t273) * qJD(1) + t212;
t28 = qJD(1) * (-rSges(5,3) * t232 + t284) + t186;
t27 = -qJD(4) * t173 + t202 * t245 + qJD(2);
t26 = -t111 * t173 + (-t112 * t163 + t113 * t165 + (-t101 * t165 - t102 * t163) * qJD(5)) * t171;
t25 = t26 * t153;
t13 = t101 * t78 + t102 * t79 - t104 * t112 + t105 * t113 + (t100 * t250 + t111 * t164) * t171;
t12 = t101 * t76 + t102 * t77 + t106 * t112 + t107 * t113 + (-t100 * t251 + t111 * t166) * t171;
t6 = qJD(5) * t188 + t30;
t5 = qJD(5) * t189 + t277;
t1 = [(t30 + ((t14 + t17 - t309) * t166 + t298 * t164) * t245) * t216 + t25 + m(3) * ((-t131 * t177 - t242) * t296 + (-t241 + (-0.2e1 * t218 - t167 + t296) * t177) * (-t131 - t290)) + (t8 + t12) * t213 + (t24 + t32) * t199 + (t23 + t31) * t198 + (t11 * (-t206 + t220 + t243) + t21 * (t155 + t210) + t10 * (t196 + t297 + t70) + t22 * (t235 + t257 + t285) + (t11 * t219 + (-t21 * qJD(4) + t11 * t274) * t171) * t164 + ((-t22 * t175 - t21 * t176) * pkin(1) + (t171 * t274 + t219) * t276 + (t21 * (-qJ(3) - t288) + t22 * (-rSges(6,3) * t171 + t219)) * t164) * qJD(1) - (-t21 + (t75 - t290) * qJD(1) + t217 + t301) * t22) * m(6) + (t29 * (t157 + t184) + t44 * (t155 + t207) + t28 * t310 + t45 * (t235 + t284) + (t29 * t233 + (-t44 * qJD(4) + t275 * t29) * t171) * t164 + ((-t45 * t175 - t44 * t176) * pkin(1) + t44 * (t171 * t275 + t233) * t166 + (-t44 * qJ(3) + t45 * (-rSges(5,3) * t171 - pkin(2) - t205)) * t164) * qJD(1) - (qJD(1) * t181 + t217 - t44) * t45) * m(5) + (t53 * (t164 * t226 + t220 + t254) + t71 * t155 + t52 * t304 + t72 * (t253 + t256) + ((-t72 * t175 - t71 * t176) * pkin(1) + t71 * (rSges(4,2) * t171 + t226) * t166 + (t71 * (-rSges(4,3) - qJ(3)) + t72 * t226) * t164) * qJD(1) - (-t71 + (-t91 - t290) * qJD(1) + t252) * t72) * m(4) + (t7 + t6 + t13) * t215 + (t5 - t277 + ((-t15 + t298 - t308) * t166 - t278) * t245) * t214; t293; 0.2e1 * (-t280 / 0.2e1 + t279 / 0.2e1) * m(6) + 0.2e1 * (t28 * t291 + t29 * t292) * m(5) + 0.2e1 * (t291 * t52 + t292 * t53) * m(4); -t173 * t293 + 0.2e1 * (m(5) * (t164 * t28 + t166 * t29) / 0.2e1 + m(6) * (t10 * t164 + t11 * t166) / 0.2e1) * t171; -t6 * t232 / 0.2e1 + (-t173 * t32 + t188) * t199 + (-t12 * t173 + t180) * t213 + (qJD(5) * t179 + t13 * t153) * t267 / 0.2e1 + (-t173 * t31 + t189) * t198 + (-t13 * t173 + t179) * t215 - t173 * (qJD(5) * t178 + t25) / 0.2e1 + t153 * (-t173 * t26 + t178) / 0.2e1 + ((t106 * t258 + t107 * t259 + t122 * t264) * t153 + (t106 * t294 + t183 * t107 + t166 * t191) * t245) * t214 + ((-t104 * t258 + t105 * t259 + t122 * t267) * t153 + (-t104 * t294 + t105 * t183 + t164 * t191) * t245) * t216 - t153 * (-t173 * t122 * t153 + ((-t163 * t258 + t165 * t259) * t153 + ((-t163 * t294 + t165 * t183) * t171 - t201 * t173) * qJD(5)) * t171) / 0.2e1 + (qJD(1) * t5 + qJD(5) * t180 + t12 * t153) * t264 / 0.2e1 + ((-t10 * t70 + t11 * t68 + t21 * t43 - t22 * t42) * t173 + (t9 * t202 + t27 * (-t250 * t70 - t251 * t68 + t203) + t204 * t114 + (-t280 + t279 + (t22 * t164 + t276) * qJD(1)) * t103) * t171 - (-t21 * t86 + t22 * t87) * t153 - (t27 * (-t164 * t87 + t166 * t86) + t204 * t190) * t245) * m(6);];
tauc = t1(:);
