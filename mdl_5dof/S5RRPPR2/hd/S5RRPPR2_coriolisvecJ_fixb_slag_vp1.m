% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:20:03
% EndTime: 2019-12-05 18:20:14
% DurationCPUTime: 5.14s
% Computational Cost: add. (11071->405), mult. (9267->529), div. (0->0), fcn. (8721->10), ass. (0->221)
t171 = qJ(1) + qJ(2);
t165 = pkin(8) + t171;
t162 = sin(t165);
t163 = cos(t165);
t176 = cos(qJ(5));
t173 = cos(pkin(9));
t174 = sin(qJ(5));
t263 = t173 * t174;
t112 = t162 * t263 + t163 * t176;
t262 = t173 * t176;
t113 = t162 * t262 - t163 * t174;
t115 = t162 * t174 + t163 * t262;
t107 = Icges(6,4) * t115;
t114 = -t162 * t176 + t163 * t263;
t172 = sin(pkin(9));
t267 = t163 * t172;
t69 = -Icges(6,2) * t114 + Icges(6,6) * t267 + t107;
t106 = Icges(6,4) * t114;
t73 = -Icges(6,1) * t115 - Icges(6,5) * t267 + t106;
t288 = t112 * t69 + t113 * t73;
t216 = -t114 * t69 - t115 * t73;
t66 = Icges(6,5) * t115 - Icges(6,6) * t114 + Icges(6,3) * t267;
t31 = -t173 * t66 + (-t174 * t69 - t176 * t73) * t172;
t170 = qJD(1) + qJD(2);
t166 = sin(t171);
t167 = cos(t171);
t220 = rSges(3,1) * t166 + rSges(3,2) * t167;
t124 = t220 * t170;
t175 = sin(qJ(1));
t277 = pkin(1) * qJD(1);
t248 = t175 * t277;
t110 = t124 + t248;
t126 = -rSges(6,3) * t173 + (rSges(6,1) * t176 - rSges(6,2) * t174) * t172;
t157 = -qJD(5) * t173 + t170;
t254 = qJD(5) * t172;
t240 = t163 * t254;
t217 = rSges(6,1) * t115 - rSges(6,2) * t114;
t249 = rSges(6,3) * t267;
t74 = t217 + t249;
t301 = t126 * t240 - t157 * t74;
t150 = rSges(5,2) * t267;
t279 = rSges(5,3) * t162;
t280 = rSges(5,1) * t173;
t207 = -t163 * t280 - t279;
t102 = -t150 - t207;
t264 = t170 * t172;
t244 = t163 * t264;
t140 = rSges(5,2) * t244;
t300 = -t170 * t102 - t140;
t289 = pkin(4) * t173;
t224 = pkin(7) * t172 + t289;
t192 = -rSges(6,3) * t172 - pkin(3) - t224;
t203 = (-rSges(6,1) * t174 - rSges(6,2) * t176) * t172;
t130 = qJD(5) * t203;
t242 = t126 * t254;
t293 = pkin(1) * qJD(1) ^ 2;
t164 = t175 * t293;
t290 = pkin(2) * t170 ^ 2;
t256 = t166 * t290 + t164;
t270 = t162 * t170;
t245 = t173 * t270;
t246 = t162 * t264;
t258 = pkin(4) * t245 + pkin(7) * t246;
t78 = -qJD(5) * t115 + t112 * t170;
t79 = -qJD(5) * t114 - t113 * t170;
t47 = rSges(6,1) * t79 + rSges(6,2) * t78 - rSges(6,3) * t246;
t149 = pkin(3) * t270;
t158 = qJD(4) * t162;
t255 = t158 - t149;
t268 = t163 * t170;
t97 = qJ(4) * t268 + t255;
t21 = t130 * t240 - t157 * t47 + (-t97 + (-qJD(4) - t242) * t162 + t258) * t170 + t256;
t271 = qJ(4) * t162;
t132 = pkin(3) * t163 + t271;
t122 = t170 * t132;
t159 = qJD(4) * t163;
t177 = cos(qJ(1));
t251 = t177 * t293;
t200 = -t167 * t290 - t251;
t190 = t200 + (0.2e1 * t159 - t122) * t170;
t206 = t224 * t170;
t241 = t162 * t254;
t80 = qJD(5) * t113 + t114 * t170;
t81 = qJD(5) * t112 - t115 * t170;
t282 = t81 * rSges(6,1) + t80 * rSges(6,2);
t48 = -rSges(6,3) * t244 + t282;
t22 = t130 * t241 + t157 * t48 + (-t206 + t242) * t268 + t190;
t120 = t170 * (-pkin(3) * t162 + qJ(4) * t163);
t266 = t166 * t170;
t253 = pkin(2) * t266;
t199 = t248 + t253;
t189 = -t120 - t158 + t199;
t261 = -t113 * rSges(6,1) + t112 * rSges(6,2);
t269 = t162 * t172;
t209 = -rSges(6,3) * t269 + t261;
t297 = -t126 * t241 - t157 * t209 + t162 * t206;
t33 = t189 + t297;
t116 = t224 * t163;
t247 = t177 * t277;
t226 = t159 - t247;
t291 = pkin(2) * t167;
t235 = -t132 - t291;
t34 = (-t116 + t235) * t170 + t226 + t301;
t299 = (t22 * (-t289 - pkin(3) + (-rSges(6,3) - pkin(7)) * t172) + (t170 * t33 - t21) * qJ(4)) * t162 + (t21 * t192 + t22 * qJ(4) + (-t34 * qJ(4) - t192 * t33) * t170) * t163;
t160 = t162 * rSges(4,2);
t281 = rSges(4,1) * t163;
t133 = -t160 + t281;
t146 = rSges(4,2) * t270;
t223 = -t281 - t291;
t298 = -t146 + (-t133 - t223) * t170;
t65 = -Icges(6,5) * t113 + Icges(6,6) * t112 - Icges(6,3) * t269;
t274 = Icges(6,4) * t113;
t68 = Icges(6,2) * t112 - Icges(6,6) * t269 - t274;
t105 = Icges(6,4) * t112;
t71 = -Icges(6,1) * t113 - Icges(6,5) * t269 + t105;
t26 = -t114 * t68 + t115 * t71 + t65 * t267;
t296 = -t170 * t116 - t282 + t301;
t194 = t162 * (Icges(6,2) * t113 + t105 + t71) - t163 * (-Icges(6,2) * t115 - t106 - t73);
t195 = t162 * (-Icges(6,1) * t112 - t274 + t68) - t163 * (Icges(6,1) * t114 + t107 + t69);
t295 = pkin(1) * t175;
t294 = pkin(1) * t177;
t292 = pkin(2) * t166;
t287 = -t112 * t68 + t113 * t71;
t278 = rSges(5,3) * t163;
t117 = -Icges(6,3) * t173 + (Icges(6,5) * t176 - Icges(6,6) * t174) * t172;
t272 = Icges(6,4) * t176;
t118 = -Icges(6,6) * t173 + (-Icges(6,2) * t174 + t272) * t172;
t273 = Icges(6,4) * t174;
t119 = -Icges(6,5) * t173 + (Icges(6,1) * t176 - t273) * t172;
t50 = -t114 * t118 + t115 * t119 + t117 * t267;
t276 = t157 * t50;
t275 = rSges(5,3) + qJ(4);
t265 = t167 * t170;
t136 = (-Icges(6,1) * t174 - t272) * t172;
t260 = t118 - t136;
t135 = (-Icges(6,2) * t176 - t273) * t172;
t259 = t119 + t135;
t257 = -rSges(4,1) * t270 - rSges(4,2) * t268;
t252 = pkin(2) * t265;
t250 = rSges(5,2) * t269;
t139 = rSges(5,1) * t245;
t243 = t139 - t255;
t238 = -pkin(3) - t280;
t237 = -t254 / 0.2e1;
t236 = t254 / 0.2e1;
t234 = t150 - t291;
t143 = rSges(3,1) * t167 - t166 * rSges(3,2);
t233 = t162 * t237;
t232 = t162 * t236;
t231 = t163 * t237;
t230 = t163 * t236;
t229 = t170 * t237;
t228 = -t158 + t253;
t227 = t159 - t252;
t125 = -rSges(3,1) * t265 + rSges(3,2) * t266;
t225 = -t292 - t295;
t219 = -rSges(4,1) * t162 - rSges(4,2) * t163;
t218 = rSges(5,2) * t172 - t280;
t215 = t162 * (Icges(6,5) * t112 + Icges(6,6) * t113) - t163 * (-Icges(6,5) * t114 - Icges(6,6) * t115);
t214 = t269 * t65 + t287;
t213 = t162 * t229;
t212 = t163 * t229;
t211 = -t120 + t228;
t210 = -t122 + t227;
t208 = t160 + t223;
t111 = -t143 * t170 - t247;
t25 = -t269 * t66 + t288;
t202 = (t162 * t214 + t163 * t25) * t172;
t27 = t267 * t66 + t216;
t201 = (-t162 * t26 + t163 * t27) * t172;
t134 = (-Icges(6,5) * t174 - Icges(6,6) * t176) * t172;
t198 = t247 + t252;
t197 = t219 - t292;
t196 = -t217 - t291;
t193 = t159 - t198;
t188 = t122 - t193;
t10 = qJD(5) * t201 + t276;
t42 = Icges(6,5) * t81 + Icges(6,6) * t80 - Icges(6,3) * t244;
t44 = Icges(6,4) * t81 + Icges(6,2) * t80 - Icges(6,6) * t244;
t46 = Icges(6,1) * t81 + Icges(6,4) * t80 - Icges(6,5) * t244;
t13 = -t173 * t42 + (-t174 * t44 + t176 * t46 + (-t174 * t71 - t176 * t68) * qJD(5)) * t172;
t41 = Icges(6,5) * t79 + Icges(6,6) * t78 - Icges(6,3) * t246;
t43 = Icges(6,4) * t79 + Icges(6,2) * t78 - Icges(6,6) * t246;
t45 = Icges(6,1) * t79 + Icges(6,4) * t78 - Icges(6,5) * t246;
t14 = -t173 * t41 + (-t174 * t43 + t176 * t45 + (t174 * t73 - t176 * t69) * qJD(5)) * t172;
t127 = qJD(5) * t134;
t128 = qJD(5) * t135;
t129 = qJD(5) * t136;
t19 = -t114 * t128 + t115 * t129 + t118 * t78 + t119 * t79 + (-t117 * t270 + t127 * t163) * t172;
t20 = t112 * t128 - t113 * t129 + t118 * t80 + t119 * t81 + (-t117 * t268 - t127 * t162) * t172;
t30 = -t173 * t65 + (-t174 * t68 + t176 * t71) * t172;
t39 = -t127 * t173 + (-t128 * t174 + t129 * t176 + (-t118 * t176 - t119 * t174) * qJD(5)) * t172;
t36 = t39 * t157;
t49 = t112 * t118 - t113 * t119 - t117 * t269;
t40 = t49 * t157;
t9 = qJD(5) * t202 + t40;
t187 = t36 + (t40 + (t288 * t163 + (t214 + t216 - t27) * t162) * t254) * t231 + (t10 - t276 + (-(-t26 - t288) * t162 + t214 * t163 + (-t216 - t287) * t163 - t25 * t162 + (-t66 * t162 ^ 2 + (-t162 * t65 - t163 * t66) * t163) * t172) * t254) * t232 + (t31 + t50) * t213 + (t30 + t49) * t212 + (t13 + t20) * t233 + (t9 + t14 + t19) * t230;
t186 = t172 * (-t163 * t261 + (-t74 + t249) * t162);
t185 = ((t170 * t214 + t112 * t43 - t113 * t45 + t69 * t80 - t73 * t81 + (-t162 * t41 - t268 * t66) * t172) * t163 + (-t170 * t25 - t112 * t44 + t113 * t46 - t68 * t80 - t71 * t81 - (-t162 * t42 - t268 * t65) * t172) * t162) * t172;
t184 = ((-t170 * t26 - t114 * t43 + t115 * t45 + t69 * t78 - t73 * t79 + (t163 * t41 - t270 * t66) * t172) * t163 + (-t170 * t27 + t114 * t44 - t115 * t46 - t68 * t78 - t71 * t79 - (t163 * t42 - t270 * t65) * t172) * t162) * t172;
t183 = ((-t170 * t30 + t14) * t163 + (-t170 * t31 - t13) * t162) * t172;
t182 = t149 + t228 - t47 + t258;
t181 = t172 * ((-t170 * t74 - t48) * t163 + (t170 * t209 - t47) * t162);
t51 = (-t158 + t139 - t97 + (-t250 - t278) * t170) * t170 + t256;
t52 = (t170 * t207 + t140) * t170 + t190;
t95 = t170 * (t162 * t218 + t278);
t60 = t189 - t95;
t61 = (-t102 + t235) * t170 + t226;
t179 = (-t51 * t275 + t52 * (-pkin(3) + t218)) * t162 + (t51 * t238 + t52 * t275) * t163 + (t61 * (-t250 + t292) - t60 * (-t271 - t279 - t291) + (-t238 * t60 - t275 * t61) * t163) * t170;
t121 = t170 * t219;
t101 = t125 * t170 - t251;
t100 = t124 * t170 + t164;
t93 = -t247 + (-t133 - t291) * t170;
t92 = -t121 + t199;
t91 = t170 * (-rSges(4,1) * t268 + t146) + t200;
t90 = -t170 * t257 + t256;
t89 = -rSges(6,1) * t114 - rSges(6,2) * t115;
t88 = rSges(6,1) * t112 + rSges(6,2) * t113;
t35 = qJD(5) * t186 + qJD(3);
t16 = qJD(5) * t181;
t1 = [m(3) * (t100 * (-t143 - t294) + t101 * (-t220 - t295) + (t111 - t125 + t247) * t110) + t187 + (t21 * (t196 - t294) + t34 * (t182 + t248) + t22 * (t225 + t261) + (-t193 - t188 - t34 + t296) * t33 + t299) * m(6) + (t51 * (t234 - t294) + t61 * (t243 + t248) + t52 * t225 + t179 + (-t226 - t188 - t61 + t300) * t60) * m(5) + (t90 * (t208 - t294) + t93 * (t199 - t257) + t91 * (t197 - t295) + (-t198 - t93 + t247 + t298) * t92) * m(4); t187 + (t21 * t196 + t22 * (t261 - t292) + (t182 - t211 - t297) * t34 + (-t227 + t210 + t296) * t33 + t299) * m(6) + (t51 * t234 - t52 * t292 + t179 + (t243 - t211 + t95) * t61 + (-t159 + t210 + t300) * t60) * m(5) + (t91 * t197 + t90 * t208 + (-t252 + t298) * t92 + (-t257 + t121) * t93) * m(4) + (-(t110 * t143 + t111 * t220) * t170 - t100 * t143 - t101 * t220 - t110 * t125 + t111 * t124) * m(3); m(6) * t16; m(5) * (t162 * t52 + t163 * t51) + m(6) * (t162 * t22 + t163 * t21); -t173 * (qJD(5) * t183 + t36) / 0.2e1 + t157 * (-t173 * t39 + t183) / 0.2e1 - (qJD(5) * t185 + t157 * t20) * t269 / 0.2e1 + (-t173 * t49 + t202) * t212 + (-t173 * t20 + t185) * t233 + (qJD(5) * t184 + t157 * t19) * t267 / 0.2e1 + (-t173 * t50 + t201) * t213 + (-t173 * t19 + t184) * t230 - t157 * (-t173 * t134 * t157 + ((-t174 * t259 - t176 * t260) * t157 + ((t174 * t194 + t176 * t195) * t172 + t215 * t173) * qJD(5)) * t172) / 0.2e1 + ((t112 * t259 + t113 * t260 - t134 * t269) * t157 + (-t112 * t194 - t113 * t195 + t215 * t269) * t254) * t232 + ((-t114 * t259 - t115 * t260 + t134 * t267) * t157 + (t114 * t194 + t115 * t195 - t215 * t267) * t254) * t231 - (t10 * t162 + t163 * t9) * t264 / 0.2e1 + (t16 * t186 + t35 * t181 + t21 * (t126 * t267 + t173 * t74) + t34 * (t173 * t47 + (-t126 * t270 + t130 * t163) * t172) + t22 * (t126 * t269 - t173 * t209) - t33 * (-t173 * t48 + (t126 * t268 + t130 * t162) * t172) - (-t33 * t88 - t34 * t89) * t157 - (t35 * (-t162 * t89 - t163 * t88) + (-t162 * t33 + t163 * t34) * t203) * t254) * m(6);];
tauc = t1(:);
