% Calculate kinetic energy for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:22:57
% EndTime: 2019-03-09 07:23:01
% DurationCPUTime: 4.16s
% Computational Cost: add. (1624->359), mult. (2164->532), div. (0->0), fcn. (2058->10), ass. (0->169)
t312 = Icges(2,4) + Icges(3,6);
t311 = Icges(2,1) + Icges(3,2);
t310 = -Icges(3,4) + Icges(2,5);
t309 = Icges(3,5) - Icges(2,6);
t308 = Icges(2,2) + Icges(3,3);
t246 = cos(qJ(1));
t307 = t312 * t246;
t243 = sin(qJ(1));
t306 = t312 * t243;
t305 = -t246 * t308 - t306;
t304 = t243 * t308 - t307;
t303 = t243 * t311 + t307;
t302 = t246 * t311 - t306;
t242 = sin(qJ(3));
t245 = cos(qJ(3));
t240 = qJ(4) + qJ(5);
t233 = cos(t240);
t277 = pkin(5) * t233;
t299 = -pkin(10) * t245 + t242 * t277;
t244 = cos(qJ(4));
t292 = t244 * pkin(4);
t298 = -pkin(9) * t245 + t242 * t292;
t289 = Icges(4,4) * t242;
t262 = Icges(4,2) * t245 + t289;
t163 = Icges(4,6) * t246 + t243 * t262;
t164 = Icges(4,6) * t243 - t246 * t262;
t288 = Icges(4,4) * t245;
t263 = Icges(4,1) * t242 + t288;
t166 = Icges(4,5) * t246 + t243 * t263;
t167 = Icges(4,5) * t243 - t246 * t263;
t200 = -Icges(4,2) * t242 + t288;
t205 = Icges(4,1) * t245 - t289;
t218 = qJD(3) * t243 + V_base(5);
t219 = qJD(3) * t246 + V_base(4);
t226 = V_base(6) + qJD(1);
t297 = (t163 * t245 + t166 * t242) * t219 + (t164 * t245 + t167 * t242) * t218 + (t200 * t245 + t205 * t242) * t226;
t294 = pkin(7) * t243;
t293 = pkin(7) * t246;
t241 = sin(qJ(4));
t285 = t241 * t243;
t284 = t241 * t246;
t283 = t242 * t243;
t282 = t242 * t246;
t281 = t243 * t244;
t280 = t243 * t245;
t279 = t244 * t246;
t278 = t245 * t246;
t275 = qJD(4) * t245;
t274 = -qJD(4) - qJD(5);
t209 = pkin(1) * t243 - qJ(2) * t246;
t273 = V_base(4) * t209 + V_base(3);
t272 = V_base(5) * pkin(6) + V_base(1);
t232 = sin(t240);
t269 = pkin(5) * t232;
t268 = -t209 - t294;
t267 = qJD(2) * t243 + t272;
t177 = t246 * t275 + t218;
t208 = qJD(4) * t242 + t226;
t266 = V_base(5) * pkin(2) + t267;
t265 = pkin(3) * t242 - pkin(8) * t245;
t264 = rSges(4,1) * t242 + rSges(4,2) * t245;
t152 = qJD(5) * t278 + t177;
t261 = Icges(4,5) * t242 + Icges(4,6) * t245;
t186 = qJD(5) * t242 + t208;
t213 = pkin(1) * t246 + qJ(2) * t243;
t257 = -qJD(2) * t246 + t226 * t213 + V_base(2);
t256 = (Icges(4,3) * t246 + t243 * t261) * t219 + (Icges(4,3) * t243 - t246 * t261) * t218 + (Icges(4,5) * t245 - Icges(4,6) * t242) * t226;
t255 = V_base(4) * t294 + (-t213 - t293) * V_base(5) + t273;
t254 = t226 * t293 + (-pkin(2) - pkin(6)) * V_base(4) + t257;
t185 = t265 * t246;
t216 = t245 * pkin(3) + t242 * pkin(8);
t253 = t218 * t216 + (t185 + t268) * t226 + t266;
t184 = t265 * t243;
t252 = -t184 * t218 - t219 * t185 + t255;
t251 = t226 * t184 - t216 * t219 + t254;
t134 = pkin(4) * t285 - t246 * t298;
t146 = pkin(9) * t242 + t245 * t292;
t250 = -t134 * t208 + t177 * t146 + t253;
t133 = pkin(4) * t284 + t243 * t298;
t178 = -t243 * t275 + t219;
t249 = -t133 * t177 + t178 * t134 + t252;
t248 = t208 * t133 - t146 * t178 + t251;
t235 = qJ(6) + t240;
t224 = cos(t235);
t223 = sin(t235);
t215 = rSges(2,1) * t246 - rSges(2,2) * t243;
t214 = -rSges(3,2) * t246 + rSges(3,3) * t243;
t212 = rSges(4,1) * t245 - rSges(4,2) * t242;
t211 = rSges(2,1) * t243 + rSges(2,2) * t246;
t210 = -rSges(3,2) * t243 - rSges(3,3) * t246;
t191 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t190 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t189 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t182 = -t242 * t279 + t285;
t181 = t241 * t282 + t281;
t180 = t242 * t281 + t284;
t179 = -t241 * t283 + t279;
t175 = qJD(6) * t242 + t186;
t174 = t232 * t243 - t233 * t282;
t173 = t232 * t282 + t233 * t243;
t172 = t232 * t246 + t233 * t283;
t171 = -t232 * t283 + t233 * t246;
t170 = rSges(4,3) * t243 - t246 * t264;
t169 = rSges(5,3) * t242 + (rSges(5,1) * t244 - rSges(5,2) * t241) * t245;
t168 = rSges(4,3) * t246 + t243 * t264;
t165 = Icges(5,5) * t242 + (Icges(5,1) * t244 - Icges(5,4) * t241) * t245;
t162 = Icges(5,6) * t242 + (Icges(5,4) * t244 - Icges(5,2) * t241) * t245;
t159 = Icges(5,3) * t242 + (Icges(5,5) * t244 - Icges(5,6) * t241) * t245;
t157 = t223 * t243 - t224 * t282;
t156 = t223 * t282 + t224 * t243;
t155 = t223 * t246 + t224 * t283;
t154 = -t223 * t283 + t224 * t246;
t153 = t274 * t280 + t219;
t151 = rSges(6,3) * t242 + (rSges(6,1) * t233 - rSges(6,2) * t232) * t245;
t149 = Icges(6,5) * t242 + (Icges(6,1) * t233 - Icges(6,4) * t232) * t245;
t148 = Icges(6,6) * t242 + (Icges(6,4) * t233 - Icges(6,2) * t232) * t245;
t147 = Icges(6,3) * t242 + (Icges(6,5) * t233 - Icges(6,6) * t232) * t245;
t145 = rSges(7,3) * t242 + (rSges(7,1) * t224 - rSges(7,2) * t223) * t245;
t144 = Icges(7,5) * t242 + (Icges(7,1) * t224 - Icges(7,4) * t223) * t245;
t143 = Icges(7,6) * t242 + (Icges(7,4) * t224 - Icges(7,2) * t223) * t245;
t142 = Icges(7,3) * t242 + (Icges(7,5) * t224 - Icges(7,6) * t223) * t245;
t141 = V_base(5) * rSges(2,3) - t211 * t226 + t272;
t140 = t215 * t226 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t139 = t211 * V_base(4) - t215 * V_base(5) + V_base(3);
t138 = (-qJD(6) + t274) * t280 + t219;
t137 = qJD(6) * t278 + t152;
t136 = pkin(10) * t242 + t245 * t277;
t132 = rSges(5,1) * t182 + rSges(5,2) * t181 + rSges(5,3) * t278;
t131 = rSges(5,1) * t180 + rSges(5,2) * t179 - rSges(5,3) * t280;
t130 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t278;
t129 = Icges(5,1) * t180 + Icges(5,4) * t179 - Icges(5,5) * t280;
t128 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t278;
t127 = Icges(5,4) * t180 + Icges(5,2) * t179 - Icges(5,6) * t280;
t126 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t278;
t125 = Icges(5,5) * t180 + Icges(5,6) * t179 - Icges(5,3) * t280;
t124 = V_base(5) * rSges(3,1) + (-t209 - t210) * t226 + t267;
t123 = t214 * t226 + (-rSges(3,1) - pkin(6)) * V_base(4) + t257;
t122 = rSges(6,1) * t174 + rSges(6,2) * t173 + rSges(6,3) * t278;
t121 = rSges(6,1) * t172 + rSges(6,2) * t171 - rSges(6,3) * t280;
t120 = Icges(6,1) * t174 + Icges(6,4) * t173 + Icges(6,5) * t278;
t119 = Icges(6,1) * t172 + Icges(6,4) * t171 - Icges(6,5) * t280;
t118 = Icges(6,4) * t174 + Icges(6,2) * t173 + Icges(6,6) * t278;
t117 = Icges(6,4) * t172 + Icges(6,2) * t171 - Icges(6,6) * t280;
t116 = Icges(6,5) * t174 + Icges(6,6) * t173 + Icges(6,3) * t278;
t115 = Icges(6,5) * t172 + Icges(6,6) * t171 - Icges(6,3) * t280;
t113 = t210 * V_base(4) + (-t213 - t214) * V_base(5) + t273;
t112 = rSges(7,1) * t157 + rSges(7,2) * t156 + rSges(7,3) * t278;
t111 = rSges(7,1) * t155 + rSges(7,2) * t154 - rSges(7,3) * t280;
t110 = Icges(7,1) * t157 + Icges(7,4) * t156 + Icges(7,5) * t278;
t109 = Icges(7,1) * t155 + Icges(7,4) * t154 - Icges(7,5) * t280;
t108 = Icges(7,4) * t157 + Icges(7,2) * t156 + Icges(7,6) * t278;
t107 = Icges(7,4) * t155 + Icges(7,2) * t154 - Icges(7,6) * t280;
t106 = Icges(7,5) * t157 + Icges(7,6) * t156 + Icges(7,3) * t278;
t105 = Icges(7,5) * t155 + Icges(7,6) * t154 - Icges(7,3) * t280;
t103 = t269 * t243 - t246 * t299;
t102 = t243 * t299 + t269 * t246;
t101 = t212 * t218 + (-t170 + t268) * t226 + t266;
t100 = t168 * t226 - t212 * t219 + t254;
t99 = -t168 * t218 + t170 * t219 + t255;
t98 = -t132 * t208 + t169 * t177 + t253;
t97 = t131 * t208 - t169 * t178 + t251;
t96 = -t131 * t177 + t132 * t178 + t252;
t95 = -t122 * t186 + t151 * t152 + t250;
t94 = t121 * t186 - t151 * t153 + t248;
t93 = -t121 * t152 + t122 * t153 + t249;
t92 = -t103 * t186 - t112 * t175 + t136 * t152 + t137 * t145 + t250;
t91 = t102 * t186 + t111 * t175 - t136 * t153 - t138 * t145 + t248;
t90 = -t102 * t152 + t103 * t153 - t111 * t137 + t112 * t138 + t249;
t1 = t178 * ((-t125 * t280 + t179 * t127 + t180 * t129) * t178 + (-t126 * t280 + t128 * t179 + t130 * t180) * t177 + (-t159 * t280 + t162 * t179 + t165 * t180) * t208) / 0.2e1 + t153 * ((-t115 * t280 + t171 * t117 + t172 * t119) * t153 + (-t116 * t280 + t118 * t171 + t120 * t172) * t152 + (-t147 * t280 + t148 * t171 + t149 * t172) * t186) / 0.2e1 + t138 * ((-t105 * t280 + t154 * t107 + t155 * t109) * t138 + (-t106 * t280 + t108 * t154 + t110 * t155) * t137 + (-t142 * t280 + t143 * t154 + t144 * t155) * t175) / 0.2e1 + m(3) * (t113 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(5) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + t219 * (t297 * t243 + t256 * t246) / 0.2e1 + t218 * (t256 * t243 - t297 * t246) / 0.2e1 + t177 * ((t125 * t278 + t127 * t181 + t129 * t182) * t178 + (t126 * t278 + t181 * t128 + t182 * t130) * t177 + (t159 * t278 + t162 * t181 + t165 * t182) * t208) / 0.2e1 + t152 * ((t115 * t278 + t117 * t173 + t119 * t174) * t153 + (t116 * t278 + t173 * t118 + t174 * t120) * t152 + (t147 * t278 + t148 * t173 + t149 * t174) * t186) / 0.2e1 + t137 * ((t105 * t278 + t107 * t156 + t109 * t157) * t138 + (t106 * t278 + t156 * t108 + t157 * t110) * t137 + (t142 * t278 + t143 * t156 + t144 * t157) * t175) / 0.2e1 + t208 * ((t125 * t178 + t126 * t177 + t159 * t208) * t242 + ((-t127 * t241 + t129 * t244) * t178 + (-t128 * t241 + t130 * t244) * t177 + (-t162 * t241 + t165 * t244) * t208) * t245) / 0.2e1 + t186 * ((t115 * t153 + t116 * t152 + t147 * t186) * t242 + ((-t117 * t232 + t119 * t233) * t153 + (-t118 * t232 + t120 * t233) * t152 + (-t148 * t232 + t149 * t233) * t186) * t245) / 0.2e1 + m(2) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(1) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + t175 * ((t105 * t138 + t106 * t137 + t142 * t175) * t242 + ((-t107 * t223 + t109 * t224) * t138 + (-t108 * t223 + t110 * t224) * t137 + (-t143 * t223 + t144 * t224) * t175) * t245) / 0.2e1 + ((-t163 * t242 + t166 * t245) * t219 + (-t164 * t242 + t167 * t245) * t218 + (-t242 * t200 + t245 * t205 + Icges(3,1) + Icges(2,3)) * t226) * t226 / 0.2e1 + ((t243 * t305 + t303 * t246 + Icges(1,4)) * V_base(5) + (t304 * t243 + t302 * t246 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t303 * t243 - t305 * t246 + Icges(1,2)) * V_base(5) + (t243 * t302 - t246 * t304 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t226 * (t243 * t310 - t246 * t309) + V_base(4) * t226 * (t243 * t309 + t310 * t246) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
