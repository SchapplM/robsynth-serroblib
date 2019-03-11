% Calculate kinetic energy for
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:06
% EndTime: 2019-03-09 05:12:10
% DurationCPUTime: 3.90s
% Computational Cost: add. (1849->319), mult. (1730->452), div. (0->0), fcn. (1512->10), ass. (0->170)
t331 = Icges(5,4) + Icges(6,6);
t330 = Icges(5,1) + Icges(6,2);
t329 = -Icges(5,2) - Icges(6,3);
t235 = pkin(10) + qJ(3);
t226 = qJ(4) + t235;
t222 = cos(t226);
t328 = t331 * t222;
t221 = sin(t226);
t327 = t331 * t221;
t326 = Icges(6,4) - Icges(5,5);
t325 = Icges(6,5) - Icges(5,6);
t324 = t329 * t221 + t328;
t323 = t330 * t222 - t327;
t322 = Icges(6,1) + Icges(5,3);
t240 = sin(qJ(1));
t242 = cos(qJ(1));
t321 = t324 * t240 + t325 * t242;
t320 = -t325 * t240 + t324 * t242;
t319 = t323 * t240 + t326 * t242;
t318 = -t326 * t240 + t323 * t242;
t317 = t329 * t222 - t327;
t316 = t330 * t221 + t328;
t315 = t325 * t221 - t326 * t222;
t195 = V_base(5) + (-qJD(3) - qJD(4)) * t242;
t219 = qJD(3) * t240 + V_base(4);
t196 = qJD(4) * t240 + t219;
t227 = V_base(6) + qJD(1);
t314 = (t317 * t221 + t316 * t222) * t227 + (-t320 * t221 + t318 * t222) * t196 + (-t321 * t221 + t319 * t222) * t195;
t313 = (-t326 * t221 - t325 * t222) * t227 + (t322 * t240 + t315 * t242) * t196 + (t315 * t240 - t322 * t242) * t195;
t236 = sin(pkin(10));
t309 = pkin(2) * t236;
t224 = sin(t235);
t308 = pkin(3) * t224;
t307 = pkin(9) * t221;
t237 = cos(pkin(10));
t306 = t237 * pkin(2);
t305 = Icges(2,4) * t240;
t304 = Icges(3,4) * t236;
t303 = Icges(3,4) * t237;
t302 = Icges(4,4) * t224;
t225 = cos(t235);
t301 = Icges(4,4) * t225;
t296 = t222 * t240;
t295 = t222 * t242;
t239 = sin(qJ(6));
t294 = t239 * t242;
t293 = t240 * t239;
t241 = cos(qJ(6));
t292 = t240 * t241;
t291 = t241 * t242;
t151 = -pkin(7) * t242 + t240 * t306;
t214 = t240 * pkin(1) - qJ(2) * t242;
t289 = -t151 - t214;
t288 = pkin(3) * t225;
t286 = qJD(5) * t221;
t285 = qJD(6) * t222;
t284 = V_base(4) * t214 + V_base(3);
t283 = V_base(5) * pkin(6) + V_base(1);
t121 = -pkin(8) * t242 + t240 * t288;
t280 = -t121 + t289;
t279 = qJD(2) * t240 + t283;
t272 = pkin(4) * t222 + qJ(5) * t221;
t163 = t272 * t240;
t278 = -t163 + t280;
t277 = V_base(5) * t309 + t279;
t276 = rSges(3,1) * t237 - rSges(3,2) * t236;
t275 = rSges(4,1) * t225 - rSges(4,2) * t224;
t274 = rSges(5,1) * t222 - rSges(5,2) * t221;
t273 = -rSges(6,2) * t222 + rSges(6,3) * t221;
t271 = Icges(3,1) * t237 - t304;
t270 = Icges(4,1) * t225 - t302;
t268 = -Icges(3,2) * t236 + t303;
t267 = -Icges(4,2) * t224 + t301;
t264 = Icges(3,5) * t237 - Icges(3,6) * t236;
t263 = Icges(4,5) * t225 - Icges(4,6) * t224;
t216 = pkin(1) * t242 + t240 * qJ(2);
t259 = -qJD(2) * t242 + t227 * t216 + V_base(2);
t218 = -qJD(3) * t242 + V_base(5);
t258 = t218 * t308 + t277;
t183 = pkin(4) * t221 - qJ(5) * t222;
t257 = t195 * t183 + t242 * t286 + t258;
t254 = (-Icges(4,3) * t242 + t240 * t263) * t218 + (Icges(4,3) * t240 + t242 * t263) * t219 + (Icges(4,5) * t224 + Icges(4,6) * t225) * t227;
t152 = pkin(7) * t240 + t242 * t306;
t253 = V_base(4) * t151 + (-t152 - t216) * V_base(5) + t284;
t252 = (-Icges(3,3) * t242 + t240 * t264) * V_base(5) + (Icges(3,3) * t240 + t242 * t264) * V_base(4) + (Icges(3,5) * t236 + Icges(3,6) * t237) * t227;
t122 = pkin(8) * t240 + t242 * t288;
t251 = t219 * t121 - t122 * t218 + t253;
t250 = t227 * t152 + (-pkin(6) - t309) * V_base(4) + t259;
t249 = -qJD(5) * t222 + t196 * t163 + t251;
t248 = t227 * t122 - t219 * t308 + t250;
t164 = t272 * t242;
t247 = t227 * t164 + t240 * t286 + t248;
t155 = -Icges(4,6) * t242 + t240 * t267;
t156 = Icges(4,6) * t240 + t242 * t267;
t157 = -Icges(4,5) * t242 + t240 * t270;
t158 = Icges(4,5) * t240 + t242 * t270;
t191 = Icges(4,2) * t225 + t302;
t192 = Icges(4,1) * t224 + t301;
t244 = (-t156 * t224 + t158 * t225) * t219 + (-t155 * t224 + t157 * t225) * t218 + (-t191 * t224 + t192 * t225) * t227;
t167 = -Icges(3,6) * t242 + t240 * t268;
t168 = Icges(3,6) * t240 + t242 * t268;
t169 = -Icges(3,5) * t242 + t240 * t271;
t170 = Icges(3,5) * t240 + t242 * t271;
t205 = Icges(3,2) * t237 + t304;
t206 = Icges(3,1) * t236 + t303;
t243 = (-t168 * t236 + t170 * t237) * V_base(4) + (-t167 * t236 + t169 * t237) * V_base(5) + (-t205 * t236 + t206 * t237) * t227;
t232 = Icges(2,4) * t242;
t217 = rSges(2,1) * t242 - t240 * rSges(2,2);
t215 = t240 * rSges(2,1) + rSges(2,2) * t242;
t213 = Icges(2,1) * t242 - t305;
t212 = Icges(2,1) * t240 + t232;
t211 = -Icges(2,2) * t240 + t232;
t210 = Icges(2,2) * t242 + t305;
t209 = Icges(2,5) * t242 - Icges(2,6) * t240;
t208 = Icges(2,5) * t240 + Icges(2,6) * t242;
t207 = rSges(3,1) * t236 + rSges(3,2) * t237;
t201 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t200 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t199 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t194 = qJD(6) * t221 + t227;
t193 = rSges(4,1) * t224 + rSges(4,2) * t225;
t187 = -pkin(5) * t242 + pkin(9) * t296;
t186 = t240 * pkin(5) + pkin(9) * t295;
t185 = rSges(5,1) * t221 + rSges(5,2) * t222;
t184 = -rSges(6,2) * t221 - rSges(6,3) * t222;
t176 = t221 * t293 - t291;
t175 = t221 * t292 + t294;
t174 = t221 * t294 + t292;
t173 = t221 * t291 - t293;
t172 = t240 * rSges(3,3) + t242 * t276;
t171 = -rSges(3,3) * t242 + t240 * t276;
t162 = t240 * rSges(4,3) + t242 * t275;
t161 = -rSges(4,3) * t242 + t240 * t275;
t160 = t242 * t285 + t196;
t159 = t240 * t285 + t195;
t150 = -rSges(6,1) * t242 + t240 * t273;
t149 = t240 * rSges(6,1) + t242 * t273;
t148 = t240 * rSges(5,3) + t242 * t274;
t147 = -rSges(5,3) * t242 + t240 * t274;
t146 = V_base(5) * rSges(2,3) - t215 * t227 + t283;
t145 = t217 * t227 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t131 = t215 * V_base(4) - t217 * V_base(5) + V_base(3);
t129 = rSges(7,3) * t221 + (-rSges(7,1) * t239 - rSges(7,2) * t241) * t222;
t127 = Icges(7,5) * t221 + (-Icges(7,1) * t239 - Icges(7,4) * t241) * t222;
t126 = Icges(7,6) * t221 + (-Icges(7,4) * t239 - Icges(7,2) * t241) * t222;
t125 = Icges(7,3) * t221 + (-Icges(7,5) * t239 - Icges(7,6) * t241) * t222;
t118 = rSges(7,1) * t176 + rSges(7,2) * t175 + rSges(7,3) * t296;
t117 = t174 * rSges(7,1) + t173 * rSges(7,2) + rSges(7,3) * t295;
t116 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t296;
t115 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t295;
t114 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t296;
t113 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t295;
t112 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t296;
t111 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t295;
t110 = t207 * V_base(5) + (-t171 - t214) * t227 + t279;
t109 = t227 * t172 + (-pkin(6) - t207) * V_base(4) + t259;
t108 = t171 * V_base(4) + (-t172 - t216) * V_base(5) + t284;
t107 = t193 * t218 + (-t161 + t289) * t227 + t277;
t106 = t227 * t162 - t219 * t193 + t250;
t105 = t161 * t219 - t162 * t218 + t253;
t104 = t185 * t195 + (-t147 + t280) * t227 + t258;
t103 = t227 * t148 - t196 * t185 + t248;
t102 = t184 * t195 + (-t150 + t278) * t227 + t257;
t101 = t227 * t149 + (-t183 - t184) * t196 + t247;
t100 = t147 * t196 - t148 * t195 + t251;
t99 = t150 * t196 + (-t149 - t164) * t195 + t249;
t98 = t195 * t307 - t118 * t194 + t129 * t159 + (-t187 + t278) * t227 + t257;
t97 = t194 * t117 - t160 * t129 + t227 * t186 + (-t183 - t307) * t196 + t247;
t96 = -t117 * t159 + t118 * t160 + t187 * t196 + (-t164 - t186) * t195 + t249;
t1 = m(2) * (t131 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + t219 * (t254 * t240 + t244 * t242) / 0.2e1 + t218 * (t244 * t240 - t254 * t242) / 0.2e1 + t194 * ((t111 * t160 + t112 * t159 + t125 * t194) * t221 + ((-t113 * t241 - t115 * t239) * t160 + (-t114 * t241 - t116 * t239) * t159 + (-t126 * t241 - t127 * t239) * t194) * t222) / 0.2e1 + m(1) * (t199 ^ 2 + t200 ^ 2 + t201 ^ 2) / 0.2e1 + t160 * ((t111 * t295 + t173 * t113 + t174 * t115) * t160 + (t112 * t295 + t173 * t114 + t174 * t116) * t159 + (t125 * t295 + t173 * t126 + t174 * t127) * t194) / 0.2e1 + t159 * ((t111 * t296 + t113 * t175 + t115 * t176) * t160 + (t112 * t296 + t114 * t175 + t116 * t176) * t159 + (t125 * t296 + t126 * t175 + t127 * t176) * t194) / 0.2e1 + (t314 * t240 - t313 * t242) * t195 / 0.2e1 + (t313 * t240 + t314 * t242) * t196 / 0.2e1 + (t209 * t227 + t252 * t240 + t243 * t242 + (-t240 * t210 + t212 * t242 + Icges(1,4)) * V_base(5) + (-t240 * t211 + t213 * t242 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t208 * t227 + t243 * t240 - t252 * t242 + (t210 * t242 + t240 * t212 + Icges(1,2)) * V_base(5) + (t211 * t242 + t240 * t213 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t156 * t225 + t158 * t224) * t219 + (t155 * t225 + t157 * t224) * t218 + (t167 * t237 + t169 * t236 + t208) * V_base(5) + (t168 * t237 + t170 * t236 + t209) * V_base(4) + (t318 * t221 + t320 * t222) * t196 + (t319 * t221 + t321 * t222) * t195 + (t191 * t225 + t192 * t224 + t205 * t237 + t206 * t236 + t316 * t221 - t317 * t222 + Icges(2,3)) * t227) * t227 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
