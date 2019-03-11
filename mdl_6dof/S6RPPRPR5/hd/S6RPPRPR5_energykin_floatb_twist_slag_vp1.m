% Calculate kinetic energy for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:27
% EndTime: 2019-03-09 01:48:30
% DurationCPUTime: 2.82s
% Computational Cost: add. (984->307), mult. (1518->420), div. (0->0), fcn. (1356->8), ass. (0->146)
t275 = Icges(2,4) + Icges(3,6) - Icges(4,6);
t274 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t273 = -Icges(3,4) + Icges(2,5) + Icges(4,5);
t272 = Icges(4,4) + Icges(3,5) - Icges(2,6);
t271 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t207 = sin(qJ(1));
t270 = t275 * t207;
t209 = cos(qJ(1));
t269 = t275 * t209;
t268 = t274 * t209 - t270;
t267 = t274 * t207 + t269;
t266 = -t271 * t209 - t270;
t265 = t271 * t207 - t269;
t262 = -pkin(2) - pkin(6);
t258 = pkin(7) * t207;
t257 = pkin(7) * t209;
t204 = cos(pkin(9));
t256 = pkin(5) * t204;
t206 = sin(qJ(4));
t254 = Icges(5,4) * t206;
t208 = cos(qJ(4));
t253 = Icges(5,4) * t208;
t249 = qJ(3) * t207;
t248 = qJ(3) * t209;
t203 = sin(pkin(9));
t247 = t203 * t209;
t246 = t206 * t209;
t202 = pkin(9) + qJ(6);
t189 = sin(t202);
t245 = t207 * t189;
t190 = cos(t202);
t244 = t207 * t190;
t243 = t207 * t203;
t242 = t207 * t204;
t241 = t207 * t208;
t240 = t208 * t209;
t238 = qJD(5) * t208;
t237 = qJD(6) * t208;
t172 = t207 * pkin(1) - qJ(2) * t209;
t236 = V_base(4) * t172 + V_base(3);
t235 = V_base(5) * pkin(6) + V_base(1);
t183 = qJD(4) * t209 + V_base(5);
t191 = V_base(6) + qJD(1);
t232 = V_base(4) * t249 + t236;
t231 = qJD(2) * t207 + t235;
t230 = -t172 - t249;
t178 = pkin(1) * t209 + t207 * qJ(2);
t229 = -t178 - t248;
t184 = -qJD(4) * t207 + V_base(4);
t228 = rSges(5,1) * t206 + rSges(5,2) * t208;
t227 = pkin(4) * t206 - qJ(5) * t208;
t226 = Icges(5,1) * t206 + t253;
t225 = Icges(5,2) * t208 + t254;
t224 = Icges(5,5) * t206 + Icges(5,6) * t208;
t223 = -qJD(2) * t209 + t191 * t178 + V_base(2);
t222 = V_base(5) * pkin(2) + qJD(3) * t209 + t231;
t221 = t230 - t257;
t220 = V_base(5) * pkin(3) + t222;
t144 = t227 * t207;
t219 = -t144 + t221;
t218 = qJD(3) * t207 + t191 * t248 + t223;
t217 = (Icges(5,3) * t209 + t207 * t224) * t183 + (-Icges(5,3) * t207 + t209 * t224) * t184 + (Icges(5,5) * t208 - Icges(5,6) * t206) * t191;
t216 = -pkin(8) * t208 + t206 * t256;
t176 = pkin(4) * t208 + qJ(5) * t206;
t215 = t183 * t176 - t209 * t238 + t220;
t214 = V_base(4) * t257 + t232 + (t229 + t258) * V_base(5);
t213 = (-pkin(3) + t262) * V_base(4) + t218;
t212 = qJD(5) * t206 + t184 * t144 + t214;
t145 = t227 * t209;
t211 = t191 * t145 - t207 * t238 + t213;
t130 = Icges(5,6) * t209 + t207 * t225;
t131 = -Icges(5,6) * t207 + t209 * t225;
t132 = Icges(5,5) * t209 + t207 * t226;
t133 = -Icges(5,5) * t207 + t209 * t226;
t161 = -Icges(5,2) * t206 + t253;
t168 = Icges(5,1) * t208 - t254;
t210 = (t131 * t208 + t133 * t206) * t184 + (t130 * t208 + t132 * t206) * t183 + (t161 * t208 + t168 * t206) * t191;
t181 = rSges(2,1) * t209 - t207 * rSges(2,2);
t180 = -rSges(3,2) * t209 + t207 * rSges(3,3);
t179 = -rSges(4,2) * t209 + t207 * rSges(4,3);
t177 = rSges(5,1) * t208 - rSges(5,2) * t206;
t175 = t207 * rSges(2,1) + rSges(2,2) * t209;
t174 = -t207 * rSges(3,2) - rSges(3,3) * t209;
t173 = t207 * rSges(4,2) + rSges(4,3) * t209;
t171 = qJD(6) * t206 + t191;
t149 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t148 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t147 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t142 = -t209 * t237 + t184;
t141 = -t207 * t237 + t183;
t140 = t204 * t246 - t243;
t139 = -t203 * t246 - t242;
t138 = t206 * t242 + t247;
t137 = t204 * t209 - t206 * t243;
t135 = -t207 * rSges(5,3) + t209 * t228;
t134 = rSges(5,3) * t209 + t207 * t228;
t127 = t190 * t246 - t245;
t126 = -t189 * t246 - t244;
t125 = t189 * t209 + t206 * t244;
t124 = t190 * t209 - t206 * t245;
t123 = rSges(6,3) * t206 + (rSges(6,1) * t204 - rSges(6,2) * t203) * t208;
t121 = Icges(6,5) * t206 + (Icges(6,1) * t204 - Icges(6,4) * t203) * t208;
t120 = Icges(6,6) * t206 + (Icges(6,4) * t204 - Icges(6,2) * t203) * t208;
t119 = Icges(6,3) * t206 + (Icges(6,5) * t204 - Icges(6,6) * t203) * t208;
t117 = rSges(7,3) * t206 + (rSges(7,1) * t190 - rSges(7,2) * t189) * t208;
t116 = Icges(7,5) * t206 + (Icges(7,1) * t190 - Icges(7,4) * t189) * t208;
t115 = Icges(7,6) * t206 + (Icges(7,4) * t190 - Icges(7,2) * t189) * t208;
t114 = Icges(7,3) * t206 + (Icges(7,5) * t190 - Icges(7,6) * t189) * t208;
t113 = pkin(8) * t206 + t208 * t256;
t112 = V_base(5) * rSges(2,3) - t175 * t191 + t235;
t111 = t181 * t191 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t110 = t175 * V_base(4) - t181 * V_base(5) + V_base(3);
t109 = -pkin(5) * t243 + t209 * t216;
t108 = pkin(5) * t247 + t207 * t216;
t107 = t140 * rSges(6,1) + t139 * rSges(6,2) - rSges(6,3) * t240;
t106 = rSges(6,1) * t138 + rSges(6,2) * t137 - rSges(6,3) * t241;
t105 = Icges(6,1) * t140 + Icges(6,4) * t139 - Icges(6,5) * t240;
t104 = Icges(6,1) * t138 + Icges(6,4) * t137 - Icges(6,5) * t241;
t103 = Icges(6,4) * t140 + Icges(6,2) * t139 - Icges(6,6) * t240;
t102 = Icges(6,4) * t138 + Icges(6,2) * t137 - Icges(6,6) * t241;
t101 = Icges(6,5) * t140 + Icges(6,6) * t139 - Icges(6,3) * t240;
t100 = Icges(6,5) * t138 + Icges(6,6) * t137 - Icges(6,3) * t241;
t99 = V_base(5) * rSges(3,1) + (-t172 - t174) * t191 + t231;
t98 = t191 * t180 + (-rSges(3,1) - pkin(6)) * V_base(4) + t223;
t97 = t127 * rSges(7,1) + t126 * rSges(7,2) - rSges(7,3) * t240;
t96 = rSges(7,1) * t125 + rSges(7,2) * t124 - rSges(7,3) * t241;
t95 = Icges(7,1) * t127 + Icges(7,4) * t126 - Icges(7,5) * t240;
t94 = Icges(7,1) * t125 + Icges(7,4) * t124 - Icges(7,5) * t241;
t93 = Icges(7,4) * t127 + Icges(7,2) * t126 - Icges(7,6) * t240;
t92 = Icges(7,4) * t125 + Icges(7,2) * t124 - Icges(7,6) * t241;
t91 = Icges(7,5) * t127 + Icges(7,6) * t126 - Icges(7,3) * t240;
t90 = Icges(7,5) * t125 + Icges(7,6) * t124 - Icges(7,3) * t241;
t89 = t174 * V_base(4) + (-t178 - t180) * V_base(5) + t236;
t88 = V_base(5) * rSges(4,1) + (-t179 + t230) * t191 + t222;
t87 = t191 * t173 + (-rSges(4,1) + t262) * V_base(4) + t218;
t86 = V_base(4) * t179 + (-t173 + t229) * V_base(5) + t232;
t85 = t183 * t177 + (-t134 + t221) * t191 + t220;
t84 = -t184 * t177 + (t135 - t258) * t191 + t213;
t83 = t184 * t134 - t183 * t135 + t214;
t82 = t183 * t123 + (-t106 + t219) * t191 + t215;
t81 = (t107 - t258) * t191 + (-t123 - t176) * t184 + t211;
t80 = t184 * t106 + (-t107 - t145) * t183 + t212;
t79 = t183 * t113 + t141 * t117 - t171 * t96 + (-t108 + t219) * t191 + t215;
t78 = -t142 * t117 + t171 * t97 + (t109 - t258) * t191 + (-t113 - t176) * t184 + t211;
t77 = t184 * t108 - t141 * t97 + t142 * t96 + (-t109 - t145) * t183 + t212;
t1 = t141 * ((t124 * t93 + t125 * t95 - t241 * t91) * t142 + (t124 * t92 + t125 * t94 - t90 * t241) * t141 + (-t114 * t241 + t115 * t124 + t116 * t125) * t171) / 0.2e1 + t142 * ((t126 * t93 + t127 * t95 - t91 * t240) * t142 + (t126 * t92 + t127 * t94 - t240 * t90) * t141 + (-t114 * t240 + t126 * t115 + t127 * t116) * t171) / 0.2e1 + m(2) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + t171 * ((t114 * t171 + t90 * t141 + t91 * t142) * t206 + ((-t189 * t93 + t190 * t95) * t142 + (-t189 * t92 + t190 * t94) * t141 + (-t115 * t189 + t116 * t190) * t171) * t208) / 0.2e1 + m(1) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + ((-t101 * t241 + t103 * t137 + t105 * t138) * t184 + (-t100 * t241 + t137 * t102 + t138 * t104) * t183 + (-t119 * t241 + t120 * t137 + t121 * t138) * t191 + t210 * t207 + t217 * t209) * t183 / 0.2e1 + ((-t101 * t240 + t139 * t103 + t140 * t105) * t184 + (-t100 * t240 + t139 * t102 + t140 * t104) * t183 + (-t119 * t240 + t139 * t120 + t140 * t121) * t191 - t217 * t207 + t210 * t209) * t184 / 0.2e1 + ((t207 * t266 + t209 * t267 + Icges(1,4)) * V_base(5) + (t265 * t207 + t268 * t209 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t207 - t266 * t209 + Icges(1,2)) * V_base(5) + (t207 * t268 - t265 * t209 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t131 * t206 + t133 * t208) * t184 + (-t130 * t206 + t132 * t208) * t183 + (t100 * t183 + t101 * t184) * t206 + ((-t103 * t203 + t105 * t204) * t184 + (-t102 * t203 + t104 * t204) * t183) * t208 + (Icges(2,3) + Icges(4,1) + Icges(3,1) + (-t120 * t203 + t121 * t204 + t168) * t208 + (-t161 + t119) * t206) * t191) * t191 / 0.2e1 + t191 * V_base(5) * (t273 * t207 - t272 * t209) + t191 * V_base(4) * (t272 * t207 + t273 * t209) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
