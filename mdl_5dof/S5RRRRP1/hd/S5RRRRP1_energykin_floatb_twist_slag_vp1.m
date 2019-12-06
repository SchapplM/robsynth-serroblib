% Calculate kinetic energy for
% S5RRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:41
% EndTime: 2019-12-05 18:44:43
% DurationCPUTime: 2.39s
% Computational Cost: add. (1271->231), mult. (1228->329), div. (0->0), fcn. (1010->8), ass. (0->129)
t266 = Icges(5,4) + Icges(6,4);
t265 = Icges(5,1) + Icges(6,1);
t264 = Icges(5,2) + Icges(6,2);
t178 = qJ(2) + qJ(3);
t174 = qJ(4) + t178;
t165 = cos(t174);
t263 = t266 * t165;
t164 = sin(t174);
t262 = t266 * t164;
t261 = Icges(5,5) + Icges(6,5);
t260 = Icges(5,6) + Icges(6,6);
t259 = -t264 * t164 + t263;
t258 = t265 * t165 - t262;
t257 = rSges(6,1) + pkin(4);
t180 = sin(qJ(1));
t182 = cos(qJ(1));
t256 = t259 * t180 - t260 * t182;
t255 = t260 * t180 + t259 * t182;
t254 = t258 * t180 - t261 * t182;
t253 = t261 * t180 + t258 * t182;
t252 = Icges(5,3) + Icges(6,3);
t251 = t264 * t165 + t262;
t250 = t265 * t164 + t263;
t249 = -t260 * t164 + t261 * t165;
t248 = rSges(6,3) + qJ(5);
t247 = -rSges(6,2) * t164 + t257 * t165;
t221 = -qJD(2) - qJD(3);
t124 = V_base(5) + (-qJD(4) + t221) * t182;
t162 = qJD(2) * t180 + V_base(4);
t140 = qJD(3) * t180 + t162;
t125 = qJD(4) * t180 + t140;
t167 = V_base(6) + qJD(1);
t246 = (-t251 * t164 + t250 * t165) * t167 + (-t255 * t164 + t253 * t165) * t125 + (-t256 * t164 + t254 * t165) * t124;
t245 = (t261 * t164 + t260 * t165) * t167 + (t252 * t180 + t249 * t182) * t125 + (t249 * t180 - t252 * t182) * t124;
t179 = sin(qJ(2));
t241 = pkin(2) * t179;
t170 = sin(t178);
t240 = pkin(3) * t170;
t181 = cos(qJ(2));
t239 = t181 * pkin(2);
t237 = t247 * t180 - t248 * t182;
t236 = t248 * t180 + t247 * t182;
t235 = Icges(2,4) * t180;
t234 = Icges(3,4) * t179;
t233 = Icges(3,4) * t181;
t232 = Icges(4,4) * t170;
t171 = cos(t178);
t231 = Icges(4,4) * t171;
t104 = -pkin(7) * t182 + t180 * t239;
t159 = t180 * pkin(1) - t182 * pkin(6);
t226 = -t104 - t159;
t224 = pkin(3) * t171;
t220 = V_base(5) * pkin(5) + V_base(1);
t81 = -pkin(8) * t182 + t180 * t224;
t217 = -t81 + t226;
t216 = rSges(6,2) * t165 + t257 * t164;
t161 = -qJD(2) * t182 + V_base(5);
t215 = t161 * t241 + t220;
t139 = t182 * t221 + V_base(5);
t214 = t139 * t240 + t215;
t213 = rSges(3,1) * t181 - rSges(3,2) * t179;
t212 = rSges(4,1) * t171 - rSges(4,2) * t170;
t211 = rSges(5,1) * t165 - rSges(5,2) * t164;
t209 = Icges(3,1) * t181 - t234;
t208 = Icges(4,1) * t171 - t232;
t205 = -Icges(3,2) * t179 + t233;
t204 = -Icges(4,2) * t170 + t231;
t201 = Icges(3,5) * t181 - Icges(3,6) * t179;
t200 = Icges(4,5) * t171 - Icges(4,6) * t170;
t160 = t182 * pkin(1) + t180 * pkin(6);
t197 = -V_base(4) * pkin(5) + t167 * t160 + V_base(2);
t196 = V_base(4) * t159 - t160 * V_base(5) + V_base(3);
t193 = (-Icges(4,3) * t182 + t180 * t200) * t139 + (Icges(4,3) * t180 + t182 * t200) * t140 + (Icges(4,5) * t170 + Icges(4,6) * t171) * t167;
t192 = (-Icges(3,3) * t182 + t180 * t201) * t161 + (Icges(3,3) * t180 + t182 * t201) * t162 + (Icges(3,5) * t179 + Icges(3,6) * t181) * t167;
t105 = pkin(7) * t180 + t182 * t239;
t191 = t162 * t104 - t105 * t161 + t196;
t190 = t167 * t105 - t162 * t241 + t197;
t82 = pkin(8) * t180 + t182 * t224;
t189 = -t139 * t82 + t140 * t81 + t191;
t188 = -t140 * t240 + t167 * t82 + t190;
t108 = -Icges(4,6) * t182 + t180 * t204;
t109 = Icges(4,6) * t180 + t182 * t204;
t110 = -Icges(4,5) * t182 + t180 * t208;
t111 = Icges(4,5) * t180 + t182 * t208;
t136 = Icges(4,2) * t171 + t232;
t137 = Icges(4,1) * t170 + t231;
t185 = (-t109 * t170 + t111 * t171) * t140 + (-t108 * t170 + t110 * t171) * t139 + (-t136 * t170 + t137 * t171) * t167;
t116 = -Icges(3,6) * t182 + t180 * t205;
t117 = Icges(3,6) * t180 + t182 * t205;
t118 = -Icges(3,5) * t182 + t180 * t209;
t119 = Icges(3,5) * t180 + t182 * t209;
t150 = Icges(3,2) * t181 + t234;
t153 = Icges(3,1) * t179 + t233;
t184 = (-t117 * t179 + t119 * t181) * t162 + (-t116 * t179 + t118 * t181) * t161 + (-t150 * t179 + t153 * t181) * t167;
t173 = Icges(2,4) * t182;
t158 = rSges(2,1) * t182 - rSges(2,2) * t180;
t157 = rSges(2,1) * t180 + rSges(2,2) * t182;
t156 = rSges(3,1) * t179 + rSges(3,2) * t181;
t155 = Icges(2,1) * t182 - t235;
t154 = Icges(2,1) * t180 + t173;
t152 = -Icges(2,2) * t180 + t173;
t151 = Icges(2,2) * t182 + t235;
t146 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t145 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t144 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t138 = rSges(4,1) * t170 + rSges(4,2) * t171;
t133 = rSges(5,1) * t164 + rSges(5,2) * t165;
t122 = rSges(3,3) * t180 + t182 * t213;
t121 = -rSges(3,3) * t182 + t180 * t213;
t113 = rSges(4,3) * t180 + t182 * t212;
t112 = -rSges(4,3) * t182 + t180 * t212;
t103 = rSges(5,3) * t180 + t182 * t211;
t101 = -rSges(5,3) * t182 + t180 * t211;
t87 = V_base(5) * rSges(2,3) - t157 * t167 + t220;
t86 = t158 * t167 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t85 = t157 * V_base(4) - t158 * V_base(5) + V_base(3);
t76 = t156 * t161 + (-t121 - t159) * t167 + t220;
t75 = t122 * t167 - t156 * t162 + t197;
t74 = t121 * t162 - t122 * t161 + t196;
t73 = t138 * t139 + (-t112 + t226) * t167 + t215;
t72 = t113 * t167 - t138 * t140 + t190;
t71 = t112 * t140 - t113 * t139 + t191;
t70 = t124 * t133 + (-t101 + t217) * t167 + t214;
t69 = t103 * t167 - t125 * t133 + t188;
t68 = qJD(5) * t180 + t216 * t124 + (t217 - t237) * t167 + t214;
t67 = -qJD(5) * t182 - t125 * t216 + t167 * t236 + t188;
t66 = t101 * t125 - t103 * t124 + t189;
t65 = -t124 * t236 + t125 * t237 + t189;
t1 = m(1) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(2) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t162 * (t192 * t180 + t184 * t182) / 0.2e1 + t161 * (t184 * t180 - t192 * t182) / 0.2e1 + m(4) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + t140 * (t193 * t180 + t185 * t182) / 0.2e1 + t139 * (t185 * t180 - t193 * t182) / 0.2e1 + m(5) * (t66 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + (t246 * t180 - t245 * t182) * t124 / 0.2e1 + (t245 * t180 + t246 * t182) * t125 / 0.2e1 + ((-t151 * t180 + t154 * t182 + Icges(1,4)) * V_base(5) + (-t152 * t180 + t155 * t182 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t151 * t182 + t154 * t180 + Icges(1,2)) * V_base(5) + (t152 * t182 + t155 * t180 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t117 * t181 + t119 * t179) * t162 + (t116 * t181 + t118 * t179) * t161 + (t109 * t171 + t111 * t170) * t140 + (t108 * t171 + t110 * t170) * t139 + (t253 * t164 + t255 * t165) * t125 + (t254 * t164 + t256 * t165) * t124 + (t171 * t136 + t170 * t137 + t181 * t150 + t179 * t153 + t250 * t164 + t251 * t165 + Icges(2,3)) * t167) * t167 / 0.2e1 + V_base(4) * t167 * (Icges(2,5) * t182 - Icges(2,6) * t180) + V_base(5) * t167 * (Icges(2,5) * t180 + Icges(2,6) * t182) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
