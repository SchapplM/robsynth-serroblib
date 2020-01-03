% Calculate kinetic energy for
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:53:57
% EndTime: 2019-12-31 19:53:59
% DurationCPUTime: 2.16s
% Computational Cost: add. (1194->229), mult. (1171->309), div. (0->0), fcn. (953->8), ass. (0->125)
t272 = Icges(5,4) - Icges(6,5);
t271 = Icges(5,1) + Icges(6,1);
t270 = Icges(5,2) + Icges(6,3);
t179 = qJ(2) + pkin(8);
t171 = qJ(4) + t179;
t166 = sin(t171);
t269 = t272 * t166;
t167 = cos(t171);
t268 = t272 * t167;
t267 = Icges(6,4) + Icges(5,5);
t266 = Icges(5,6) - Icges(6,6);
t265 = t166 * t270 - t268;
t264 = t167 * t271 - t269;
t263 = rSges(6,1) + pkin(4);
t262 = rSges(6,3) + qJ(5);
t182 = sin(qJ(1));
t184 = cos(qJ(1));
t261 = t182 * t265 + t184 * t266;
t260 = -t182 * t266 + t184 * t265;
t259 = t182 * t264 - t184 * t267;
t258 = t182 * t267 + t184 * t264;
t257 = Icges(6,2) + Icges(5,3);
t256 = Icges(3,3) + Icges(4,3);
t255 = -t167 * t270 - t269;
t254 = t166 * t271 + t268;
t253 = -t166 * t266 + t167 * t267;
t169 = sin(t179);
t170 = cos(t179);
t181 = sin(qJ(2));
t183 = cos(qJ(2));
t252 = Icges(3,5) * t183 + Icges(4,5) * t170 - Icges(3,6) * t181 - Icges(4,6) * t169;
t251 = t166 * t262 + t167 * t263;
t233 = Icges(4,4) * t170;
t206 = -Icges(4,2) * t169 + t233;
t111 = -Icges(4,6) * t184 + t182 * t206;
t112 = Icges(4,6) * t182 + t184 * t206;
t234 = Icges(4,4) * t169;
t210 = Icges(4,1) * t170 - t234;
t113 = -Icges(4,5) * t184 + t182 * t210;
t114 = Icges(4,5) * t182 + t184 * t210;
t235 = Icges(3,4) * t183;
t207 = -Icges(3,2) * t181 + t235;
t121 = -Icges(3,6) * t184 + t182 * t207;
t122 = Icges(3,6) * t182 + t184 * t207;
t236 = Icges(3,4) * t181;
t211 = Icges(3,1) * t183 - t236;
t123 = -Icges(3,5) * t184 + t182 * t211;
t124 = Icges(3,5) * t182 + t184 * t211;
t139 = Icges(4,2) * t170 + t234;
t140 = Icges(4,1) * t169 + t233;
t153 = Icges(3,2) * t183 + t236;
t156 = Icges(3,1) * t181 + t235;
t164 = -qJD(2) * t184 + V_base(5);
t165 = qJD(2) * t182 + V_base(4);
t172 = V_base(6) + qJD(1);
t250 = (-t139 * t169 + t140 * t170 - t153 * t181 + t156 * t183) * t172 + (-t112 * t169 + t114 * t170 - t122 * t181 + t124 * t183) * t165 + (-t111 * t169 + t113 * t170 - t121 * t181 + t123 * t183) * t164;
t142 = V_base(5) + (-qJD(2) - qJD(4)) * t184;
t143 = qJD(4) * t182 + t165;
t249 = (t166 * t255 + t167 * t254) * t172 + (t166 * t260 + t167 * t258) * t143 + (t166 * t261 + t167 * t259) * t142;
t248 = (Icges(3,5) * t181 + Icges(4,5) * t169 + Icges(3,6) * t183 + Icges(4,6) * t170) * t172 + (t182 * t256 + t184 * t252) * t165 + (t182 * t252 - t184 * t256) * t164;
t247 = (t166 * t267 + t167 * t266) * t172 + (t182 * t257 + t184 * t253) * t143 + (t182 * t253 - t184 * t257) * t142;
t241 = pkin(2) * t181;
t240 = pkin(3) * t169;
t239 = t183 * pkin(2);
t237 = Icges(2,4) * t182;
t228 = -rSges(6,2) * t184 + t182 * t251;
t227 = t182 * rSges(6,2) + t184 * t251;
t107 = -qJ(3) * t184 + t182 * t239;
t162 = pkin(1) * t182 - pkin(6) * t184;
t226 = -t107 - t162;
t225 = t166 * t263 - t167 * t262;
t224 = pkin(3) * t170;
t222 = qJD(5) * t166;
t221 = V_base(5) * pkin(5) + V_base(1);
t84 = -pkin(7) * t184 + t182 * t224;
t218 = -t84 + t226;
t217 = qJD(3) * t182 + t164 * t241 + t221;
t216 = rSges(3,1) * t183 - rSges(3,2) * t181;
t215 = rSges(4,1) * t170 - rSges(4,2) * t169;
t214 = rSges(5,1) * t167 - rSges(5,2) * t166;
t199 = t164 * t240 + t217;
t163 = pkin(1) * t184 + t182 * pkin(6);
t198 = -V_base(4) * pkin(5) + t163 * t172 + V_base(2);
t197 = t162 * V_base(4) - t163 * V_base(5) + V_base(3);
t196 = t107 * t165 + t197;
t108 = qJ(3) * t182 + t184 * t239;
t191 = -qJD(3) * t184 + t108 * t172 + t198;
t85 = pkin(7) * t182 + t184 * t224;
t190 = t165 * t84 + (-t108 - t85) * t164 + t196;
t189 = t172 * t85 + (-t240 - t241) * t165 + t191;
t175 = Icges(2,4) * t184;
t161 = rSges(2,1) * t184 - rSges(2,2) * t182;
t160 = rSges(2,1) * t182 + rSges(2,2) * t184;
t159 = rSges(3,1) * t181 + rSges(3,2) * t183;
t158 = Icges(2,1) * t184 - t237;
t157 = Icges(2,1) * t182 + t175;
t155 = -Icges(2,2) * t182 + t175;
t154 = Icges(2,2) * t184 + t237;
t149 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t148 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t147 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t141 = rSges(4,1) * t169 + rSges(4,2) * t170;
t135 = rSges(5,1) * t166 + rSges(5,2) * t167;
t126 = rSges(3,3) * t182 + t184 * t216;
t125 = -rSges(3,3) * t184 + t182 * t216;
t116 = rSges(4,3) * t182 + t184 * t215;
t115 = -rSges(4,3) * t184 + t182 * t215;
t106 = t182 * rSges(5,3) + t184 * t214;
t104 = -rSges(5,3) * t184 + t182 * t214;
t102 = V_base(5) * rSges(2,3) - t160 * t172 + t221;
t101 = t161 * t172 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t88 = t160 * V_base(4) - t161 * V_base(5) + V_base(3);
t81 = t159 * t164 + (-t125 - t162) * t172 + t221;
t80 = t126 * t172 - t159 * t165 + t198;
t79 = t125 * t165 - t126 * t164 + t197;
t78 = t141 * t164 + (-t115 + t226) * t172 + t217;
t77 = t172 * t116 + (-t141 - t241) * t165 + t191;
t76 = t115 * t165 + (-t108 - t116) * t164 + t196;
t75 = t135 * t142 + (-t104 + t218) * t172 + t199;
t74 = t106 * t172 - t135 * t143 + t189;
t73 = t184 * t222 + t225 * t142 + (t218 - t228) * t172 + t199;
t72 = -t143 * t225 + t172 * t227 + t182 * t222 + t189;
t71 = t104 * t143 - t106 * t142 + t190;
t70 = -qJD(5) * t167 - t142 * t227 + t143 * t228 + t190;
t1 = m(1) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(2) * (t101 ^ 2 + t102 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + (t249 * t182 - t247 * t184) * t142 / 0.2e1 + (t247 * t182 + t249 * t184) * t143 / 0.2e1 + (t250 * t182 - t248 * t184) * t164 / 0.2e1 + (t248 * t182 + t250 * t184) * t165 / 0.2e1 + ((-t182 * t154 + t157 * t184 + Icges(1,4)) * V_base(5) + (-t182 * t155 + t184 * t158 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t184 * t154 + t182 * t157 + Icges(1,2)) * V_base(5) + (t155 * t184 + t182 * t158 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t112 * t170 + t114 * t169 + t122 * t183 + t124 * t181) * t165 + (t111 * t170 + t113 * t169 + t121 * t183 + t123 * t181) * t164 + (t166 * t258 - t167 * t260) * t143 + (t166 * t259 - t167 * t261) * t142 + (t170 * t139 + t169 * t140 + t183 * t153 + t181 * t156 + t254 * t166 - t255 * t167 + Icges(2,3)) * t172) * t172 / 0.2e1 + t172 * V_base(4) * (Icges(2,5) * t184 - Icges(2,6) * t182) + V_base(5) * t172 * (Icges(2,5) * t182 + Icges(2,6) * t184) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
