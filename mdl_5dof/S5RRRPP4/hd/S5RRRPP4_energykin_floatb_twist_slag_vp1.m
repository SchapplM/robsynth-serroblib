% Calculate kinetic energy for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:48
% EndTime: 2019-12-31 20:54:50
% DurationCPUTime: 2.22s
% Computational Cost: add. (1212->228), mult. (1189->309), div. (0->0), fcn. (971->8), ass. (0->124)
t266 = Icges(5,4) - Icges(6,5);
t265 = Icges(5,1) + Icges(6,1);
t264 = Icges(5,2) + Icges(6,3);
t179 = qJ(2) + qJ(3);
t169 = pkin(8) + t179;
t166 = sin(t169);
t263 = t266 * t166;
t167 = cos(t169);
t262 = t266 * t167;
t261 = Icges(6,4) + Icges(5,5);
t260 = Icges(5,6) - Icges(6,6);
t259 = t264 * t166 - t262;
t258 = t265 * t167 - t263;
t257 = rSges(6,1) + pkin(4);
t256 = rSges(6,3) + qJ(5);
t181 = sin(qJ(1));
t183 = cos(qJ(1));
t255 = t259 * t181 + t260 * t183;
t254 = -t260 * t181 + t259 * t183;
t253 = t258 * t181 - t261 * t183;
t252 = t261 * t181 + t258 * t183;
t251 = -t264 * t167 - t263;
t250 = t265 * t166 + t262;
t249 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t173 = sin(t179);
t174 = cos(t179);
t248 = Icges(4,5) * t174 - Icges(4,6) * t173 - t260 * t166 + t261 * t167;
t247 = t256 * t166 + t257 * t167;
t233 = Icges(4,4) * t174;
t206 = -Icges(4,2) * t173 + t233;
t111 = -Icges(4,6) * t183 + t206 * t181;
t112 = Icges(4,6) * t181 + t206 * t183;
t234 = Icges(4,4) * t173;
t210 = Icges(4,1) * t174 - t234;
t113 = -Icges(4,5) * t183 + t210 * t181;
t114 = Icges(4,5) * t181 + t210 * t183;
t139 = Icges(4,2) * t174 + t234;
t140 = Icges(4,1) * t173 + t233;
t142 = V_base(5) + (-qJD(2) - qJD(3)) * t183;
t165 = qJD(2) * t181 + V_base(4);
t143 = qJD(3) * t181 + t165;
t170 = V_base(6) + qJD(1);
t246 = (-t139 * t173 + t140 * t174 + t251 * t166 + t250 * t167) * t170 + (-t112 * t173 + t114 * t174 + t254 * t166 + t252 * t167) * t143 + (-t111 * t173 + t113 * t174 + t255 * t166 + t253 * t167) * t142;
t245 = (Icges(4,5) * t173 + Icges(4,6) * t174 + t261 * t166 + t260 * t167) * t170 + (t249 * t181 + t248 * t183) * t143 + (t248 * t181 - t249 * t183) * t142;
t180 = sin(qJ(2));
t241 = pkin(2) * t180;
t240 = pkin(3) * t173;
t182 = cos(qJ(2));
t239 = t182 * pkin(2);
t237 = Icges(2,4) * t181;
t236 = Icges(3,4) * t180;
t235 = Icges(3,4) * t182;
t228 = -rSges(6,2) * t183 + t247 * t181;
t227 = rSges(6,2) * t181 + t247 * t183;
t107 = -pkin(7) * t183 + t239 * t181;
t162 = t181 * pkin(1) - t183 * pkin(6);
t226 = -t107 - t162;
t225 = t257 * t166 - t256 * t167;
t224 = pkin(3) * t174;
t222 = qJD(5) * t166;
t221 = V_base(5) * pkin(5) + V_base(1);
t84 = -qJ(4) * t183 + t224 * t181;
t218 = -t84 + t226;
t164 = -qJD(2) * t183 + V_base(5);
t217 = t164 * t241 + t221;
t216 = rSges(3,1) * t182 - rSges(3,2) * t180;
t215 = rSges(4,1) * t174 - rSges(4,2) * t173;
t214 = rSges(5,1) * t167 - rSges(5,2) * t166;
t211 = Icges(3,1) * t182 - t236;
t207 = -Icges(3,2) * t180 + t235;
t203 = Icges(3,5) * t182 - Icges(3,6) * t180;
t199 = qJD(4) * t181 + t142 * t240 + t217;
t163 = t183 * pkin(1) + t181 * pkin(6);
t198 = -V_base(4) * pkin(5) + t170 * t163 + V_base(2);
t197 = V_base(4) * t162 - t163 * V_base(5) + V_base(3);
t193 = (-Icges(3,3) * t183 + t203 * t181) * t164 + (Icges(3,3) * t181 + t203 * t183) * t165 + (Icges(3,5) * t180 + Icges(3,6) * t182) * t170;
t108 = pkin(7) * t181 + t239 * t183;
t192 = t165 * t107 - t108 * t164 + t197;
t191 = t170 * t108 - t165 * t241 + t198;
t190 = t143 * t84 + t192;
t85 = qJ(4) * t181 + t224 * t183;
t189 = -qJD(4) * t183 + t170 * t85 + t191;
t121 = -Icges(3,6) * t183 + t207 * t181;
t122 = Icges(3,6) * t181 + t207 * t183;
t123 = -Icges(3,5) * t183 + t211 * t181;
t124 = Icges(3,5) * t181 + t211 * t183;
t153 = Icges(3,2) * t182 + t236;
t156 = Icges(3,1) * t180 + t235;
t185 = (-t122 * t180 + t124 * t182) * t165 + (-t121 * t180 + t123 * t182) * t164 + (-t153 * t180 + t156 * t182) * t170;
t175 = Icges(2,4) * t183;
t161 = rSges(2,1) * t183 - rSges(2,2) * t181;
t160 = rSges(2,1) * t181 + rSges(2,2) * t183;
t159 = rSges(3,1) * t180 + rSges(3,2) * t182;
t158 = Icges(2,1) * t183 - t237;
t157 = Icges(2,1) * t181 + t175;
t155 = -Icges(2,2) * t181 + t175;
t154 = Icges(2,2) * t183 + t237;
t149 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t148 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t147 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t141 = rSges(4,1) * t173 + rSges(4,2) * t174;
t136 = rSges(5,1) * t166 + rSges(5,2) * t167;
t126 = rSges(3,3) * t181 + t216 * t183;
t125 = -rSges(3,3) * t183 + t216 * t181;
t116 = rSges(4,3) * t181 + t215 * t183;
t115 = -rSges(4,3) * t183 + t215 * t181;
t106 = rSges(5,3) * t181 + t214 * t183;
t104 = -rSges(5,3) * t183 + t214 * t181;
t102 = V_base(5) * rSges(2,3) - t160 * t170 + t221;
t101 = t161 * t170 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t88 = t160 * V_base(4) - t161 * V_base(5) + V_base(3);
t81 = t159 * t164 + (-t125 - t162) * t170 + t221;
t80 = t126 * t170 - t159 * t165 + t198;
t79 = t125 * t165 - t126 * t164 + t197;
t78 = t141 * t142 + (-t115 + t226) * t170 + t217;
t77 = t116 * t170 - t141 * t143 + t191;
t76 = t115 * t143 - t116 * t142 + t192;
t75 = t136 * t142 + (-t104 + t218) * t170 + t199;
t74 = t106 * t170 + (-t136 - t240) * t143 + t189;
t73 = t183 * t222 + t225 * t142 + (t218 - t228) * t170 + t199;
t72 = t181 * t222 + t227 * t170 + (-t225 - t240) * t143 + t189;
t71 = t104 * t143 + (-t106 - t85) * t142 + t190;
t70 = -qJD(5) * t167 + t228 * t143 + (-t85 - t227) * t142 + t190;
t1 = m(1) * (t147 ^ 2 + t148 ^ 2 + t149 ^ 2) / 0.2e1 + m(2) * (t101 ^ 2 + t102 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + t165 * (t193 * t181 + t185 * t183) / 0.2e1 + t164 * (t185 * t181 - t193 * t183) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + ((-t154 * t181 + t157 * t183 + Icges(1,4)) * V_base(5) + (-t155 * t181 + t158 * t183 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t154 * t183 + t157 * t181 + Icges(1,2)) * V_base(5) + (t155 * t183 + t158 * t181 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t246 * t181 - t245 * t183) * t142 / 0.2e1 + (t245 * t181 + t246 * t183) * t143 / 0.2e1 + ((t122 * t182 + t124 * t180) * t165 + (t121 * t182 + t123 * t180) * t164 + (t112 * t174 + t114 * t173 + t252 * t166 - t254 * t167) * t143 + (t111 * t174 + t113 * t173 + t253 * t166 - t255 * t167) * t142 + (t174 * t139 + t173 * t140 + t182 * t153 + t180 * t156 + t250 * t166 - t251 * t167 + Icges(2,3)) * t170) * t170 / 0.2e1 + V_base(4) * t170 * (Icges(2,5) * t183 - Icges(2,6) * t181) + V_base(5) * t170 * (Icges(2,5) * t181 + Icges(2,6) * t183) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
