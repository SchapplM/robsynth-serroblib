% Calculate kinetic energy for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP5_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:31
% EndTime: 2019-12-31 20:57:33
% DurationCPUTime: 2.04s
% Computational Cost: add. (989->203), mult. (1166->279), div. (0->0), fcn. (948->6), ass. (0->111)
t254 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t253 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t252 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t169 = qJ(2) + qJ(3);
t165 = sin(t169);
t251 = t254 * t165;
t166 = cos(t169);
t250 = t254 * t166;
t249 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t248 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t247 = t252 * t165 - t250;
t246 = t253 * t166 - t251;
t245 = rSges(6,1) + pkin(4);
t171 = sin(qJ(1));
t173 = cos(qJ(1));
t244 = t247 * t171 + t248 * t173;
t243 = -t248 * t171 + t247 * t173;
t242 = t246 * t171 - t249 * t173;
t241 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t240 = t249 * t171 + t246 * t173;
t239 = -t252 * t166 - t251;
t238 = t253 * t165 + t250;
t237 = -t248 * t165 + t249 * t166;
t236 = rSges(6,3) + qJ(5);
t235 = rSges(6,2) * t165 + t245 * t166;
t137 = V_base(5) + (-qJD(2) - qJD(3)) * t173;
t161 = qJD(2) * t171 + V_base(4);
t138 = qJD(3) * t171 + t161;
t163 = V_base(6) + qJD(1);
t234 = (t239 * t165 + t238 * t166) * t163 + (t243 * t165 + t240 * t166) * t138 + (t244 * t165 + t242 * t166) * t137;
t233 = (t249 * t165 + t248 * t166) * t163 + (t241 * t171 + t237 * t173) * t138 + (t237 * t171 - t241 * t173) * t137;
t170 = sin(qJ(2));
t229 = pkin(2) * t170;
t172 = cos(qJ(2));
t227 = pkin(2) * t172;
t158 = t171 * pkin(1) - t173 * pkin(6);
t84 = -pkin(7) * t173 + t227 * t171;
t225 = -t158 - t84;
t224 = Icges(2,4) * t171;
t223 = Icges(3,4) * t170;
t222 = Icges(3,4) * t172;
t215 = t235 * t171 + t236 * t173;
t214 = -t236 * t171 + t235 * t173;
t213 = qJD(4) * t165;
t212 = V_base(5) * pkin(5) + V_base(1);
t202 = pkin(3) * t166 + qJ(4) * t165;
t119 = t202 * t171;
t209 = -t119 + t225;
t208 = -rSges(6,2) * t166 + t245 * t165;
t160 = -qJD(2) * t173 + V_base(5);
t207 = t160 * t229 + t212;
t206 = rSges(3,1) * t172 - rSges(3,2) * t170;
t205 = rSges(4,1) * t166 - rSges(4,2) * t165;
t204 = rSges(5,1) * t166 + rSges(5,3) * t165;
t201 = Icges(3,1) * t172 - t223;
t197 = -Icges(3,2) * t170 + t222;
t193 = Icges(3,5) * t172 - Icges(3,6) * t170;
t133 = pkin(3) * t165 - qJ(4) * t166;
t189 = t137 * t133 + t173 * t213 + t207;
t159 = t173 * pkin(1) + t171 * pkin(6);
t188 = -V_base(4) * pkin(5) + t163 * t159 + V_base(2);
t187 = V_base(4) * t158 - t159 * V_base(5) + V_base(3);
t183 = (-Icges(3,3) * t173 + t193 * t171) * t160 + (Icges(3,3) * t171 + t193 * t173) * t161 + (Icges(3,5) * t170 + Icges(3,6) * t172) * t163;
t85 = pkin(7) * t171 + t227 * t173;
t182 = -t160 * t85 + t161 * t84 + t187;
t181 = -t161 * t229 + t163 * t85 + t188;
t120 = t202 * t173;
t180 = t163 * t120 + t171 * t213 + t181;
t179 = -qJD(4) * t166 + t138 * t119 + t182;
t113 = -Icges(3,6) * t173 + t197 * t171;
t114 = Icges(3,6) * t171 + t197 * t173;
t115 = -Icges(3,5) * t173 + t201 * t171;
t116 = Icges(3,5) * t171 + t201 * t173;
t147 = Icges(3,2) * t172 + t223;
t150 = Icges(3,1) * t170 + t222;
t175 = (-t114 * t170 + t116 * t172) * t161 + (-t113 * t170 + t115 * t172) * t160 + (-t147 * t170 + t150 * t172) * t163;
t167 = Icges(2,4) * t173;
t155 = rSges(2,1) * t173 - rSges(2,2) * t171;
t154 = rSges(2,1) * t171 + rSges(2,2) * t173;
t153 = rSges(3,1) * t170 + rSges(3,2) * t172;
t152 = Icges(2,1) * t173 - t224;
t151 = Icges(2,1) * t171 + t167;
t149 = -Icges(2,2) * t171 + t167;
t148 = Icges(2,2) * t173 + t224;
t143 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t142 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t141 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t136 = rSges(4,1) * t165 + rSges(4,2) * t166;
t135 = rSges(5,1) * t165 - rSges(5,3) * t166;
t118 = rSges(3,3) * t171 + t206 * t173;
t117 = -rSges(3,3) * t173 + t206 * t171;
t110 = rSges(4,3) * t171 + t205 * t173;
t109 = rSges(5,2) * t171 + t204 * t173;
t107 = -rSges(4,3) * t173 + t205 * t171;
t106 = -rSges(5,2) * t173 + t204 * t171;
t83 = V_base(5) * rSges(2,3) - t154 * t163 + t212;
t82 = t155 * t163 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t80 = t154 * V_base(4) - t155 * V_base(5) + V_base(3);
t76 = t153 * t160 + (-t117 - t158) * t163 + t212;
t75 = t118 * t163 - t153 * t161 + t188;
t74 = t117 * t161 - t118 * t160 + t187;
t73 = t136 * t137 + (-t107 + t225) * t163 + t207;
t72 = t110 * t163 - t136 * t138 + t181;
t71 = t107 * t138 - t110 * t137 + t182;
t70 = t135 * t137 + (-t106 + t209) * t163 + t189;
t69 = t109 * t163 + (-t133 - t135) * t138 + t180;
t68 = -qJD(5) * t171 + t208 * t137 + (t209 - t215) * t163 + t189;
t67 = qJD(5) * t173 + t214 * t163 + (-t133 - t208) * t138 + t180;
t66 = t106 * t138 + (-t109 - t120) * t137 + t179;
t65 = t215 * t138 + (-t120 - t214) * t137 + t179;
t1 = m(1) * (t141 ^ 2 + t142 ^ 2 + t143 ^ 2) / 0.2e1 + m(2) * (t80 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(3) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t161 * (t183 * t171 + t175 * t173) / 0.2e1 + t160 * (t175 * t171 - t183 * t173) / 0.2e1 + m(4) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + ((-t148 * t171 + t151 * t173 + Icges(1,4)) * V_base(5) + (-t149 * t171 + t152 * t173 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t173 * t148 + t151 * t171 + Icges(1,2)) * V_base(5) + (t149 * t173 + t152 * t171 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t234 * t171 - t233 * t173) * t137 / 0.2e1 + (t233 * t171 + t234 * t173) * t138 / 0.2e1 + ((t114 * t172 + t116 * t170) * t161 + (t113 * t172 + t115 * t170) * t160 + (t240 * t165 - t243 * t166) * t138 + (t242 * t165 - t244 * t166) * t137 + (t147 * t172 + t150 * t170 + t238 * t165 - t239 * t166 + Icges(2,3)) * t163) * t163 / 0.2e1 + t163 * V_base(4) * (Icges(2,5) * t173 - Icges(2,6) * t171) + V_base(5) * t163 * (Icges(2,5) * t171 + Icges(2,6) * t173) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
