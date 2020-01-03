% Calculate kinetic energy for
% S5RRRPP8
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:26
% EndTime: 2019-12-31 21:07:28
% DurationCPUTime: 2.03s
% Computational Cost: add. (955->230), mult. (1888->327), div. (0->0), fcn. (1864->6), ass. (0->111)
t251 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t250 = -Icges(4,4) - Icges(5,6) + Icges(6,6);
t249 = -Icges(4,5) - Icges(6,5) + Icges(5,4);
t248 = Icges(4,2) + Icges(6,2) + Icges(5,3);
t247 = Icges(4,6) - Icges(5,5) - Icges(6,4);
t246 = -Icges(4,3) - Icges(6,1) - Icges(5,1);
t245 = rSges(6,1) + pkin(4);
t244 = rSges(6,3) + qJ(5);
t194 = sin(qJ(3));
t197 = cos(qJ(3));
t199 = cos(qJ(1));
t223 = t197 * t199;
t196 = sin(qJ(1));
t198 = cos(qJ(2));
t224 = t196 * t198;
t156 = t194 * t224 + t223;
t222 = t199 * t194;
t157 = t197 * t224 - t222;
t195 = sin(qJ(2));
t227 = t195 * t196;
t243 = t248 * t156 + t250 * t157 - t247 * t227;
t242 = t250 * t156 + t251 * t157 - t249 * t227;
t158 = -t196 * t197 + t198 * t222;
t159 = t196 * t194 + t198 * t223;
t225 = t195 * t199;
t241 = t250 * t158 + t251 * t159 - t249 * t225;
t240 = t248 * t158 + t250 * t159 - t247 * t225;
t239 = -t247 * t156 - t249 * t157 - t246 * t227;
t238 = -t247 * t158 - t249 * t159 - t246 * t225;
t237 = t246 * t198 + (-t247 * t194 - t249 * t197) * t195;
t236 = t247 * t198 + (t248 * t194 + t250 * t197) * t195;
t235 = t249 * t198 + (t250 * t194 + t251 * t197) * t195;
t230 = Icges(2,4) * t196;
t229 = Icges(3,4) * t195;
t228 = Icges(3,4) * t198;
t226 = t195 * t197;
t221 = rSges(6,2) * t156 + t244 * t157 + t245 * t227;
t220 = t158 * rSges(6,2) + t244 * t159 + t245 * t225;
t219 = (rSges(6,2) * t194 + rSges(6,3) * t197) * t195 + qJ(5) * t226 - t245 * t198;
t218 = qJD(3) * t195;
t217 = V_base(5) * pkin(5) + V_base(1);
t186 = qJD(2) * t196 + V_base(4);
t190 = V_base(6) + qJD(1);
t214 = pkin(2) * t198 + pkin(7) * t195;
t185 = -qJD(2) * t199 + V_base(5);
t213 = rSges(3,1) * t198 - rSges(3,2) * t195;
t212 = Icges(3,1) * t198 - t229;
t211 = -Icges(3,2) * t195 + t228;
t210 = Icges(3,5) * t198 - Icges(3,6) * t195;
t184 = pkin(1) * t199 + t196 * pkin(6);
t209 = -V_base(4) * pkin(5) + t190 * t184 + V_base(2);
t183 = t196 * pkin(1) - pkin(6) * t199;
t208 = V_base(4) * t183 - t184 * V_base(5) + V_base(3);
t162 = t214 * t196;
t182 = pkin(2) * t195 - pkin(7) * t198;
t207 = t185 * t182 + (-t162 - t183) * t190 + t217;
t206 = (-Icges(3,3) * t199 + t196 * t210) * t185 + (Icges(3,3) * t196 + t199 * t210) * t186 + (Icges(3,5) * t195 + Icges(3,6) * t198) * t190;
t163 = t214 * t199;
t205 = t190 * t163 - t182 * t186 + t209;
t154 = t196 * t218 + t185;
t160 = (pkin(3) * t197 + qJ(4) * t194) * t195;
t204 = qJD(4) * t158 + t154 * t160 + t207;
t203 = t186 * t162 - t163 * t185 + t208;
t122 = pkin(3) * t159 + qJ(4) * t158;
t178 = -qJD(3) * t198 + t190;
t202 = qJD(4) * t156 + t178 * t122 + t205;
t121 = pkin(3) * t157 + qJ(4) * t156;
t155 = t199 * t218 + t186;
t201 = qJD(4) * t195 * t194 + t155 * t121 + t203;
t135 = -Icges(3,6) * t199 + t196 * t211;
t136 = Icges(3,6) * t196 + t199 * t211;
t138 = -Icges(3,5) * t199 + t196 * t212;
t139 = Icges(3,5) * t196 + t199 * t212;
t172 = Icges(3,2) * t198 + t229;
t175 = Icges(3,1) * t195 + t228;
t200 = (-t136 * t195 + t139 * t198) * t186 + (-t135 * t195 + t138 * t198) * t185 + (-t172 * t195 + t175 * t198) * t190;
t192 = Icges(2,4) * t199;
t181 = rSges(2,1) * t199 - t196 * rSges(2,2);
t180 = t196 * rSges(2,1) + rSges(2,2) * t199;
t179 = rSges(3,1) * t195 + rSges(3,2) * t198;
t177 = Icges(2,1) * t199 - t230;
t176 = Icges(2,1) * t196 + t192;
t174 = -Icges(2,2) * t196 + t192;
t173 = Icges(2,2) * t199 + t230;
t168 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t167 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t166 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t150 = -rSges(5,1) * t198 + (-rSges(5,2) * t197 + rSges(5,3) * t194) * t195;
t148 = t196 * rSges(3,3) + t199 * t213;
t147 = -rSges(3,3) * t199 + t196 * t213;
t146 = -rSges(4,3) * t198 + (rSges(4,1) * t197 - rSges(4,2) * t194) * t195;
t126 = V_base(5) * rSges(2,3) - t180 * t190 + t217;
t125 = t181 * t190 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t124 = t180 * V_base(4) - t181 * V_base(5) + V_base(3);
t120 = t159 * rSges(4,1) - t158 * rSges(4,2) + rSges(4,3) * t225;
t119 = rSges(4,1) * t157 - rSges(4,2) * t156 + rSges(4,3) * t227;
t118 = rSges(5,1) * t225 - t159 * rSges(5,2) + t158 * rSges(5,3);
t116 = rSges(5,1) * t227 - rSges(5,2) * t157 + rSges(5,3) * t156;
t94 = t179 * t185 + (-t147 - t183) * t190 + t217;
t93 = t148 * t190 - t179 * t186 + t209;
t92 = t147 * t186 - t148 * t185 + t208;
t91 = -t119 * t178 + t146 * t154 + t207;
t90 = t120 * t178 - t146 * t155 + t205;
t89 = t119 * t155 - t120 * t154 + t203;
t88 = t150 * t154 + (-t116 - t121) * t178 + t204;
t87 = t118 * t178 + (-t150 - t160) * t155 + t202;
t86 = t116 * t155 + (-t118 - t122) * t154 + t201;
t85 = qJD(5) * t159 + t219 * t154 + (-t121 - t221) * t178 + t204;
t84 = qJD(5) * t157 + t220 * t178 + (-t160 - t219) * t155 + t202;
t83 = qJD(5) * t226 + t221 * t155 + (-t122 - t220) * t154 + t201;
t1 = m(1) * (t166 ^ 2 + t167 ^ 2 + t168 ^ 2) / 0.2e1 + m(2) * (t124 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + t186 * (t206 * t196 + t200 * t199) / 0.2e1 + t185 * (t200 * t196 - t206 * t199) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + ((t136 * t198 + t139 * t195) * t186 + (t135 * t198 + t138 * t195) * t185 + (t172 * t198 + t175 * t195 + Icges(2,3)) * t190) * t190 / 0.2e1 + ((-t196 * t173 + t176 * t199 + Icges(1,4)) * V_base(5) + (-t196 * t174 + t177 * t199 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t173 * t199 + t196 * t176 + Icges(1,2)) * V_base(5) + (t174 * t199 + t196 * t177 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t236 * t156 + t235 * t157 + t237 * t227) * t178 + (t240 * t156 + t241 * t157 + t238 * t227) * t155 + (t243 * t156 + t242 * t157 + t239 * t227) * t154) * t154 / 0.2e1 + ((t236 * t158 + t235 * t159 + t237 * t225) * t178 + (t240 * t158 + t241 * t159 + t238 * t225) * t155 + (t243 * t158 + t242 * t159 + t239 * t225) * t154) * t155 / 0.2e1 + ((-t239 * t154 - t238 * t155 - t237 * t178) * t198 + ((t236 * t194 + t235 * t197) * t178 + (t240 * t194 + t241 * t197) * t155 + (t243 * t194 + t242 * t197) * t154) * t195) * t178 / 0.2e1 + V_base(4) * t190 * (Icges(2,5) * t199 - Icges(2,6) * t196) + V_base(5) * t190 * (Icges(2,5) * t196 + Icges(2,6) * t199) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
