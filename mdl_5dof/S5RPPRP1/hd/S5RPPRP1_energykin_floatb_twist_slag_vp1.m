% Calculate kinetic energy for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:32
% EndTime: 2019-12-05 17:35:34
% DurationCPUTime: 2.12s
% Computational Cost: add. (1213->248), mult. (1371->326), div. (0->0), fcn. (1263->8), ass. (0->121)
t253 = Icges(5,1) + Icges(6,1);
t252 = Icges(5,4) + Icges(6,4);
t251 = -Icges(6,5) - Icges(5,5);
t250 = Icges(5,2) + Icges(6,2);
t249 = -Icges(6,6) - Icges(5,6);
t248 = -Icges(6,3) - Icges(5,3);
t179 = qJ(1) + pkin(7);
t174 = sin(t179);
t175 = cos(t179);
t185 = cos(qJ(4));
t181 = cos(pkin(8));
t183 = sin(qJ(4));
t214 = t181 * t183;
t131 = t174 * t214 + t175 * t185;
t213 = t181 * t185;
t215 = t175 * t183;
t132 = -t174 * t213 + t215;
t180 = sin(pkin(8));
t218 = t174 * t180;
t247 = -t249 * t131 - t251 * t132 + t248 * t218;
t133 = t174 * t185 - t175 * t214;
t217 = t174 * t183;
t134 = t175 * t213 + t217;
t216 = t175 * t180;
t246 = -t249 * t133 - t251 * t134 - t248 * t216;
t245 = t250 * t131 + t252 * t132 + t249 * t218;
t244 = t250 * t133 + t252 * t134 - t249 * t216;
t243 = t252 * t131 + t253 * t132 + t251 * t218;
t242 = t252 * t133 + t253 * t134 - t251 * t216;
t241 = t248 * t181 + (t249 * t183 - t251 * t185) * t180;
t240 = t249 * t181 + (-t250 * t183 + t252 * t185) * t180;
t239 = t251 * t181 + (-t252 * t183 + t253 * t185) * t180;
t184 = sin(qJ(1));
t186 = cos(qJ(1));
t238 = -Icges(2,5) * t184 - Icges(3,5) * t174 - Icges(2,6) * t186 - Icges(3,6) * t175;
t237 = Icges(2,5) * t186 + Icges(3,5) * t175 - Icges(2,6) * t184 - Icges(3,6) * t174;
t219 = Icges(4,4) * t181;
t196 = -Icges(4,2) * t180 + t219;
t112 = Icges(4,6) * t175 - t174 * t196;
t220 = Icges(4,4) * t180;
t197 = Icges(4,1) * t181 - t220;
t114 = Icges(4,5) * t175 - t174 * t197;
t221 = Icges(3,4) * t175;
t236 = -Icges(3,1) * t174 - t112 * t180 + t114 * t181 - t221;
t113 = Icges(4,6) * t174 + t175 * t196;
t115 = Icges(4,5) * t174 + t175 * t197;
t222 = Icges(3,4) * t174;
t235 = Icges(3,1) * t175 - t113 * t180 + t115 * t181 - t222;
t227 = pkin(4) * t185;
t234 = qJ(5) * t180 + t181 * t227;
t229 = pkin(1) * t184;
t228 = pkin(1) * t186;
t226 = -pkin(5) - qJ(2);
t224 = Icges(2,4) * t184;
t223 = Icges(2,4) * t186;
t212 = rSges(6,1) * t132 + rSges(6,2) * t131 - rSges(6,3) * t218 + pkin(4) * t215 - t234 * t174;
t211 = rSges(6,1) * t134 + rSges(6,2) * t133 + rSges(6,3) * t216 + pkin(4) * t217 + t234 * t175;
t210 = (-qJ(5) - rSges(6,3)) * t181 + (rSges(6,1) * t185 - rSges(6,2) * t183 + t227) * t180;
t209 = qJD(4) * t180;
t208 = qJD(5) * t180;
t207 = V_base(6) * pkin(5) + V_base(2);
t176 = V_base(4) + qJD(1);
t147 = pkin(2) * t175 + qJ(3) * t174;
t204 = -t147 - t228;
t145 = -pkin(2) * t174 + qJ(3) * t175;
t203 = qJD(3) * t174 + t176 * t145 + V_base(3);
t202 = V_base(6) * qJ(2) + t207;
t201 = qJD(3) * t175 + t202;
t200 = pkin(3) * t181 + pkin(6) * t180;
t199 = V_base(5) * t228 + V_base(6) * t229 + qJD(2) + V_base(1);
t198 = rSges(4,1) * t181 - rSges(4,2) * t180;
t195 = Icges(4,5) * t181 - Icges(4,6) * t180;
t156 = Icges(4,2) * t181 + t220;
t157 = Icges(4,1) * t180 + t219;
t192 = t156 * t180 - t157 * t181;
t191 = V_base(5) * t147 + t199;
t190 = (Icges(4,3) * t175 - t174 * t195) * V_base(5) + (Icges(4,3) * t174 + t175 * t195) * V_base(6) + (Icges(4,5) * t180 + Icges(4,6) * t181) * t176;
t136 = t200 * t175;
t160 = pkin(3) * t180 - pkin(6) * t181;
t189 = V_base(6) * t160 + (-t136 + t204) * t176 + t201;
t135 = t200 * t174;
t188 = V_base(5) * t136 + (t135 - t145) * V_base(6) + t191;
t187 = (-t160 + t226) * V_base(5) + t203 + (-t135 - t229) * t176;
t168 = rSges(2,1) * t186 - t184 * rSges(2,2);
t167 = -t184 * rSges(2,1) - rSges(2,2) * t186;
t166 = Icges(2,1) * t186 - t224;
t165 = -Icges(2,1) * t184 - t223;
t164 = -Icges(2,2) * t184 + t223;
t163 = -Icges(2,2) * t186 - t224;
t159 = -qJD(4) * t181 + t176;
t158 = rSges(4,1) * t180 + rSges(4,2) * t181;
t154 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t153 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t152 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t150 = -t174 * t209 + V_base(5);
t149 = t175 * t209 + V_base(6);
t148 = rSges(3,1) * t175 - rSges(3,2) * t174;
t146 = -rSges(3,1) * t174 - rSges(3,2) * t175;
t142 = -Icges(3,2) * t174 + t221;
t141 = -Icges(3,2) * t175 - t222;
t130 = -rSges(5,3) * t181 + (rSges(5,1) * t185 - rSges(5,2) * t183) * t180;
t119 = V_base(6) * rSges(2,3) - t168 * t176 + t207;
t118 = t167 * t176 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t117 = rSges(4,3) * t174 + t175 * t198;
t116 = rSges(4,3) * t175 - t174 * t198;
t109 = -t167 * V_base(6) + t168 * V_base(5) + V_base(1);
t108 = V_base(6) * rSges(3,3) + (-t148 - t228) * t176 + t202;
t107 = V_base(3) + (t146 - t229) * t176 + (-rSges(3,3) + t226) * V_base(5);
t106 = -t146 * V_base(6) + t148 * V_base(5) + t199;
t105 = rSges(5,1) * t134 + rSges(5,2) * t133 + rSges(5,3) * t216;
t103 = rSges(5,1) * t132 + rSges(5,2) * t131 - rSges(5,3) * t218;
t87 = V_base(6) * t158 + (-t117 + t204) * t176 + t201;
t86 = (t116 - t229) * t176 + (-t158 + t226) * V_base(5) + t203;
t85 = t117 * V_base(5) + (-t116 - t145) * V_base(6) + t191;
t84 = -t159 * t105 + t149 * t130 + t189;
t83 = t103 * t159 - t130 * t150 + t187;
t82 = -t103 * t149 + t105 * t150 + t188;
t81 = t149 * t210 - t159 * t211 - t174 * t208 + t189;
t80 = -t150 * t210 + t159 * t212 + t175 * t208 + t187;
t79 = -qJD(5) * t181 - t149 * t212 + t150 * t211 + t188;
t1 = m(1) * (t152 ^ 2 + t153 ^ 2 + t154 ^ 2) / 0.2e1 + m(2) * (t109 ^ 2 + t118 ^ 2 + t119 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(4) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(5) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + ((t240 * t133 + t239 * t134 + t241 * t216) * t159 + (t245 * t133 + t243 * t134 + t247 * t216) * t150 + (t244 * t133 + t242 * t134 + t246 * t216) * t149) * t149 / 0.2e1 + ((t240 * t131 + t239 * t132 - t241 * t218) * t159 + (t245 * t131 + t243 * t132 - t247 * t218) * t150 + (t244 * t131 + t242 * t132 - t246 * t218) * t149) * t150 / 0.2e1 + ((-t246 * t149 - t247 * t150 - t241 * t159) * t181 + ((-t240 * t183 + t239 * t185) * t159 + (-t245 * t183 + t243 * t185) * t150 + (-t244 * t183 + t242 * t185) * t149) * t180) * t159 / 0.2e1 + ((t113 * t181 + t115 * t180 + t237) * V_base(6) + (t112 * t181 + t114 * t180 + t238) * V_base(5) + (t181 * t156 + t180 * t157 + Icges(2,3) + Icges(3,3)) * t176) * t176 / 0.2e1 + (t190 * t175 + (t192 * t174 + t238) * t176 + (-t142 * t175 - t164 * t186 - t184 * t166 - t235 * t174 + Icges(1,6)) * V_base(6) + (-t175 * t141 - t186 * t163 - t184 * t165 - t236 * t174 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t190 * t174 + (-t192 * t175 + t237) * t176 + (-t174 * t142 - t184 * t164 + t186 * t166 + t235 * t175 + Icges(1,3)) * V_base(6) + (-t141 * t174 - t184 * t163 + t165 * t186 + t236 * t175 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
