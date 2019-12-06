% Calculate kinetic energy for
% S5PRPRP2
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:27
% EndTime: 2019-12-05 15:30:29
% DurationCPUTime: 2.31s
% Computational Cost: add. (1200->246), mult. (1371->332), div. (0->0), fcn. (1263->8), ass. (0->123)
t254 = Icges(5,1) + Icges(6,1);
t253 = Icges(5,4) + Icges(6,4);
t252 = -Icges(6,5) - Icges(5,5);
t251 = Icges(5,2) + Icges(6,2);
t250 = -Icges(6,6) - Icges(5,6);
t249 = -Icges(6,3) - Icges(5,3);
t183 = pkin(7) + qJ(2);
t177 = sin(t183);
t178 = cos(t183);
t190 = cos(qJ(4));
t186 = cos(pkin(8));
t189 = sin(qJ(4));
t221 = t186 * t189;
t134 = -t177 * t221 - t178 * t190;
t220 = t186 * t190;
t222 = t178 * t189;
t135 = t177 * t220 - t222;
t184 = sin(pkin(8));
t225 = t177 * t184;
t248 = -t250 * t134 - t252 * t135 - t249 * t225;
t136 = t177 * t190 - t178 * t221;
t224 = t177 * t189;
t137 = t178 * t220 + t224;
t223 = t178 * t184;
t247 = -t250 * t136 - t252 * t137 - t249 * t223;
t246 = t251 * t134 + t253 * t135 - t250 * t225;
t245 = t251 * t136 + t253 * t137 - t250 * t223;
t244 = t253 * t134 + t254 * t135 - t252 * t225;
t241 = t253 * t136 + t254 * t137 - t252 * t223;
t240 = t249 * t186 + (t250 * t189 - t252 * t190) * t184;
t239 = t250 * t186 + (-t251 * t189 + t253 * t190) * t184;
t238 = t252 * t186 + (-t253 * t189 + t254 * t190) * t184;
t187 = cos(pkin(7));
t233 = pkin(1) * t187;
t232 = pkin(4) * t190;
t231 = -pkin(5) - qJ(1);
t185 = sin(pkin(7));
t229 = Icges(2,4) * t185;
t228 = Icges(3,4) * t177;
t227 = Icges(4,4) * t184;
t226 = Icges(4,4) * t186;
t198 = qJ(5) * t184 + t186 * t232;
t219 = rSges(6,1) * t135 + rSges(6,2) * t134 + rSges(6,3) * t225 - pkin(4) * t222 + t177 * t198;
t218 = rSges(6,1) * t137 + rSges(6,2) * t136 + rSges(6,3) * t223 + pkin(4) * t224 + t178 * t198;
t217 = (-qJ(5) - rSges(6,3)) * t186 + (rSges(6,1) * t190 - rSges(6,2) * t189 + t232) * t184;
t216 = qJD(4) * t184;
t215 = qJD(5) * t184;
t208 = pkin(1) * V_base(6);
t214 = t187 * t208 + V_base(2);
t213 = V_base(5) * qJ(1) + V_base(1);
t209 = qJD(1) + V_base(3);
t180 = V_base(6) + qJD(2);
t150 = pkin(2) * t178 + qJ(3) * t177;
t207 = -t150 - t233;
t206 = V_base(4) * t185 * pkin(1) + t209;
t205 = pkin(3) * t186 + pkin(6) * t184;
t148 = pkin(2) * t177 - qJ(3) * t178;
t204 = V_base(4) * t148 + t206;
t203 = rSges(4,1) * t186 - rSges(4,2) * t184;
t202 = Icges(4,1) * t186 - t227;
t201 = -Icges(4,2) * t184 + t226;
t200 = Icges(4,5) * t186 - Icges(4,6) * t184;
t199 = -qJD(3) * t178 + t180 * t150 + t214;
t197 = V_base(5) * pkin(5) - t185 * t208 + t213;
t196 = qJD(3) * t177 + t197;
t195 = (-Icges(4,3) * t178 + t177 * t200) * V_base(5) + (Icges(4,3) * t177 + t178 * t200) * V_base(4) + (Icges(4,5) * t184 + Icges(4,6) * t186) * t180;
t138 = t205 * t177;
t139 = t205 * t178;
t194 = V_base(4) * t138 + (-t139 + t207) * V_base(5) + t204;
t171 = pkin(3) * t184 - pkin(6) * t186;
t193 = t180 * t139 + (-t171 + t231) * V_base(4) + t199;
t192 = V_base(5) * t171 + (-t138 - t148) * t180 + t196;
t115 = -Icges(4,6) * t178 + t177 * t201;
t116 = Icges(4,6) * t177 + t178 * t201;
t117 = -Icges(4,5) * t178 + t177 * t202;
t118 = Icges(4,5) * t177 + t178 * t202;
t161 = Icges(4,2) * t186 + t227;
t164 = Icges(4,1) * t184 + t226;
t191 = (-t116 * t184 + t118 * t186) * V_base(4) + (-t115 * t184 + t117 * t186) * V_base(5) + (-t161 * t184 + t164 * t186) * t180;
t179 = Icges(2,4) * t187;
t175 = Icges(3,4) * t178;
t170 = -qJD(4) * t186 + t180;
t169 = rSges(2,1) * t187 - rSges(2,2) * t185;
t168 = rSges(2,1) * t185 + rSges(2,2) * t187;
t167 = rSges(4,1) * t184 + rSges(4,2) * t186;
t166 = Icges(2,1) * t187 - t229;
t165 = Icges(2,1) * t185 + t179;
t163 = -Icges(2,2) * t185 + t179;
t162 = Icges(2,2) * t187 + t229;
t157 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t156 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t155 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t153 = t178 * t216 + V_base(4);
t152 = t177 * t216 + V_base(5);
t151 = rSges(3,1) * t178 - rSges(3,2) * t177;
t149 = rSges(3,1) * t177 + rSges(3,2) * t178;
t147 = Icges(3,1) * t178 - t228;
t146 = Icges(3,1) * t177 + t175;
t145 = -Icges(3,2) * t177 + t175;
t144 = Icges(3,2) * t178 + t228;
t143 = Icges(3,5) * t178 - Icges(3,6) * t177;
t142 = Icges(3,5) * t177 + Icges(3,6) * t178;
t133 = -t186 * rSges(5,3) + (rSges(5,1) * t190 - rSges(5,2) * t189) * t184;
t123 = V_base(5) * rSges(2,3) - t168 * V_base(6) + t213;
t122 = t169 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t120 = rSges(4,3) * t177 + t178 * t203;
t119 = -rSges(4,3) * t178 + t177 * t203;
t112 = t168 * V_base(4) - t169 * V_base(5) + t209;
t111 = V_base(5) * rSges(3,3) - t149 * t180 + t197;
t110 = t151 * t180 + (-rSges(3,3) + t231) * V_base(4) + t214;
t109 = t149 * V_base(4) + (-t151 - t233) * V_base(5) + t206;
t108 = rSges(5,1) * t137 + rSges(5,2) * t136 + rSges(5,3) * t223;
t106 = rSges(5,1) * t135 + rSges(5,2) * t134 + rSges(5,3) * t225;
t90 = t167 * V_base(5) + (-t119 - t148) * t180 + t196;
t89 = t120 * t180 + (-t167 + t231) * V_base(4) + t199;
t88 = t119 * V_base(4) + (-t120 + t207) * V_base(5) + t204;
t87 = -t106 * t170 + t133 * t152 + t192;
t86 = t108 * t170 - t133 * t153 + t193;
t85 = t106 * t153 - t108 * t152 + t194;
t84 = t152 * t217 - t170 * t219 + t178 * t215 + t192;
t83 = -t153 * t217 + t170 * t218 + t177 * t215 + t193;
t82 = -qJD(5) * t186 - t152 * t218 + t153 * t219 + t194;
t1 = m(1) * (t155 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + m(2) * (t112 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(3) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(4) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + ((t239 * t134 + t238 * t135 + t240 * t225) * t170 + (t245 * t134 + t241 * t135 + t247 * t225) * t153 + (t246 * t134 + t244 * t135 + t248 * t225) * t152) * t152 / 0.2e1 + ((t239 * t136 + t238 * t137 + t240 * t223) * t170 + (t245 * t136 + t241 * t137 + t247 * t223) * t153 + (t246 * t136 + t244 * t137 + t248 * t223) * t152) * t153 / 0.2e1 + ((-t248 * t152 - t247 * t153 - t240 * t170) * t186 + ((-t239 * t189 + t238 * t190) * t170 + (-t245 * t189 + t241 * t190) * t153 + (-t246 * t189 + t244 * t190) * t152) * t184) * t170 / 0.2e1 + ((t115 * t186 + t117 * t184 + t142) * V_base(5) + (t116 * t186 + t118 * t184 + t143) * V_base(4) + (t186 * t161 + t184 * t164 + Icges(3,3)) * t180) * t180 / 0.2e1 + (t143 * t180 + t177 * t195 + t178 * t191 + (-t144 * t177 + t146 * t178 - t162 * t185 + t165 * t187 + Icges(1,4)) * V_base(5) + (-t177 * t145 + t178 * t147 - t185 * t163 + t187 * t166 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t142 * t180 + t177 * t191 - t178 * t195 + (t178 * t144 + t177 * t146 + t187 * t162 + t185 * t165 + Icges(1,2)) * V_base(5) + (t145 * t178 + t147 * t177 + t163 * t187 + t166 * t185 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t185 + Icges(2,6) * t187 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t187 - Icges(2,6) * t185 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
