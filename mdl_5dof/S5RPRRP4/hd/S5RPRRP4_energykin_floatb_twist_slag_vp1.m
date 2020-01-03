% Calculate kinetic energy for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:09
% EndTime: 2020-01-03 11:49:13
% DurationCPUTime: 3.58s
% Computational Cost: add. (1307->287), mult. (1945->410), div. (0->0), fcn. (1895->8), ass. (0->138)
t273 = Icges(5,1) + Icges(6,1);
t272 = Icges(5,4) + Icges(6,4);
t271 = -Icges(6,5) - Icges(5,5);
t270 = Icges(5,2) + Icges(6,2);
t269 = -Icges(6,6) - Icges(5,6);
t268 = -Icges(6,3) - Icges(5,3);
t199 = qJ(3) + qJ(4);
t193 = sin(t199);
t194 = cos(t199);
t205 = cos(qJ(1));
t201 = cos(pkin(8));
t203 = sin(qJ(1));
t240 = t201 * t203;
t150 = -t193 * t240 - t194 * t205;
t151 = -t193 * t205 + t194 * t240;
t200 = sin(pkin(8));
t242 = t200 * t203;
t267 = -t269 * t150 - t271 * t151 - t268 * t242;
t239 = t201 * t205;
t152 = t193 * t239 - t194 * t203;
t153 = -t193 * t203 - t194 * t239;
t241 = t200 * t205;
t266 = -t269 * t152 - t271 * t153 + t268 * t241;
t265 = t270 * t150 + t272 * t151 - t269 * t242;
t264 = t270 * t152 + t272 * t153 + t269 * t241;
t263 = t272 * t150 + t273 * t151 - t271 * t242;
t262 = t272 * t152 + t273 * t153 + t271 * t241;
t261 = t268 * t201 + (t269 * t193 - t271 * t194) * t200;
t260 = t269 * t201 + (-t270 * t193 + t272 * t194) * t200;
t259 = t271 * t201 + (-t272 * t193 + t273 * t194) * t200;
t243 = Icges(3,4) * t201;
t220 = -Icges(3,2) * t200 + t243;
t143 = -Icges(3,6) * t205 + t203 * t220;
t244 = Icges(3,4) * t200;
t221 = Icges(3,1) * t201 - t244;
t145 = -Icges(3,5) * t205 + t203 * t221;
t245 = Icges(2,4) * t205;
t258 = Icges(2,1) * t203 - t143 * t200 + t145 * t201 + t245;
t144 = -Icges(3,6) * t203 - t205 * t220;
t146 = -Icges(3,5) * t203 - t205 * t221;
t195 = Icges(2,4) * t203;
t257 = -Icges(2,1) * t205 - t144 * t200 + t146 * t201 + t195;
t233 = pkin(4) * t194;
t256 = qJ(5) * t200 + t201 * t233;
t204 = cos(qJ(3));
t249 = t204 * pkin(3);
t255 = pkin(7) * t200 + t201 * t249;
t224 = pkin(4) * t193;
t247 = rSges(6,1) * t151 + rSges(6,2) * t150 + rSges(6,3) * t242 + t203 * t256 - t224 * t205;
t246 = rSges(6,1) * t153 + rSges(6,2) * t152 - rSges(6,3) * t241 - t224 * t203 - t205 * t256;
t202 = sin(qJ(3));
t238 = t202 * t203;
t237 = t202 * t205;
t236 = t203 * t204;
t235 = t204 * t205;
t234 = (-qJ(5) - rSges(6,3)) * t201 + (rSges(6,1) * t194 - rSges(6,2) * t193 + t233) * t200;
t231 = qJD(3) * t200;
t230 = qJD(5) * t200;
t229 = -qJD(3) - qJD(4);
t188 = -pkin(1) * t205 - qJ(2) * t203;
t228 = V_base(5) * t188 + V_base(1);
t227 = V_base(6) * pkin(5) + V_base(2);
t173 = t203 * t231 + V_base(5);
t192 = V_base(4) + qJD(1);
t223 = pkin(2) * t201 + pkin(6) * t200;
t222 = rSges(3,1) * t201 - rSges(3,2) * t200;
t219 = Icges(3,5) * t201 - Icges(3,6) * t200;
t175 = Icges(3,2) * t201 + t244;
t176 = Icges(3,1) * t200 + t243;
t216 = t175 * t200 - t176 * t201;
t186 = pkin(1) * t203 - qJ(2) * t205;
t215 = -qJD(2) * t203 + t192 * t186 + V_base(3);
t214 = -qJD(2) * t205 + t227;
t161 = t223 * t203;
t162 = t223 * t205;
t213 = -V_base(5) * t162 + (-t161 - t186) * V_base(6) + t228;
t212 = -(-Icges(3,3) * t205 + t203 * t219) * V_base(5) - (-Icges(3,3) * t203 - t205 * t219) * V_base(6) - (Icges(3,5) * t200 + Icges(3,6) * t201) * t192;
t179 = t200 * pkin(2) - t201 * pkin(6);
t211 = V_base(6) * t179 + (t162 - t188) * t192 + t214;
t210 = t192 * t161 + (-pkin(5) - t179) * V_base(5) + t215;
t119 = -pkin(3) * t237 + t203 * t255;
t120 = -pkin(3) * t238 - t205 * t255;
t172 = -t205 * t231 + V_base(6);
t209 = -t119 * t172 + t173 * t120 + t213;
t128 = -pkin(7) * t201 + t200 * t249;
t178 = -qJD(3) * t201 + t192;
t208 = -t120 * t178 + t172 * t128 + t211;
t207 = t178 * t119 - t128 * t173 + t210;
t189 = -rSges(2,1) * t205 + rSges(2,2) * t203;
t187 = rSges(2,1) * t203 + rSges(2,2) * t205;
t183 = Icges(2,2) * t203 - t245;
t182 = Icges(2,2) * t205 + t195;
t181 = -Icges(2,5) * t205 + Icges(2,6) * t203;
t180 = Icges(2,5) * t203 + Icges(2,6) * t205;
t177 = rSges(3,1) * t200 + rSges(3,2) * t201;
t170 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t169 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t168 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t163 = t201 * t229 + t192;
t160 = -t201 * t235 - t238;
t159 = t201 * t237 - t236;
t158 = t201 * t236 - t237;
t157 = -t201 * t238 - t235;
t155 = qJD(4) * t242 + t173;
t154 = t229 * t241 + V_base(6);
t149 = -rSges(3,3) * t203 - t205 * t222;
t148 = -rSges(3,3) * t205 + t203 * t222;
t147 = -rSges(4,3) * t201 + (rSges(4,1) * t204 - rSges(4,2) * t202) * t200;
t139 = -Icges(4,5) * t201 + (Icges(4,1) * t204 - Icges(4,4) * t202) * t200;
t138 = -Icges(4,6) * t201 + (Icges(4,4) * t204 - Icges(4,2) * t202) * t200;
t137 = -Icges(4,3) * t201 + (Icges(4,5) * t204 - Icges(4,6) * t202) * t200;
t136 = -rSges(5,3) * t201 + (rSges(5,1) * t194 - rSges(5,2) * t193) * t200;
t127 = V_base(6) * rSges(2,3) - t189 * t192 + t227;
t126 = t187 * t192 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t125 = -t187 * V_base(6) + t189 * V_base(5) + V_base(1);
t122 = rSges(4,1) * t160 + rSges(4,2) * t159 - rSges(4,3) * t241;
t121 = rSges(4,1) * t158 + rSges(4,2) * t157 + rSges(4,3) * t242;
t118 = Icges(4,1) * t160 + Icges(4,4) * t159 - Icges(4,5) * t241;
t117 = Icges(4,1) * t158 + Icges(4,4) * t157 + Icges(4,5) * t242;
t116 = Icges(4,4) * t160 + Icges(4,2) * t159 - Icges(4,6) * t241;
t115 = Icges(4,4) * t158 + Icges(4,2) * t157 + Icges(4,6) * t242;
t114 = Icges(4,5) * t160 + Icges(4,6) * t159 - Icges(4,3) * t241;
t113 = Icges(4,5) * t158 + Icges(4,6) * t157 + Icges(4,3) * t242;
t112 = rSges(5,1) * t153 + rSges(5,2) * t152 - rSges(5,3) * t241;
t110 = rSges(5,1) * t151 + rSges(5,2) * t150 + rSges(5,3) * t242;
t94 = t177 * V_base(6) + (-t149 - t188) * t192 + t214;
t93 = t148 * t192 + (-pkin(5) - t177) * V_base(5) + t215;
t90 = t149 * V_base(5) + (-t148 - t186) * V_base(6) + t228;
t89 = -t122 * t178 + t147 * t172 + t211;
t88 = t121 * t178 - t147 * t173 + t210;
t87 = -t121 * t172 + t122 * t173 + t213;
t86 = -t112 * t163 + t136 * t154 + t208;
t85 = t110 * t163 - t136 * t155 + t207;
t84 = -t110 * t154 + t112 * t155 + t209;
t83 = t154 * t234 - t163 * t246 + t203 * t230 + t208;
t82 = -t155 * t234 + t163 * t247 - t205 * t230 + t207;
t81 = -qJD(5) * t201 - t154 * t247 + t155 * t246 + t209;
t1 = m(1) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(2) * (t125 ^ 2 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + m(3) * (t90 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t178 * ((-t113 * t173 - t114 * t172 - t137 * t178) * t201 + ((-t138 * t202 + t139 * t204) * t178 + (-t115 * t202 + t117 * t204) * t173 + (-t116 * t202 + t118 * t204) * t172) * t200) / 0.2e1 + t173 * ((t137 * t242 + t138 * t157 + t139 * t158) * t178 + (t113 * t242 + t157 * t115 + t158 * t117) * t173 + (t114 * t242 + t116 * t157 + t118 * t158) * t172) / 0.2e1 + t172 * ((-t137 * t241 + t138 * t159 + t139 * t160) * t178 + (-t113 * t241 + t115 * t159 + t117 * t160) * t173 + (-t114 * t241 + t159 * t116 + t160 * t118) * t172) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + ((t152 * t260 + t153 * t259 - t241 * t261) * t163 + (t265 * t152 + t263 * t153 - t241 * t267) * t155 + (t264 * t152 + t262 * t153 - t266 * t241) * t154) * t154 / 0.2e1 + ((t150 * t260 + t151 * t259 + t242 * t261) * t163 + (t265 * t150 + t263 * t151 + t267 * t242) * t155 + (t150 * t264 + t151 * t262 + t242 * t266) * t154) * t155 / 0.2e1 + ((-t266 * t154 - t155 * t267 - t261 * t163) * t201 + ((-t193 * t260 + t194 * t259) * t163 + (-t193 * t265 + t263 * t194) * t155 + (-t193 * t264 + t262 * t194) * t154) * t200) * t163 / 0.2e1 + ((t144 * t201 + t146 * t200 + t181) * V_base(6) + (t143 * t201 + t145 * t200 + t180) * V_base(5) + (t201 * t175 + t200 * t176 + Icges(2,3)) * t192) * t192 / 0.2e1 + (t212 * t205 + (-t216 * t203 + t180) * t192 + (t183 * t205 + t203 * t257 + Icges(1,6)) * V_base(6) + (t205 * t182 + t203 * t258 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t212 * t203 + (t216 * t205 + t181) * t192 + (t203 * t183 - t205 * t257 + Icges(1,3)) * V_base(6) + (t182 * t203 - t258 * t205 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
