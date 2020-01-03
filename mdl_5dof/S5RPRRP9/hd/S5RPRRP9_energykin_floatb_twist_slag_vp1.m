% Calculate kinetic energy for
% S5RPRRP9
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:32
% EndTime: 2019-12-31 18:48:34
% DurationCPUTime: 2.50s
% Computational Cost: add. (1172->235), mult. (1149->323), div. (0->0), fcn. (931->8), ass. (0->130)
t267 = Icges(5,4) - Icges(6,5);
t266 = Icges(5,1) + Icges(6,1);
t265 = Icges(5,2) + Icges(6,3);
t179 = pkin(8) + qJ(3);
t171 = qJ(4) + t179;
t167 = cos(t171);
t264 = t267 * t167;
t166 = sin(t171);
t263 = t267 * t166;
t262 = Icges(6,4) + Icges(5,5);
t261 = Icges(5,6) - Icges(6,6);
t260 = t265 * t166 - t264;
t259 = t266 * t167 - t263;
t258 = rSges(6,1) + pkin(4);
t257 = rSges(6,3) + qJ(5);
t183 = sin(qJ(1));
t184 = cos(qJ(1));
t256 = t260 * t183 + t261 * t184;
t255 = -t261 * t183 + t260 * t184;
t254 = t259 * t183 - t262 * t184;
t253 = t262 * t183 + t259 * t184;
t252 = Icges(6,2) + Icges(5,3);
t251 = -t265 * t167 - t263;
t250 = t266 * t166 + t264;
t249 = -t261 * t166 + t262 * t167;
t248 = t257 * t166 + t258 * t167;
t142 = V_base(5) + (-qJD(3) - qJD(4)) * t184;
t164 = qJD(3) * t183 + V_base(4);
t143 = qJD(4) * t183 + t164;
t172 = V_base(6) + qJD(1);
t247 = (t251 * t166 + t250 * t167) * t172 + (t255 * t166 + t253 * t167) * t143 + (t256 * t166 + t254 * t167) * t142;
t246 = (t262 * t166 + t261 * t167) * t172 + (t252 * t183 + t249 * t184) * t143 + (t249 * t183 - t252 * t184) * t142;
t180 = sin(pkin(8));
t242 = pkin(2) * t180;
t169 = sin(t179);
t241 = pkin(3) * t169;
t181 = cos(pkin(8));
t240 = t181 * pkin(2);
t239 = Icges(2,4) * t183;
t238 = Icges(3,4) * t180;
t237 = Icges(3,4) * t181;
t236 = Icges(4,4) * t169;
t170 = cos(t179);
t235 = Icges(4,4) * t170;
t229 = -rSges(6,2) * t184 + t248 * t183;
t228 = t183 * rSges(6,2) + t248 * t184;
t107 = -pkin(6) * t184 + t240 * t183;
t159 = t183 * pkin(1) - qJ(2) * t184;
t227 = -t107 - t159;
t226 = t258 * t166 - t257 * t167;
t225 = pkin(3) * t170;
t223 = qJD(5) * t166;
t222 = V_base(4) * t159 + V_base(3);
t221 = V_base(5) * pkin(5) + V_base(1);
t84 = -pkin(7) * t184 + t183 * t225;
t218 = -t84 + t227;
t217 = qJD(2) * t183 + t221;
t216 = V_base(5) * t242 + t217;
t215 = rSges(3,1) * t181 - rSges(3,2) * t180;
t214 = rSges(4,1) * t170 - rSges(4,2) * t169;
t213 = rSges(5,1) * t167 - rSges(5,2) * t166;
t210 = Icges(3,1) * t181 - t238;
t209 = Icges(4,1) * t170 - t236;
t206 = -Icges(3,2) * t180 + t237;
t205 = -Icges(4,2) * t169 + t235;
t202 = Icges(3,5) * t181 - Icges(3,6) * t180;
t201 = Icges(4,5) * t170 - Icges(4,6) * t169;
t161 = pkin(1) * t184 + t183 * qJ(2);
t198 = -qJD(2) * t184 + t172 * t161 + V_base(2);
t163 = -qJD(3) * t184 + V_base(5);
t197 = t163 * t241 + t216;
t194 = (-Icges(4,3) * t184 + t183 * t201) * t163 + (Icges(4,3) * t183 + t184 * t201) * t164 + (Icges(4,5) * t169 + Icges(4,6) * t170) * t172;
t108 = pkin(6) * t183 + t240 * t184;
t193 = V_base(4) * t107 + (-t108 - t161) * V_base(5) + t222;
t192 = (-Icges(3,3) * t184 + t183 * t202) * V_base(5) + (Icges(3,3) * t183 + t184 * t202) * V_base(4) + (Icges(3,5) * t180 + Icges(3,6) * t181) * t172;
t85 = pkin(7) * t183 + t184 * t225;
t191 = -t163 * t85 + t164 * t84 + t193;
t190 = t172 * t108 + (-pkin(5) - t242) * V_base(4) + t198;
t189 = -t164 * t241 + t172 * t85 + t190;
t111 = -Icges(4,6) * t184 + t183 * t205;
t112 = Icges(4,6) * t183 + t184 * t205;
t113 = -Icges(4,5) * t184 + t183 * t209;
t114 = Icges(4,5) * t183 + t184 * t209;
t139 = Icges(4,2) * t170 + t236;
t140 = Icges(4,1) * t169 + t235;
t186 = (-t112 * t169 + t114 * t170) * t164 + (-t111 * t169 + t113 * t170) * t163 + (-t139 * t169 + t140 * t170) * t172;
t121 = -Icges(3,6) * t184 + t183 * t206;
t122 = Icges(3,6) * t183 + t184 * t206;
t123 = -Icges(3,5) * t184 + t183 * t210;
t124 = Icges(3,5) * t183 + t184 * t210;
t150 = Icges(3,2) * t181 + t238;
t151 = Icges(3,1) * t180 + t237;
t185 = (-t122 * t180 + t124 * t181) * V_base(4) + (-t121 * t180 + t123 * t181) * V_base(5) + (-t150 * t180 + t151 * t181) * t172;
t176 = Icges(2,4) * t184;
t162 = rSges(2,1) * t184 - t183 * rSges(2,2);
t160 = t183 * rSges(2,1) + rSges(2,2) * t184;
t158 = Icges(2,1) * t184 - t239;
t157 = Icges(2,1) * t183 + t176;
t156 = -Icges(2,2) * t183 + t176;
t155 = Icges(2,2) * t184 + t239;
t154 = Icges(2,5) * t184 - Icges(2,6) * t183;
t153 = Icges(2,5) * t183 + Icges(2,6) * t184;
t152 = rSges(3,1) * t180 + rSges(3,2) * t181;
t148 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t147 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t146 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t141 = rSges(4,1) * t169 + rSges(4,2) * t170;
t135 = rSges(5,1) * t166 + rSges(5,2) * t167;
t126 = t183 * rSges(3,3) + t184 * t215;
t125 = -rSges(3,3) * t184 + t183 * t215;
t116 = t183 * rSges(4,3) + t184 * t214;
t115 = -rSges(4,3) * t184 + t183 * t214;
t106 = t183 * rSges(5,3) + t184 * t213;
t104 = -rSges(5,3) * t184 + t183 * t213;
t102 = V_base(5) * rSges(2,3) - t160 * t172 + t221;
t101 = t162 * t172 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t88 = t160 * V_base(4) - t162 * V_base(5) + V_base(3);
t81 = t152 * V_base(5) + (-t125 - t159) * t172 + t217;
t80 = t172 * t126 + (-pkin(5) - t152) * V_base(4) + t198;
t79 = t125 * V_base(4) + (-t126 - t161) * V_base(5) + t222;
t78 = t141 * t163 + (-t115 + t227) * t172 + t216;
t77 = t172 * t116 - t164 * t141 + t190;
t76 = t115 * t164 - t116 * t163 + t193;
t75 = t135 * t142 + (-t104 + t218) * t172 + t197;
t74 = t172 * t106 - t143 * t135 + t189;
t73 = t184 * t223 + t226 * t142 + (t218 - t229) * t172 + t197;
t72 = -t226 * t143 + t172 * t228 + t183 * t223 + t189;
t71 = t104 * t143 - t106 * t142 + t191;
t70 = -qJD(5) * t167 - t228 * t142 + t229 * t143 + t191;
t1 = m(1) * (t146 ^ 2 + t147 ^ 2 + t148 ^ 2) / 0.2e1 + m(2) * (t101 ^ 2 + t102 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + t164 * (t183 * t194 + t184 * t186) / 0.2e1 + t163 * (t183 * t186 - t184 * t194) / 0.2e1 + m(5) * (t71 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + (t247 * t183 - t246 * t184) * t142 / 0.2e1 + (t246 * t183 + t247 * t184) * t143 / 0.2e1 + (t154 * t172 + t183 * t192 + t184 * t185 + (-t183 * t155 + t157 * t184 + Icges(1,4)) * V_base(5) + (-t183 * t156 + t158 * t184 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t153 * t172 + t183 * t185 - t184 * t192 + (t155 * t184 + t183 * t157 + Icges(1,2)) * V_base(5) + (t156 * t184 + t183 * t158 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t112 * t170 + t114 * t169) * t164 + (t111 * t170 + t113 * t169) * t163 + (t121 * t181 + t123 * t180 + t153) * V_base(5) + (t122 * t181 + t124 * t180 + t154) * V_base(4) + (t253 * t166 - t255 * t167) * t143 + (t254 * t166 - t256 * t167) * t142 + (t170 * t139 + t169 * t140 + t181 * t150 + t180 * t151 + t250 * t166 - t251 * t167 + Icges(2,3)) * t172) * t172 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
