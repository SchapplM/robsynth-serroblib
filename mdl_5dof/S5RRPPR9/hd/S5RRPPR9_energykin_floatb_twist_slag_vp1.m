% Calculate kinetic energy for
% S5RRPPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR9_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR9_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:36
% EndTime: 2019-12-31 19:40:39
% DurationCPUTime: 2.53s
% Computational Cost: add. (725->232), mult. (1337->326), div. (0->0), fcn. (1175->6), ass. (0->120)
t259 = Icges(3,4) + Icges(5,4) - Icges(4,5);
t258 = Icges(3,1) + Icges(4,1) + Icges(5,2);
t257 = Icges(5,1) + Icges(3,2) + Icges(4,3);
t185 = cos(qJ(2));
t256 = t259 * t185;
t182 = sin(qJ(2));
t255 = t259 * t182;
t254 = Icges(4,4) + Icges(3,5) + Icges(5,6);
t253 = Icges(5,5) + Icges(3,6) - Icges(4,6);
t252 = t257 * t182 - t256;
t251 = t258 * t185 - t255;
t250 = Icges(4,2) + Icges(3,3) + Icges(5,3);
t183 = sin(qJ(1));
t186 = cos(qJ(1));
t249 = t252 * t183 + t253 * t186;
t248 = -t253 * t183 + t252 * t186;
t247 = t251 * t183 - t254 * t186;
t246 = t254 * t183 + t251 * t186;
t245 = -t257 * t185 - t255;
t244 = t258 * t182 + t256;
t243 = -t253 * t182 + t254 * t185;
t172 = -qJD(2) * t186 + V_base(5);
t173 = qJD(2) * t183 + V_base(4);
t176 = V_base(6) + qJD(1);
t242 = (t245 * t182 + t244 * t185) * t176 + (t248 * t182 + t246 * t185) * t173 + (t249 * t182 + t247 * t185) * t172;
t241 = (t254 * t182 + t253 * t185) * t176 + (t250 * t183 + t243 * t186) * t173 + (t243 * t183 - t250 * t186) * t172;
t237 = pkin(3) * t182;
t236 = Icges(2,4) * t183;
t181 = sin(qJ(5));
t229 = t181 * t186;
t228 = t183 * t181;
t184 = cos(qJ(5));
t227 = t183 * t184;
t226 = t183 * t185;
t225 = t184 * t186;
t224 = t185 * t186;
t209 = pkin(2) * t185 + qJ(3) * t182;
t135 = t209 * t183;
t169 = t183 * pkin(1) - pkin(6) * t186;
t223 = -t135 - t169;
t136 = t209 * t186;
t143 = pkin(3) * t224 - t183 * qJ(4);
t222 = -t136 - t143;
t221 = qJD(3) * t182;
t220 = qJD(5) * t185;
t219 = V_base(5) * pkin(5) + V_base(1);
t142 = pkin(3) * t226 + qJ(4) * t186;
t216 = -t142 + t223;
t163 = pkin(2) * t182 - qJ(3) * t185;
t215 = -t163 - t237;
t214 = t172 * t163 + t186 * t221 + t219;
t213 = pkin(4) * t182 + pkin(7) * t185;
t212 = rSges(3,1) * t185 - rSges(3,2) * t182;
t211 = rSges(4,1) * t185 + rSges(4,3) * t182;
t210 = rSges(5,1) * t182 - rSges(5,2) * t185;
t171 = pkin(1) * t186 + t183 * pkin(6);
t199 = -V_base(4) * pkin(5) + t176 * t171 + V_base(2);
t198 = V_base(4) * t169 - t171 * V_base(5) + V_base(3);
t194 = t176 * t136 + t183 * t221 + t199;
t193 = -qJD(4) * t183 + t172 * t237 + t214;
t192 = -qJD(3) * t185 + t173 * t135 + t198;
t191 = qJD(4) * t186 + t176 * t143 + t194;
t190 = t173 * t142 + t192;
t179 = Icges(2,4) * t186;
t170 = -pkin(4) * t185 + pkin(7) * t182;
t168 = rSges(2,1) * t186 - t183 * rSges(2,2);
t167 = -rSges(5,1) * t185 - rSges(5,2) * t182;
t166 = t183 * rSges(2,1) + rSges(2,2) * t186;
t165 = rSges(3,1) * t182 + rSges(3,2) * t185;
t164 = rSges(4,1) * t182 - rSges(4,3) * t185;
t162 = qJD(5) * t182 + t176;
t161 = Icges(2,1) * t186 - t236;
t160 = Icges(2,1) * t183 + t179;
t156 = -Icges(2,2) * t183 + t179;
t155 = Icges(2,2) * t186 + t236;
t146 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t145 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t144 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t139 = t213 * t186;
t138 = t213 * t183;
t134 = t182 * t225 - t228;
t133 = -t182 * t229 - t227;
t132 = t182 * t227 + t229;
t131 = -t182 * t228 + t225;
t130 = t186 * t220 + t173;
t129 = t183 * t220 + t172;
t126 = t183 * rSges(3,3) + t186 * t212;
t125 = t183 * rSges(4,2) + t186 * t211;
t124 = -t183 * rSges(5,3) + t186 * t210;
t123 = rSges(6,3) * t182 + (-rSges(6,1) * t184 + rSges(6,2) * t181) * t185;
t122 = -rSges(3,3) * t186 + t183 * t212;
t121 = -rSges(4,2) * t186 + t183 * t211;
t120 = rSges(5,3) * t186 + t183 * t210;
t113 = Icges(6,5) * t182 + (-Icges(6,1) * t184 + Icges(6,4) * t181) * t185;
t106 = Icges(6,6) * t182 + (-Icges(6,4) * t184 + Icges(6,2) * t181) * t185;
t99 = Icges(6,3) * t182 + (-Icges(6,5) * t184 + Icges(6,6) * t181) * t185;
t95 = V_base(5) * rSges(2,3) - t166 * t176 + t219;
t94 = t168 * t176 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t93 = t166 * V_base(4) - t168 * V_base(5) + V_base(3);
t92 = t134 * rSges(6,1) + t133 * rSges(6,2) + rSges(6,3) * t224;
t91 = rSges(6,1) * t132 + rSges(6,2) * t131 + rSges(6,3) * t226;
t90 = Icges(6,1) * t134 + Icges(6,4) * t133 + Icges(6,5) * t224;
t89 = Icges(6,1) * t132 + Icges(6,4) * t131 + Icges(6,5) * t226;
t88 = Icges(6,4) * t134 + Icges(6,2) * t133 + Icges(6,6) * t224;
t87 = Icges(6,4) * t132 + Icges(6,2) * t131 + Icges(6,6) * t226;
t86 = Icges(6,5) * t134 + Icges(6,6) * t133 + Icges(6,3) * t224;
t85 = Icges(6,5) * t132 + Icges(6,6) * t131 + Icges(6,3) * t226;
t84 = t165 * t172 + (-t122 - t169) * t176 + t219;
t83 = t126 * t176 - t165 * t173 + t199;
t82 = t122 * t173 - t126 * t172 + t198;
t81 = t164 * t172 + (-t121 + t223) * t176 + t214;
t80 = t125 * t176 + (-t163 - t164) * t173 + t194;
t79 = t121 * t173 + (-t125 - t136) * t172 + t192;
t78 = t167 * t172 + (-t120 + t216) * t176 + t193;
t77 = t124 * t176 + (-t167 + t215) * t173 + t191;
t76 = t120 * t173 + (-t124 + t222) * t172 + t190;
t75 = t123 * t129 - t162 * t91 + t170 * t172 + (-t138 + t216) * t176 + t193;
t74 = -t123 * t130 + t139 * t176 + t162 * t92 + (-t170 + t215) * t173 + t191;
t73 = -t129 * t92 + t130 * t91 + t138 * t173 + (-t139 + t222) * t172 + t190;
t1 = m(1) * (t144 ^ 2 + t145 ^ 2 + t146 ^ 2) / 0.2e1 + m(2) * (t93 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(3) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(5) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + t130 * ((t133 * t88 + t134 * t90 + t86 * t224) * t130 + (t133 * t87 + t134 * t89 + t224 * t85) * t129 + (t133 * t106 + t134 * t113 + t224 * t99) * t162) / 0.2e1 + t129 * ((t131 * t88 + t132 * t90 + t226 * t86) * t130 + (t131 * t87 + t132 * t89 + t85 * t226) * t129 + (t106 * t131 + t113 * t132 + t226 * t99) * t162) / 0.2e1 + t162 * ((t85 * t129 + t86 * t130 + t99 * t162) * t182 + ((t181 * t88 - t184 * t90) * t130 + (t181 * t87 - t184 * t89) * t129 + (t106 * t181 - t113 * t184) * t162) * t185) / 0.2e1 + ((-t183 * t155 + t160 * t186 + Icges(1,4)) * V_base(5) + (-t183 * t156 + t186 * t161 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t186 * t155 + t183 * t160 + Icges(1,2)) * V_base(5) + (t156 * t186 + t183 * t161 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (t242 * t183 - t241 * t186) * t172 / 0.2e1 + (t241 * t183 + t242 * t186) * t173 / 0.2e1 + ((t246 * t182 - t248 * t185) * t173 + (t247 * t182 - t249 * t185) * t172 + (t244 * t182 - t245 * t185 + Icges(2,3)) * t176) * t176 / 0.2e1 + V_base(4) * t176 * (Icges(2,5) * t186 - Icges(2,6) * t183) + V_base(5) * t176 * (Icges(2,5) * t183 + Icges(2,6) * t186) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
