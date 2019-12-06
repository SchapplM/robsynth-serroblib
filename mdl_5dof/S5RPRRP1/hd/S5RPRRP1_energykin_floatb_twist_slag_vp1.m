% Calculate kinetic energy for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:10
% EndTime: 2019-12-05 17:59:12
% DurationCPUTime: 2.11s
% Computational Cost: add. (770->203), mult. (962->269), div. (0->0), fcn. (744->6), ass. (0->112)
t257 = Icges(5,4) + Icges(6,4);
t256 = Icges(5,1) + Icges(6,1);
t255 = Icges(5,2) + Icges(6,2);
t153 = qJ(3) + qJ(4);
t148 = cos(t153);
t254 = t257 * t148;
t147 = sin(t153);
t253 = t257 * t147;
t252 = Icges(5,5) + Icges(6,5);
t251 = Icges(5,6) + Icges(6,6);
t250 = t148 * t255 + t253;
t249 = t147 * t256 + t254;
t248 = rSges(6,1) + pkin(4);
t247 = Icges(2,4) + Icges(3,6);
t155 = sin(qJ(1));
t157 = cos(qJ(1));
t246 = t155 * t250 + t157 * t251;
t245 = t155 * t251 - t157 * t250;
t244 = t155 * t249 + t157 * t252;
t243 = t155 * t252 - t157 * t249;
t242 = Icges(2,1) + Icges(3,2);
t241 = -Icges(3,4) + Icges(2,5);
t240 = Icges(3,5) - Icges(2,6);
t239 = Icges(2,2) + Icges(3,3);
t238 = Icges(5,3) + Icges(6,3);
t237 = -t147 * t255 + t254;
t236 = t148 * t256 - t253;
t235 = t147 * t252 + t148 * t251;
t234 = rSges(6,3) + qJ(5);
t233 = t247 * t157;
t232 = rSges(6,2) * t148 + t147 * t248;
t231 = t247 * t155;
t140 = qJD(3) * t155 + V_base(5);
t109 = qJD(4) * t155 + t140;
t141 = qJD(3) * t157 + V_base(4);
t110 = qJD(4) * t157 + t141;
t143 = V_base(6) + qJD(1);
t230 = t109 * (t147 * t243 + t148 * t245) + t110 * (t147 * t244 + t148 * t246) + t143 * (t147 * t236 + t148 * t237);
t229 = -t157 * t239 - t231;
t228 = t155 * t239 - t233;
t227 = t155 * t242 + t233;
t226 = t157 * t242 - t231;
t223 = (-t147 * t251 + t148 * t252) * t143 + (t155 * t235 + t157 * t238) * t110 + (t155 * t238 - t157 * t235) * t109;
t154 = sin(qJ(3));
t156 = cos(qJ(3));
t206 = Icges(4,4) * t156;
t124 = -Icges(4,2) * t154 + t206;
t207 = Icges(4,4) * t154;
t129 = Icges(4,1) * t156 - t207;
t176 = Icges(4,2) * t156 + t207;
t92 = Icges(4,6) * t157 + t155 * t176;
t93 = Icges(4,6) * t155 - t157 * t176;
t179 = Icges(4,1) * t154 + t206;
t94 = Icges(4,5) * t157 + t155 * t179;
t95 = Icges(4,5) * t155 - t157 * t179;
t219 = (t154 * t94 + t156 * t92) * t141 + (t154 * t95 + t156 * t93) * t140 + (t124 * t156 + t129 * t154) * t143;
t215 = pkin(3) * t154;
t214 = pkin(3) * t156;
t213 = t155 * pkin(6);
t212 = t157 * pkin(6);
t210 = t155 * t234 - t157 * t232;
t209 = t155 * t232 + t157 * t234;
t132 = pkin(1) * t155 - qJ(2) * t157;
t198 = t132 * V_base(4) + V_base(3);
t197 = V_base(5) * pkin(5) + V_base(1);
t193 = -rSges(6,2) * t147 + t148 * t248;
t192 = -t132 - t213;
t191 = qJD(2) * t155 + t197;
t98 = pkin(7) * t155 - t157 * t215;
t190 = t192 - t98;
t189 = V_base(5) * pkin(2) + t191;
t188 = rSges(4,1) * t154 + rSges(4,2) * t156;
t187 = rSges(5,1) * t147 + rSges(5,2) * t148;
t173 = Icges(4,5) * t154 + Icges(4,6) * t156;
t136 = pkin(1) * t157 + qJ(2) * t155;
t167 = -qJD(2) * t157 + t136 * t143 + V_base(2);
t166 = t140 * t214 + t189;
t163 = (Icges(4,5) * t156 - Icges(4,6) * t154) * t143 + (Icges(4,3) * t155 - t157 * t173) * t140 + (Icges(4,3) * t157 + t155 * t173) * t141;
t162 = V_base(4) * t213 + (-t136 - t212) * V_base(5) + t198;
t161 = t143 * t212 + (-pkin(2) - pkin(5)) * V_base(4) + t167;
t99 = pkin(7) * t157 + t155 * t215;
t160 = -t140 * t99 + t141 * t98 + t162;
t159 = -t141 * t214 + t143 * t99 + t161;
t138 = rSges(2,1) * t157 - rSges(2,2) * t155;
t137 = -rSges(3,2) * t157 + rSges(3,3) * t155;
t135 = rSges(4,1) * t156 - rSges(4,2) * t154;
t134 = rSges(2,1) * t155 + rSges(2,2) * t157;
t133 = -rSges(3,2) * t155 - rSges(3,3) * t157;
t115 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t114 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t113 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t108 = rSges(5,1) * t148 - rSges(5,2) * t147;
t97 = rSges(4,3) * t155 - t157 * t188;
t96 = rSges(4,3) * t157 + t155 * t188;
t88 = rSges(5,3) * t155 - t157 * t187;
t86 = rSges(5,3) * t157 + t155 * t187;
t72 = V_base(5) * rSges(2,3) - t134 * t143 + t197;
t71 = t138 * t143 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t69 = t134 * V_base(4) - t138 * V_base(5) + V_base(3);
t66 = V_base(5) * rSges(3,1) + (-t132 - t133) * t143 + t191;
t65 = t137 * t143 + (-rSges(3,1) - pkin(5)) * V_base(4) + t167;
t64 = t133 * V_base(4) + (-t136 - t137) * V_base(5) + t198;
t63 = t135 * t140 + (t192 - t97) * t143 + t189;
t62 = -t135 * t141 + t143 * t96 + t161;
t61 = -t140 * t96 + t141 * t97 + t162;
t60 = t108 * t109 + (t190 - t88) * t143 + t166;
t59 = -t108 * t110 + t143 * t86 + t159;
t58 = -t109 * t86 + t110 * t88 + t160;
t57 = qJD(5) * t157 + t193 * t109 + (t190 - t210) * t143 + t166;
t56 = qJD(5) * t155 - t110 * t193 + t143 * t209 + t159;
t55 = -t109 * t209 + t110 * t210 + t160;
t1 = m(1) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(2) * (t69 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(3) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + t141 * (t219 * t155 + t163 * t157) / 0.2e1 + t140 * (t163 * t155 - t219 * t157) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + (t223 * t155 - t157 * t230) * t109 / 0.2e1 + (t155 * t230 + t223 * t157) * t110 / 0.2e1 + ((t155 * t229 + t157 * t227 + Icges(1,4)) * V_base(5) + (t155 * t228 + t157 * t226 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t227 * t155 - t157 * t229 + Icges(1,2)) * V_base(5) + (t155 * t226 - t157 * t228 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t154 * t92 + t156 * t94) * t141 + (-t154 * t93 + t156 * t95) * t140 + (-t147 * t246 + t148 * t244) * t110 + (-t147 * t245 + t148 * t243) * t109 + (-t154 * t124 + t156 * t129 - t237 * t147 + t236 * t148 + Icges(3,1) + Icges(2,3)) * t143) * t143 / 0.2e1 + t143 * V_base(5) * (t155 * t241 - t157 * t240) + t143 * V_base(4) * (t155 * t240 + t157 * t241) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
