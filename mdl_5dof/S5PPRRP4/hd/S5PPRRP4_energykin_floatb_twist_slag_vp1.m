% Calculate kinetic energy for
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:23
% EndTime: 2019-12-31 17:34:25
% DurationCPUTime: 1.93s
% Computational Cost: add. (724->189), mult. (1282->236), div. (0->0), fcn. (1326->6), ass. (0->100)
t237 = Icges(5,4) + Icges(6,4);
t236 = Icges(5,1) + Icges(6,1);
t235 = Icges(5,2) + Icges(6,2);
t154 = cos(qJ(4));
t234 = t237 * t154;
t153 = sin(qJ(4));
t233 = t237 * t153;
t232 = Icges(5,5) + Icges(6,5);
t231 = Icges(5,6) + Icges(6,6);
t230 = t235 * t153 - t234;
t229 = t236 * t154 - t233;
t228 = rSges(6,1) + pkin(4);
t227 = Icges(2,4) - Icges(3,5);
t192 = sin(pkin(7));
t193 = cos(pkin(7));
t198 = sin(qJ(3));
t199 = cos(qJ(3));
t109 = -t192 * t198 - t193 * t199;
t110 = -t192 * t199 + t193 * t198;
t226 = -t231 * t109 + t230 * t110;
t225 = t230 * t109 + t231 * t110;
t224 = t232 * t109 + t229 * t110;
t223 = t229 * t109 - t232 * t110;
t222 = Icges(2,1) + Icges(3,1);
t221 = Icges(3,4) + Icges(2,5);
t220 = Icges(2,2) + Icges(3,3);
t219 = Icges(2,6) - Icges(3,6);
t218 = Icges(5,3) + Icges(6,3);
t217 = -t235 * t154 - t233;
t216 = t236 * t153 + t234;
t215 = t231 * t153 - t232 * t154;
t214 = rSges(6,3) + qJ(5);
t213 = t227 * t192;
t212 = t227 * t193;
t211 = rSges(6,2) * t153 - t228 * t154;
t210 = -t220 * t193 - t213;
t209 = t220 * t192 - t212;
t208 = t222 * t192 + t212;
t207 = t222 * t193 - t213;
t104 = -qJD(4) * t109 + V_base(5);
t105 = qJD(4) * t110 + V_base(4);
t149 = V_base(6) - qJD(3);
t206 = (t217 * t153 + t216 * t154) * t149 + (t225 * t153 + t223 * t154) * t105 + (t226 * t153 + t224 * t154) * t104;
t205 = (-t232 * t153 - t231 * t154) * t149 + (t215 * t109 + t218 * t110) * t105 + (-t218 * t109 + t215 * t110) * t104;
t195 = -t109 * t214 + t211 * t110;
t194 = t211 * t109 + t110 * t214;
t191 = Icges(4,4) * t109;
t186 = V_base(5) * qJ(1) + V_base(1);
t182 = qJD(1) + V_base(3);
t181 = t193 * pkin(2);
t180 = t192 * pkin(2);
t179 = rSges(6,2) * t154 + t228 * t153;
t176 = qJD(2) * t192 + t186;
t126 = pkin(1) * t192 - qJ(2) * t193;
t175 = V_base(4) * t126 + t182;
t129 = pkin(1) * t193 + qJ(2) * t192;
t174 = -t129 - t181;
t173 = V_base(4) * t180 + t175;
t172 = -rSges(5,1) * t154 + rSges(5,2) * t153;
t164 = -qJD(2) * t193 + V_base(6) * t129 + V_base(2);
t161 = V_base(4) * pkin(5) + V_base(6) * t181 + t164;
t160 = (-t180 - t126) * V_base(6) + t176;
t101 = -pkin(3) * t109 + pkin(6) * t110;
t159 = -V_base(4) * qJ(1) + t149 * t101 + t161;
t100 = -pkin(3) * t110 - pkin(6) * t109;
t158 = V_base(4) * t100 + (-t101 + t174) * V_base(5) + t173;
t157 = -V_base(5) * pkin(5) + t160;
t139 = -t153 * rSges(5,1) - rSges(5,2) * t154;
t131 = rSges(2,1) * t193 - rSges(2,2) * t192;
t130 = rSges(3,1) * t193 + rSges(3,3) * t192;
t128 = rSges(2,1) * t192 + rSges(2,2) * t193;
t127 = rSges(3,1) * t192 - rSges(3,3) * t193;
t113 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t112 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t111 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t106 = Icges(4,4) * t110;
t103 = V_base(5) * rSges(2,3) - t128 * V_base(6) + t186;
t102 = t131 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t99 = -rSges(4,1) * t109 - rSges(4,2) * t110;
t98 = -rSges(4,1) * t110 + rSges(4,2) * t109;
t97 = t128 * V_base(4) - t131 * V_base(5) + t182;
t96 = -Icges(4,1) * t109 - t106;
t95 = -Icges(4,1) * t110 + t191;
t94 = -Icges(4,2) * t110 - t191;
t93 = Icges(4,2) * t109 - t106;
t88 = V_base(5) * rSges(3,2) + (-t126 - t127) * V_base(6) + t176;
t87 = V_base(6) * t130 + (-qJ(1) - rSges(3,2)) * V_base(4) + t164;
t86 = t110 * rSges(5,3) + t109 * t172;
t84 = -t109 * rSges(5,3) + t110 * t172;
t70 = t127 * V_base(4) + (-t129 - t130) * V_base(5) + t175;
t67 = -t149 * t98 + (-pkin(5) - rSges(4,3)) * V_base(5) + t160;
t66 = t149 * t99 + (rSges(4,3) - qJ(1)) * V_base(4) + t161;
t65 = V_base(4) * t98 + (-t99 + t174) * V_base(5) + t173;
t64 = t104 * t139 + (-t100 - t84) * t149 + t157;
t63 = -t105 * t139 + t149 * t86 + t159;
t62 = -t104 * t86 + t105 * t84 + t158;
t61 = qJD(5) * t110 - t179 * t104 + (-t100 - t195) * t149 + t157;
t60 = -qJD(5) * t109 + t105 * t179 + t149 * t194 + t159;
t59 = -t104 * t194 + t105 * t195 + t158;
t1 = m(1) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(2) * (t102 ^ 2 + t103 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t70 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + (-t205 * t109 + t206 * t110) * t104 / 0.2e1 + (t206 * t109 + t205 * t110) * t105 / 0.2e1 + ((t223 * t153 - t225 * t154) * t105 + (t224 * t153 - t226 * t154) * t104 + (t216 * t153 - t217 * t154 + Icges(4,3)) * t149) * t149 / 0.2e1 + ((-t109 * t95 - t110 * t93 + t210 * t192 + t208 * t193 + Icges(1,4)) * V_base(5) + (-t109 * t96 - t110 * t94 + t209 * t192 + t207 * t193 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t109 * t93 - t110 * t95 + t208 * t192 - t210 * t193 + Icges(1,2)) * V_base(5) + (t109 * t94 - t110 * t96 + t207 * t192 - t209 * t193 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t149 * (Icges(4,5) * t109 + Icges(4,6) * t110) + V_base(5) * t149 * (Icges(4,5) * t110 - Icges(4,6) * t109) + (Icges(1,6) * V_base(5) + Icges(1,5) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) + (t219 * V_base(5) + t221 * V_base(4)) * t193 + (-t219 * V_base(4) + t221 * V_base(5)) * t192) * V_base(6);
T = t1;
