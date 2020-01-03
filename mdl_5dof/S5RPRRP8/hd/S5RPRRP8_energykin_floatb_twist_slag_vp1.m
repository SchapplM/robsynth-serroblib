% Calculate kinetic energy for
% S5RPRRP8
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:04
% EndTime: 2019-12-31 18:47:06
% DurationCPUTime: 2.05s
% Computational Cost: add. (748->190), mult. (1267->240), div. (0->0), fcn. (1303->6), ass. (0->101)
t238 = Icges(5,4) - Icges(6,5);
t237 = Icges(5,1) + Icges(6,1);
t236 = Icges(5,2) + Icges(6,3);
t156 = cos(qJ(4));
t235 = t238 * t156;
t155 = sin(qJ(4));
t234 = t238 * t155;
t233 = Icges(6,4) + Icges(5,5);
t232 = Icges(5,6) - Icges(6,6);
t231 = -t236 * t155 + t235;
t230 = t237 * t156 - t234;
t229 = -rSges(6,1) - pkin(4);
t228 = rSges(6,3) + qJ(5);
t227 = Icges(2,4) - Icges(3,5);
t196 = sin(qJ(3));
t197 = sin(qJ(1));
t198 = cos(qJ(3));
t199 = cos(qJ(1));
t111 = -t197 * t196 - t199 * t198;
t112 = t199 * t196 - t197 * t198;
t226 = t232 * t111 + t231 * t112;
t225 = t231 * t111 - t232 * t112;
t224 = t233 * t111 + t230 * t112;
t223 = t230 * t111 - t233 * t112;
t222 = Icges(2,1) + Icges(3,1);
t221 = Icges(3,4) + Icges(2,5);
t220 = Icges(2,2) + Icges(3,3);
t219 = Icges(6,2) + Icges(5,3);
t218 = Icges(2,6) - Icges(3,6);
t217 = t236 * t156 + t234;
t216 = t237 * t155 + t235;
t215 = t232 * t155 - t233 * t156;
t214 = t227 * t197;
t213 = t227 * t199;
t212 = -t228 * t155 + t229 * t156;
t211 = -t220 * t199 - t214;
t210 = t220 * t197 - t213;
t207 = t222 * t197 + t213;
t206 = t222 * t199 - t214;
t106 = -qJD(4) * t111 + V_base(5);
t107 = qJD(4) * t112 + V_base(4);
t149 = V_base(6) + qJD(1);
t145 = -qJD(3) + t149;
t205 = (-t155 * t217 + t156 * t216) * t145 + (-t225 * t155 + t223 * t156) * t107 + (-t226 * t155 + t224 * t156) * t106;
t204 = (-t233 * t155 - t232 * t156) * t145 + (t215 * t111 + t112 * t219) * t107 + (-t111 * t219 + t215 * t112) * t106;
t195 = -t111 * rSges(6,2) + t212 * t112;
t194 = t112 * rSges(6,2) + t212 * t111;
t193 = Icges(4,4) * t111;
t188 = t229 * t155 + t228 * t156;
t187 = qJD(5) * t155;
t137 = t197 * pkin(1) - t199 * qJ(2);
t186 = V_base(4) * t137 + V_base(3);
t185 = V_base(5) * pkin(5) + V_base(1);
t182 = t199 * pkin(2);
t181 = t197 * pkin(2);
t178 = V_base(4) * t181 + t186;
t177 = qJD(2) * t197 + t185;
t140 = t199 * pkin(1) + t197 * qJ(2);
t176 = -t140 - t182;
t175 = -rSges(5,1) * t156 + rSges(5,2) * t155;
t166 = -qJD(2) * t199 + t149 * t140 + V_base(2);
t163 = V_base(4) * pkin(6) + t149 * t182 + t166;
t162 = (-t181 - t137) * t149 + t177;
t102 = -pkin(3) * t111 + pkin(7) * t112;
t161 = -V_base(4) * pkin(5) + t145 * t102 + t163;
t101 = -pkin(3) * t112 - pkin(7) * t111;
t160 = V_base(4) * t101 + (-t102 + t176) * V_base(5) + t178;
t159 = -V_base(5) * pkin(6) + t162;
t142 = t199 * rSges(2,1) - t197 * rSges(2,2);
t141 = t199 * rSges(3,1) + t197 * rSges(3,3);
t139 = t197 * rSges(2,1) + t199 * rSges(2,2);
t138 = t197 * rSges(3,1) - t199 * rSges(3,3);
t136 = -t155 * rSges(5,1) - rSges(5,2) * t156;
t115 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t114 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t113 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t109 = Icges(4,4) * t112;
t105 = V_base(5) * rSges(2,3) - t139 * t149 + t185;
t104 = t142 * t149 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t103 = t139 * V_base(4) - t142 * V_base(5) + V_base(3);
t100 = -rSges(4,1) * t111 - rSges(4,2) * t112;
t99 = -rSges(4,1) * t112 + rSges(4,2) * t111;
t98 = -Icges(4,1) * t111 - t109;
t97 = -Icges(4,1) * t112 + t193;
t96 = -Icges(4,2) * t112 - t193;
t95 = Icges(4,2) * t111 - t109;
t88 = V_base(5) * rSges(3,2) + (-t137 - t138) * t149 + t177;
t87 = t149 * t141 + (-pkin(5) - rSges(3,2)) * V_base(4) + t166;
t86 = t112 * rSges(5,3) + t111 * t175;
t84 = -t111 * rSges(5,3) + t112 * t175;
t70 = t138 * V_base(4) + (-t140 - t141) * V_base(5) + t186;
t69 = -t145 * t99 + (-pkin(6) - rSges(4,3)) * V_base(5) + t162;
t68 = t145 * t100 + (rSges(4,3) - pkin(5)) * V_base(4) + t163;
t67 = V_base(4) * t99 + (-t100 + t176) * V_base(5) + t178;
t66 = t106 * t136 + (-t101 - t84) * t145 + t159;
t65 = -t107 * t136 + t145 * t86 + t161;
t64 = -t106 * t86 + t107 * t84 + t160;
t63 = -t111 * t187 + t188 * t106 + (-t101 - t195) * t145 + t159;
t62 = -t188 * t107 - t112 * t187 + t194 * t145 + t161;
t61 = qJD(5) * t156 - t194 * t106 + t195 * t107 + t160;
t1 = m(1) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(2) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(3) * (t70 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (-t204 * t111 + t205 * t112) * t106 / 0.2e1 + (t205 * t111 + t204 * t112) * t107 / 0.2e1 + ((t223 * t155 + t225 * t156) * t107 + (t224 * t155 + t226 * t156) * t106 + (t155 * t216 + t156 * t217 + Icges(4,3)) * t145) * t145 / 0.2e1 + ((-t111 * t97 - t112 * t95 + t211 * t197 + t207 * t199 + Icges(1,4)) * V_base(5) + (-t111 * t98 - t112 * t96 + t210 * t197 + t206 * t199 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t111 * t95 - t112 * t97 + t207 * t197 - t211 * t199 + Icges(1,2)) * V_base(5) + (t111 * t96 - t112 * t98 + t206 * t197 - t210 * t199 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t145 * (Icges(4,5) * t111 + Icges(4,6) * t112) + V_base(5) * t145 * (Icges(4,5) * t112 - Icges(4,6) * t111) + ((Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * t149 + (t218 * V_base(5) + t221 * V_base(4)) * t199 + (-t218 * V_base(4) + t221 * V_base(5)) * t197) * t149 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
