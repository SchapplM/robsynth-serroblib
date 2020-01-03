% Calculate kinetic energy for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR10_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:55
% EndTime: 2019-12-31 18:03:58
% DurationCPUTime: 3.00s
% Computational Cost: add. (940->263), mult. (1531->360), div. (0->0), fcn. (1469->8), ass. (0->134)
t260 = Icges(3,4) - Icges(4,5);
t259 = Icges(3,1) + Icges(4,1);
t258 = Icges(3,2) + Icges(4,3);
t191 = cos(pkin(8));
t257 = t260 * t191;
t190 = sin(pkin(8));
t256 = t260 * t190;
t255 = Icges(4,4) + Icges(3,5);
t254 = Icges(3,6) - Icges(4,6);
t253 = t258 * t190 - t257;
t252 = t259 * t191 - t256;
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t251 = t253 * t193 + t254 * t195;
t250 = -t254 * t193 + t253 * t195;
t249 = t252 * t193 - t255 * t195;
t248 = t255 * t193 + t252 * t195;
t247 = -t258 * t191 - t256;
t246 = t259 * t190 + t257;
t245 = t254 * t190 - t255 * t191;
t244 = Icges(2,2) + Icges(4,2) + Icges(3,3);
t182 = V_base(6) + qJD(1);
t187 = Icges(2,4) * t195;
t235 = Icges(2,4) * t193;
t243 = (-t255 * t190 - t254 * t191) * t182 + (t193 * t245 + t195 * t244 + t235) * V_base(5) + (-t193 * t244 + t195 * t245 + t187) * V_base(4);
t242 = (t190 * t247 + t191 * t246) * t182 + (Icges(2,1) * t193 + t190 * t251 + t249 * t191 + t187) * V_base(5) + (Icges(2,1) * t195 + t190 * t250 + t191 * t248 - t235) * V_base(4);
t240 = pkin(3) * t190;
t239 = pkin(3) * t191;
t194 = cos(qJ(4));
t238 = pkin(4) * t194;
t163 = pkin(2) * t190 - qJ(3) * t191;
t237 = -pkin(5) - t163;
t192 = sin(qJ(4));
t230 = t190 * t192;
t229 = t191 * t192;
t217 = pkin(2) * t191 + qJ(3) * t190;
t142 = t217 * t193;
t172 = pkin(1) * t193 - qJ(2) * t195;
t228 = -t142 - t172;
t143 = t217 * t195;
t174 = pkin(1) * t195 + qJ(2) * t193;
t227 = -t143 - t174;
t226 = qJD(3) * t190;
t225 = V_base(4) * t172 + V_base(3);
t224 = V_base(5) * pkin(5) + V_base(1);
t152 = t195 * pkin(6) + t193 * t239;
t221 = -t152 + t228;
t176 = qJD(4) * t195 + V_base(5);
t220 = qJD(2) * t193 + t224;
t219 = rSges(3,1) * t191 - rSges(3,2) * t190;
t218 = rSges(4,1) * t191 + rSges(4,3) * t190;
t189 = qJ(4) + qJ(5);
t185 = sin(t189);
t186 = cos(t189);
t141 = -t185 * t191 + t186 * t190;
t210 = t185 * t190 + t186 * t191;
t148 = t190 * t194 - t229;
t209 = t191 * t194 + t230;
t208 = -qJD(2) * t195 + t182 * t174 + V_base(2);
t207 = V_base(5) * t163 + t195 * t226 + t220;
t206 = -qJD(3) * t191 + V_base(4) * t142 + t225;
t205 = V_base(5) * t240 + t207;
t204 = pkin(4) * t230 + t191 * t238;
t203 = t182 * t143 + t193 * t226 + t208;
t153 = -t193 * pkin(6) + t195 * t239;
t200 = V_base(4) * t152 + (-t153 + t227) * V_base(5) + t206;
t199 = t182 * t153 + (t237 - t240) * V_base(4) + t203;
t177 = -qJD(4) * t193 + V_base(4);
t175 = rSges(2,1) * t195 - rSges(2,2) * t193;
t173 = rSges(2,1) * t193 + rSges(2,2) * t195;
t167 = Icges(2,5) * t195 - Icges(2,6) * t193;
t166 = Icges(2,5) * t193 + Icges(2,6) * t195;
t165 = rSges(3,1) * t190 + rSges(3,2) * t191;
t164 = rSges(4,1) * t190 - rSges(4,3) * t191;
t156 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t155 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t154 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t150 = V_base(4) + (-qJD(4) - qJD(5)) * t193;
t149 = qJD(5) * t195 + t176;
t138 = t209 * t195;
t137 = t148 * t195;
t136 = t209 * t193;
t135 = t148 * t193;
t133 = rSges(3,3) * t193 + t195 * t219;
t132 = rSges(4,2) * t193 + t195 * t218;
t131 = -rSges(3,3) * t195 + t193 * t219;
t130 = -rSges(4,2) * t195 + t193 * t218;
t116 = -pkin(4) * t229 + t190 * t238;
t115 = t210 * t195;
t114 = t141 * t195;
t113 = t210 * t193;
t112 = t141 * t193;
t111 = V_base(5) * rSges(2,3) - t173 * t182 + t224;
t110 = t175 * t182 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t109 = t173 * V_base(4) - t175 * V_base(5) + V_base(3);
t108 = rSges(5,1) * t148 - rSges(5,2) * t209;
t107 = Icges(5,1) * t148 - Icges(5,4) * t209;
t106 = Icges(5,4) * t148 - Icges(5,2) * t209;
t105 = Icges(5,5) * t148 - Icges(5,6) * t209;
t104 = rSges(6,1) * t141 - rSges(6,2) * t210;
t103 = Icges(6,1) * t141 - Icges(6,4) * t210;
t102 = Icges(6,4) * t141 - Icges(6,2) * t210;
t101 = Icges(6,5) * t141 - Icges(6,6) * t210;
t100 = -pkin(7) * t193 + t195 * t204;
t99 = pkin(7) * t195 + t193 * t204;
t98 = rSges(5,1) * t138 + rSges(5,2) * t137 - rSges(5,3) * t193;
t97 = rSges(5,1) * t136 + rSges(5,2) * t135 + rSges(5,3) * t195;
t96 = Icges(5,1) * t138 + Icges(5,4) * t137 - Icges(5,5) * t193;
t95 = Icges(5,1) * t136 + Icges(5,4) * t135 + Icges(5,5) * t195;
t94 = Icges(5,4) * t138 + Icges(5,2) * t137 - Icges(5,6) * t193;
t93 = Icges(5,4) * t136 + Icges(5,2) * t135 + Icges(5,6) * t195;
t92 = Icges(5,5) * t138 + Icges(5,6) * t137 - Icges(5,3) * t193;
t91 = Icges(5,5) * t136 + Icges(5,6) * t135 + Icges(5,3) * t195;
t90 = rSges(6,1) * t115 + rSges(6,2) * t114 - rSges(6,3) * t193;
t89 = rSges(6,1) * t113 + rSges(6,2) * t112 + rSges(6,3) * t195;
t88 = Icges(6,1) * t115 + Icges(6,4) * t114 - Icges(6,5) * t193;
t87 = Icges(6,1) * t113 + Icges(6,4) * t112 + Icges(6,5) * t195;
t86 = Icges(6,4) * t115 + Icges(6,2) * t114 - Icges(6,6) * t193;
t85 = Icges(6,4) * t113 + Icges(6,2) * t112 + Icges(6,6) * t195;
t84 = Icges(6,5) * t115 + Icges(6,6) * t114 - Icges(6,3) * t193;
t83 = Icges(6,5) * t113 + Icges(6,6) * t112 + Icges(6,3) * t195;
t82 = t165 * V_base(5) + (-t131 - t172) * t182 + t220;
t81 = t133 * t182 + (-pkin(5) - t165) * V_base(4) + t208;
t80 = t131 * V_base(4) + (-t133 - t174) * V_base(5) + t225;
t79 = t164 * V_base(5) + (-t130 + t228) * t182 + t207;
t78 = t132 * t182 + (-t164 + t237) * V_base(4) + t203;
t77 = t130 * V_base(4) + (-t132 + t227) * V_base(5) + t206;
t76 = t108 * t176 + (-t97 + t221) * t182 + t205;
t75 = -t108 * t177 + t182 * t98 + t199;
t74 = -t176 * t98 + t177 * t97 + t200;
t73 = t104 * t149 + t116 * t176 + (-t89 - t99 + t221) * t182 + t205;
t72 = -t104 * t150 - t116 * t177 + (t100 + t90) * t182 + t199;
t71 = -t100 * t176 - t149 * t90 + t150 * t89 + t177 * t99 + t200;
t1 = m(1) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(2) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(3) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + t177 * ((t137 * t94 + t138 * t96 - t193 * t92) * t177 + (t137 * t93 + t138 * t95 - t193 * t91) * t176 + (-t105 * t193 + t106 * t137 + t107 * t138) * t182) / 0.2e1 + t176 * ((t135 * t94 + t136 * t96 + t195 * t92) * t177 + (t135 * t93 + t136 * t95 + t195 * t91) * t176 + (t105 * t195 + t106 * t135 + t107 * t136) * t182) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + t150 * ((t114 * t86 + t115 * t88 - t193 * t84) * t150 + (t114 * t85 + t115 * t87 - t193 * t83) * t149 + (-t101 * t193 + t102 * t114 + t103 * t115) * t182) / 0.2e1 + t149 * ((t112 * t86 + t113 * t88 + t195 * t84) * t150 + (t112 * t85 + t113 * t87 + t195 * t83) * t149 + (t101 * t195 + t102 * t112 + t103 * t113) * t182) / 0.2e1 + (Icges(1,1) * V_base(4) + t167 * t182 - t243 * t193 + t242 * t195) * V_base(4) / 0.2e1 + (Icges(1,2) * V_base(5) + t166 * t182 + t242 * t193 + t243 * t195) * V_base(5) / 0.2e1 + ((t148 * t96 - t209 * t94) * t177 + (t148 * t95 - t209 * t93) * t176 + (t141 * t88 - t210 * t86) * t150 + (t141 * t87 - t210 * t85) * t149 + (t249 * t190 - t191 * t251 + t166) * V_base(5) + (t190 * t248 - t191 * t250 + t167) * V_base(4) + (-t210 * t102 + t141 * t103 - t209 * t106 + t148 * t107 + t246 * t190 - t247 * t191 + Icges(2,3)) * t182) * t182 / 0.2e1 + V_base(5) * V_base(4) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
