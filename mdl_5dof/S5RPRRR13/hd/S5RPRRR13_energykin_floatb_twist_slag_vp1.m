% Calculate kinetic energy for
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR13_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:28
% EndTime: 2019-12-31 19:14:31
% DurationCPUTime: 2.81s
% Computational Cost: add. (922->280), mult. (1448->410), div. (0->0), fcn. (1342->8), ass. (0->134)
t248 = Icges(2,4) + Icges(3,6);
t247 = Icges(2,1) + Icges(3,2);
t246 = -Icges(3,4) + Icges(2,5);
t245 = Icges(3,5) - Icges(2,6);
t244 = Icges(2,2) + Icges(3,3);
t191 = cos(qJ(1));
t243 = t248 * t191;
t188 = sin(qJ(1));
t242 = t248 * t188;
t241 = -t244 * t191 - t242;
t240 = t244 * t188 - t243;
t239 = t247 * t188 + t243;
t238 = t247 * t191 - t242;
t187 = sin(qJ(3));
t190 = cos(qJ(3));
t189 = cos(qJ(4));
t230 = pkin(4) * t189;
t235 = -pkin(8) * t190 + t187 * t230;
t227 = Icges(4,4) * t187;
t204 = Icges(4,2) * t190 + t227;
t118 = Icges(4,6) * t191 + t188 * t204;
t119 = Icges(4,6) * t188 - t191 * t204;
t226 = Icges(4,4) * t190;
t205 = Icges(4,1) * t187 + t226;
t121 = Icges(4,5) * t191 + t188 * t205;
t122 = Icges(4,5) * t188 - t191 * t205;
t152 = -Icges(4,2) * t187 + t226;
t157 = Icges(4,1) * t190 - t227;
t170 = qJD(3) * t188 + V_base(5);
t171 = qJD(3) * t191 + V_base(4);
t175 = V_base(6) + qJD(1);
t234 = (t118 * t190 + t121 * t187) * t171 + (t119 * t190 + t122 * t187) * t170 + (t152 * t190 + t157 * t187) * t175;
t232 = pkin(6) * t188;
t231 = pkin(6) * t191;
t186 = sin(qJ(4));
t223 = t186 * t188;
t222 = t186 * t191;
t221 = t187 * t188;
t220 = t187 * t191;
t219 = t188 * t189;
t218 = t188 * t190;
t217 = t189 * t191;
t216 = t190 * t191;
t215 = qJD(4) * t190;
t161 = pkin(1) * t188 - qJ(2) * t191;
t214 = V_base(4) * t161 + V_base(3);
t213 = V_base(5) * pkin(5) + V_base(1);
t210 = -t161 - t232;
t209 = qJD(2) * t188 + t213;
t131 = t191 * t215 + t170;
t160 = qJD(4) * t187 + t175;
t208 = V_base(5) * pkin(2) + t209;
t207 = pkin(3) * t187 - pkin(7) * t190;
t206 = rSges(4,1) * t187 + rSges(4,2) * t190;
t203 = Icges(4,5) * t187 + Icges(4,6) * t190;
t165 = pkin(1) * t191 + qJ(2) * t188;
t199 = -qJD(2) * t191 + t175 * t165 + V_base(2);
t198 = (Icges(4,3) * t191 + t188 * t203) * t171 + (Icges(4,3) * t188 - t191 * t203) * t170 + (Icges(4,5) * t190 - Icges(4,6) * t187) * t175;
t197 = V_base(4) * t232 + (-t165 - t231) * V_base(5) + t214;
t196 = t175 * t231 + (-pkin(2) - pkin(5)) * V_base(4) + t199;
t139 = t207 * t191;
t168 = t190 * pkin(3) + t187 * pkin(7);
t195 = t170 * t168 + (t139 + t210) * t175 + t208;
t138 = t207 * t188;
t194 = -t138 * t170 - t171 * t139 + t197;
t193 = t175 * t138 - t168 * t171 + t196;
t185 = qJ(4) + qJ(5);
t181 = cos(t185);
t180 = sin(t185);
t167 = rSges(2,1) * t191 - rSges(2,2) * t188;
t166 = -rSges(3,2) * t191 + rSges(3,3) * t188;
t164 = rSges(4,1) * t190 - rSges(4,2) * t187;
t163 = rSges(2,1) * t188 + rSges(2,2) * t191;
t162 = -rSges(3,2) * t188 - rSges(3,3) * t191;
t144 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t143 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t142 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t140 = qJD(5) * t187 + t160;
t136 = -t187 * t217 + t223;
t135 = t186 * t220 + t219;
t134 = t187 * t219 + t222;
t133 = -t186 * t221 + t217;
t132 = -t188 * t215 + t171;
t129 = t180 * t188 - t181 * t220;
t128 = t180 * t220 + t181 * t188;
t127 = t180 * t191 + t181 * t221;
t126 = -t180 * t221 + t181 * t191;
t125 = rSges(4,3) * t188 - t191 * t206;
t124 = rSges(5,3) * t187 + (rSges(5,1) * t189 - rSges(5,2) * t186) * t190;
t123 = rSges(4,3) * t191 + t188 * t206;
t120 = Icges(5,5) * t187 + (Icges(5,1) * t189 - Icges(5,4) * t186) * t190;
t117 = Icges(5,6) * t187 + (Icges(5,4) * t189 - Icges(5,2) * t186) * t190;
t114 = Icges(5,3) * t187 + (Icges(5,5) * t189 - Icges(5,6) * t186) * t190;
t112 = (-qJD(4) - qJD(5)) * t218 + t171;
t111 = qJD(5) * t216 + t131;
t110 = rSges(6,3) * t187 + (rSges(6,1) * t181 - rSges(6,2) * t180) * t190;
t108 = Icges(6,5) * t187 + (Icges(6,1) * t181 - Icges(6,4) * t180) * t190;
t107 = Icges(6,6) * t187 + (Icges(6,4) * t181 - Icges(6,2) * t180) * t190;
t106 = Icges(6,3) * t187 + (Icges(6,5) * t181 - Icges(6,6) * t180) * t190;
t105 = pkin(8) * t187 + t190 * t230;
t104 = V_base(5) * rSges(2,3) - t163 * t175 + t213;
t103 = t167 * t175 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t102 = t163 * V_base(4) - t167 * V_base(5) + V_base(3);
t101 = pkin(4) * t223 - t191 * t235;
t100 = pkin(4) * t222 + t188 * t235;
t99 = rSges(5,1) * t136 + rSges(5,2) * t135 + rSges(5,3) * t216;
t98 = rSges(5,1) * t134 + rSges(5,2) * t133 - rSges(5,3) * t218;
t97 = Icges(5,1) * t136 + Icges(5,4) * t135 + Icges(5,5) * t216;
t96 = Icges(5,1) * t134 + Icges(5,4) * t133 - Icges(5,5) * t218;
t95 = Icges(5,4) * t136 + Icges(5,2) * t135 + Icges(5,6) * t216;
t94 = Icges(5,4) * t134 + Icges(5,2) * t133 - Icges(5,6) * t218;
t93 = Icges(5,5) * t136 + Icges(5,6) * t135 + Icges(5,3) * t216;
t92 = Icges(5,5) * t134 + Icges(5,6) * t133 - Icges(5,3) * t218;
t91 = V_base(5) * rSges(3,1) + (-t161 - t162) * t175 + t209;
t90 = t166 * t175 + (-rSges(3,1) - pkin(5)) * V_base(4) + t199;
t89 = rSges(6,1) * t129 + rSges(6,2) * t128 + rSges(6,3) * t216;
t88 = rSges(6,1) * t127 + rSges(6,2) * t126 - rSges(6,3) * t218;
t87 = Icges(6,1) * t129 + Icges(6,4) * t128 + Icges(6,5) * t216;
t86 = Icges(6,1) * t127 + Icges(6,4) * t126 - Icges(6,5) * t218;
t85 = Icges(6,4) * t129 + Icges(6,2) * t128 + Icges(6,6) * t216;
t84 = Icges(6,4) * t127 + Icges(6,2) * t126 - Icges(6,6) * t218;
t83 = Icges(6,5) * t129 + Icges(6,6) * t128 + Icges(6,3) * t216;
t82 = Icges(6,5) * t127 + Icges(6,6) * t126 - Icges(6,3) * t218;
t81 = t162 * V_base(4) + (-t165 - t166) * V_base(5) + t214;
t80 = t164 * t170 + (-t125 + t210) * t175 + t208;
t79 = t123 * t175 - t164 * t171 + t196;
t78 = -t123 * t170 + t125 * t171 + t197;
t77 = t124 * t131 - t160 * t99 + t195;
t76 = -t124 * t132 + t160 * t98 + t193;
t75 = -t131 * t98 + t132 * t99 + t194;
t74 = -t101 * t160 + t105 * t131 + t110 * t111 - t140 * t89 + t195;
t73 = t100 * t160 - t105 * t132 - t110 * t112 + t140 * t88 + t193;
t72 = -t100 * t131 + t101 * t132 - t111 * t88 + t112 * t89 + t194;
t1 = m(1) * (t142 ^ 2 + t143 ^ 2 + t144 ^ 2) / 0.2e1 + m(2) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(3) * (t81 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t171 * (t188 * t234 + t198 * t191) / 0.2e1 + t170 * (t198 * t188 - t191 * t234) / 0.2e1 + m(5) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + t132 * ((t133 * t94 + t134 * t96 - t92 * t218) * t132 + (t133 * t95 + t134 * t97 - t218 * t93) * t131 + (-t114 * t218 + t117 * t133 + t120 * t134) * t160) / 0.2e1 + t131 * ((t135 * t94 + t136 * t96 + t216 * t92) * t132 + (t135 * t95 + t136 * t97 + t93 * t216) * t131 + (t114 * t216 + t117 * t135 + t120 * t136) * t160) / 0.2e1 + t160 * ((t114 * t160 + t131 * t93 + t132 * t92) * t187 + ((-t186 * t94 + t189 * t96) * t132 + (-t186 * t95 + t189 * t97) * t131 + (-t117 * t186 + t120 * t189) * t160) * t190) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + t112 * ((t126 * t84 + t127 * t86 - t82 * t218) * t112 + (t126 * t85 + t127 * t87 - t218 * t83) * t111 + (-t106 * t218 + t107 * t126 + t108 * t127) * t140) / 0.2e1 + t111 * ((t128 * t84 + t129 * t86 + t216 * t82) * t112 + (t128 * t85 + t129 * t87 + t83 * t216) * t111 + (t106 * t216 + t107 * t128 + t108 * t129) * t140) / 0.2e1 + t140 * ((t106 * t140 + t111 * t83 + t112 * t82) * t187 + ((-t180 * t84 + t181 * t86) * t112 + (-t180 * t85 + t181 * t87) * t111 + (-t107 * t180 + t108 * t181) * t140) * t190) / 0.2e1 + ((-t118 * t187 + t121 * t190) * t171 + (-t119 * t187 + t122 * t190) * t170 + (-t187 * t152 + t190 * t157 + Icges(3,1) + Icges(2,3)) * t175) * t175 / 0.2e1 + ((t188 * t241 + t239 * t191 + Icges(1,4)) * V_base(5) + (t240 * t188 + t238 * t191 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t239 * t188 - t241 * t191 + Icges(1,2)) * V_base(5) + (t188 * t238 - t191 * t240 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t175 * (t246 * t188 - t245 * t191) + V_base(4) * t175 * (t245 * t188 + t246 * t191) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
