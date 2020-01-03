% Calculate joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:53
% EndTime: 2019-12-31 20:54:59
% DurationCPUTime: 2.34s
% Computational Cost: add. (2767->227), mult. (2608->324), div. (0->0), fcn. (2371->8), ass. (0->126)
t228 = Icges(5,4) - Icges(6,5);
t227 = Icges(5,1) + Icges(6,1);
t226 = Icges(5,2) + Icges(6,3);
t128 = qJ(2) + qJ(3);
t116 = pkin(8) + t128;
t114 = cos(t116);
t225 = t228 * t114;
t113 = sin(t116);
t224 = t228 * t113;
t223 = Icges(6,4) + Icges(5,5);
t222 = Icges(5,6) - Icges(6,6);
t221 = -t226 * t113 + t225;
t220 = t227 * t114 - t224;
t130 = sin(qJ(1));
t132 = cos(qJ(1));
t219 = t221 * t130 - t222 * t132;
t218 = t222 * t130 + t221 * t132;
t217 = t220 * t130 - t223 * t132;
t216 = t223 * t130 + t220 * t132;
t215 = Icges(6,2) + Icges(4,3) + Icges(5,3);
t117 = sin(t128);
t118 = cos(t128);
t214 = Icges(4,5) * t118 - Icges(4,6) * t117 - t222 * t113 + t223 * t114;
t210 = rSges(6,1) + pkin(4);
t213 = t226 * t114 + t224;
t212 = t227 * t113 + t225;
t209 = rSges(6,3) + qJ(5);
t211 = t209 * t113 + t210 * t114;
t208 = -t214 * t130 + t215 * t132;
t207 = t215 * t130 + t214 * t132;
t182 = Icges(4,4) * t118;
t146 = -Icges(4,2) * t117 + t182;
t70 = -Icges(4,6) * t132 + t146 * t130;
t183 = Icges(4,4) * t117;
t150 = Icges(4,1) * t118 - t183;
t72 = -Icges(4,5) * t132 + t150 * t130;
t206 = t219 * t113 - t217 * t114 + t117 * t70 - t118 * t72;
t71 = Icges(4,6) * t130 + t146 * t132;
t73 = Icges(4,5) * t130 + t150 * t132;
t205 = t218 * t113 - t216 * t114 + t117 * t71 - t118 * t73;
t126 = t130 ^ 2;
t204 = t130 * pkin(6);
t203 = t130 * rSges(6,2);
t176 = t114 * t132;
t177 = t113 * t132;
t202 = t210 * t176 + t209 * t177 + t203;
t164 = rSges(4,1) * t118 - rSges(4,2) * t117;
t201 = Icges(4,5) * t117 + Icges(4,6) * t118 + t223 * t113 + t222 * t114;
t92 = Icges(4,2) * t118 + t183;
t93 = Icges(4,1) * t117 + t182;
t200 = -t213 * t113 + t212 * t114 - t117 * t92 + t118 * t93;
t199 = -t210 * t113 + t209 * t114;
t127 = t132 ^ 2;
t198 = t208 * t127 + (t205 * t130 + (-t206 + t207) * t132) * t130;
t133 = -pkin(7) - pkin(6);
t197 = t130 / 0.2e1;
t196 = -t132 / 0.2e1;
t129 = sin(qJ(2));
t195 = pkin(2) * t129;
t194 = pkin(3) * t117;
t131 = cos(qJ(2));
t115 = t131 * pkin(2) + pkin(1);
t108 = t132 * t115;
t95 = pkin(3) * t118 + t115;
t89 = t132 * t95;
t193 = t132 * (-t108 + t89) + (-t115 + t95) * t126;
t124 = t132 * pkin(6);
t192 = t130 * (t124 + (-pkin(1) + t115) * t130) + t132 * (-t132 * pkin(1) + t108 - t204);
t138 = t130 * rSges(4,3) + t164 * t132;
t28 = t130 * (-t132 * rSges(4,3) + t164 * t130) + t132 * t138;
t191 = rSges(3,1) * t131;
t189 = rSges(3,2) * t129;
t187 = t132 * rSges(3,3);
t185 = Icges(3,4) * t129;
t184 = Icges(3,4) * t131;
t175 = t130 * rSges(3,3) + t132 * t191;
t173 = t126 + t127;
t172 = (t207 * t126 + ((-t205 + t208) * t130 + t206 * t132) * t132) * t130;
t94 = t117 * rSges(4,1) + t118 * rSges(4,2);
t170 = -t94 - t195;
t88 = t113 * rSges(5,1) + t114 * rSges(5,2);
t169 = -t88 - t194;
t139 = rSges(5,1) * t176 - rSges(5,2) * t177 + t130 * rSges(5,3);
t163 = rSges(5,1) * t114 - rSges(5,2) * t113;
t10 = t130 * (-t132 * rSges(5,3) + t163 * t130) + t132 * t139 + t193;
t125 = -qJ(4) + t133;
t168 = -t130 * t125 + t89;
t167 = -t194 + t199;
t166 = -t194 - t195;
t165 = -t189 + t191;
t7 = t193 + (t202 - t203) * t132 + t211 * t126;
t151 = Icges(3,1) * t131 - t185;
t147 = -Icges(3,2) * t129 + t184;
t143 = Icges(3,5) * t131 - Icges(3,6) * t129;
t137 = t166 - t88;
t136 = t166 + t199;
t135 = t198 * t132 + t172;
t134 = (t216 * t113 + t218 * t114 + t117 * t73 + t118 * t71 + t201 * t130 + t200 * t132) * t197 + (t217 * t113 + t219 * t114 + t117 * t72 + t118 * t70 + t200 * t130 - t201 * t132) * t196;
t107 = t132 * rSges(2,1) - t130 * rSges(2,2);
t106 = -t130 * rSges(2,1) - t132 * rSges(2,2);
t105 = t129 * rSges(3,1) + t131 * rSges(3,2);
t75 = Icges(3,3) * t130 + t143 * t132;
t74 = -Icges(3,3) * t132 + t143 * t130;
t67 = t170 * t132;
t66 = t170 * t130;
t47 = t204 + (pkin(1) - t189) * t132 + t175;
t46 = t187 + t124 + (-pkin(1) - t165) * t130;
t45 = t169 * t132;
t44 = t169 * t130;
t39 = t137 * t132;
t38 = t137 * t130;
t37 = -t130 * t133 + t108 + t138;
t36 = (rSges(4,3) - t133) * t132 + (-t115 - t164) * t130;
t33 = t132 * (-t132 * t189 + t175) + (t165 * t130 - t187) * t130;
t32 = t139 + t168;
t31 = (rSges(5,3) - t125) * t132 + (-t163 - t95) * t130;
t30 = t167 * t132;
t29 = t167 * t130;
t27 = t136 * t132;
t26 = t136 * t130;
t21 = t168 + t202;
t20 = (rSges(6,2) - t125) * t132 + (-t211 - t95) * t130;
t11 = t28 + t192;
t6 = t10 + t192;
t1 = t7 + t192;
t2 = [t131 * (Icges(3,2) * t131 + t185) + t129 * (Icges(3,1) * t129 + t184) + t117 * t93 + t118 * t92 + Icges(2,3) + t213 * t114 + t212 * t113 + m(5) * (t31 ^ 2 + t32 ^ 2) + m(6) * (t20 ^ 2 + t21 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t46 ^ 2 + t47 ^ 2) + m(2) * (t106 ^ 2 + t107 ^ 2); (t127 / 0.2e1 + t126 / 0.2e1) * (Icges(3,5) * t129 + Icges(3,6) * t131) + m(3) * (-t130 * t47 - t132 * t46) * t105 + m(5) * (t39 * t31 + t38 * t32) + m(6) * (t27 * t20 + t26 * t21) + m(4) * (t67 * t36 + t66 * t37) + t134 + (t129 * (-Icges(3,5) * t132 + t151 * t130) + t131 * (-Icges(3,6) * t132 + t147 * t130)) * t196 + (t129 * (Icges(3,5) * t130 + t151 * t132) + t131 * (Icges(3,6) * t130 + t147 * t132)) * t197; m(6) * (t1 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(5) * (t38 ^ 2 + t39 ^ 2 + t6 ^ 2) + m(4) * (t11 ^ 2 + t66 ^ 2 + t67 ^ 2) + t130 * t126 * t75 + m(3) * (t173 * t105 ^ 2 + t33 ^ 2) + t172 + (-t127 * t74 + (-t130 * t74 + t132 * t75) * t130 + t198) * t132; m(5) * (t45 * t31 + t44 * t32) + m(6) * (t30 * t20 + t29 * t21) + m(4) * (-t130 * t37 - t132 * t36) * t94 + t134; m(6) * (t7 * t1 + t29 * t26 + t30 * t27) + m(5) * (t10 * t6 + t44 * t38 + t45 * t39) + m(4) * (t28 * t11 + (-t130 * t66 - t132 * t67) * t94) + t135; m(5) * (t10 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t29 ^ 2 + t30 ^ 2 + t7 ^ 2) + m(4) * (t173 * t94 ^ 2 + t28 ^ 2) + t135; m(5) * (t130 * t31 - t132 * t32) + m(6) * (t130 * t20 - t132 * t21); m(6) * (t130 * t27 - t132 * t26) + m(5) * (t130 * t39 - t132 * t38); m(5) * (t130 * t45 - t132 * t44) + m(6) * (t130 * t30 - t132 * t29); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t173; m(6) * (t130 * t21 + t132 * t20) * t113; m(6) * (-t114 * t1 + (t130 * t26 + t132 * t27) * t113); m(6) * (-t114 * t7 + (t130 * t29 + t132 * t30) * t113); 0; m(6) * (t173 * t113 ^ 2 + t114 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
