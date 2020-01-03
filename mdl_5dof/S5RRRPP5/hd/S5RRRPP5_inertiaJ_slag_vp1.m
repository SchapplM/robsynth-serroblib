% Calculate joint inertia matrix for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:36
% EndTime: 2019-12-31 20:57:41
% DurationCPUTime: 2.23s
% Computational Cost: add. (2156->207), mult. (2668->306), div. (0->0), fcn. (2426->6), ass. (0->107)
t215 = Icges(4,4) - Icges(6,4) - Icges(5,5);
t214 = Icges(4,1) + Icges(5,1) + Icges(6,1);
t213 = -Icges(4,2) - Icges(6,2) - Icges(5,3);
t122 = qJ(2) + qJ(3);
t115 = cos(t122);
t212 = t215 * t115;
t114 = sin(t122);
t211 = t215 * t114;
t210 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t209 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t208 = t213 * t114 + t212;
t207 = t214 * t115 - t211;
t124 = sin(qJ(1));
t126 = cos(qJ(1));
t206 = t208 * t124 - t209 * t126;
t205 = t209 * t124 + t208 * t126;
t204 = t207 * t124 - t210 * t126;
t203 = t210 * t124 + t207 * t126;
t202 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t201 = t209 * t114 - t210 * t115;
t194 = rSges(6,1) + pkin(4);
t200 = t213 * t115 - t211;
t199 = t214 * t114 + t212;
t198 = t201 * t124 + t202 * t126;
t197 = t202 * t124 - t201 * t126;
t196 = t206 * t114 - t204 * t115;
t195 = t205 * t114 - t203 * t115;
t120 = t124 ^ 2;
t193 = t124 * pkin(6);
t168 = t115 * t126;
t169 = t114 * t126;
t192 = rSges(6,2) * t169 + t194 * t168;
t191 = -t210 * t114 - t209 * t115;
t190 = t200 * t114 + t199 * t115;
t121 = t126 ^ 2;
t189 = t198 * t121 + (t195 * t124 + (-t196 + t197) * t126) * t124;
t188 = m(5) / 0.2e1;
t187 = m(6) / 0.2e1;
t186 = t124 / 0.2e1;
t185 = -t126 / 0.2e1;
t123 = sin(qJ(2));
t184 = pkin(2) * t123;
t125 = cos(qJ(2));
t112 = t125 * pkin(2) + pkin(1);
t103 = t126 * t112;
t119 = t126 * pkin(6);
t183 = t124 * (t119 + (-pkin(1) + t112) * t124) + t126 * (-t126 * pkin(1) + t103 - t193);
t131 = rSges(4,1) * t168 - rSges(4,2) * t169 + t124 * rSges(4,3);
t155 = rSges(4,1) * t115 - rSges(4,2) * t114;
t28 = t124 * (-t126 * rSges(4,3) + t155 * t124) + t126 * t131;
t166 = pkin(3) * t168 + qJ(4) * t169;
t182 = t120 * (pkin(3) * t115 + qJ(4) * t114) + t126 * t166;
t91 = t114 * pkin(3) - t115 * qJ(4);
t181 = -t114 * rSges(5,1) + t115 * rSges(5,3) - t91;
t180 = rSges(3,1) * t125;
t179 = rSges(3,2) * t123;
t178 = t126 * rSges(3,3);
t177 = Icges(3,4) * t123;
t176 = Icges(3,4) * t125;
t165 = t124 * rSges(3,3) + t126 * t180;
t164 = t120 + t121;
t163 = (t197 * t120 + ((-t195 + t198) * t124 + t196 * t126) * t126) * t124;
t127 = -pkin(7) - pkin(6);
t162 = -rSges(6,3) - qJ(5) - t127;
t161 = rSges(5,1) * t168 + t124 * rSges(5,2) + rSges(5,3) * t169;
t94 = t114 * rSges(4,1) + t115 * rSges(4,2);
t160 = -t94 - t184;
t11 = t124 * (-t126 * rSges(5,2) + (rSges(5,1) * t115 + rSges(5,3) * t114) * t124) + t126 * t161 + t182;
t159 = -t124 * t127 + t103;
t158 = t181 - t184;
t157 = t115 * rSges(6,2) - t194 * t114 - t91;
t156 = -t179 + t180;
t9 = t182 + t192 * t126 + (rSges(6,2) * t114 + t194 * t115) * t120;
t143 = Icges(3,1) * t125 - t177;
t139 = -Icges(3,2) * t123 + t176;
t135 = Icges(3,5) * t125 - Icges(3,6) * t123;
t130 = t157 - t184;
t129 = t189 * t126 + t163;
t128 = (t203 * t114 + t205 * t115 - t191 * t124 + t190 * t126) * t186 + (t204 * t114 + t206 * t115 + t190 * t124 + t191 * t126) * t185;
t101 = t126 * rSges(2,1) - t124 * rSges(2,2);
t100 = -t124 * rSges(2,1) - t126 * rSges(2,2);
t99 = t123 * rSges(3,1) + t125 * rSges(3,2);
t74 = Icges(3,3) * t124 + t135 * t126;
t73 = -Icges(3,3) * t126 + t135 * t124;
t51 = t160 * t126;
t50 = t160 * t124;
t41 = t193 + (pkin(1) - t179) * t126 + t165;
t40 = t178 + t119 + (-pkin(1) - t156) * t124;
t39 = t181 * t126;
t38 = t181 * t124;
t37 = t131 + t159;
t36 = (rSges(4,3) - t127) * t126 + (-t112 - t155) * t124;
t35 = t158 * t126;
t34 = t158 * t124;
t33 = t157 * t126;
t32 = t157 * t124;
t31 = t126 * (-t126 * t179 + t165) + (t156 * t124 - t178) * t124;
t30 = t130 * t126;
t29 = t130 * t124;
t27 = t159 + t161 + t166;
t26 = (rSges(5,2) - t127) * t126 + (-t112 + (-rSges(5,1) - pkin(3)) * t115 + (-rSges(5,3) - qJ(4)) * t114) * t124;
t13 = t162 * t124 + t103 + t166 + t192;
t12 = t162 * t126 + (-t112 + (-rSges(6,2) - qJ(4)) * t114 + (-pkin(3) - t194) * t115) * t124;
t10 = t28 + t183;
t8 = t11 + t183;
t7 = t9 + t183;
t1 = [t123 * (Icges(3,1) * t123 + t176) + t125 * (Icges(3,2) * t125 + t177) + Icges(2,3) - t200 * t115 + t199 * t114 + m(6) * (t12 ^ 2 + t13 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2) + m(3) * (t40 ^ 2 + t41 ^ 2) + m(2) * (t100 ^ 2 + t101 ^ 2); (t120 / 0.2e1 + t121 / 0.2e1) * (Icges(3,5) * t123 + Icges(3,6) * t125) + m(3) * (-t124 * t41 - t126 * t40) * t99 + m(6) * (t30 * t12 + t29 * t13) + m(5) * (t35 * t26 + t34 * t27) + m(4) * (t51 * t36 + t50 * t37) + t128 + (t123 * (-Icges(3,5) * t126 + t143 * t124) + t125 * (-Icges(3,6) * t126 + t139 * t124)) * t185 + (t123 * (Icges(3,5) * t124 + t143 * t126) + t125 * (Icges(3,6) * t124 + t139 * t126)) * t186; m(6) * (t29 ^ 2 + t30 ^ 2 + t7 ^ 2) + m(5) * (t34 ^ 2 + t35 ^ 2 + t8 ^ 2) + m(4) * (t10 ^ 2 + t50 ^ 2 + t51 ^ 2) + t124 * t120 * t74 + m(3) * (t164 * t99 ^ 2 + t31 ^ 2) + t163 + (-t121 * t73 + (-t124 * t73 + t126 * t74) * t124 + t189) * t126; m(6) * (t33 * t12 + t32 * t13) + m(5) * (t39 * t26 + t38 * t27) + m(4) * (-t124 * t37 - t126 * t36) * t94 + t128; m(6) * (t32 * t29 + t33 * t30 + t9 * t7) + m(5) * (t11 * t8 + t38 * t34 + t39 * t35) + m(4) * (t28 * t10 + (-t124 * t50 - t126 * t51) * t94) + t129; m(6) * (t32 ^ 2 + t33 ^ 2 + t9 ^ 2) + m(5) * (t11 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(4) * (t164 * t94 ^ 2 + t28 ^ 2) + t129; 0.2e1 * ((t12 * t126 + t124 * t13) * t187 + (t124 * t27 + t126 * t26) * t188) * t114; m(6) * (-t115 * t7 + (t124 * t29 + t126 * t30) * t114) + m(5) * (-t115 * t8 + (t124 * t34 + t126 * t35) * t114); m(6) * (-t115 * t9 + (t124 * t32 + t126 * t33) * t114) + m(5) * (-t115 * t11 + (t124 * t38 + t126 * t39) * t114); 0.2e1 * (t188 + t187) * (t164 * t114 ^ 2 + t115 ^ 2); m(6) * (-t124 * t12 + t126 * t13); m(6) * (-t124 * t30 + t126 * t29); m(6) * (-t124 * t33 + t126 * t32); 0; m(6) * t164;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
