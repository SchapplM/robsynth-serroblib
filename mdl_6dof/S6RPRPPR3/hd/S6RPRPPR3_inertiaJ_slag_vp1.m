% Calculate joint inertia matrix for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:12
% EndTime: 2019-03-09 02:44:15
% DurationCPUTime: 1.80s
% Computational Cost: add. (2909->269), mult. (3368->391), div. (0->0), fcn. (3378->8), ass. (0->132)
t206 = Icges(4,1) + Icges(5,1);
t205 = Icges(6,1) + Icges(5,3);
t202 = Icges(4,6) + Icges(6,5);
t199 = Icges(4,5) + Icges(5,4) + Icges(6,6);
t124 = cos(qJ(3));
t204 = (Icges(6,4) - Icges(5,5)) * t124;
t121 = sin(qJ(3));
t203 = (Icges(4,4) - Icges(5,5)) * t121;
t201 = t121 * t205 - t204;
t200 = t124 * t206 - t203;
t198 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t197 = -t199 * t124 + (-Icges(5,6) + t202) * t121;
t118 = qJ(1) + pkin(9);
t115 = sin(t118);
t196 = -t115 / 0.2e1;
t183 = t115 / 0.2e1;
t116 = cos(t118);
t195 = -t116 / 0.2e1;
t194 = t116 / 0.2e1;
t193 = t121 / 0.2e1;
t163 = t116 * t124;
t164 = t116 * t121;
t123 = cos(qJ(6));
t120 = sin(qJ(6));
t162 = t120 * t121;
t72 = -t115 * t123 - t116 * t162;
t161 = t121 * t123;
t73 = -t115 * t120 + t116 * t161;
t36 = t73 * rSges(7,1) + t72 * rSges(7,2) + rSges(7,3) * t163;
t191 = pkin(5) * t164 + pkin(8) * t163 + t36;
t190 = t197 * t115 + t198 * t116;
t189 = t198 * t115 - t197 * t116;
t113 = t115 ^ 2;
t114 = t116 ^ 2;
t188 = m(5) / 0.2e1;
t187 = m(6) / 0.2e1;
t186 = m(7) / 0.2e1;
t185 = -pkin(3) - pkin(4);
t93 = rSges(4,1) * t121 + rSges(4,2) * t124;
t184 = m(4) * t93;
t122 = sin(qJ(1));
t182 = pkin(1) * t122;
t74 = Icges(7,3) * t121 + (-Icges(7,5) * t123 + Icges(7,6) * t120) * t124;
t75 = Icges(7,6) * t121 + (-Icges(7,4) * t123 + Icges(7,2) * t120) * t124;
t181 = t124 * t120 * t75 + t121 * t74;
t174 = pkin(3) * t163 + qJ(4) * t164;
t180 = t113 * (pkin(3) * t124 + qJ(4) * t121) + t116 * t174;
t91 = pkin(3) * t121 - qJ(4) * t124;
t179 = -rSges(5,1) * t121 + rSges(5,3) * t124 - t91;
t178 = t116 * rSges(5,2);
t177 = t116 * rSges(4,3);
t76 = Icges(7,5) * t121 + (-Icges(7,1) * t123 + Icges(7,4) * t120) * t124;
t176 = t123 * t76;
t175 = -rSges(6,3) - qJ(5);
t172 = Icges(4,4) * t124;
t171 = Icges(6,4) * t121;
t167 = qJ(5) * t115;
t166 = qJ(5) * t116;
t165 = t115 * t124;
t159 = t113 + t114;
t158 = -m(5) - m(6) - m(7);
t157 = t187 + t186;
t70 = -t115 * t162 + t116 * t123;
t71 = t115 * t161 + t116 * t120;
t29 = Icges(7,5) * t71 + Icges(7,6) * t70 + Icges(7,3) * t165;
t31 = Icges(7,4) * t71 + Icges(7,2) * t70 + Icges(7,6) * t165;
t33 = Icges(7,1) * t71 + Icges(7,4) * t70 + Icges(7,5) * t165;
t10 = t121 * t29 + (t120 * t31 - t123 * t33) * t124;
t17 = t74 * t165 + t70 * t75 + t71 * t76;
t156 = t17 / 0.2e1 + t10 / 0.2e1;
t30 = Icges(7,5) * t73 + Icges(7,6) * t72 + Icges(7,3) * t163;
t32 = Icges(7,4) * t73 + Icges(7,2) * t72 + Icges(7,6) * t163;
t34 = Icges(7,1) * t73 + Icges(7,4) * t72 + Icges(7,5) * t163;
t11 = t121 * t30 + (t120 * t32 - t123 * t34) * t124;
t18 = t74 * t163 + t72 * t75 + t73 * t76;
t155 = t18 / 0.2e1 + t11 / 0.2e1;
t154 = rSges(5,1) * t163 + t115 * rSges(5,2) + rSges(5,3) * t164;
t125 = cos(qJ(1));
t117 = t125 * pkin(1);
t153 = t116 * pkin(2) + t115 * pkin(7) + t117;
t152 = -pkin(4) * t121 - t91;
t151 = t116 * pkin(7) - t182;
t106 = pkin(4) * t163;
t150 = t115 * (pkin(4) * t165 + t166) + t116 * (t106 - t167) + t180;
t149 = rSges(6,1) * t124 + rSges(6,2) * t121 + t152;
t148 = t199 * t193 + (-Icges(5,6) / 0.2e1 + t202 / 0.2e1) * t124;
t147 = rSges(6,1) * t164 - rSges(6,2) * t163;
t146 = -rSges(7,1) * t71 - rSges(7,2) * t70;
t145 = t153 + t174;
t144 = rSges(4,1) * t124 - rSges(4,2) * t121;
t77 = rSges(7,3) * t121 + (-rSges(7,1) * t123 + rSges(7,2) * t120) * t124;
t137 = pkin(5) * t124 - pkin(8) * t121 + t152 - t77;
t133 = -Icges(4,2) * t121 + t172;
t131 = -Icges(6,2) * t124 + t171;
t127 = rSges(4,1) * t163 - rSges(4,2) * t164 + t115 * rSges(4,3);
t126 = t106 + t145;
t96 = rSges(2,1) * t125 - t122 * rSges(2,2);
t94 = -t122 * rSges(2,1) - rSges(2,2) * t125;
t79 = rSges(3,1) * t116 - rSges(3,2) * t115 + t117;
t78 = -rSges(3,1) * t115 - rSges(3,2) * t116 - t182;
t42 = t179 * t116;
t41 = t179 * t115;
t40 = t149 * t116;
t39 = t149 * t115;
t38 = t127 + t153;
t37 = t177 + (-pkin(2) - t144) * t115 + t151;
t35 = rSges(7,3) * t165 - t146;
t28 = t145 + t154;
t27 = t178 + (-pkin(2) + (-rSges(5,1) - pkin(3)) * t124 + (-rSges(5,3) - qJ(4)) * t121) * t115 + t151;
t26 = t116 * t127 + (t144 * t115 - t177) * t115;
t25 = t137 * t116;
t24 = t137 * t115;
t23 = t175 * t115 + t126 + t147;
t22 = t175 * t116 + (-pkin(2) + (-rSges(6,1) - qJ(4)) * t121 + (rSges(6,2) + t185) * t124) * t115 + t151;
t21 = (-t124 * t176 + t181) * t121;
t20 = t121 * t36 - t77 * t163;
t19 = -t121 * t35 + t77 * t165;
t16 = t116 * t154 + (-t178 + (rSges(5,1) * t124 + rSges(5,3) * t121) * t115) * t115 + t180;
t15 = t126 - t167 + t191;
t14 = -t166 + (-pkin(2) + (-pkin(5) - qJ(4)) * t121 + (-rSges(7,3) - pkin(8) + t185) * t124) * t115 + t146 + t151;
t13 = (-t115 * t36 + t116 * t35) * t124;
t12 = t116 * t147 + (rSges(6,1) * t121 - rSges(6,2) * t124) * t113 + t150;
t9 = t30 * t163 + t32 * t72 + t34 * t73;
t8 = t29 * t163 + t31 * t72 + t33 * t73;
t7 = t30 * t165 + t32 * t70 + t34 * t71;
t6 = t29 * t165 + t31 * t70 + t33 * t71;
t5 = t191 * t116 + (t35 + (pkin(5) * t121 + pkin(8) * t124) * t115) * t115 + t150;
t4 = t115 * t9 - t116 * t8;
t3 = t115 * t7 - t116 * t6;
t2 = t121 * t18 + (t115 * t8 + t116 * t9) * t124;
t1 = t121 * t17 + (t115 * t6 + t116 * t7) * t124;
t43 = [Icges(2,3) + Icges(3,3) + m(7) * (t14 ^ 2 + t15 ^ 2) + m(6) * (t22 ^ 2 + t23 ^ 2) + m(5) * (t27 ^ 2 + t28 ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2) + m(3) * (t78 ^ 2 + t79 ^ 2) + m(2) * (t94 ^ 2 + t96 ^ 2) + t181 + (t171 - t176 + (Icges(4,2) + t205) * t124 + t203) * t124 + (t172 + (Icges(6,2) + t206) * t121 + t204) * t121; 0; m(3) + m(4) - t158; m(7) * (t14 * t25 + t15 * t24) + m(6) * (t22 * t40 + t23 * t39) + m(5) * (t27 * t42 + t28 * t41) + (-t37 * t184 + t148 * t116 + (Icges(5,6) * t195 + t133 * t196 + t201 * t183 + t202 * t194) * t124 - t156) * t116 + (-t38 * t184 + t148 * t115 + (Icges(5,6) * t196 + t133 * t194 + t202 * t183 + t201 * t195) * t124 + t155) * t115 + (t131 * t195 * t115 + t200 * t116 * t196 + (t199 * t115 + t131 * t116) * t183 + (t200 * t115 + t199 * t116) * t194) * t121; m(4) * t26 + m(5) * t16 + m(6) * t12 + m(7) * t5; m(7) * (t24 ^ 2 + t25 ^ 2 + t5 ^ 2) + m(6) * (t12 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t16 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(4) * (t159 * t93 ^ 2 + t26 ^ 2) + (t190 * t114 - t3) * t116 + (t4 + t189 * t113 + (t190 * t115 + t189 * t116) * t116) * t115; 0.2e1 * ((t115 * t15 + t116 * t14) * t186 + (t115 * t23 + t116 * t22) * t187 + (t115 * t28 + t116 * t27) * t188) * t121; t158 * t124; m(7) * (-t124 * t5 + (t115 * t24 + t116 * t25) * t121) + m(6) * (-t12 * t124 + (t115 * t39 + t116 * t40) * t121) + m(5) * (-t124 * t16 + (t115 * t41 + t116 * t42) * t121); 0.2e1 * (t188 + t157) * (t159 * t121 ^ 2 + t124 ^ 2); m(7) * (-t115 * t14 + t116 * t15) + m(6) * (-t115 * t22 + t116 * t23); 0; m(7) * (-t115 * t25 + t116 * t24) + m(6) * (-t115 * t40 + t116 * t39); 0; 0.2e1 * t157 * t159; m(7) * (t14 * t19 + t15 * t20) + t21 + (t156 * t115 + t155 * t116) * t124; m(7) * t13; m(7) * (t13 * t5 + t19 * t25 + t20 * t24) + (-t10 * t116 + t11 * t115) * t193 + t1 * t195 + t2 * t183 + (t3 * t183 + t194 * t4) * t124; m(7) * (-t124 * t13 + (t115 * t20 + t116 * t19) * t121); m(7) * (-t115 * t19 + t116 * t20); t121 * t21 + m(7) * (t13 ^ 2 + t19 ^ 2 + t20 ^ 2) + (t116 * t2 + t115 * t1 + t121 * (t10 * t115 + t11 * t116)) * t124;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t43(1) t43(2) t43(4) t43(7) t43(11) t43(16); t43(2) t43(3) t43(5) t43(8) t43(12) t43(17); t43(4) t43(5) t43(6) t43(9) t43(13) t43(18); t43(7) t43(8) t43(9) t43(10) t43(14) t43(19); t43(11) t43(12) t43(13) t43(14) t43(15) t43(20); t43(16) t43(17) t43(18) t43(19) t43(20) t43(21);];
Mq  = res;
