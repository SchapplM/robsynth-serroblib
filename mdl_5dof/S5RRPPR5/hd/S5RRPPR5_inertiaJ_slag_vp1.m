% Calculate joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:54
% EndTime: 2019-12-31 19:28:59
% DurationCPUTime: 1.63s
% Computational Cost: add. (2076->224), mult. (2694->335), div. (0->0), fcn. (2753->8), ass. (0->116)
t193 = Icges(4,1) + Icges(5,1);
t187 = Icges(5,4) + Icges(4,5);
t192 = Icges(4,6) - Icges(5,6);
t102 = qJ(2) + pkin(8);
t96 = sin(t102);
t191 = (Icges(4,4) - Icges(5,5)) * t96;
t107 = sin(qJ(2));
t110 = cos(qJ(2));
t82 = rSges(3,1) * t107 + rSges(3,2) * t110;
t190 = m(3) * t82;
t189 = Icges(3,5) * t107;
t188 = Icges(3,6) * t110;
t97 = cos(t102);
t186 = t193 * t97 - t191;
t185 = Icges(5,2) + Icges(3,3) + Icges(4,3);
t184 = Icges(3,5) * t110 - Icges(3,6) * t107 + t187 * t97 - t192 * t96;
t183 = t189 / 0.2e1 + t188 / 0.2e1;
t106 = sin(qJ(5));
t109 = cos(qJ(5));
t68 = -t106 * t97 + t109 * t96;
t182 = t68 / 0.2e1;
t108 = sin(qJ(1));
t181 = -t108 / 0.2e1;
t180 = t108 / 0.2e1;
t111 = cos(qJ(1));
t179 = -t111 / 0.2e1;
t178 = t111 / 0.2e1;
t123 = t106 * t96 + t109 * t97;
t177 = -t123 / 0.2e1;
t176 = t108 * pkin(6);
t175 = Icges(6,5) * t182 + Icges(6,6) * t177;
t174 = -t184 * t108 + t185 * t111;
t173 = t185 * t108 + t184 * t111;
t103 = t108 ^ 2;
t104 = t111 ^ 2;
t172 = m(5) / 0.2e1;
t171 = pkin(2) * t107;
t101 = t111 * pkin(6);
t94 = pkin(2) * t110 + pkin(1);
t91 = t111 * t94;
t170 = t108 * (t101 + (-pkin(1) + t94) * t108) + t111 * (-pkin(1) * t111 - t176 + t91);
t59 = t68 * t111;
t60 = t123 * t111;
t169 = t60 * rSges(6,1) + t59 * rSges(6,2);
t153 = t111 * t97;
t154 = t111 * t96;
t168 = pkin(3) * t153 + qJ(4) * t154;
t166 = rSges(3,1) * t110;
t167 = t108 * rSges(3,3) + t111 * t166;
t165 = rSges(3,2) * t107;
t163 = Icges(4,4) * t97;
t58 = t123 * t108;
t162 = Icges(6,4) * t58;
t160 = Icges(5,5) * t97;
t159 = Icges(6,5) * t58;
t57 = t68 * t108;
t158 = Icges(6,2) * t57;
t157 = Icges(6,6) * t57;
t156 = qJ(4) * t96;
t155 = t111 * rSges(3,3);
t150 = Icges(6,3) * t111;
t149 = t103 + t104;
t148 = t172 + m(6) / 0.2e1;
t112 = Icges(6,6) * t111 + t158 + t162;
t113 = Icges(6,1) * t58 + Icges(6,4) * t57 + Icges(6,5) * t111;
t31 = Icges(6,4) * t68 - Icges(6,2) * t123;
t32 = Icges(6,1) * t68 - Icges(6,4) * t123;
t147 = t111 * t175 + t57 * t31 / 0.2e1 + t58 * t32 / 0.2e1 + t112 * t177 + t113 * t182;
t105 = -qJ(3) - pkin(6);
t146 = -rSges(6,3) - pkin(7) - t105;
t20 = Icges(6,4) * t60 + Icges(6,2) * t59 - Icges(6,6) * t108;
t21 = Icges(6,1) * t60 + Icges(6,4) * t59 - Icges(6,5) * t108;
t145 = t108 * t175 - t31 * t59 / 0.2e1 - t32 * t60 / 0.2e1 + t20 * t123 / 0.2e1 - t21 * t68 / 0.2e1;
t144 = rSges(5,1) * t153 + t108 * rSges(5,2) + rSges(5,3) * t154;
t143 = -pkin(3) * t96 + qJ(4) * t97 - t171;
t142 = -rSges(4,1) * t96 - rSges(4,2) * t97 - t171;
t141 = t103 * (pkin(3) * t97 + t156) + t111 * t168 + t170;
t140 = -t108 * t105 + t91;
t139 = -rSges(5,1) * t96 + rSges(5,3) * t97 + t143;
t137 = rSges(4,1) * t97 - rSges(4,2) * t96;
t136 = -t58 * rSges(6,1) - t57 * rSges(6,2);
t19 = Icges(6,5) * t60 + Icges(6,6) * t59 - Icges(6,3) * t108;
t135 = t111 * ((t111 * t19 + t57 * t20 + t58 * t21) * t108 - (Icges(6,1) * t58 ^ 2 + (t158 + 0.2e1 * t162) * t57 + (t150 + 0.2e1 * t157 + 0.2e1 * t159) * t111) * t111) - t108 * ((-t108 * t19 + t20 * t59 + t21 * t60) * t108 - (t60 * t113 + t59 * t112 - t108 * (t150 + t157 + t159)) * t111);
t130 = -t165 + t166;
t127 = -Icges(4,2) * t96 + t163;
t124 = Icges(5,3) * t96 + t160;
t33 = rSges(6,1) * t68 - rSges(6,2) * t123;
t114 = -pkin(4) * t96 + t143 - t33;
t15 = t114 * t108;
t16 = t114 * t111;
t120 = t108 * t15 + t111 * t16;
t119 = rSges(4,1) * t153 - rSges(4,2) * t154 + t108 * rSges(4,3);
t13 = t146 * t111 + (-t156 - t94 + (-pkin(3) - pkin(4)) * t97) * t108 + t136;
t89 = pkin(4) * t153;
t14 = t108 * t146 + t168 + t169 + t89 + t91;
t115 = m(6) * (t108 * t14 + t111 * t13);
t84 = rSges(2,1) * t111 - t108 * rSges(2,2);
t83 = -t108 * rSges(2,1) - rSges(2,2) * t111;
t39 = t142 * t111;
t38 = t142 * t108;
t37 = t176 + (pkin(1) - t165) * t111 + t167;
t36 = t155 + t101 + (-pkin(1) - t130) * t108;
t29 = t119 + t140;
t28 = (rSges(4,3) - t105) * t111 + (-t137 - t94) * t108;
t26 = t139 * t111;
t25 = t139 * t108;
t24 = t111 * (-t111 * t165 + t167) + (t108 * t130 - t155) * t108;
t23 = -rSges(6,3) * t108 + t169;
t22 = rSges(6,3) * t111 - t136;
t18 = t140 + t144 + t168;
t17 = (rSges(5,2) - t105) * t111 + (-t94 + (-rSges(5,1) - pkin(3)) * t97 + (-rSges(5,3) - qJ(4)) * t96) * t108;
t12 = t111 * t119 + (-t111 * rSges(4,3) + t108 * t137) * t108 + t170;
t11 = -t108 * t22 - t111 * t23;
t8 = t111 * t144 + (-t111 * rSges(5,2) + (rSges(5,1) * t97 + rSges(5,3) * t96) * t108) * t108 + t141;
t3 = (t23 + t89) * t111 + (t108 * t97 * pkin(4) + t22) * t108 + t141;
t1 = [Icges(3,2) * t110 ^ 2 - t123 * t31 + t68 * t32 + Icges(2,3) + m(6) * (t13 ^ 2 + t14 ^ 2) + m(5) * (t17 ^ 2 + t18 ^ 2) + m(4) * (t28 ^ 2 + t29 ^ 2) + m(3) * (t36 ^ 2 + t37 ^ 2) + m(2) * (t83 ^ 2 + t84 ^ 2) + ((Icges(4,2) + Icges(5,3)) * t97 + t191) * t97 + (t193 * t96 - t160 + t163) * t96 + (Icges(3,1) * t107 + 0.2e1 * Icges(3,4) * t110) * t107; m(6) * (t13 * t16 + t14 * t15) + m(5) * (t17 * t26 + t18 * t25) + m(4) * (t28 * t39 + t29 * t38) + (t187 * t96 + t192 * t97 + t188 + t189) * (t104 / 0.2e1 + t103 / 0.2e1) + ((Icges(4,6) * t178 + Icges(5,6) * t179 + t124 * t180 + t127 * t181) * t97 + (t187 * t178 + t186 * t181) * t96 - t147 - t36 * t190 + t183 * t111) * t111 + ((Icges(4,6) * t180 + Icges(5,6) * t181 + t124 * t179 + t127 * t178) * t97 + (t186 * t178 + t187 * t180) * t96 - t145 - t37 * t190 + t183 * t108) * t108; m(6) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2 + t8 ^ 2) + m(4) * (t12 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(3) * (t149 * t82 ^ 2 + t24 ^ 2) - t135 + t173 * t108 * t103 + (t174 * t104 + (t174 * t108 + t173 * t111) * t108) * t111; m(6) * (t108 * t13 - t111 * t14) + m(5) * (t108 * t17 - t111 * t18) + m(4) * (t108 * t28 - t111 * t29); m(6) * (t108 * t16 - t111 * t15) + m(5) * (t108 * t26 - t111 * t25) + m(4) * (t108 * t39 - t111 * t38); 0.2e1 * (m(4) / 0.2e1 + t148) * t149; 0.2e1 * (t115 / 0.2e1 + (t108 * t18 + t111 * t17) * t172) * t96; m(6) * (t120 * t96 - t97 * t3) + m(5) * (-t97 * t8 + (t108 * t25 + t111 * t26) * t96); 0; 0.2e1 * t148 * (t149 * t96 ^ 2 + t97 ^ 2); t108 * t145 + t111 * t147 + t115 * t33; m(6) * (t11 * t3 + t120 * t33) + t135; 0; m(6) * (t149 * t33 * t96 - t11 * t97); m(6) * (t149 * t33 ^ 2 + t11 ^ 2) - t135;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
