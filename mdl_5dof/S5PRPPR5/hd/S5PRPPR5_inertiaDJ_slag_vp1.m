% Calculate time derivative of joint inertia matrix for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:38:03
% EndTime: 2019-12-31 17:38:09
% DurationCPUTime: 3.60s
% Computational Cost: add. (4279->321), mult. (12686->511), div. (0->0), fcn. (13371->8), ass. (0->150)
t147 = sin(qJ(2));
t149 = cos(qJ(2));
t222 = ((-Icges(3,6) + Icges(4,6)) * t149 + (-Icges(4,4) - Icges(3,5)) * t147) * qJD(2);
t144 = sin(pkin(7));
t141 = t144 ^ 2;
t145 = cos(pkin(7));
t142 = t145 ^ 2;
t215 = t141 + t142;
t221 = t222 * t144;
t220 = t222 * t145;
t143 = sin(pkin(8));
t205 = cos(pkin(8));
t184 = t147 * t205;
t132 = -t143 * t149 + t184;
t213 = t215 * qJD(2);
t210 = 2 * m(6);
t209 = m(4) / 0.2e1;
t208 = m(5) / 0.2e1;
t207 = m(6) / 0.2e1;
t135 = pkin(2) * t147 - qJ(3) * t149;
t206 = t215 * (-qJD(2) * t135 + qJD(3) * t147);
t131 = t143 * t147 + t149 * t205;
t130 = t131 * qJD(2);
t146 = sin(qJ(5));
t200 = t130 * t146;
t148 = cos(qJ(5));
t199 = t130 * t148;
t197 = t144 * t147;
t196 = t144 * t149;
t195 = t145 * t147;
t194 = t145 * t149;
t177 = pkin(2) * t149 + qJ(3) * t147;
t192 = t215 * t177;
t126 = qJD(2) * t177 - qJD(3) * t149;
t178 = rSges(4,1) * t149 + rSges(4,3) * t147;
t191 = -qJD(2) * t178 - t126;
t136 = rSges(4,1) * t147 - rSges(4,3) * t149;
t190 = -t135 - t136;
t189 = qJD(2) * t147;
t188 = qJD(2) * t149;
t187 = qJD(5) * t132;
t185 = -pkin(3) * t147 - t135;
t183 = -pkin(3) * t189 * t215 + t206;
t182 = t192 + (t144 * t196 + t145 * t194) * pkin(3);
t181 = -rSges(5,1) * t132 + rSges(5,2) * t131 + t185;
t180 = -pkin(3) * t188 - t126;
t137 = rSges(3,1) * t147 + rSges(3,2) * t149;
t106 = t132 * t144;
t107 = t131 * t144;
t85 = -t107 * t146 + t145 * t148;
t86 = t107 * t148 + t145 * t146;
t50 = Icges(6,4) * t86 + Icges(6,2) * t85 - Icges(6,6) * t106;
t52 = Icges(6,1) * t86 + Icges(6,4) * t85 - Icges(6,5) * t106;
t176 = -t146 * t50 + t148 * t52;
t108 = t132 * t145;
t109 = t131 * t145;
t87 = -t109 * t146 - t144 * t148;
t88 = t109 * t148 - t144 * t146;
t51 = Icges(6,4) * t88 + Icges(6,2) * t87 - Icges(6,6) * t108;
t53 = Icges(6,1) * t88 + Icges(6,4) * t87 - Icges(6,5) * t108;
t175 = -t146 * t51 + t148 * t53;
t67 = rSges(6,3) * t131 + (rSges(6,1) * t148 - rSges(6,2) * t146) * t132;
t173 = -pkin(4) * t132 - pkin(6) * t131 + t185 - t67;
t129 = -qJD(2) * t184 + t143 * t188;
t163 = -rSges(5,1) * t130 + rSges(5,2) * t129 + t180;
t160 = -t146 * t187 + t199;
t161 = -t148 * t187 - t200;
t46 = rSges(6,1) * t160 + rSges(6,2) * t161 + rSges(6,3) * t129;
t162 = -pkin(4) * t130 - pkin(6) * t129 + t180 - t46;
t151 = qJD(2) * t132;
t97 = t145 * t151;
t96 = t145 * t130;
t95 = t144 * t151;
t94 = t144 * t130;
t93 = t190 * t145;
t92 = t190 * t144;
t80 = t191 * t145;
t79 = t191 * t144;
t78 = t137 * t213;
t77 = Icges(5,1) * t109 + Icges(5,4) * t108 - Icges(5,5) * t144;
t76 = Icges(5,1) * t107 + Icges(5,4) * t106 + Icges(5,5) * t145;
t75 = Icges(5,4) * t109 + Icges(5,2) * t108 - Icges(5,6) * t144;
t74 = Icges(5,4) * t107 + Icges(5,2) * t106 + Icges(5,6) * t145;
t73 = -Icges(5,1) * t97 + Icges(5,4) * t96;
t72 = -Icges(5,1) * t95 + Icges(5,4) * t94;
t71 = -Icges(5,4) * t97 + Icges(5,2) * t96;
t70 = -Icges(5,4) * t95 + Icges(5,2) * t94;
t69 = -Icges(5,5) * t97 + Icges(5,6) * t96;
t68 = -Icges(5,5) * t95 + Icges(5,6) * t94;
t66 = Icges(6,5) * t131 + (Icges(6,1) * t148 - Icges(6,4) * t146) * t132;
t65 = Icges(6,6) * t131 + (Icges(6,4) * t148 - Icges(6,2) * t146) * t132;
t64 = Icges(6,3) * t131 + (Icges(6,5) * t148 - Icges(6,6) * t146) * t132;
t63 = t181 * t145;
t62 = t181 * t144;
t61 = qJD(5) * t87 - t148 * t97;
t60 = -qJD(5) * t88 + t146 * t97;
t59 = qJD(5) * t85 - t148 * t95;
t58 = -qJD(5) * t86 + t146 * t95;
t57 = t163 * t145;
t56 = t163 * t144;
t55 = rSges(6,1) * t88 + rSges(6,2) * t87 - rSges(6,3) * t108;
t54 = rSges(6,1) * t86 + rSges(6,2) * t85 - rSges(6,3) * t106;
t49 = Icges(6,5) * t88 + Icges(6,6) * t87 - Icges(6,3) * t108;
t48 = Icges(6,5) * t86 + Icges(6,6) * t85 - Icges(6,3) * t106;
t47 = t178 * t215 + t192;
t45 = Icges(6,1) * t160 + Icges(6,4) * t161 + Icges(6,5) * t129;
t44 = Icges(6,4) * t160 + Icges(6,2) * t161 + Icges(6,6) * t129;
t43 = Icges(6,5) * t160 + Icges(6,6) * t161 + Icges(6,3) * t129;
t42 = -t136 * t213 + t206;
t41 = t173 * t145;
t40 = t173 * t144;
t39 = rSges(6,1) * t61 + rSges(6,2) * t60 - rSges(6,3) * t96;
t38 = rSges(6,1) * t59 + rSges(6,2) * t58 - rSges(6,3) * t94;
t37 = Icges(6,1) * t61 + Icges(6,4) * t60 - Icges(6,5) * t96;
t36 = Icges(6,1) * t59 + Icges(6,4) * t58 - Icges(6,5) * t94;
t35 = Icges(6,4) * t61 + Icges(6,2) * t60 - Icges(6,6) * t96;
t34 = Icges(6,4) * t59 + Icges(6,2) * t58 - Icges(6,6) * t94;
t33 = Icges(6,5) * t61 + Icges(6,6) * t60 - Icges(6,3) * t96;
t32 = Icges(6,5) * t59 + Icges(6,6) * t58 - Icges(6,3) * t94;
t31 = t108 * t67 + t131 * t55;
t30 = -t106 * t67 - t131 * t54;
t29 = t162 * t145;
t28 = t162 * t144;
t27 = t144 * (rSges(5,1) * t107 + rSges(5,2) * t106) + t145 * (rSges(5,1) * t109 + rSges(5,2) * t108) + t182;
t26 = t144 * (-rSges(5,1) * t95 + rSges(5,2) * t94) + t145 * (-rSges(5,1) * t97 + rSges(5,2) * t96) + t183;
t24 = t106 * t55 - t108 * t54;
t23 = -t108 * t64 + t65 * t87 + t66 * t88;
t22 = -t106 * t64 + t65 * t85 + t66 * t86;
t21 = t131 * t49 + t132 * t175;
t20 = t131 * t48 + t132 * t176;
t19 = -t108 * t49 + t51 * t87 + t53 * t88;
t18 = -t108 * t48 + t50 * t87 + t52 * t88;
t17 = -t106 * t49 + t51 * t85 + t53 * t86;
t16 = -t106 * t48 + t50 * t85 + t52 * t86;
t15 = (pkin(4) * t109 - pkin(6) * t108 + t55) * t145 + (pkin(4) * t107 - pkin(6) * t106 + t54) * t144 + t182;
t14 = t108 * t46 + t129 * t55 + t131 * t39 + t67 * t96;
t13 = -t106 * t46 - t129 * t54 - t131 * t38 - t67 * t94;
t12 = (-pkin(4) * t97 - pkin(6) * t96 + t39) * t145 + (-pkin(4) * t95 - pkin(6) * t94 + t38) * t144 + t183;
t11 = t106 * t39 - t108 * t38 - t54 * t96 + t55 * t94;
t10 = -t108 * t33 + t35 * t87 + t37 * t88 - t49 * t96 + t51 * t60 + t53 * t61;
t9 = -t108 * t32 + t34 * t87 + t36 * t88 - t48 * t96 + t50 * t60 + t52 * t61;
t8 = -t106 * t33 + t35 * t85 + t37 * t86 - t49 * t94 + t51 * t58 + t53 * t59;
t7 = -t106 * t32 + t34 * t85 + t36 * t86 - t48 * t94 + t50 * t58 + t52 * t59;
t6 = t129 * t49 + t131 * t33 + t175 * t130 + (-t146 * t35 + t148 * t37 + (-t146 * t53 - t148 * t51) * qJD(5)) * t132;
t5 = t129 * t48 + t131 * t32 + t176 * t130 + (-t146 * t34 + t148 * t36 + (-t146 * t52 - t148 * t50) * qJD(5)) * t132;
t4 = t10 * t144 - t145 * t9;
t3 = t144 * t8 - t145 * t7;
t2 = -t10 * t108 - t19 * t96 - t9 * t106 - t18 * t94 + (-t108 * t43 + t44 * t87 + t45 * t88 + t60 * t65 + t61 * t66 - t64 * t96) * t131 + t23 * t129;
t1 = -t8 * t108 - t17 * t96 - t7 * t106 - t16 * t94 + (-t106 * t43 + t44 * t85 + t45 * t86 + t58 * t65 + t59 * t66 - t64 * t94) * t131 + t22 * t129;
t25 = [0; -m(3) * t78 + m(4) * t42 + m(5) * t26 + m(6) * t12; (t12 * t15 + t28 * t40 + t29 * t41) * t210 + 0.2e1 * m(5) * (t26 * t27 + t56 * t62 + t57 * t63) + 0.2e1 * m(4) * (t42 * t47 + t79 * t92 + t80 * t93) + 0.2e1 * m(3) * (qJD(2) * t137 - t78) * t215 * (rSges(3,1) * t149 - rSges(3,2) * t147) + (-t3 + (t106 * t70 + t107 * t72 + t145 * t68 + t74 * t94 - t76 * t95) * t145 - t221 * t142) * t145 + (t4 + (t108 * t71 + t109 * t73 - t144 * t69 + t75 * t96 - t77 * t97) * t144 + t220 * t141 + (-t108 * t70 - t109 * t72 - t74 * t96 + t76 * t97 - t106 * t71 - t107 * t73 - t75 * t94 + t77 * t95 + (t68 - t221) * t144 + (-t69 + t220) * t145) * t145) * t144; (m(4) + m(5) + m(6)) * t189; m(6) * (-t12 * t149 + t195 * t29 + t197 * t28) + m(5) * (-t149 * t26 + t195 * t57 + t197 * t56) + m(4) * (-t149 * t42 + t195 * t80 + t197 * t79) + 0.2e1 * ((t147 * t15 + t194 * t41 + t196 * t40) * t207 + (t147 * t27 + t194 * t63 + t196 * t62) * t208 + (t147 * t47 + t194 * t93 + t196 * t92) * t209) * qJD(2); 0.4e1 * (t209 + t208 + t207) * (-0.1e1 + t215) * t147 * t188; 0; m(6) * (-t144 * t29 + t145 * t28) + m(5) * (-t144 * t57 + t145 * t56); 0; 0; m(6) * t11; t144 * t2 / 0.2e1 - t145 * t1 / 0.2e1 - t96 * (t144 * t19 - t145 * t18) / 0.2e1 - t108 * t4 / 0.2e1 - t94 * (t144 * t17 - t145 * t16) / 0.2e1 - t106 * t3 / 0.2e1 + t129 * (t144 * t21 - t145 * t20) / 0.2e1 + t131 * (t144 * t6 - t145 * t5) / 0.2e1 + m(6) * (t11 * t15 + t12 * t24 + t13 * t41 + t14 * t40 + t28 * t31 + t29 * t30); m(6) * (-t11 * t149 + (t13 * t145 + t14 * t144) * t147 + (t147 * t24 + (t144 * t31 + t145 * t30) * t149) * qJD(2)); m(6) * (-t13 * t144 + t14 * t145); (t11 * t24 + t13 * t30 + t14 * t31) * t210 - t96 * (-t106 * t18 - t108 * t19) - t108 * t2 - t94 * (-t106 * t16 - t108 * t17) - t106 * t1 + t129 * (-t20 * t106 - t21 * t108) + (-t5 * t106 - t6 * t108 + (t199 * t66 - t200 * t65 + (-t146 * t44 + t148 * t45 + (-t146 * t66 - t148 * t65) * qJD(5)) * t132 + t131 * t43) * t131 + (-t21 - t23) * t96 + (-t20 - t22) * t94 + (0.3e1 * t131 * t64 + 0.2e1 * (-t146 * t65 + t148 * t66) * t132) * t129) * t131;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t25(1), t25(2), t25(4), t25(7), t25(11); t25(2), t25(3), t25(5), t25(8), t25(12); t25(4), t25(5), t25(6), t25(9), t25(13); t25(7), t25(8), t25(9), t25(10), t25(14); t25(11), t25(12), t25(13), t25(14), t25(15);];
Mq = res;
