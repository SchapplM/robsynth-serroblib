% Calculate time derivative of joint inertia matrix for
% S5RPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR8_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR8_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR8_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR8_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:01:04
% EndTime: 2019-12-31 18:01:08
% DurationCPUTime: 1.57s
% Computational Cost: add. (4095->200), mult. (4982->293), div. (0->0), fcn. (4840->8), ass. (0->103)
t99 = sin(qJ(5));
t151 = qJD(5) * t99;
t164 = cos(qJ(1));
t97 = sin(pkin(8));
t144 = t164 * t97;
t163 = sin(qJ(1));
t147 = t163 * pkin(1);
t98 = cos(pkin(8));
t90 = pkin(3) * t98 + pkin(2);
t180 = pkin(3) * t144 - t163 * t90 - t147;
t157 = t164 * pkin(1) + t163 * qJ(2);
t61 = rSges(3,1) * t164 + rSges(3,3) * t163 + t157;
t100 = cos(qJ(5));
t178 = qJD(1) - qJD(4);
t149 = pkin(8) + qJ(4);
t135 = sin(t149);
t136 = cos(t149);
t62 = -t163 * t135 - t164 * t136;
t50 = t178 * t62;
t63 = t164 * t135 - t163 * t136;
t124 = -t100 * t50 + t63 * t151;
t51 = t178 * t63;
t122 = t100 * t51 + t62 * t151;
t179 = -t164 * pkin(2) - t157;
t143 = t163 * t97;
t140 = pkin(3) * t143 + t164 * t90 + t157;
t150 = qJD(5) * t100;
t121 = t150 * t62 - t51 * t99;
t123 = t150 * t63 + t50 * t99;
t169 = 2 * m(6);
t156 = rSges(6,1) * t100;
t162 = rSges(6,2) * t99;
t73 = (-t156 + t162) * qJD(5);
t82 = -t99 * rSges(6,1) - rSges(6,2) * t100;
t177 = -t62 * (Icges(6,5) * t124 + Icges(6,6) * t123 + Icges(6,3) * t51) + t63 * (Icges(6,5) * t122 + Icges(6,6) * t121 + Icges(6,3) * t50) + t82 * t73 * t169;
t152 = Icges(6,4) * t100;
t126 = Icges(6,2) * t99 - t152;
t36 = -Icges(6,6) * t62 + t126 * t63;
t155 = Icges(6,4) * t99;
t127 = -Icges(6,1) * t100 + t155;
t38 = -Icges(6,5) * t62 + t127 * t63;
t132 = t100 * t38 - t36 * t99;
t173 = t132 * t62;
t37 = Icges(6,6) * t63 + t126 * t62;
t39 = Icges(6,5) * t63 + t127 * t62;
t131 = t100 * t39 - t37 * t99;
t174 = t131 * t63;
t125 = -Icges(6,5) * t100 + Icges(6,6) * t99;
t34 = -Icges(6,3) * t62 + t125 * t63;
t35 = Icges(6,3) * t63 + t125 * t62;
t176 = t63 * t34 - t62 * t35 - t173 - t174;
t119 = t122 * rSges(6,1) + t50 * rSges(6,3);
t9 = -rSges(6,2) * t121 - t51 * pkin(4) - t50 * pkin(7) - t119;
t120 = t124 * rSges(6,1) + t51 * rSges(6,3);
t8 = rSges(6,2) * t123 - t50 * pkin(4) + t51 * pkin(7) + t120;
t142 = -pkin(4) + t162;
t159 = t63 * rSges(6,3) - t62 * t156;
t27 = -pkin(7) * t63 - t142 * t62 - t159;
t160 = t62 * rSges(6,3) + t63 * t156;
t26 = -t62 * pkin(7) + t142 * t63 - t160;
t141 = (Icges(6,2) * t100 + t127 + t155) * t151;
t71 = t126 * qJD(5);
t77 = -Icges(6,1) * t99 - t152;
t172 = -(qJD(5) * t77 + t71) * t100 - t141;
t170 = 2 * m(5);
t95 = t164 * qJ(2);
t158 = qJD(1) * t95 + qJD(2) * t163;
t31 = -t51 * rSges(5,1) + t50 * rSges(5,2);
t30 = -t50 * rSges(5,1) - t51 * rSges(5,2);
t52 = -t63 * rSges(5,1) + t62 * rSges(5,2);
t53 = rSges(5,1) * t62 + rSges(5,2) * t63;
t114 = -pkin(2) * t163 - t147;
t110 = t164 * t98 + t143;
t109 = -t163 * t98 + t144;
t108 = t95 + t180;
t105 = -rSges(3,1) * t163 + rSges(3,3) * t164 - t147;
t103 = t180 * qJD(1) + t158;
t93 = qJD(2) * t164;
t102 = -t140 * qJD(1) + t93;
t70 = t125 * qJD(5);
t101 = (-t100 * t37 - t99 * t39) * t50 / 0.2e1 + (-t100 * t36 - t99 * t38) * t51 / 0.2e1 - (-qJD(5) * t132 - t100 * (Icges(6,4) * t124 + Icges(6,2) * t123 + Icges(6,6) * t51) - t99 * (Icges(6,1) * t124 + Icges(6,4) * t123 + Icges(6,5) * t51) - t62 * t70) * t62 / 0.2e1 + (-qJD(5) * t131 - t100 * (Icges(6,4) * t122 + Icges(6,2) * t121 + Icges(6,6) * t50) - t99 * (Icges(6,1) * t122 + Icges(6,4) * t121 + Icges(6,5) * t50) + t63 * t70) * t63 / 0.2e1 + (t63 * t50 - t62 * t51) * (-Icges(6,5) * t99 - Icges(6,6) * t100);
t74 = t82 ^ 2;
t66 = t109 * qJD(1);
t65 = t110 * qJD(1);
t60 = t95 + t105;
t55 = -t61 * qJD(1) + t93;
t54 = qJD(1) * t105 + t158;
t45 = rSges(4,1) * t110 - rSges(4,2) * t109 - t179;
t44 = rSges(4,1) * t109 + rSges(4,2) * t110 + t114 + t95;
t43 = -t65 * rSges(4,1) + t66 * rSges(4,2) + t179 * qJD(1) + t93;
t42 = t66 * rSges(4,1) + t65 * rSges(4,2) + qJD(1) * t114 + t158;
t41 = t162 * t62 + t159;
t40 = t162 * t63 - t160;
t33 = -t53 + t140;
t32 = t108 - t52;
t25 = t102 - t30;
t24 = t103 - t31;
t23 = t140 - t27;
t22 = t108 - t26;
t7 = t102 - t8;
t6 = t103 - t9;
t1 = t50 * t40 + t63 * t120 - t51 * t41 + t62 * t119 + (t121 * t62 + t123 * t63) * rSges(6,2);
t2 = [-t77 * t150 - t100 * t71 + (t22 * t7 + t23 * t6) * t169 + (t24 * t33 + t25 * t32) * t170 + 0.2e1 * m(3) * (t54 * t61 + t55 * t60) + 0.2e1 * m(4) * (t42 * t45 + t43 * t44) - t141; m(6) * (t163 * t7 - t164 * t6 + (t163 * t23 + t164 * t22) * qJD(1)) + m(5) * (t163 * t25 - t164 * t24 + (t163 * t33 + t164 * t32) * qJD(1)) + m(3) * (t163 * t55 - t164 * t54 + (t163 * t61 + t164 * t60) * qJD(1)) + m(4) * (t163 * t43 - t164 * t42 + (t163 * t45 + t164 * t44) * qJD(1)); 0; 0; 0; 0; m(6) * (t22 * t8 + t23 * t9 + t26 * t7 + t27 * t6) + m(5) * (t24 * t53 + t25 * t52 + t30 * t32 + t31 * t33) - t172; m(5) * (t30 * t163 - t31 * t164 + (t163 * t53 + t164 * t52) * qJD(1)) + m(6) * (t8 * t163 - t9 * t164 + (t163 * t27 + t164 * t26) * qJD(1)); 0; (t30 * t52 + t31 * t53) * t170 + (t26 * t8 + t27 * t9) * t169 + t172; m(6) * ((-t22 * t62 - t23 * t63) * t73 + (t22 * t51 - t23 * t50 - t6 * t63 - t62 * t7) * t82) + t101; m(6) * ((-t163 * t62 + t164 * t63) * t73 + (t163 * t51 + t164 * t50 + (-t163 * t63 - t164 * t62) * qJD(1)) * t82); -m(6) * t1; m(6) * ((-t26 * t62 - t27 * t63) * t73 + (t26 * t51 - t27 * t50 - t62 * t8 - t63 * t9) * t82) - t101; ((t40 * t1 + t74 * t50) * t169 + (t174 + t176) * t51 + (0.3e1 * t50 * t35 + t177) * t63) * t63 + ((t132 - t35) * t51 * t63 + (0.3e1 * t51 * t34 + t177) * t62 - (t173 + t176 + (t34 + t131) * t63) * t50 + (t41 * t1 - t74 * t51) * t169) * t62;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
