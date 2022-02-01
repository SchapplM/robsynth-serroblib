% Calculate time derivative of joint inertia matrix for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:18:48
% EndTime: 2022-01-23 09:18:51
% DurationCPUTime: 1.31s
% Computational Cost: add. (3646->175), mult. (2294->251), div. (0->0), fcn. (1646->10), ass. (0->106)
t99 = pkin(9) + qJ(5);
t93 = sin(t99);
t141 = qJD(5) * t93;
t95 = cos(t99);
t140 = qJD(5) * t95;
t100 = qJD(1) + qJD(3);
t172 = t100 * t93;
t171 = t100 * t95;
t101 = qJ(1) + pkin(8);
t97 = qJ(3) + t101;
t90 = sin(t97);
t170 = (rSges(5,3) + qJ(4)) * t90;
t103 = cos(pkin(9));
t134 = -rSges(5,1) * t103 - pkin(3);
t91 = cos(t97);
t168 = t91 * qJ(4) + t134 * t90;
t156 = rSges(6,2) * t93;
t157 = rSges(6,1) * t95;
t167 = -t156 + t157;
t166 = -pkin(2) * cos(t101) - cos(qJ(1)) * pkin(1);
t149 = Icges(6,4) * t95;
t119 = -Icges(6,2) * t93 + t149;
t41 = Icges(6,6) * t90 + t119 * t91;
t150 = Icges(6,4) * t93;
t120 = Icges(6,1) * t95 - t150;
t43 = Icges(6,5) * t90 + t120 * t91;
t122 = t41 * t93 - t43 * t95;
t165 = t122 * t90;
t40 = -Icges(6,6) * t91 + t119 * t90;
t42 = -Icges(6,5) * t91 + t120 * t90;
t124 = t40 * t93 - t42 * t95;
t164 = t124 * t91;
t118 = Icges(6,5) * t95 - Icges(6,6) * t93;
t62 = Icges(6,2) * t95 + t150;
t63 = Icges(6,1) * t93 + t149;
t121 = t62 * t93 - t63 * t95;
t163 = t118 * qJD(5) + t100 * t121;
t38 = -Icges(6,3) * t91 + t118 * t90;
t162 = 2 * m(4);
t161 = 2 * m(5);
t160 = 2 * m(6);
t83 = t90 * rSges(6,3);
t154 = t91 * rSges(6,3) + t90 * t156;
t151 = rSges(5,2) * sin(pkin(9));
t39 = Icges(6,3) * t90 + t118 * t91;
t145 = t100 * t39;
t144 = t100 * t90;
t143 = t100 * t91;
t104 = -pkin(7) - qJ(4);
t142 = t104 * t90;
t138 = t100 * t156;
t139 = -t91 * t138 + (-rSges(6,1) * t141 - rSges(6,2) * t140) * t90;
t137 = t100 * t151;
t92 = pkin(4) * t103 + pkin(3);
t133 = -t92 - t157;
t58 = t91 * rSges(4,1) - rSges(4,2) * t90;
t49 = -rSges(4,1) * t143 + rSges(4,2) * t144;
t132 = -pkin(2) * sin(t101) - sin(qJ(1)) * pkin(1);
t57 = -rSges(4,1) * t90 - rSges(4,2) * t91;
t68 = rSges(6,1) * t93 + rSges(6,2) * t95;
t109 = -t104 * t91 + t133 * t90;
t32 = t109 + t154;
t28 = t132 + t32;
t45 = t167 * t91 + t83;
t33 = t91 * t92 - t142 + t45;
t29 = t33 - t166;
t127 = t28 * t91 + t29 * t90;
t126 = t32 * t91 + t33 * t90;
t125 = t40 * t95 + t42 * t93;
t123 = t41 * t95 + t43 * t93;
t61 = Icges(6,5) * t93 + Icges(6,6) * t95;
t117 = (t120 - t62) * t141 + (t119 + t63) * t140;
t48 = t57 * t100;
t116 = t132 * qJD(1);
t115 = t166 * qJD(1);
t114 = (-qJD(5) * t122 + t163 * t90 - t40 * t171 - t42 * t172) * t90 / 0.2e1 - (-qJD(5) * t124 - t163 * t91 + t41 * t171 + t43 * t172) * t91 / 0.2e1 + (-t121 * t90 - t61 * t91 + t125) * t144 / 0.2e1 + (-t121 * t91 + t61 * t90 + t123) * t143 / 0.2e1;
t111 = qJD(5) * t61;
t35 = t170 + (-t151 - t134) * t91;
t34 = t91 * rSges(5,3) + t90 * t151 + t168;
t108 = -qJD(5) * t68 * t91 + rSges(6,3) * t143 + t90 * t138;
t79 = qJD(4) * t90;
t26 = rSges(5,3) * t143 + t168 * t100 + t90 * t137 + t79;
t80 = qJD(4) * t91;
t27 = t91 * t137 + t80 + (t134 * t91 - t170) * t100;
t13 = t80 - t139 + (t133 * t91 + t142 - t83) * t100;
t12 = t100 * t109 + t108 + t79;
t56 = t167 * qJD(5);
t47 = t58 - t166;
t46 = t132 + t57;
t44 = t90 * t157 - t154;
t37 = t115 + t49;
t36 = t48 + t116;
t31 = t35 - t166;
t30 = t132 + t34;
t21 = -t111 * t90 + t145;
t20 = -t100 * t38 - t91 * t111;
t19 = t115 + t27;
t18 = t116 + t26;
t11 = t115 + t13;
t10 = t116 + t12;
t9 = -t122 * t91 + t39 * t90;
t8 = t38 * t90 - t164;
t7 = -t39 * t91 - t165;
t6 = -t124 * t90 - t38 * t91;
t3 = ((-t45 + t83) * t100 + t139) * t90 + (t100 * t44 + t108) * t91;
t1 = [(t10 * t29 + t11 * t28) * t160 + (t18 * t31 + t19 * t30) * t161 + (t36 * t47 + t37 * t46) * t162 + t117; 0; 0; m(6) * (t10 * t33 + t11 * t32 + t12 * t29 + t13 * t28) + m(5) * (t18 * t35 + t19 * t34 + t26 * t31 + t27 * t30) + m(4) * (t36 * t58 + t37 * t57 + t49 * t46 + t48 * t47) + t117; 0; (t12 * t33 + t13 * t32) * t160 + (t26 * t35 + t27 * t34) * t161 + (t48 * t58 + t49 * t57) * t162 + t117; m(6) * (-t10 * t91 + t100 * t127 + t11 * t90) + m(5) * (-t18 * t91 + t19 * t90 + (t30 * t91 + t31 * t90) * t100); 0; m(6) * (t100 * t126 - t12 * t91 + t13 * t90) + m(5) * (-t26 * t91 + t27 * t90 + (t34 * t91 + t35 * t90) * t100); 0; m(6) * (-t127 * t56 + (-t10 * t90 - t11 * t91 + (t28 * t90 - t29 * t91) * t100) * t68) + t114; m(6) * t3; m(6) * (-t126 * t56 + (-t12 * t90 - t13 * t91 + (t32 * t90 - t33 * t91) * t100) * t68) + t114; 0; ((t44 * t90 + t45 * t91) * t3 + (t90 ^ 2 + t91 ^ 2) * t68 * t56) * t160 + (-t8 * t91 + t9 * t90) * t143 + t90 * ((t90 * t20 + (t8 + t165) * t100) * t90 + (t9 * t100 + (t140 * t40 + t141 * t42) * t91 + (-t21 - t123 * qJD(5) + (-t124 + t39) * t100) * t90) * t91) + (-t6 * t91 + t7 * t90) * t144 - t91 * ((t91 * t21 + (t7 + t164) * t100) * t91 + (t6 * t100 + (-t140 * t41 - t141 * t43 + t145) * t90 + (t125 * qJD(5) - t122 * t100 - t20) * t91) * t90);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
