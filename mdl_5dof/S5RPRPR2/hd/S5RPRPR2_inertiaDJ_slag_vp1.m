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
% Datum: 2020-01-03 11:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:33:38
% EndTime: 2020-01-03 11:33:45
% DurationCPUTime: 1.78s
% Computational Cost: add. (3646->187), mult. (2294->267), div. (0->0), fcn. (1646->10), ass. (0->105)
t108 = qJ(1) + pkin(8);
t143 = cos(qJ(1)) * pkin(1) + pkin(2) * cos(t108);
t106 = pkin(9) + qJ(5);
t99 = sin(t106);
t138 = qJD(5) * t99;
t101 = cos(t106);
t137 = qJD(5) * t101;
t107 = qJD(1) + qJD(3);
t154 = Icges(6,4) * t99;
t124 = Icges(6,1) * t101 - t154;
t103 = qJ(3) + t108;
t95 = sin(t103);
t96 = cos(t103);
t167 = Icges(6,5) * t95 + t124 * t96;
t176 = t167 * t107;
t142 = Icges(6,4) * t101;
t123 = -Icges(6,2) * t99 + t142;
t168 = Icges(6,6) * t95 + t123 * t96;
t175 = t168 * t107;
t155 = rSges(5,2) * sin(pkin(9));
t174 = pkin(3) - t155;
t173 = rSges(5,3) + qJ(4);
t156 = rSges(6,1) * t101;
t160 = rSges(6,2) * t99;
t172 = t156 - t160;
t57 = t95 * rSges(4,1) + t96 * rSges(4,2);
t129 = -t101 * t167 + t168 * t99;
t171 = t129 * t95;
t144 = pkin(2) * sin(t108) + sin(qJ(1)) * pkin(1);
t122 = Icges(6,5) * t101 - Icges(6,6) * t99;
t169 = Icges(6,3) * t95 + t122 * t96;
t165 = 2 * m(4);
t164 = 2 * m(5);
t163 = 2 * m(6);
t146 = t107 * t95;
t71 = t96 * t156;
t159 = rSges(6,3) * t146 + t107 * t71;
t158 = t143 * qJD(1);
t110 = cos(pkin(9));
t157 = rSges(5,1) * t110;
t40 = -Icges(6,6) * t96 + t123 * t95;
t148 = t107 * t40;
t42 = -Icges(6,5) * t96 + t124 * t95;
t147 = t107 * t42;
t145 = t107 * t96;
t111 = -pkin(7) - qJ(4);
t82 = t96 * t111;
t141 = qJ(4) * t107;
t140 = qJD(5) * t95;
t139 = qJD(5) * t96;
t136 = t107 * t160;
t81 = t96 * t157;
t135 = t107 * t155;
t58 = t96 * rSges(4,1) - rSges(4,2) * t95;
t49 = rSges(4,1) * t145 - rSges(4,2) * t146;
t44 = -rSges(6,3) * t96 + t172 * t95;
t97 = pkin(4) * t110 + pkin(3);
t32 = t95 * t97 + t44 + t82;
t28 = t32 + t144;
t45 = -t95 * rSges(6,3) + t96 * t160 - t71;
t33 = -t111 * t95 + t96 * t97 - t45;
t29 = t33 + t143;
t132 = t28 * t96 - t29 * t95;
t131 = t32 * t96 - t33 * t95;
t67 = rSges(6,1) * t99 + rSges(6,2) * t101;
t130 = t101 * t42 - t40 * t99;
t63 = Icges(6,2) * t101 + t154;
t64 = Icges(6,1) * t99 + t142;
t128 = t101 * t64 - t63 * t99;
t62 = Icges(6,5) * t99 + Icges(6,6) * t101;
t121 = (t124 - t63) * t138 + (t123 + t64) * t137;
t48 = t57 * t107;
t120 = t130 * t96;
t119 = t144 * qJD(1);
t115 = -t122 * qJD(5) + t128 * t107;
t118 = -(t129 * qJD(5) + t101 * (t63 * t139 + t148) + (t64 * t139 + t147) * t99 + t115 * t95) * t95 / 0.2e1 - (t130 * qJD(5) + t101 * (-t63 * t140 + t175) + (-t64 * t140 + t176) * t99 + t115 * t96) * t96 / 0.2e1 + (t101 * t40 + t128 * t95 + t42 * t99 - t62 * t96) * t146 / 0.2e1 - (-t101 * t168 - t128 * t96 - t167 * t99 - t62 * t95) * t145 / 0.2e1;
t117 = t67 * qJD(5);
t35 = t173 * t95 + t174 * t96 + t81;
t34 = -t173 * t96 + (t157 + t174) * t95;
t38 = -Icges(6,3) * t96 + t122 * t95;
t114 = rSges(6,3) * t145 - t96 * t117 + t95 * t136;
t83 = qJD(4) * t95;
t26 = t95 * t135 + t96 * t141 + rSges(5,3) * t145 + t83 + (-pkin(3) - t157) * t146;
t27 = t107 * t81 + t95 * t141 + rSges(5,3) * t146 + pkin(3) * t145 + (-qJD(4) - t135) * t96;
t12 = t83 + (-t82 + (-t97 - t156) * t95) * t107 + t114;
t13 = t97 * t145 + (-qJD(4) - t136) * t96 + (-t107 * t111 - t117) * t95 + t159;
t56 = t172 * qJD(5);
t47 = t58 + t143;
t46 = t144 + t57;
t37 = t49 + t158;
t36 = -t48 - t119;
t31 = t35 + t143;
t30 = t34 + t144;
t21 = t169 * t107 - t62 * t140;
t20 = t107 * t38 + t62 * t139;
t19 = t27 + t158;
t18 = -t119 + t26;
t11 = t13 + t158;
t10 = -t119 + t12;
t9 = -t129 * t96 + t169 * t95;
t8 = -t38 * t95 - t120;
t7 = t169 * t96 + t171;
t6 = t130 * t95 - t38 * t96;
t3 = (t107 * t44 + t114) * t96 + (-t95 * t117 + (t45 + (-t156 - t160) * t96) * t107 + t159) * t95;
t1 = [(t36 * t47 + t37 * t46) * t165 + (t18 * t31 + t19 * t30) * t164 + (t10 * t29 + t11 * t28) * t163 + t121; 0; 0; m(4) * (t36 * t58 + t37 * t57 + t49 * t46 - t48 * t47) + m(5) * (t18 * t35 + t19 * t34 + t26 * t31 + t27 * t30) + m(6) * (t10 * t33 + t11 * t32 + t12 * t29 + t13 * t28) + t121; 0; (t12 * t33 + t13 * t32) * t163 + (t26 * t35 + t27 * t34) * t164 + (-t48 * t58 + t49 * t57) * t165 + t121; m(5) * (-t18 * t96 - t19 * t95 + (-t30 * t96 + t31 * t95) * t107) + m(6) * (-t10 * t96 - t132 * t107 - t11 * t95); 0; m(6) * (-t131 * t107 - t12 * t96 - t13 * t95) + m(5) * (-t26 * t96 - t27 * t95 + (-t34 * t96 + t35 * t95) * t107); 0; m(6) * (t132 * t56 + (-t10 * t95 + t11 * t96 + (-t28 * t95 - t29 * t96) * t107) * t67) + t118; m(6) * t3; m(6) * (t131 * t56 + (-t12 * t95 + t13 * t96 + (-t32 * t95 - t33 * t96) * t107) * t67) + t118; 0; ((t44 * t95 - t45 * t96) * t3 + (t95 ^ 2 + t96 ^ 2) * t67 * t56) * t163 + (-t6 * t96 - t7 * t95) * t146 - t96 * ((t96 * t21 + (-t7 - t120) * t107) * t96 + (t6 * t107 + (-t137 * t168 - t138 * t167) * t95 + (t20 + (qJD(5) * t42 - t175) * t99 + (qJD(5) * t40 + t176) * t101) * t96) * t95) - (-t8 * t96 - t9 * t95) * t145 - t95 * ((t95 * t20 + (t8 - t171) * t107) * t95 + (-t9 * t107 + (-t40 * t137 - t42 * t138) * t96 + (t21 + (qJD(5) * t167 + t148) * t99 + (qJD(5) * t168 - t147) * t101) * t95) * t96);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
