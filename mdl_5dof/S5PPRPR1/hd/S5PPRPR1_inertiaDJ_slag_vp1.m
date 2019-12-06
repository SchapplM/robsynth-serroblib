% Calculate time derivative of joint inertia matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:53
% EndTime: 2019-12-05 15:00:59
% DurationCPUTime: 2.63s
% Computational Cost: add. (6374->275), mult. (6983->466), div. (0->0), fcn. (6517->8), ass. (0->154)
t118 = pkin(8) + qJ(3);
t112 = sin(t118);
t114 = cos(t118);
t120 = sin(pkin(7));
t115 = t120 ^ 2;
t122 = cos(pkin(7));
t116 = t122 ^ 2;
t203 = t115 + t116;
t201 = t203 * qJD(3);
t207 = m(4) * (rSges(4,1) * t112 + rSges(4,2) * t114) * t201;
t206 = 0.1e1 - t203;
t119 = sin(pkin(9));
t121 = cos(pkin(9));
t156 = rSges(5,1) * t121 - rSges(5,2) * t119;
t202 = t114 * rSges(5,3) - t112 * t156;
t123 = -pkin(6) - qJ(4);
t169 = qJ(4) + t123;
t110 = pkin(4) * t121 + pkin(3);
t192 = pkin(3) - t110;
t200 = t112 * t192 - t114 * t169;
t197 = 2 * m(6);
t196 = m(5) / 0.2e1;
t195 = m(6) / 0.2e1;
t194 = t120 / 0.2e1;
t193 = t122 / 0.2e1;
t108 = t112 * pkin(3) - t114 * qJ(4);
t191 = t203 * (-qJD(3) * t108 + qJD(4) * t112);
t154 = pkin(3) * t114 + qJ(4) * t112;
t90 = qJD(3) * t154 - qJD(4) * t114;
t190 = -(rSges(5,3) * t112 + t114 * t156) * qJD(3) - t90;
t189 = t203 * t154;
t117 = pkin(9) + qJ(5);
t111 = sin(t117);
t113 = cos(t117);
t140 = Icges(6,5) * t113 - Icges(6,6) * t111;
t164 = qJD(5) * t112;
t187 = t114 * ((-Icges(6,5) * t111 - Icges(6,6) * t113) * t164 + (Icges(6,3) * t112 + t114 * t140) * qJD(3));
t80 = -Icges(6,3) * t114 + t112 * t140;
t186 = t114 * t80;
t174 = t114 * t122;
t101 = -t111 * t174 + t113 * t120;
t102 = t111 * t120 + t113 * t174;
t176 = t112 * t122;
t175 = t114 * t120;
t100 = -t111 * t122 + t113 * t175;
t177 = t112 * t120;
t99 = -t111 * t175 - t113 * t122;
t48 = Icges(6,5) * t100 + Icges(6,6) * t99 + Icges(6,3) * t177;
t50 = Icges(6,4) * t100 + Icges(6,2) * t99 + Icges(6,6) * t177;
t52 = Icges(6,1) * t100 + Icges(6,4) * t99 + Icges(6,5) * t177;
t18 = t101 * t50 + t102 * t52 + t176 * t48;
t185 = t120 * t18;
t49 = Icges(6,5) * t102 + Icges(6,6) * t101 + Icges(6,3) * t176;
t51 = Icges(6,4) * t102 + Icges(6,2) * t101 + Icges(6,6) * t176;
t53 = Icges(6,1) * t102 + Icges(6,4) * t101 + Icges(6,5) * t176;
t17 = t100 * t53 + t177 * t49 + t51 * t99;
t184 = t122 * t17;
t183 = -t108 + t202;
t180 = Icges(6,4) * t111;
t179 = Icges(6,4) * t113;
t173 = t119 * t120;
t172 = t119 * t122;
t171 = t120 * t121;
t170 = t121 * t122;
t168 = qJD(3) * t112;
t167 = qJD(3) * t114;
t166 = qJD(3) * t120;
t165 = qJD(3) * t122;
t155 = rSges(6,1) * t113 - rSges(6,2) * t111;
t47 = (-rSges(6,1) * t111 - rSges(6,2) * t113) * t164 + (rSges(6,3) * t112 + t114 * t155) * qJD(3);
t163 = -t47 - (-t112 * t169 - t114 * t192) * qJD(3) - t90;
t84 = -t114 * rSges(6,3) + t112 * t155;
t162 = -t108 + t200 - t84;
t161 = t112 * t166;
t160 = t112 * t165;
t159 = t114 * t166;
t158 = t114 * t165;
t153 = -t111 * t50 + t113 * t52;
t152 = -t111 * t51 + t113 * t53;
t141 = -Icges(6,2) * t111 + t179;
t81 = -Icges(6,6) * t114 + t112 * t141;
t143 = Icges(6,1) * t113 - t180;
t82 = -Icges(6,5) * t114 + t112 * t143;
t151 = t111 * t81 - t113 * t82;
t20 = t112 * t153 - t114 * t48;
t21 = t112 * t152 - t114 * t49;
t146 = t20 * t120 + t21 * t122;
t54 = rSges(6,1) * t100 + rSges(6,2) * t99 + rSges(6,3) * t177;
t55 = rSges(6,1) * t102 + rSges(6,2) * t101 + rSges(6,3) * t176;
t145 = -t120 * t55 + t122 * t54;
t66 = -qJD(5) * t100 + t111 * t161;
t67 = qJD(5) * t99 - t113 * t161;
t32 = Icges(6,5) * t67 + Icges(6,6) * t66 + Icges(6,3) * t159;
t139 = t112 * t32 + t167 * t48;
t68 = -qJD(5) * t102 + t111 * t160;
t69 = qJD(5) * t101 - t113 * t160;
t33 = Icges(6,5) * t69 + Icges(6,6) * t68 + Icges(6,3) * t158;
t138 = t112 * t33 + t167 * t49;
t131 = qJD(3) * (-Icges(4,5) * t112 - Icges(4,6) * t114);
t127 = t110 * t114 - t112 * t123 - t154;
t126 = qJD(3) * (Icges(5,5) * t114 + (-Icges(5,1) * t121 + Icges(5,4) * t119) * t112);
t125 = qJD(3) * (Icges(5,6) * t114 + (-Icges(5,4) * t121 + Icges(5,2) * t119) * t112);
t106 = t114 * t170 + t173;
t105 = -t114 * t172 + t171;
t104 = t114 * t171 - t172;
t103 = -t114 * t173 - t170;
t94 = t122 * t131;
t93 = t120 * t131;
t75 = t122 * t126;
t74 = t120 * t126;
t73 = t122 * t125;
t72 = t120 * t125;
t64 = t183 * t122;
t63 = t183 * t120;
t46 = (-Icges(6,1) * t111 - t179) * t164 + (Icges(6,5) * t112 + t114 * t143) * qJD(3);
t45 = (-Icges(6,2) * t113 - t180) * t164 + (Icges(6,6) * t112 + t114 * t141) * qJD(3);
t43 = t190 * t122;
t42 = t190 * t120;
t41 = t162 * t122;
t40 = t162 * t120;
t39 = rSges(6,1) * t69 + rSges(6,2) * t68 + rSges(6,3) * t158;
t38 = rSges(6,1) * t67 + rSges(6,2) * t66 + rSges(6,3) * t159;
t37 = Icges(6,1) * t69 + Icges(6,4) * t68 + Icges(6,5) * t158;
t36 = Icges(6,1) * t67 + Icges(6,4) * t66 + Icges(6,5) * t159;
t35 = Icges(6,4) * t69 + Icges(6,2) * t68 + Icges(6,6) * t158;
t34 = Icges(6,4) * t67 + Icges(6,2) * t66 + Icges(6,6) * t159;
t31 = -t114 * t55 - t176 * t84;
t30 = t114 * t54 + t177 * t84;
t29 = -t112 * t151 - t186;
t28 = t163 * t122;
t27 = t163 * t120;
t26 = t145 * t112;
t25 = t101 * t81 + t102 * t82 + t176 * t80;
t24 = t100 * t82 + t177 * t80 + t81 * t99;
t23 = t120 * (rSges(5,1) * t104 + rSges(5,2) * t103 + rSges(5,3) * t177) + t122 * (rSges(5,1) * t106 + rSges(5,2) * t105 + rSges(5,3) * t176) + t189;
t22 = t202 * t201 + t191;
t19 = t101 * t51 + t102 * t53 + t176 * t49;
t16 = t100 * t52 + t177 * t48 + t50 * t99;
t15 = -t47 * t176 - t114 * t39 + (t112 * t55 - t174 * t84) * qJD(3);
t14 = t47 * t177 + t114 * t38 + (-t112 * t54 + t175 * t84) * qJD(3);
t13 = (t122 * t127 + t55) * t122 + (t120 * t127 + t54) * t120 + t189;
t12 = (-t120 * t39 + t122 * t38) * t112 + t145 * t167;
t11 = t120 * t38 + t122 * t39 + t200 * t201 + t191;
t10 = t101 * t35 + t102 * t37 + t122 * t138 + t51 * t68 + t53 * t69;
t9 = t101 * t34 + t102 * t36 + t122 * t139 + t50 * t68 + t52 * t69;
t8 = t100 * t37 + t120 * t138 + t35 * t99 + t51 * t66 + t53 * t67;
t7 = t100 * t36 + t120 * t139 + t34 * t99 + t50 * t66 + t52 * t67;
t6 = (qJD(3) * t152 - t33) * t114 + (qJD(3) * t49 - t111 * t35 + t113 * t37 + (-t111 * t53 - t113 * t51) * qJD(5)) * t112;
t5 = (qJD(3) * t153 - t32) * t114 + (qJD(3) * t48 - t111 * t34 + t113 * t36 + (-t111 * t52 - t113 * t50) * qJD(5)) * t112;
t4 = t10 * t120 - t122 * t9;
t3 = t120 * t8 - t122 * t7;
t2 = -(t101 * t45 + t102 * t46 + t68 * t81 + t69 * t82) * t114 + (t9 * t120 + (t10 - t187) * t122) * t112 + (t25 * t112 + (t185 + (t19 - t186) * t122) * t114) * qJD(3);
t1 = -(t100 * t46 + t45 * t99 + t66 * t81 + t67 * t82) * t114 + (t8 * t122 + (t7 - t187) * t120) * t112 + (t24 * t112 + (t184 + (t16 - t186) * t120) * t114) * qJD(3);
t44 = [0; 0; 0; m(5) * t22 + m(6) * t11 - t207; m(5) * (t120 * t43 - t122 * t42) + m(6) * (t120 * t28 - t122 * t27); (t11 * t13 + t27 * t40 + t28 * t41) * t197 + 0.2e1 * m(5) * (t22 * t23 + t42 * t63 + t43 * t64) + (-t116 * t93 - t3 + (t103 * t72 + t104 * t74) * t122) * t122 + (t4 + t115 * t94 + (t105 * t73 + t106 * t75) * t120 + (-t103 * t73 - t104 * t75 - t105 * t72 - t106 * t74 - t120 * t93 + t122 * t94) * t122) * t120 + 0.2e1 * t206 * (rSges(4,1) * t114 - rSges(4,2) * t112) * t207; (m(5) + m(6)) * t168; 0; m(6) * (-t11 * t114 + t176 * t28 + t177 * t27) + m(5) * (-t114 * t22 + t176 * t43 + t177 * t42) + 0.2e1 * ((t112 * t13 + t174 * t41 + t175 * t40) * t195 + (t112 * t23 + t174 * t64 + t175 * t63) * t196) * qJD(3); -0.4e1 * (t196 + t195) * t206 * t112 * t167; m(6) * t12; m(6) * (t120 * t14 - t122 * t15); -t122 * t1 / 0.2e1 + m(6) * (t11 * t26 + t12 * t13 + t14 * t41 + t15 * t40 + t27 * t31 + t28 * t30) - t114 * (t6 * t120 - t5 * t122) / 0.2e1 + t2 * t194 + (t193 * t4 + t194 * t3) * t112 + (t112 * (t120 * t21 - t122 * t20) / 0.2e1 + ((t19 * t120 - t18 * t122) * t193 + (t17 * t120 - t16 * t122) * t194) * t114) * qJD(3); m(6) * (-t114 * t12 + (t120 * t15 + t122 * t14) * t112 + (t112 * t26 + (t120 * t31 + t122 * t30) * t114) * qJD(3)); (t12 * t26 + t14 * t30 + t15 * t31) * t197 + (-t25 * t114 + (t122 * t19 + t185) * t112) * t158 + t2 * t176 + (-t24 * t114 + (t120 * t16 + t184) * t112) * t159 + t1 * t177 + (t112 * t146 - t29 * t114) * t168 - t114 * ((t187 + (t114 * t151 + t146) * qJD(3)) * t114 + (t6 * t122 + t5 * t120 - (qJD(3) * t80 - t111 * t45 + t113 * t46 + (-t111 * t82 - t113 * t81) * qJD(5)) * t114 + t29 * qJD(3)) * t112);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t44(1), t44(2), t44(4), t44(7), t44(11); t44(2), t44(3), t44(5), t44(8), t44(12); t44(4), t44(5), t44(6), t44(9), t44(13); t44(7), t44(8), t44(9), t44(10), t44(14); t44(11), t44(12), t44(13), t44(14), t44(15);];
Mq = res;
