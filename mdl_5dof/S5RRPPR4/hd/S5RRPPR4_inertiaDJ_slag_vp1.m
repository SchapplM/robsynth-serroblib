% Calculate time derivative of joint inertia matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:40
% EndTime: 2019-12-31 19:27:43
% DurationCPUTime: 1.67s
% Computational Cost: add. (4309->215), mult. (4274->303), div. (0->0), fcn. (4012->8), ass. (0->116)
t110 = sin(qJ(5));
t156 = qJD(5) * t110;
t158 = qJ(1) + qJ(2);
t146 = cos(t158);
t145 = sin(t158);
t157 = t146 * pkin(2) + t145 * qJ(3);
t149 = t146 * pkin(3) + t157;
t112 = cos(qJ(5));
t155 = qJD(5) * t112;
t109 = qJD(1) + qJD(2);
t161 = sin(pkin(8));
t162 = cos(pkin(8));
t76 = -t145 * t161 - t146 * t162;
t61 = t76 * t109;
t77 = -t145 * t162 + t146 * t161;
t128 = t110 * t61 + t77 * t155;
t53 = t146 * rSges(4,1) + t145 * rSges(4,3) + t157;
t168 = rSges(6,2) * t110;
t169 = rSges(6,1) * t112;
t188 = t168 - t169;
t152 = t76 * t156;
t62 = t77 * t109;
t163 = t112 * t62;
t125 = t152 + t163;
t165 = t110 * t62;
t126 = t76 * t155 - t165;
t151 = t77 * t156;
t164 = t112 * t61;
t127 = t151 - t164;
t136 = rSges(6,1) * t110 + rSges(6,2) * t112;
t178 = 2 * m(6);
t85 = t188 * qJD(5);
t187 = -(t127 * Icges(6,5) + t128 * Icges(6,6) + Icges(6,3) * t62) * t76 + (t125 * Icges(6,5) + t126 * Icges(6,6) + Icges(6,3) * t61) * t77 - t136 * t85 * t178;
t159 = Icges(6,4) * t112;
t130 = Icges(6,2) * t110 - t159;
t34 = -Icges(6,6) * t76 + t130 * t77;
t160 = Icges(6,4) * t110;
t131 = -Icges(6,1) * t112 + t160;
t36 = -Icges(6,5) * t76 + t131 * t77;
t135 = t110 * t34 - t112 * t36;
t184 = t135 * t76;
t35 = Icges(6,6) * t77 + t130 * t76;
t37 = Icges(6,5) * t77 + t131 * t76;
t134 = t110 * t35 - t112 * t37;
t185 = t134 * t77;
t129 = -Icges(6,5) * t112 + Icges(6,6) * t110;
t32 = -Icges(6,3) * t76 + t129 * t77;
t33 = Icges(6,3) * t77 + t129 * t76;
t186 = t32 * t77 - t33 * t76 + t184 + t185;
t181 = 2 * m(3);
t180 = 2 * m(4);
t179 = 2 * m(5);
t111 = sin(qJ(1));
t173 = t111 * pkin(1);
t171 = rSges(6,1) * t163 + t61 * rSges(6,3);
t101 = t146 * qJ(3);
t170 = qJD(3) * t145 + t109 * t101;
t167 = pkin(1) * qJD(1);
t154 = t111 * t167;
t113 = cos(qJ(1));
t153 = t113 * t167;
t148 = pkin(4) + t169;
t147 = t76 * rSges(6,3) - t77 * t168;
t144 = (-Icges(6,2) * t112 - t131 - t160) * t156;
t141 = t145 * pkin(2);
t140 = t145 * rSges(3,2);
t139 = t109 * t146;
t43 = -t76 * rSges(5,1) - t77 * rSges(5,2) + t149;
t132 = -rSges(6,1) * t151 - rSges(6,2) * t128 - t62 * rSges(6,3);
t39 = t77 * rSges(6,3) + t188 * t76;
t81 = t146 * rSges(3,1) - t140;
t82 = t129 * qJD(5);
t122 = (-t110 * t37 - t112 * t35) * t61 / 0.2e1 + (-t110 * t36 - t112 * t34) * t62 / 0.2e1 - (t135 * qJD(5) - t110 * (t127 * Icges(6,1) + t128 * Icges(6,4) + Icges(6,5) * t62) - t112 * (t127 * Icges(6,4) + t128 * Icges(6,2) + Icges(6,6) * t62) - t76 * t82) * t76 / 0.2e1 + (t134 * qJD(5) - t110 * (t125 * Icges(6,1) + t126 * Icges(6,4) + Icges(6,5) * t61) - t112 * (t125 * Icges(6,4) + t126 * Icges(6,2) + Icges(6,6) * t61) + t77 * t82) * t77 / 0.2e1 + (t61 * t77 - t76 * t62) * (-Icges(6,5) * t110 - Icges(6,6) * t112);
t68 = -rSges(3,1) * t139 + t109 * t140;
t83 = t130 * qJD(5);
t91 = -Icges(6,1) * t110 - t159;
t121 = (-qJD(5) * t91 - t83) * t112 + t144;
t120 = -t145 * pkin(3) - t141;
t119 = -t145 * rSges(4,1) - t141;
t80 = -t145 * rSges(3,1) - t146 * rSges(3,2);
t117 = t101 + t120;
t25 = -t76 * pkin(4) + pkin(7) * t77 + t149 + t39;
t67 = t80 * t109;
t52 = t146 * rSges(4,3) + t101 + t119;
t115 = t120 * t109 + t170;
t46 = rSges(4,3) * t139 + t119 * t109 + t170;
t30 = t62 * rSges(5,1) - t61 * rSges(5,2) + t115;
t42 = t77 * rSges(5,1) - t76 * rSges(5,2) + t117;
t99 = qJD(3) * t146;
t114 = -t109 * t149 + t99;
t24 = t76 * pkin(7) + t148 * t77 + t117 + t147;
t31 = t61 * rSges(5,1) + t62 * rSges(5,2) + t114;
t47 = -t109 * t53 + t99;
t13 = -t62 * pkin(7) + t148 * t61 + t114 + t132;
t12 = rSges(6,1) * t152 + t126 * rSges(6,2) + t62 * pkin(4) + t61 * pkin(7) + t115 + t171;
t108 = t113 * pkin(1);
t87 = t136 ^ 2;
t75 = t108 + t81;
t74 = t80 - t173;
t55 = t68 - t153;
t54 = t67 - t154;
t49 = t108 + t53;
t48 = t52 - t173;
t45 = t47 - t153;
t44 = t46 - t154;
t41 = t108 + t43;
t40 = t42 - t173;
t38 = -t77 * t169 - t147;
t29 = t31 - t153;
t28 = t30 - t154;
t23 = t108 + t25;
t22 = t24 - t173;
t11 = t13 - t153;
t10 = t12 - t154;
t1 = t61 * t38 + t77 * (-rSges(6,1) * t164 - t132) - t62 * t39 + (t136 * t76 * qJD(5) - rSges(6,2) * t165 + t171) * t76;
t2 = [(t10 * t23 + t11 * t22) * t178 - t91 * t155 - t112 * t83 + (t44 * t49 + t45 * t48) * t180 + (t28 * t41 + t29 * t40) * t179 + (t54 * t75 + t55 * t74) * t181 + t144; m(6) * (t10 * t25 + t11 * t24 + t12 * t23 + t13 * t22) + m(4) * (t44 * t53 + t45 * t52 + t46 * t49 + t47 * t48) + m(5) * (t28 * t43 + t29 * t42 + t30 * t41 + t31 * t40) + m(3) * (t54 * t81 + t55 * t80 + t67 * t75 + t68 * t74) + t121; (t12 * t25 + t13 * t24) * t178 + (t46 * t53 + t47 * t52) * t180 + (t30 * t43 + t31 * t42) * t179 + (t67 * t81 + t68 * t80) * t181 + t121; m(6) * (t145 * t11 - t146 * t10 + (t145 * t23 + t146 * t22) * t109) + m(4) * (t145 * t45 - t146 * t44 + (t145 * t49 + t146 * t48) * t109) + m(5) * (t145 * t29 - t146 * t28 + (t145 * t41 + t146 * t40) * t109); m(6) * (t145 * t13 - t146 * t12 + (t145 * t25 + t146 * t24) * t109) + m(4) * (t145 * t47 - t146 * t46 + (t145 * t53 + t146 * t52) * t109) + m(5) * (t145 * t31 - t146 * t30 + (t145 * t43 + t146 * t42) * t109); 0; 0; 0; 0; 0; m(6) * ((-t22 * t76 - t23 * t77) * t85 - (-t10 * t77 - t11 * t76 + t22 * t62 - t23 * t61) * t136) + t122; m(6) * ((-t24 * t76 - t25 * t77) * t85 - (-t12 * t77 - t13 * t76 + t24 * t62 - t25 * t61) * t136) + t122; m(6) * ((-t145 * t76 + t146 * t77) * t85 - (t145 * t62 + t146 * t61 + (-t145 * t77 - t146 * t76) * t109) * t136); -m(6) * t1; ((t38 * t1 + t87 * t61) * t178 + (-t185 + t186) * t62 + (0.3e1 * t33 * t61 + t187) * t77) * t77 + ((-t135 - t33) * t62 * t77 + (0.3e1 * t32 * t62 + t187) * t76 + (t184 - t186 + (-t32 + t134) * t77) * t61 + (t39 * t1 - t87 * t62) * t178) * t76;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
