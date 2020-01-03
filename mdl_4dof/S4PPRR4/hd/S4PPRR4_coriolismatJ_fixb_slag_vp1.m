% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:27
% EndTime: 2019-12-31 16:18:33
% DurationCPUTime: 2.84s
% Computational Cost: add. (7268->255), mult. (9903->420), div. (0->0), fcn. (10859->6), ass. (0->160)
t122 = cos(pkin(6));
t195 = t122 ^ 2;
t121 = sin(pkin(6));
t196 = t121 ^ 2;
t202 = t195 + t196;
t120 = pkin(7) + qJ(3);
t118 = sin(t120);
t119 = cos(t120);
t209 = Icges(4,5) * t118 + Icges(4,6) * t119;
t123 = sin(qJ(4));
t124 = cos(qJ(4));
t148 = rSges(5,1) * t124 - rSges(5,2) * t123;
t91 = -t119 * rSges(5,3) + t118 * t148;
t83 = t121 * t91;
t84 = t122 * t91;
t203 = -m(5) / 0.2e1;
t167 = Icges(5,4) * t123;
t140 = Icges(5,1) * t124 - t167;
t89 = -Icges(5,5) * t119 + t118 * t140;
t175 = (-Icges(5,2) * t124 - t167) * t118 + t89;
t166 = Icges(5,4) * t124;
t137 = -Icges(5,2) * t123 + t166;
t87 = -Icges(5,6) * t119 + t118 * t137;
t174 = (-Icges(5,1) * t123 - t166) * t118 - t87;
t159 = t119 * (-Icges(5,5) * t123 - Icges(5,6) * t124) * t118;
t135 = Icges(5,5) * t124 - Icges(5,6) * t123;
t85 = -Icges(5,3) * t119 + t118 * t135;
t200 = t209 * t202;
t197 = 0.2e1 * m(5);
t194 = m(5) / 0.2e1;
t193 = -t119 / 0.2e1;
t192 = t121 / 0.2e1;
t191 = -t122 / 0.2e1;
t190 = t122 / 0.2e1;
t189 = qJD(3) / 0.2e1;
t155 = t122 * t124;
t158 = t121 * t123;
t109 = -t119 * t158 - t155;
t156 = t122 * t123;
t157 = t121 * t124;
t110 = t119 * t157 - t156;
t169 = Icges(5,4) * t110;
t161 = t118 * t121;
t58 = Icges(5,2) * t109 + Icges(5,6) * t161 + t169;
t188 = Icges(5,1) * t109 - t169 - t58;
t111 = -t119 * t156 + t157;
t112 = t119 * t155 + t158;
t168 = Icges(5,4) * t112;
t160 = t118 * t122;
t59 = Icges(5,2) * t111 + Icges(5,6) * t160 + t168;
t187 = Icges(5,1) * t111 - t168 - t59;
t103 = Icges(5,4) * t109;
t60 = Icges(5,1) * t110 + Icges(5,5) * t161 + t103;
t186 = -Icges(5,2) * t110 + t103 + t60;
t104 = Icges(5,4) * t111;
t61 = Icges(5,1) * t112 + Icges(5,5) * t160 + t104;
t185 = -Icges(5,2) * t112 + t104 + t61;
t183 = m(5) * qJD(4);
t56 = Icges(5,5) * t110 + Icges(5,6) * t109 + Icges(5,3) * t161;
t181 = t119 * t56;
t57 = Icges(5,5) * t112 + Icges(5,6) * t111 + Icges(5,3) * t160;
t180 = t119 * t57;
t179 = t119 * t85;
t36 = t109 * t59 + t110 * t61 + t161 * t57;
t178 = t36 * t122;
t37 = t111 * t58 + t112 * t60 + t160 * t56;
t177 = t37 * t121;
t115 = pkin(3) * t118 - pkin(5) * t119;
t173 = -t115 - t91;
t116 = pkin(3) * t119 + pkin(5) * t118;
t92 = t118 * rSges(5,3) + t119 * t148;
t172 = -t116 - t92;
t62 = rSges(5,1) * t110 + rSges(5,2) * t109 + rSges(5,3) * t161;
t41 = (t121 * t92 - t62) * t118;
t63 = rSges(5,1) * t112 + rSges(5,2) * t111 + rSges(5,3) * t160;
t42 = (-t122 * t92 + t63) * t118;
t18 = t121 * t41 - t122 * t42;
t154 = t18 * qJD(2);
t153 = t18 * qJD(4);
t146 = -t121 * t63 + t122 * t62;
t34 = t146 * t119 + (t121 * t84 - t122 * t83) * t118;
t75 = rSges(5,1) * t109 - rSges(5,2) * t110;
t76 = rSges(5,1) * t111 - rSges(5,2) * t112;
t53 = t121 * t75 + t122 * t76;
t27 = (t34 / 0.4e1 - t53 / 0.4e1) * t197;
t152 = t27 * qJD(1);
t151 = qJD(4) * t118;
t113 = rSges(4,1) * t118 + rSges(4,2) * t119;
t145 = -t123 * t58 + t124 * t60;
t39 = t118 * t145 - t181;
t144 = -t123 * t59 + t124 * t61;
t40 = t118 * t144 - t180;
t147 = t39 * t121 + t40 * t122;
t143 = -t123 * t87 + t124 * t89;
t134 = -t121 * t85 - t145;
t133 = -t122 * t85 - t144;
t69 = Icges(5,5) * t109 - Icges(5,6) * t110;
t23 = t109 * t186 + t110 * t188 + t161 * t69;
t70 = Icges(5,5) * t111 - Icges(5,6) * t112;
t24 = t109 * t185 + t110 * t187 + t161 * t70;
t11 = t121 * t24 - t122 * t23;
t25 = t111 * t186 + t112 * t188 + t160 * t69;
t26 = t111 * t185 + t112 * t187 + t160 * t70;
t12 = t121 * t26 - t122 * t25;
t132 = t11 * t191 + t12 * t192;
t131 = (Icges(5,3) * t118 + t119 * t135 - t143) * t119;
t126 = t118 * t134 + t181;
t125 = t118 * t133 + t180;
t108 = (-rSges(5,1) * t123 - rSges(5,2) * t124) * t118;
t98 = t209 * t122;
t97 = t209 * t121;
t90 = Icges(5,5) * t118 + t119 * t140;
t88 = Icges(5,6) * t118 + t119 * t137;
t82 = t89 * t122;
t81 = t89 * t121;
t80 = t87 * t122;
t79 = t87 * t121;
t68 = t172 * t122;
t67 = t173 * t122;
t66 = t172 * t121;
t65 = t173 * t121;
t64 = t202 * t113;
t55 = -t108 * t160 - t119 * t76;
t54 = t108 * t161 + t119 * t75;
t52 = -t119 * t63 - t160 * t91;
t51 = t119 * t62 + t161 * t91;
t49 = (-t121 * t76 + t122 * t75) * t118;
t48 = t118 * t143 - t179;
t47 = t146 * t118;
t46 = (-t115 * t122 - t84) * t122 + (-t115 * t121 - t83) * t121;
t45 = t111 * t87 + t112 * t89 + t160 * t85;
t44 = t109 * t87 + t110 * t89 + t161 * t85;
t43 = (t116 * t122 + t63) * t122 + (t116 * t121 + t62) * t121;
t38 = t111 * t59 + t112 * t61 + t160 * t57;
t35 = t109 * t58 + t110 * t60 + t161 * t56;
t32 = -t119 * t70 + (-t123 * t185 + t124 * t187) * t118;
t31 = -t119 * t69 + (-t123 * t186 + t124 * t188) * t118;
t30 = -t133 * t119 + (t123 * t80 - t124 * t82 + t57) * t118;
t29 = -t134 * t119 + (t123 * t79 - t124 * t81 + t56) * t118;
t28 = (t34 + t53) * t194;
t22 = -t111 * t80 - t112 * t82 + t122 * t125;
t21 = -t111 * t79 - t112 * t81 + t122 * t126;
t20 = -t109 * t80 - t110 * t82 + t121 * t125;
t19 = -t109 * t79 - t110 * t81 + t121 * t126;
t17 = m(5) * t18 * t189;
t16 = t43 * t53 + (-t121 * t65 - t122 * t67) * t108;
t15 = t118 * t147 - t119 * t48;
t14 = -t119 * t45 + (t38 * t122 + t177) * t118;
t13 = -t119 * t44 + (t35 * t121 + t178) * t118;
t10 = t121 * t22 - t122 * t21;
t9 = t121 * t20 - t122 * t19;
t8 = -(t111 * t175 + t112 * t174) * t119 + (t25 * t121 + (t26 - t159) * t122) * t118;
t7 = -(t109 * t175 + t110 * t174) * t119 + (t24 * t122 + (t23 - t159) * t121) * t118;
t6 = t34 * t47 + t41 * t51 + t42 * t52;
t5 = (t131 + t147) * t119 + (t30 * t122 + t29 * t121 - (-t123 * t88 + t124 * t90 + t85) * t119 + t48) * t118;
t4 = (-t111 * t88 - t112 * t90 + t177 + (t38 - t179) * t122) * t119 + (t21 * t121 + t45 + (t22 - t131) * t122) * t118;
t3 = (-t109 * t88 - t110 * t90 + t178 + (t35 - t179) * t121) * t119 + (t20 * t122 + t44 + (t19 - t131) * t121) * t118;
t2 = m(5) * t16 + t132;
t1 = m(5) * t6 + (t14 * t190 + t13 * t192 - t5 / 0.2e1) * t119 + (t4 * t190 + t3 * t192 + t15 / 0.2e1) * t118;
t33 = [0, 0, t28 * qJD(4) + 0.2e1 * (-m(4) * t64 / 0.2e1 + t46 * t194) * qJD(3), t28 * qJD(3) + t183 * t49; 0, 0, ((t121 * t68 - t122 * t66) * t189 + t153 / 0.4e1) * t197, t17 + (t121 * t54 - t122 * t55) * t183; -t27 * qJD(4), t153 * t203, (m(5) * (t43 * t46 + t65 * t66 + t67 * t68) + m(4) * (t113 - t64) * t202 * (rSges(4,1) * t119 - rSges(4,2) * t118) + (-t196 * t98 + (t121 * t97 - t200) * t122 + t10) * t192 + (-t195 * t97 + (t122 * t98 - t200) * t121 + t9) * t191) * qJD(3) + t2 * qJD(4), -t152 + t2 * qJD(3) + t154 * t203 + (-t15 / 0.2e1 + (t12 / 0.2e1 - t4 / 0.2e1) * t122 + (t11 / 0.2e1 - t3 / 0.2e1) * t121) * t151 + (t191 * t7 + t192 * t8 + (t5 / 0.2e1 + (t31 / 0.2e1 - t14 / 0.2e1) * t122 + (-t32 / 0.2e1 - t13 / 0.2e1) * t121) * t119 + (t49 * t43 + t47 * t53 + t54 * t67 + t55 * t65 + (-t121 * t52 - t122 * t51) * t108 - t6) * m(5)) * qJD(4); t27 * qJD(3), t17, t152 + (t4 * t192 + t10 * t160 / 0.2e1 + t3 * t191 + t9 * t161 / 0.2e1 + t118 * (t121 * t40 - t122 * t39) / 0.2e1 + (t121 * t30 - t122 * t29) * t193 - t132 + ((t121 * t38 - t122 * t37) * t190 + (t121 * t36 - t122 * t35) * t192) * t119) * qJD(3) + t1 * qJD(4) + (t154 / 0.2e1 + (t34 * t43 + t41 * t67 + t42 * t65 + t46 * t47 + t51 * t68 + t52 * t66 - t16) * qJD(3)) * m(5), t1 * qJD(3) + (m(5) * (t47 * t49 + t51 * t54 + t52 * t55) - t119 ^ 2 * t159 / 0.2e1) * qJD(4) + (t8 * t190 + t7 * t192 + (t32 * t122 + t31 * t121 - (-t175 * t123 + t174 * t124) * t119) * t193) * t151;];
Cq = t33;
