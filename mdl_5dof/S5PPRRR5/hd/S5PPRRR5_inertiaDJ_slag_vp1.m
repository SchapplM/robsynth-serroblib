% Calculate time derivative of joint inertia matrix for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:38
% EndTime: 2019-12-31 17:35:41
% DurationCPUTime: 1.33s
% Computational Cost: add. (3357->159), mult. (3694->253), div. (0->0), fcn. (3622->8), ass. (0->96)
t91 = cos(qJ(5));
t134 = qJD(5) * t91;
t156 = -qJD(5) / 0.2e1;
t89 = sin(qJ(5));
t135 = qJD(5) * t89;
t110 = -Icges(6,5) * t91 + Icges(6,6) * t89;
t155 = t110 * t156;
t92 = cos(qJ(3));
t154 = pkin(3) * t92;
t133 = qJ(3) + qJ(4);
t125 = sin(t133);
t132 = qJD(3) + qJD(4);
t100 = t132 * t125;
t126 = cos(t133);
t101 = t132 * t126;
t136 = sin(pkin(8));
t137 = cos(pkin(8));
t53 = t137 * t100 - t136 * t101;
t64 = -t136 * t125 - t137 * t126;
t108 = t64 * t134 + t53 * t89;
t146 = rSges(6,2) * t89;
t147 = rSges(6,1) * t91;
t153 = t146 - t147;
t65 = t137 * t125 - t136 * t126;
t129 = t65 * t135;
t54 = t136 * t100 + t137 * t101;
t142 = t54 * t91;
t105 = t129 - t142;
t143 = t54 * t89;
t106 = t65 * t134 + t143;
t131 = t64 * t135;
t144 = t53 * t91;
t107 = t131 - t144;
t118 = rSges(6,1) * t89 + rSges(6,2) * t91;
t148 = 2 * m(6);
t75 = t153 * qJD(5);
t152 = -t64 * (t105 * Icges(6,5) + t106 * Icges(6,6) - Icges(6,3) * t53) - t118 * t75 * t148 + t65 * (t107 * Icges(6,5) + t108 * Icges(6,6) + Icges(6,3) * t54);
t138 = Icges(6,4) * t91;
t111 = Icges(6,2) * t89 - t138;
t33 = Icges(6,6) * t65 + t111 * t64;
t139 = Icges(6,4) * t89;
t112 = -Icges(6,1) * t91 + t139;
t35 = Icges(6,5) * t65 + t112 * t64;
t116 = t33 * t89 - t35 * t91;
t103 = t116 * t65;
t32 = -Icges(6,6) * t64 + t111 * t65;
t34 = -Icges(6,5) * t64 + t112 * t65;
t117 = t32 * t89 - t34 * t91;
t104 = t117 * t64;
t30 = -Icges(6,3) * t64 + t110 * t65;
t31 = Icges(6,3) * t65 + t110 * t64;
t151 = -t30 * t65 + t31 * t64 - t103 - t104;
t149 = 2 * m(5);
t140 = -rSges(6,1) * t142 - t53 * rSges(6,3);
t41 = -t54 * rSges(5,1) + t53 * rSges(5,2);
t42 = -t65 * rSges(5,1) + t64 * rSges(5,2);
t128 = pkin(4) + t147;
t127 = -t65 * rSges(6,3) - t64 * t146;
t90 = sin(qJ(3));
t124 = t137 * t90;
t123 = t136 * t90;
t122 = t136 * t92;
t119 = pkin(3) * t124;
t40 = t53 * rSges(5,1) + t54 * rSges(5,2);
t43 = t64 * rSges(5,1) + t65 * rSges(5,2);
t36 = -t64 * rSges(6,3) + t153 * t65;
t109 = -rSges(6,1) * t131 - t108 * rSges(6,2) - t54 * rSges(6,3);
t78 = -Icges(6,5) * t89 - Icges(6,6) * t91;
t102 = (-t32 * t91 - t34 * t89) * t53 / 0.2e1 - (-t33 * t91 - t35 * t89) * t54 / 0.2e1 + (-t54 * t78 + t116 * t156 + (t107 * Icges(6,4) + t108 * Icges(6,2) + Icges(6,6) * t54) * t91 / 0.2e1 + (t107 * Icges(6,1) + t108 * Icges(6,4) + Icges(6,5) * t54) * t89 / 0.2e1 + t65 * t155) * t65 + (t117 * qJD(5) / 0.2e1 - (t105 * Icges(6,4) + t106 * Icges(6,2) - Icges(6,6) * t53) * t91 / 0.2e1 - (t105 * Icges(6,1) + t106 * Icges(6,4) - Icges(6,5) * t53) * t89 / 0.2e1 - t53 * t78 + t64 * t155) * t64;
t99 = (-pkin(3) * t122 + t119) * qJD(3);
t70 = -t137 * t92 - t123;
t98 = -t122 + t124;
t26 = -t65 * pkin(4) - t64 * pkin(7) + t36;
t97 = (Icges(6,1) * t89 - t111 + t138) * t134 + (-Icges(6,2) * t91 - t112 - t139) * t135;
t27 = -t65 * pkin(7) + t128 * t64 + t127;
t96 = t154 * t136 - t119;
t95 = -pkin(3) * t123 - t154 * t137;
t68 = t70 * qJD(3);
t94 = pkin(3) * t68;
t12 = -t54 * pkin(7) + t128 * t53 + t109;
t13 = rSges(6,1) * t129 + t106 * rSges(6,2) - t54 * pkin(4) - t53 * pkin(7) + t140;
t77 = t118 ^ 2;
t67 = t98 * qJD(3);
t45 = rSges(4,1) * t68 + rSges(4,2) * t67;
t44 = rSges(4,1) * t67 - rSges(4,2) * t68;
t39 = t43 + t95;
t38 = t96 + t42;
t37 = -t64 * t147 - t127;
t29 = t40 + t99;
t28 = t94 + t41;
t23 = t95 + t27;
t22 = t26 + t96;
t11 = t12 + t99;
t10 = t94 + t13;
t1 = t54 * t36 + t53 * t37 + t64 * (-rSges(6,1) * t144 - t109) + (t118 * t65 * qJD(5) + rSges(6,2) * t143 + t140) * t65;
t2 = [0; 0; 0; 0; m(4) * (t45 * t136 - t44 * t137) + m(5) * (t28 * t136 - t29 * t137) + m(6) * (t10 * t136 - t11 * t137); (t10 * t22 + t11 * t23) * t148 + (t28 * t38 + t29 * t39) * t149 + 0.2e1 * m(4) * ((-rSges(4,1) * t98 + rSges(4,2) * t70) * t45 + (rSges(4,1) * t70 + rSges(4,2) * t98) * t44) + t97; 0; m(5) * (t41 * t136 - t40 * t137) + m(6) * (-t12 * t137 + t13 * t136); m(6) * (t10 * t26 + t11 * t27 + t12 * t23 + t13 * t22) + m(5) * (t28 * t42 + t29 * t43 + t38 * t41 + t39 * t40) + t97; (t40 * t43 + t41 * t42) * t149 + (t12 * t27 + t13 * t26) * t148 + t97; m(6) * t1; m(6) * (-(-t136 * t53 + t137 * t54) * t118 + (-t136 * t64 + t137 * t65) * t75); m(6) * ((-t22 * t64 - t23 * t65) * t75 - (-t10 * t64 - t11 * t65 - t22 * t53 - t23 * t54) * t118) + t102; m(6) * ((-t26 * t64 - t27 * t65) * t75 - (-t12 * t65 - t13 * t64 - t26 * t53 - t27 * t54) * t118) + t102; ((t36 * t1 + t77 * t54) * t148 + (t103 + t151) * t53 + (0.3e1 * t54 * t31 + t152) * t65) * t65 + ((t117 + t31) * t53 * t65 + (-0.3e1 * t53 * t30 + t152) * t64 + (t104 + t151 + (-t30 + t116) * t65) * t54 + (t37 * t1 + t77 * t53) * t148) * t64;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
