% Calculate time derivative of joint inertia matrix for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:32
% DurationCPUTime: 1.47s
% Computational Cost: add. (4199->177), mult. (4836->264), div. (0->0), fcn. (4704->8), ass. (0->88)
t135 = qJ(1) + pkin(8);
t125 = cos(t135);
t124 = sin(t135);
t134 = t125 * pkin(2) + t124 * qJ(3) + cos(qJ(1)) * pkin(1);
t129 = t125 * pkin(3) + t134;
t53 = rSges(4,1) * t125 + rSges(4,3) * t124 + t134;
t86 = sin(qJ(5));
t137 = qJD(5) * t86;
t163 = qJD(1) - qJD(4);
t148 = sin(qJ(4));
t149 = cos(qJ(4));
t58 = -t124 * t148 - t125 * t149;
t48 = t163 * t58;
t59 = -t124 * t149 + t125 * t148;
t88 = cos(qJ(5));
t112 = t59 * t137 - t48 * t88;
t49 = t163 * t59;
t110 = t58 * t137 + t49 * t88;
t136 = qJD(5) * t88;
t111 = t136 * t58 - t49 * t86;
t113 = t136 * t59 + t48 * t86;
t155 = 2 * m(6);
t146 = rSges(6,2) * t86;
t147 = rSges(6,1) * t88;
t65 = (t146 - t147) * qJD(5);
t74 = -rSges(6,1) * t86 - rSges(6,2) * t88;
t162 = -t58 * (Icges(6,5) * t112 + Icges(6,6) * t113 + Icges(6,3) * t49) + t74 * t65 * t155 + t59 * (Icges(6,5) * t110 + Icges(6,6) * t111 + Icges(6,3) * t48);
t138 = Icges(6,4) * t88;
t115 = Icges(6,2) * t86 - t138;
t34 = -Icges(6,6) * t58 + t115 * t59;
t139 = Icges(6,4) * t86;
t116 = -Icges(6,1) * t88 + t139;
t36 = -Icges(6,5) * t58 + t116 * t59;
t121 = t34 * t86 - t36 * t88;
t159 = t121 * t58;
t35 = Icges(6,6) * t59 + t115 * t58;
t37 = Icges(6,5) * t59 + t116 * t58;
t120 = t35 * t86 - t37 * t88;
t160 = t120 * t59;
t114 = -Icges(6,5) * t88 + Icges(6,6) * t86;
t32 = -Icges(6,3) * t58 + t114 * t59;
t33 = Icges(6,3) * t59 + t114 * t58;
t161 = t59 * t32 - t58 * t33 + t159 + t160;
t108 = t110 * rSges(6,1) + t48 * rSges(6,3);
t9 = -rSges(6,2) * t111 - t49 * pkin(4) - t48 * pkin(7) - t108;
t109 = t112 * rSges(6,1) + t49 * rSges(6,3);
t8 = rSges(6,2) * t113 - t48 * pkin(4) + t49 * pkin(7) + t109;
t131 = -pkin(4) + t146;
t141 = t59 * rSges(6,3) - t58 * t147;
t27 = -t59 * pkin(7) - t131 * t58 - t141;
t142 = t58 * rSges(6,3) + t59 * t147;
t26 = -t58 * pkin(7) + t131 * t59 - t142;
t130 = (Icges(6,2) * t88 + t116 + t139) * t137;
t63 = t115 * qJD(5);
t69 = -Icges(6,1) * t86 - t138;
t158 = -(qJD(5) * t69 + t63) * t88 - t130;
t156 = 2 * m(5);
t82 = t125 * qJ(3);
t140 = qJD(1) * t82 + qJD(3) * t124;
t31 = -t49 * rSges(5,1) + t48 * rSges(5,2);
t30 = -t48 * rSges(5,1) - t49 * rSges(5,2);
t50 = -t59 * rSges(5,1) + t58 * rSges(5,2);
t51 = t58 * rSges(5,1) + t59 * rSges(5,2);
t103 = -pkin(2) * t124 - sin(qJ(1)) * pkin(1);
t97 = -pkin(3) * t124 + t103;
t94 = t82 + t97;
t93 = qJD(1) * t97 + t140;
t92 = -rSges(4,1) * t124 + rSges(4,3) * t125 + t103;
t80 = qJD(3) * t125;
t91 = -t129 * qJD(1) + t80;
t62 = t114 * qJD(5);
t90 = (-t35 * t88 - t37 * t86) * t48 / 0.2e1 + (-t34 * t88 - t36 * t86) * t49 / 0.2e1 - (qJD(5) * t121 - (Icges(6,4) * t112 + Icges(6,2) * t113 + Icges(6,6) * t49) * t88 - (Icges(6,1) * t112 + Icges(6,4) * t113 + Icges(6,5) * t49) * t86 - t58 * t62) * t58 / 0.2e1 + (qJD(5) * t120 - (Icges(6,4) * t110 + Icges(6,2) * t111 + Icges(6,6) * t48) * t88 - (Icges(6,1) * t110 + Icges(6,4) * t111 + Icges(6,5) * t48) * t86 + t59 * t62) * t59 / 0.2e1 + (t48 * t59 - t49 * t58) * (-Icges(6,5) * t86 - Icges(6,6) * t88);
t66 = t74 ^ 2;
t52 = t82 + t92;
t43 = -t53 * qJD(1) + t80;
t42 = qJD(1) * t92 + t140;
t41 = -t51 + t129;
t40 = -t50 + t94;
t39 = t146 * t58 + t141;
t38 = t146 * t59 - t142;
t25 = -t30 + t91;
t24 = -t31 + t93;
t23 = t129 - t27;
t22 = -t26 + t94;
t7 = -t8 + t91;
t6 = -t9 + t93;
t1 = t48 * t38 + t59 * t109 - t49 * t39 + t58 * t108 + (t111 * t58 + t113 * t59) * rSges(6,2);
t2 = [(t22 * t7 + t23 * t6) * t155 - t69 * t136 - t88 * t63 + (t24 * t41 + t25 * t40) * t156 + 0.2e1 * m(4) * (t42 * t53 + t43 * t52) - t130; 0; 0; m(6) * (t124 * t7 - t125 * t6 + (t124 * t23 + t125 * t22) * qJD(1)) + m(5) * (t124 * t25 - t125 * t24 + (t124 * t41 + t125 * t40) * qJD(1)) + m(4) * (t124 * t43 - t125 * t42 + (t124 * t53 + t125 * t52) * qJD(1)); 0; 0; m(6) * (t22 * t8 + t23 * t9 + t26 * t7 + t27 * t6) + m(5) * (t24 * t51 + t25 * t50 + t30 * t40 + t31 * t41) - t158; 0; m(5) * (t30 * t124 - t31 * t125 + (t124 * t51 + t125 * t50) * qJD(1)) + m(6) * (t8 * t124 - t9 * t125 + (t124 * t27 + t125 * t26) * qJD(1)); (t30 * t50 + t31 * t51) * t156 + (t26 * t8 + t27 * t9) * t155 + t158; m(6) * ((-t22 * t58 - t23 * t59) * t65 + (t22 * t49 - t23 * t48 - t58 * t7 - t59 * t6) * t74) + t90; m(6) * t1; m(6) * ((-t124 * t58 + t125 * t59) * t65 + (t124 * t49 + t125 * t48 + (-t124 * t59 - t125 * t58) * qJD(1)) * t74); m(6) * ((-t26 * t58 - t27 * t59) * t65 + (t26 * t49 - t27 * t48 - t58 * t8 - t59 * t9) * t74) - t90; ((t38 * t1 + t66 * t48) * t155 + (-t160 + t161) * t49 + (0.3e1 * t48 * t33 + t162) * t59) * t59 + ((-t121 - t33) * t49 * t59 + (0.3e1 * t49 * t32 + t162) * t58 - (-t159 + t161 + (t32 - t120) * t59) * t48 + (t39 * t1 - t66 * t49) * t155) * t58;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
