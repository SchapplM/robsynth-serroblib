% Calculate time derivative of joint inertia matrix for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR9_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR9_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR9_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:40
% DurationCPUTime: 1.44s
% Computational Cost: add. (4151->175), mult. (4756->262), div. (0->0), fcn. (4654->6), ass. (0->88)
t85 = sin(qJ(5));
t131 = qJD(5) * t85;
t129 = pkin(8) + qJ(2);
t120 = cos(t129);
t119 = sin(t129);
t134 = t120 * pkin(2) + t119 * qJ(3);
t128 = t120 * pkin(3) + t134;
t53 = t120 * rSges(4,1) + t119 * rSges(4,3) + t134;
t157 = qJD(2) - qJD(4);
t143 = sin(qJ(4));
t144 = cos(qJ(4));
t58 = -t119 * t143 - t120 * t144;
t48 = t157 * t58;
t59 = -t119 * t144 + t120 * t143;
t86 = cos(qJ(5));
t106 = t59 * t131 - t48 * t86;
t49 = t157 * t59;
t104 = t58 * t131 + t49 * t86;
t130 = qJD(5) * t86;
t105 = t58 * t130 - t49 * t85;
t107 = t59 * t130 + t48 * t85;
t149 = 2 * m(6);
t141 = rSges(6,2) * t85;
t142 = rSges(6,1) * t86;
t65 = (t141 - t142) * qJD(5);
t74 = -t85 * rSges(6,1) - rSges(6,2) * t86;
t156 = -t58 * (t106 * Icges(6,5) + t107 * Icges(6,6) + Icges(6,3) * t49) + t74 * t65 * t149 + t59 * (t104 * Icges(6,5) + t105 * Icges(6,6) + Icges(6,3) * t48);
t132 = Icges(6,4) * t86;
t109 = Icges(6,2) * t85 - t132;
t34 = -Icges(6,6) * t58 + t109 * t59;
t133 = Icges(6,4) * t85;
t110 = -Icges(6,1) * t86 + t133;
t36 = -Icges(6,5) * t58 + t110 * t59;
t116 = t34 * t85 - t36 * t86;
t153 = t116 * t58;
t35 = Icges(6,6) * t59 + t109 * t58;
t37 = Icges(6,5) * t59 + t110 * t58;
t115 = t35 * t85 - t37 * t86;
t154 = t115 * t59;
t108 = -Icges(6,5) * t86 + Icges(6,6) * t85;
t32 = -Icges(6,3) * t58 + t108 * t59;
t33 = Icges(6,3) * t59 + t108 * t58;
t155 = t59 * t32 - t58 * t33 + t153 + t154;
t102 = rSges(6,1) * t104 + t48 * rSges(6,3);
t9 = -t105 * rSges(6,2) - t49 * pkin(4) - t48 * pkin(7) - t102;
t103 = rSges(6,1) * t106 + t49 * rSges(6,3);
t8 = t107 * rSges(6,2) - t48 * pkin(4) + t49 * pkin(7) + t103;
t125 = -pkin(4) + t141;
t136 = t59 * rSges(6,3) - t58 * t142;
t27 = -t59 * pkin(7) - t125 * t58 - t136;
t137 = t58 * rSges(6,3) + t59 * t142;
t26 = -t58 * pkin(7) + t125 * t59 - t137;
t124 = (Icges(6,2) * t86 + t110 + t133) * t131;
t63 = t109 * qJD(5);
t69 = -Icges(6,1) * t85 - t132;
t152 = -(qJD(5) * t69 + t63) * t86 - t124;
t150 = 2 * m(5);
t82 = t120 * qJ(3);
t135 = qJD(2) * t82 + qJD(3) * t119;
t31 = -t49 * rSges(5,1) + t48 * rSges(5,2);
t30 = -t48 * rSges(5,1) - t49 * rSges(5,2);
t50 = -t59 * rSges(5,1) + t58 * rSges(5,2);
t51 = t58 * rSges(5,1) + t59 * rSges(5,2);
t113 = t119 * pkin(2);
t95 = -t119 * pkin(3) - t113;
t92 = t82 + t95;
t90 = t95 * qJD(2) + t135;
t89 = -t119 * rSges(4,1) + t120 * rSges(4,3) - t113;
t80 = qJD(3) * t120;
t88 = -qJD(2) * t128 + t80;
t62 = t108 * qJD(5);
t87 = -(-t35 * t86 - t85 * t37) * t48 / 0.2e1 - (-t34 * t86 - t85 * t36) * t49 / 0.2e1 + (t116 * qJD(5) - (t106 * Icges(6,4) + t107 * Icges(6,2) + Icges(6,6) * t49) * t86 - t85 * (t106 * Icges(6,1) + t107 * Icges(6,4) + Icges(6,5) * t49) - t58 * t62) * t58 / 0.2e1 - (t115 * qJD(5) - (t104 * Icges(6,4) + t105 * Icges(6,2) + Icges(6,6) * t48) * t86 - t85 * (t104 * Icges(6,1) + t105 * Icges(6,4) + Icges(6,5) * t48) + t59 * t62) * t59 / 0.2e1 + (-t59 * t48 + t58 * t49) * (-Icges(6,5) * t85 - Icges(6,6) * t86);
t66 = t74 ^ 2;
t52 = t82 + t89;
t47 = -qJD(2) * t53 + t80;
t46 = t89 * qJD(2) + t135;
t41 = -t51 + t128;
t40 = -t50 + t92;
t39 = t58 * t141 + t136;
t38 = t59 * t141 - t137;
t25 = -t30 + t88;
t24 = -t31 + t90;
t23 = t128 - t27;
t22 = -t26 + t92;
t7 = -t8 + t88;
t6 = -t9 + t90;
t1 = t48 * t38 + t59 * t103 - t49 * t39 + t58 * t102 + (t58 * t105 + t59 * t107) * rSges(6,2);
t2 = [0; 0; -t69 * t130 - t86 * t63 + (t22 * t7 + t23 * t6) * t149 + (t24 * t41 + t25 * t40) * t150 + 0.2e1 * m(4) * (t46 * t53 + t47 * t52) - t124; 0; m(6) * (t119 * t7 - t120 * t6 + (t119 * t23 + t120 * t22) * qJD(2)) + m(5) * (t119 * t25 - t120 * t24 + (t119 * t41 + t120 * t40) * qJD(2)) + m(4) * (t119 * t47 - t120 * t46 + (t119 * t53 + t120 * t52) * qJD(2)); 0; 0; m(6) * (t22 * t8 + t23 * t9 + t26 * t7 + t27 * t6) + m(5) * (t24 * t51 + t25 * t50 + t30 * t40 + t31 * t41) - t152; m(5) * (t30 * t119 - t31 * t120 + (t119 * t51 + t120 * t50) * qJD(2)) + m(6) * (t8 * t119 - t9 * t120 + (t119 * t27 + t120 * t26) * qJD(2)); (t30 * t50 + t31 * t51) * t150 + (t26 * t8 + t27 * t9) * t149 + t152; m(6) * t1; m(6) * ((-t22 * t58 - t23 * t59) * t65 + (t22 * t49 - t23 * t48 - t58 * t7 - t59 * t6) * t74) - t87; m(6) * ((-t119 * t58 + t120 * t59) * t65 + (t119 * t49 + t120 * t48 + (-t119 * t59 - t120 * t58) * qJD(2)) * t74); m(6) * ((-t26 * t58 - t27 * t59) * t65 + (t26 * t49 - t27 * t48 - t58 * t8 - t59 * t9) * t74) + t87; ((t38 * t1 + t66 * t48) * t149 + (-t154 + t155) * t49 + (0.3e1 * t48 * t33 + t156) * t59) * t59 + ((-t116 - t33) * t49 * t59 + (0.3e1 * t49 * t32 + t156) * t58 - (-t153 + t155 + (t32 - t115) * t59) * t48 + (t1 * t39 - t49 * t66) * t149) * t58;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
