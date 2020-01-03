% Calculate time derivative of joint inertia matrix for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:58
% DurationCPUTime: 0.66s
% Computational Cost: add. (1025->120), mult. (1260->193), div. (0->0), fcn. (934->4), ass. (0->73)
t53 = sin(qJ(4));
t54 = cos(qJ(4));
t70 = rSges(5,1) * t53 + rSges(5,2) * t54;
t52 = pkin(6) + qJ(2);
t51 = cos(t52);
t58 = t51 * t70;
t50 = sin(t52);
t83 = Icges(5,4) * t53;
t62 = Icges(5,2) * t54 + t83;
t22 = Icges(5,6) * t51 + t62 * t50;
t82 = Icges(5,4) * t54;
t64 = Icges(5,1) * t53 + t82;
t24 = Icges(5,5) * t51 + t64 * t50;
t67 = t22 * t54 + t24 * t53;
t103 = t67 * t51;
t74 = qJD(4) * t54;
t77 = qJD(2) * t51;
t72 = t50 * rSges(5,1) * t74 + t70 * t77;
t73 = -rSges(5,3) - pkin(2) - pkin(5);
t75 = qJD(4) * t53;
t85 = qJ(3) * t77 + qJD(3) * t50;
t6 = (-rSges(5,2) * t75 + t73 * qJD(2)) * t50 + t72 + t85;
t90 = t53 * rSges(5,2);
t38 = t54 * rSges(5,1) - t90;
t43 = qJD(3) * t51;
t76 = qJD(4) * t51;
t7 = t43 + t38 * t76 + (t73 * t51 + (-qJ(3) - t70) * t50) * qJD(2);
t102 = t50 * t7 - t51 * t6;
t60 = Icges(5,5) * t53 + Icges(5,6) * t54;
t101 = -Icges(5,3) * t50 + t60 * t51;
t100 = -Icges(5,6) * t50 + t62 * t51;
t99 = -Icges(5,5) * t50 + t64 * t51;
t98 = 2 * m(5);
t48 = t50 ^ 2;
t49 = t51 ^ 2;
t97 = m(5) * t38;
t94 = rSges(4,2) - pkin(2);
t91 = t50 * rSges(5,3);
t46 = t51 * rSges(5,3);
t89 = t53 * t22;
t88 = t53 * t100;
t87 = t54 * t24;
t86 = t54 * t99;
t84 = t51 * pkin(2) + t50 * qJ(3);
t20 = Icges(5,3) * t51 + t60 * t50;
t78 = qJD(2) * t20;
t26 = t50 * t70 + t46;
t31 = t70 * qJD(4);
t71 = (t48 + t49) * t31;
t66 = -t100 * t54 - t53 * t99;
t65 = Icges(5,1) * t54 - t83;
t63 = -Icges(5,2) * t53 + t82;
t61 = Icges(5,5) * t54 - Icges(5,6) * t53;
t59 = t51 * rSges(4,3) + t94 * t50;
t57 = t66 * t50;
t56 = qJD(4) * t65;
t55 = qJD(4) * t63;
t45 = t51 * qJ(3);
t27 = t91 - t58;
t19 = -t51 * rSges(4,2) + t50 * rSges(4,3) + t84;
t18 = t45 + t59;
t17 = t43 + (t94 * t51 + (-rSges(4,3) - qJ(3)) * t50) * qJD(2);
t16 = t59 * qJD(2) + t85;
t15 = t51 * pkin(5) + t26 + t84;
t14 = t73 * t50 + t45 + t58;
t9 = t61 * t50 * qJD(4) + t101 * qJD(2);
t8 = -t61 * t76 + t78;
t5 = -t101 * t50 - t66 * t51;
t4 = t50 * t20 - t103;
t3 = -t101 * t51 + t57;
t2 = t51 * t20 + t67 * t50;
t1 = -t50 * t72 + (-t38 * t49 + t48 * t90) * qJD(4) + ((-t26 + t46) * t51 + (-t27 + t58 + t91) * t50) * qJD(2);
t10 = [0; 0; 0.2e1 * m(4) * (t19 * t16 + t18 * t17) + (t14 * t7 + t15 * t6) * t98 - t53 * t56 - t64 * t74 - t54 * t55 + t62 * t75; 0; m(4) * (-t51 * t16 + t50 * t17 + (t18 * t51 + t19 * t50) * qJD(2)) + m(5) * ((t14 * t51 + t15 * t50) * qJD(2) + t102); 0; m(5) * t1; m(5) * (t102 * t38 - (t14 * t50 - t15 * t51) * t31) + (-t67 * qJD(4) - t53 * (t100 * qJD(2) + t50 * t55) + t54 * (t99 * qJD(2) + t50 * t56)) * t51 / 0.2e1 + (-t66 * qJD(4) - t53 * (t22 * qJD(2) - t63 * t76) + t54 * (t24 * qJD(2) - t65 * t76)) * t50 / 0.2e1 - (t49 / 0.2e1 + t48 / 0.2e1) * t60 * qJD(4) + ((t15 * t97 + t89 / 0.2e1 - t87 / 0.2e1) * t50 + (t14 * t97 + t88 / 0.2e1 - t86 / 0.2e1) * t51) * qJD(2); -m(5) * t71; ((-t50 * t26 + t51 * t27) * t1 - t38 * t71) * t98 - qJD(2) * t50 * (t2 * t51 + t3 * t50) + t51 * ((t51 * t9 + (t3 + t103) * qJD(2)) * t51 + (-t2 * qJD(2) + (t100 * t75 - t74 * t99) * t50 + (t8 + (t87 - t89) * qJD(4) + (-t20 + t66) * qJD(2)) * t51) * t50) + (t4 * t51 + t5 * t50) * t77 + t50 * ((t50 * t8 + (-t4 + t57) * qJD(2)) * t50 + (t5 * qJD(2) + (t22 * t75 - t24 * t74 + t78) * t51 + (t9 + (t86 - t88) * qJD(4) + t67 * qJD(2)) * t50) * t51);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
