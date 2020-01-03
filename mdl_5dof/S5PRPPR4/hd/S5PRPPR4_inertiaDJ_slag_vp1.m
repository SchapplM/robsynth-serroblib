% Calculate time derivative of joint inertia matrix for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:47
% EndTime: 2019-12-31 17:36:48
% DurationCPUTime: 0.87s
% Computational Cost: add. (2102->188), mult. (2910->284), div. (0->0), fcn. (2733->6), ass. (0->90)
t69 = pkin(7) + qJ(2);
t67 = sin(t69);
t68 = cos(t69);
t116 = t67 ^ 2 + t68 ^ 2;
t70 = sin(pkin(8));
t71 = cos(pkin(8));
t72 = sin(qJ(5));
t73 = cos(qJ(5));
t55 = t70 * t73 - t71 * t72;
t115 = qJD(5) * t55;
t63 = t68 * qJ(3);
t106 = -rSges(6,3) - pkin(6);
t96 = qJ(4) * t70;
t78 = -t96 - pkin(2) + (-pkin(3) - pkin(4)) * t71;
t75 = t106 * t68 + t78 * t67;
t48 = t55 * t67;
t81 = t70 * t72 + t71 * t73;
t49 = t81 * t67;
t83 = -rSges(6,1) * t49 - rSges(6,2) * t48;
t14 = t63 + t75 + t83;
t101 = t68 * t71;
t97 = t68 * pkin(2) + t67 * qJ(3);
t87 = pkin(3) * t101 + t68 * t96 + t97;
t50 = t55 * t68;
t51 = t81 * t68;
t99 = t51 * rSges(6,1) + t50 * rSges(6,2);
t15 = pkin(4) * t101 + t106 * t67 + t87 + t99;
t114 = t14 * t68 + t15 * t67;
t52 = t81 * qJD(5);
t95 = qJD(2) * t67;
t32 = -t52 * t68 - t55 * t95;
t33 = t115 * t68 - t81 * t95;
t100 = t33 * rSges(6,1) + t32 * rSges(6,2);
t93 = qJD(4) * t70;
t94 = qJD(2) * t68;
t98 = qJ(3) * t94 + qJD(3) * t67;
t89 = t68 * t93 + t98;
t6 = t75 * qJD(2) + t100 + t89;
t34 = qJD(2) * t50 - t52 * t67;
t35 = qJD(2) * t51 + t115 * t67;
t84 = rSges(6,1) * t35 + rSges(6,2) * t34;
t61 = qJD(3) * t68;
t86 = -t67 * t93 + t61;
t7 = ((-qJ(3) - t106) * t67 + t78 * t68) * qJD(2) - t84 + t86;
t113 = t6 * t67 + t68 * t7;
t112 = 2 * m(6);
t92 = Icges(6,5) * qJD(2);
t91 = Icges(6,6) * qJD(2);
t90 = Icges(6,3) * qJD(2);
t41 = -rSges(6,1) * t52 - rSges(6,2) * t115;
t88 = t41 * t116;
t85 = rSges(4,1) * t71 - rSges(4,2) * t70;
t80 = -pkin(2) - t85;
t77 = -pkin(2) + (-rSges(5,1) - pkin(3)) * t71 + (-rSges(5,3) - qJ(4)) * t70;
t76 = rSges(4,3) * t68 + t80 * t67;
t74 = rSges(5,2) * t68 + t77 * t67;
t45 = rSges(6,1) * t55 - rSges(6,2) * t81;
t44 = Icges(6,1) * t55 - Icges(6,4) * t81;
t43 = Icges(6,4) * t55 - Icges(6,2) * t81;
t40 = -Icges(6,1) * t52 - Icges(6,4) * t115;
t39 = -Icges(6,4) * t52 - Icges(6,2) * t115;
t38 = -Icges(6,5) * t52 - Icges(6,6) * t115;
t37 = rSges(4,3) * t67 + t85 * t68 + t97;
t36 = t63 + t76;
t29 = t61 + ((-rSges(4,3) - qJ(3)) * t67 + t80 * t68) * qJD(2);
t28 = t76 * qJD(2) + t98;
t27 = -rSges(6,3) * t67 + t99;
t26 = rSges(6,3) * t68 - t83;
t25 = Icges(6,1) * t51 + Icges(6,4) * t50 - Icges(6,5) * t67;
t24 = Icges(6,1) * t49 + Icges(6,4) * t48 + Icges(6,5) * t68;
t23 = Icges(6,4) * t51 + Icges(6,2) * t50 - Icges(6,6) * t67;
t22 = Icges(6,4) * t49 + Icges(6,2) * t48 + Icges(6,6) * t68;
t21 = Icges(6,5) * t51 + Icges(6,6) * t50 - Icges(6,3) * t67;
t20 = Icges(6,5) * t49 + Icges(6,6) * t48 + Icges(6,3) * t68;
t19 = rSges(5,2) * t67 + (rSges(5,1) * t71 + rSges(5,3) * t70) * t68 + t87;
t18 = t63 + t74;
t17 = ((-rSges(5,2) - qJ(3)) * t67 + t77 * t68) * qJD(2) + t86;
t16 = t74 * qJD(2) + t89;
t13 = Icges(6,1) * t35 + Icges(6,4) * t34 - t67 * t92;
t12 = Icges(6,1) * t33 + Icges(6,4) * t32 - t68 * t92;
t11 = Icges(6,4) * t35 + Icges(6,2) * t34 - t67 * t91;
t10 = Icges(6,4) * t33 + Icges(6,2) * t32 - t68 * t91;
t9 = Icges(6,5) * t35 + Icges(6,6) * t34 - t67 * t90;
t8 = Icges(6,5) * t33 + Icges(6,6) * t32 - t68 * t90;
t5 = -t21 * t67 + t23 * t50 + t25 * t51;
t4 = -t20 * t67 + t22 * t50 + t24 * t51;
t3 = t21 * t68 + t23 * t48 + t25 * t49;
t2 = t20 * t68 + t22 * t48 + t24 * t49;
t1 = -t67 * t84 - t68 * t100 + (t116 * rSges(6,3) - t68 * t26 + t67 * t27) * qJD(2);
t30 = [0; 0; -t52 * t44 + t55 * t40 - t115 * t43 - t81 * t39 + 0.2e1 * m(4) * (t28 * t37 + t29 * t36) + 0.2e1 * m(5) * (t16 * t19 + t17 * t18) + (t14 * t7 + t15 * t6) * t112; 0; m(4) * (-t28 * t68 + t29 * t67 + (t36 * t68 + t37 * t67) * qJD(2)) + m(5) * (-t16 * t68 + t17 * t67 + (t18 * t68 + t19 * t67) * qJD(2)) + m(6) * (t114 * qJD(2) - t6 * t68 + t67 * t7); 0; 0; 0.2e1 * (m(5) * (t16 * t67 + t17 * t68 - t18 * t95 + t19 * t94) / 0.2e1 + m(6) * (-t14 * t95 + t15 * t94 + t113) / 0.2e1) * t70; 0; 0; m(6) * t1; m(6) * (t113 * t45 + t114 * t41) + (m(6) * (-t14 * t67 + t15 * t68) * t45 - (-t23 * t81 + t25 * t55 + t43 * t50 + t44 * t51) * t68 / 0.2e1) * qJD(2) - (-t81 * t10 + t55 * t12 - t115 * t23 - t52 * t25 + t32 * t43 + t33 * t44 - t67 * t38 + t50 * t39 + t51 * t40 + (-t22 * t81 + t24 * t55 + t43 * t48 + t44 * t49) * qJD(2)) * t67 / 0.2e1 + (-t11 * t81 - t115 * t22 + t55 * t13 - t52 * t24 + t34 * t43 + t35 * t44 + t68 * t38 + t48 * t39 + t49 * t40) * t68 / 0.2e1; 0; m(6) * (-t1 * t71 + t70 * t88); ((-t26 * t67 - t27 * t68) * t1 + t45 * t88) * t112 - (t4 * t68 - t5 * t67) * t94 - t67 * (-(t10 * t50 + t12 * t51 + t23 * t32 + t25 * t33 - t67 * t8) * t67 + (t11 * t50 + t13 * t51 + t22 * t32 + t24 * t33 - t67 * t9) * t68 + (-t4 * t67 - t5 * t68) * qJD(2)) - (t2 * t68 - t3 * t67) * t95 + t68 * (-(t10 * t48 + t12 * t49 + t23 * t34 + t25 * t35 + t68 * t8) * t67 + (t11 * t48 + t13 * t49 + t22 * t34 + t24 * t35 + t68 * t9) * t68 + (-t2 * t67 - t3 * t68) * qJD(2));];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t30(1), t30(2), t30(4), t30(7), t30(11); t30(2), t30(3), t30(5), t30(8), t30(12); t30(4), t30(5), t30(6), t30(9), t30(13); t30(7), t30(8), t30(9), t30(10), t30(14); t30(11), t30(12), t30(13), t30(14), t30(15);];
Mq = res;
