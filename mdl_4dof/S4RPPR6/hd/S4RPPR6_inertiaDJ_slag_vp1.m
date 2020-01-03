% Calculate time derivative of joint inertia matrix for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:39
% DurationCPUTime: 0.82s
% Computational Cost: add. (1021->187), mult. (2820->283), div. (0->0), fcn. (2653->6), ass. (0->89)
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t115 = t70 ^ 2 + t72 ^ 2;
t67 = sin(pkin(6));
t68 = cos(pkin(6));
t69 = sin(qJ(4));
t71 = cos(qJ(4));
t55 = t67 * t71 - t68 * t69;
t114 = qJD(4) * t55;
t63 = t72 * qJ(2);
t105 = -rSges(5,3) - pkin(5);
t95 = qJ(3) * t67;
t77 = -t95 - pkin(1) + (-pkin(2) - pkin(3)) * t68;
t74 = t105 * t72 + t77 * t70;
t48 = t55 * t70;
t80 = t67 * t69 + t68 * t71;
t49 = t80 * t70;
t82 = -t49 * rSges(5,1) - t48 * rSges(5,2);
t14 = t63 + t74 + t82;
t100 = t68 * t72;
t96 = t72 * pkin(1) + t70 * qJ(2);
t86 = pkin(2) * t100 + t72 * t95 + t96;
t50 = t55 * t72;
t51 = t80 * t72;
t98 = t51 * rSges(5,1) + t50 * rSges(5,2);
t15 = pkin(3) * t100 + t105 * t70 + t86 + t98;
t113 = t14 * t72 + t15 * t70;
t92 = qJD(3) * t67;
t93 = qJD(1) * t72;
t97 = qJ(2) * t93 + qJD(2) * t70;
t88 = t72 * t92 + t97;
t52 = t80 * qJD(4);
t94 = qJD(1) * t70;
t30 = -t52 * t72 - t55 * t94;
t31 = t114 * t72 - t80 * t94;
t99 = t31 * rSges(5,1) + t30 * rSges(5,2);
t6 = t74 * qJD(1) + t88 + t99;
t32 = qJD(1) * t50 - t52 * t70;
t33 = qJD(1) * t51 + t114 * t70;
t83 = t33 * rSges(5,1) + t32 * rSges(5,2);
t61 = qJD(2) * t72;
t85 = -t70 * t92 + t61;
t7 = ((-qJ(2) - t105) * t70 + t77 * t72) * qJD(1) - t83 + t85;
t112 = t6 * t70 + t7 * t72;
t111 = 2 * m(5);
t91 = Icges(5,5) * qJD(1);
t90 = Icges(5,6) * qJD(1);
t89 = Icges(5,3) * qJD(1);
t39 = -rSges(5,1) * t52 - rSges(5,2) * t114;
t87 = t39 * t115;
t84 = rSges(3,1) * t68 - rSges(3,2) * t67;
t79 = -pkin(1) - t84;
t76 = -pkin(1) + (-rSges(4,1) - pkin(2)) * t68 + (-rSges(4,3) - qJ(3)) * t67;
t75 = rSges(3,3) * t72 + t79 * t70;
t73 = rSges(4,2) * t72 + t76 * t70;
t45 = t70 * rSges(3,3) + t84 * t72 + t96;
t44 = t63 + t75;
t43 = rSges(5,1) * t55 - rSges(5,2) * t80;
t42 = Icges(5,1) * t55 - Icges(5,4) * t80;
t41 = Icges(5,4) * t55 - Icges(5,2) * t80;
t38 = -Icges(5,1) * t52 - Icges(5,4) * t114;
t37 = -Icges(5,4) * t52 - Icges(5,2) * t114;
t36 = -Icges(5,5) * t52 - Icges(5,6) * t114;
t35 = t61 + ((-rSges(3,3) - qJ(2)) * t70 + t79 * t72) * qJD(1);
t34 = t75 * qJD(1) + t97;
t29 = t70 * rSges(4,2) + (rSges(4,1) * t68 + rSges(4,3) * t67) * t72 + t86;
t28 = t63 + t73;
t25 = -t70 * rSges(5,3) + t98;
t24 = rSges(5,3) * t72 - t82;
t23 = Icges(5,1) * t51 + Icges(5,4) * t50 - Icges(5,5) * t70;
t22 = Icges(5,1) * t49 + Icges(5,4) * t48 + Icges(5,5) * t72;
t21 = Icges(5,4) * t51 + Icges(5,2) * t50 - Icges(5,6) * t70;
t20 = Icges(5,4) * t49 + Icges(5,2) * t48 + Icges(5,6) * t72;
t19 = Icges(5,5) * t51 + Icges(5,6) * t50 - Icges(5,3) * t70;
t18 = Icges(5,5) * t49 + Icges(5,6) * t48 + Icges(5,3) * t72;
t17 = ((-rSges(4,2) - qJ(2)) * t70 + t76 * t72) * qJD(1) + t85;
t16 = t73 * qJD(1) + t88;
t13 = Icges(5,1) * t33 + Icges(5,4) * t32 - t70 * t91;
t12 = Icges(5,1) * t31 + Icges(5,4) * t30 - t72 * t91;
t11 = Icges(5,4) * t33 + Icges(5,2) * t32 - t70 * t90;
t10 = Icges(5,4) * t31 + Icges(5,2) * t30 - t72 * t90;
t9 = Icges(5,5) * t33 + Icges(5,6) * t32 - t70 * t89;
t8 = Icges(5,5) * t31 + Icges(5,6) * t30 - t72 * t89;
t5 = -t19 * t70 + t21 * t50 + t23 * t51;
t4 = -t18 * t70 + t20 * t50 + t22 * t51;
t3 = t19 * t72 + t48 * t21 + t49 * t23;
t2 = t18 * t72 + t48 * t20 + t49 * t22;
t1 = -t70 * t83 - t72 * t99 + (t115 * rSges(5,3) - t72 * t24 + t70 * t25) * qJD(1);
t26 = [0.2e1 * m(3) * (t34 * t45 + t35 * t44) + 0.2e1 * m(4) * (t16 * t29 + t17 * t28) - t52 * t42 + t55 * t38 - t114 * t41 - t80 * t37 + (t14 * t7 + t15 * t6) * t111; m(3) * (-t72 * t34 + t70 * t35 + (t44 * t72 + t45 * t70) * qJD(1)) + m(4) * (-t72 * t16 + t70 * t17 + (t28 * t72 + t29 * t70) * qJD(1)) + m(5) * (t113 * qJD(1) - t72 * t6 + t70 * t7); 0; 0.2e1 * (m(4) * (t16 * t70 + t17 * t72 - t28 * t94 + t29 * t93) / 0.2e1 + m(5) * (-t14 * t94 + t15 * t93 + t112) / 0.2e1) * t67; 0; 0; m(5) * (t112 * t43 + t113 * t39) + (m(5) * (-t14 * t70 + t15 * t72) * t43 - (-t21 * t80 + t55 * t23 + t50 * t41 + t51 * t42) * t72 / 0.2e1) * qJD(1) - (-t80 * t10 + t55 * t12 - t114 * t21 - t52 * t23 + t30 * t41 + t31 * t42 - t70 * t36 + t50 * t37 + t51 * t38 + (-t20 * t80 + t55 * t22 + t48 * t41 + t49 * t42) * qJD(1)) * t70 / 0.2e1 + (-t11 * t80 - t114 * t20 + t55 * t13 - t52 * t22 + t32 * t41 + t33 * t42 + t72 * t36 + t48 * t37 + t49 * t38) * t72 / 0.2e1; 0; m(5) * (-t1 * t68 + t67 * t87); ((-t70 * t24 - t25 * t72) * t1 + t43 * t87) * t111 - (t4 * t72 - t5 * t70) * t93 - t70 * (-(t50 * t10 + t51 * t12 + t30 * t21 + t31 * t23 - t70 * t8) * t70 + (t50 * t11 + t51 * t13 + t30 * t20 + t31 * t22 - t70 * t9) * t72 + (-t4 * t70 - t5 * t72) * qJD(1)) - (t2 * t72 - t3 * t70) * t94 + t72 * (-(t48 * t10 + t49 * t12 + t32 * t21 + t33 * t23 + t72 * t8) * t70 + (t48 * t11 + t49 * t13 + t32 * t20 + t33 * t22 + t72 * t9) * t72 + (-t2 * t70 - t3 * t72) * qJD(1));];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t26(1), t26(2), t26(4), t26(7); t26(2), t26(3), t26(5), t26(8); t26(4), t26(5), t26(6), t26(9); t26(7), t26(8), t26(9), t26(10);];
Mq = res;
