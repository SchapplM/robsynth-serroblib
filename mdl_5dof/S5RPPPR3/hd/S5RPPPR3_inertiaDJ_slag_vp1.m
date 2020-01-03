% Calculate time derivative of joint inertia matrix for
% S5RPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:43:52
% EndTime: 2019-12-31 17:43:54
% DurationCPUTime: 0.89s
% Computational Cost: add. (2150->195), mult. (2990->286), div. (0->0), fcn. (2783->8), ass. (0->92)
t70 = qJ(1) + pkin(7);
t67 = sin(t70);
t68 = cos(t70);
t122 = t67 ^ 2 + t68 ^ 2;
t71 = sin(pkin(8));
t72 = cos(pkin(8));
t73 = sin(qJ(5));
t75 = cos(qJ(5));
t55 = t71 * t75 - t72 * t73;
t121 = qJD(5) * t55;
t63 = t68 * qJ(3);
t111 = -rSges(6,3) - pkin(6);
t114 = sin(qJ(1)) * pkin(1);
t102 = qJ(4) * t71;
t82 = -t102 - pkin(2) + (-pkin(3) - pkin(4)) * t72;
t78 = t111 * t68 + t67 * t82 - t114;
t48 = t55 * t67;
t85 = t71 * t73 + t72 * t75;
t49 = t85 * t67;
t87 = -rSges(6,1) * t49 - rSges(6,2) * t48;
t14 = t63 + t78 + t87;
t50 = t55 * t68;
t51 = t85 * t68;
t104 = t51 * rSges(6,1) + t50 * rSges(6,2);
t106 = t68 * t72;
t69 = cos(qJ(1)) * pkin(1);
t94 = t68 * pkin(2) + t67 * qJ(3) + t69;
t92 = pkin(3) * t106 + t68 * t102 + t94;
t15 = pkin(4) * t106 + t111 * t67 + t104 + t92;
t120 = t14 * t68 + t15 * t67;
t101 = qJD(1) * t67;
t52 = t85 * qJD(5);
t32 = -t101 * t55 - t52 * t68;
t33 = -t101 * t85 + t121 * t68;
t105 = t33 * rSges(6,1) + t32 * rSges(6,2);
t100 = qJD(1) * t68;
t103 = qJ(3) * t100 + qJD(3) * t67;
t99 = qJD(4) * t71;
t95 = t68 * t99 + t103;
t6 = qJD(1) * t78 + t105 + t95;
t34 = qJD(1) * t50 - t52 * t67;
t35 = qJD(1) * t51 + t121 * t67;
t88 = t35 * rSges(6,1) + t34 * rSges(6,2);
t61 = qJD(3) * t68;
t91 = -t67 * t99 + t61;
t7 = (-t69 + (-qJ(3) - t111) * t67 + t82 * t68) * qJD(1) - t88 + t91;
t119 = t6 * t67 + t68 * t7;
t118 = 2 * m(6);
t98 = Icges(6,5) * qJD(1);
t97 = Icges(6,6) * qJD(1);
t96 = Icges(6,3) * qJD(1);
t41 = -rSges(6,1) * t52 - rSges(6,2) * t121;
t93 = t41 * t122;
t89 = rSges(4,1) * t72 - rSges(4,2) * t71;
t84 = -pkin(2) - t89;
t80 = -pkin(2) + (-rSges(5,1) - pkin(3)) * t72 + (-rSges(5,3) - qJ(4)) * t71;
t79 = rSges(4,3) * t68 + t67 * t84 - t114;
t77 = rSges(5,2) * t68 + t67 * t80 - t114;
t45 = rSges(6,1) * t55 - rSges(6,2) * t85;
t44 = Icges(6,1) * t55 - Icges(6,4) * t85;
t43 = Icges(6,4) * t55 - Icges(6,2) * t85;
t40 = -Icges(6,1) * t52 - Icges(6,4) * t121;
t39 = -Icges(6,4) * t52 - Icges(6,2) * t121;
t38 = -Icges(6,5) * t52 - Icges(6,6) * t121;
t37 = rSges(4,3) * t67 + t68 * t89 + t94;
t36 = t63 + t79;
t29 = -rSges(6,3) * t67 + t104;
t28 = rSges(6,3) * t68 - t87;
t27 = Icges(6,1) * t51 + Icges(6,4) * t50 - Icges(6,5) * t67;
t26 = Icges(6,1) * t49 + Icges(6,4) * t48 + Icges(6,5) * t68;
t25 = Icges(6,4) * t51 + Icges(6,2) * t50 - Icges(6,6) * t67;
t24 = Icges(6,4) * t49 + Icges(6,2) * t48 + Icges(6,6) * t68;
t23 = Icges(6,5) * t51 + Icges(6,6) * t50 - Icges(6,3) * t67;
t22 = Icges(6,5) * t49 + Icges(6,6) * t48 + Icges(6,3) * t68;
t21 = t61 + (-t69 + (-rSges(4,3) - qJ(3)) * t67 + t84 * t68) * qJD(1);
t20 = qJD(1) * t79 + t103;
t19 = rSges(5,2) * t67 + (rSges(5,1) * t72 + rSges(5,3) * t71) * t68 + t92;
t18 = t63 + t77;
t17 = (-t69 + (-rSges(5,2) - qJ(3)) * t67 + t80 * t68) * qJD(1) + t91;
t16 = qJD(1) * t77 + t95;
t13 = Icges(6,1) * t35 + Icges(6,4) * t34 - t67 * t98;
t12 = Icges(6,1) * t33 + Icges(6,4) * t32 - t68 * t98;
t11 = Icges(6,4) * t35 + Icges(6,2) * t34 - t67 * t97;
t10 = Icges(6,4) * t33 + Icges(6,2) * t32 - t68 * t97;
t9 = Icges(6,5) * t35 + Icges(6,6) * t34 - t67 * t96;
t8 = Icges(6,5) * t33 + Icges(6,6) * t32 - t68 * t96;
t5 = -t23 * t67 + t25 * t50 + t27 * t51;
t4 = -t22 * t67 + t24 * t50 + t26 * t51;
t3 = t23 * t68 + t25 * t48 + t27 * t49;
t2 = t22 * t68 + t24 * t48 + t26 * t49;
t1 = -t67 * t88 - t68 * t105 + (t122 * rSges(6,3) - t68 * t28 + t67 * t29) * qJD(1);
t30 = [0.2e1 * m(4) * (t20 * t37 + t21 * t36) + 0.2e1 * m(5) * (t16 * t19 + t17 * t18) - t52 * t44 + t55 * t40 - t121 * t43 - t85 * t39 + (t14 * t7 + t15 * t6) * t118; 0; 0; m(4) * (-t20 * t68 + t21 * t67 + (t36 * t68 + t37 * t67) * qJD(1)) + m(5) * (-t16 * t68 + t17 * t67 + (t18 * t68 + t19 * t67) * qJD(1)) + m(6) * (t120 * qJD(1) - t6 * t68 + t67 * t7); 0; 0; 0.2e1 * (m(5) * (t100 * t19 - t101 * t18 + t16 * t67 + t17 * t68) / 0.2e1 + m(6) * (t100 * t15 - t101 * t14 + t119) / 0.2e1) * t71; 0; 0; 0; m(6) * (t119 * t45 + t120 * t41) + (m(6) * (-t14 * t67 + t15 * t68) * t45 - (-t25 * t85 + t27 * t55 + t43 * t50 + t44 * t51) * t68 / 0.2e1) * qJD(1) - (-t85 * t10 + t55 * t12 - t121 * t25 - t52 * t27 + t32 * t43 + t33 * t44 - t67 * t38 + t50 * t39 + t51 * t40 + (-t24 * t85 + t26 * t55 + t43 * t48 + t44 * t49) * qJD(1)) * t67 / 0.2e1 + (-t11 * t85 - t121 * t24 + t55 * t13 - t52 * t26 + t34 * t43 + t35 * t44 + t68 * t38 + t48 * t39 + t49 * t40) * t68 / 0.2e1; m(6) * t1; 0; m(6) * (-t1 * t72 + t71 * t93); ((-t28 * t67 - t29 * t68) * t1 + t45 * t93) * t118 - (t4 * t68 - t5 * t67) * t100 - t67 * (-(t10 * t50 + t12 * t51 + t25 * t32 + t27 * t33 - t67 * t8) * t67 + (t11 * t50 + t13 * t51 + t24 * t32 + t26 * t33 - t67 * t9) * t68 + (-t4 * t67 - t5 * t68) * qJD(1)) - (t2 * t68 - t3 * t67) * t101 + t68 * (-(t10 * t48 + t12 * t49 + t25 * t34 + t27 * t35 + t68 * t8) * t67 + (t11 * t48 + t13 * t49 + t24 * t34 + t26 * t35 + t68 * t9) * t68 + (-t2 * t67 - t3 * t68) * qJD(1));];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t30(1), t30(2), t30(4), t30(7), t30(11); t30(2), t30(3), t30(5), t30(8), t30(12); t30(4), t30(5), t30(6), t30(9), t30(13); t30(7), t30(8), t30(9), t30(10), t30(14); t30(11), t30(12), t30(13), t30(14), t30(15);];
Mq = res;
