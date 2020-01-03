% Calculate time derivative of joint inertia matrix for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR5_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR5_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:45
% DurationCPUTime: 1.05s
% Computational Cost: add. (1059->150), mult. (2560->225), div. (0->0), fcn. (2542->6), ass. (0->78)
t110 = sin(qJ(1));
t111 = cos(qJ(1));
t100 = t111 * pkin(1) + t110 * qJ(2);
t93 = t111 * pkin(2) + t100;
t32 = t111 * rSges(3,1) + t110 * rSges(3,3) + t100;
t96 = sin(pkin(6));
t97 = cos(pkin(6));
t39 = -t110 * t97 + t111 * t96;
t36 = t39 * qJD(1);
t62 = cos(qJ(4));
t38 = -t110 * t96 - t111 * t97;
t61 = sin(qJ(4));
t95 = qJD(4) * t61;
t91 = t38 * t95;
t71 = t62 * t36 + t91;
t117 = 2 * m(5);
t35 = t38 * qJD(1);
t108 = rSges(5,2) * t61;
t109 = rSges(5,1) * t62;
t83 = t108 - t109;
t44 = t83 * qJD(4);
t49 = -t61 * rSges(5,1) - rSges(5,2) * t62;
t94 = qJD(4) * t62;
t89 = t38 * t94;
t72 = -t61 * t36 + t89;
t90 = t39 * t95;
t73 = -t62 * t35 + t90;
t88 = t39 * t94;
t74 = t61 * t35 + t88;
t123 = -t38 * (t73 * Icges(5,5) + t74 * Icges(5,6) + Icges(5,3) * t36) + t49 * t44 * t117 + t39 * (t71 * Icges(5,5) + t72 * Icges(5,6) + Icges(5,3) * t35);
t98 = Icges(5,4) * t62;
t77 = Icges(5,2) * t61 - t98;
t20 = -Icges(5,6) * t38 + t77 * t39;
t99 = Icges(5,4) * t61;
t78 = -Icges(5,1) * t62 + t99;
t22 = -Icges(5,5) * t38 + t78 * t39;
t82 = t20 * t61 - t22 * t62;
t120 = t82 * t38;
t21 = Icges(5,6) * t39 + t77 * t38;
t23 = Icges(5,5) * t39 + t78 * t38;
t81 = t21 * t61 - t23 * t62;
t121 = t81 * t39;
t76 = -Icges(5,5) * t62 + Icges(5,6) * t61;
t18 = -Icges(5,3) * t38 + t76 * t39;
t19 = Icges(5,3) * t39 + t76 * t38;
t122 = t39 * t18 - t38 * t19 + t120 + t121;
t112 = rSges(5,3) + pkin(5);
t102 = t39 * rSges(5,3) - t38 * t109;
t58 = t111 * qJ(2);
t101 = qJD(1) * t58 + qJD(2) * t110;
t92 = t110 * pkin(1);
t75 = pkin(3) - t83;
t70 = t71 * rSges(5,1) + t35 * rSges(5,3);
t69 = -t110 * pkin(2) - t92;
t66 = t58 + t69;
t65 = -t110 * rSges(3,1) + t111 * rSges(3,3) - t92;
t64 = t69 * qJD(1) + t101;
t56 = qJD(2) * t111;
t63 = -t93 * qJD(1) + t56;
t48 = -Icges(5,1) * t61 - t98;
t47 = -Icges(5,2) * t62 - t99;
t45 = t49 ^ 2;
t41 = t76 * qJD(4);
t31 = t58 + t65;
t29 = -t32 * qJD(1) + t56;
t28 = qJD(1) * t65 + t101;
t27 = -rSges(4,1) * t38 - rSges(4,2) * t39 + t93;
t26 = t39 * rSges(4,1) - t38 * rSges(4,2) + t66;
t25 = t38 * t108 + t102;
t24 = -t38 * rSges(5,3) + t83 * t39;
t17 = t35 * rSges(4,1) + t36 * rSges(4,2) + t63;
t16 = t36 * rSges(4,1) - t35 * rSges(4,2) + t64;
t15 = pkin(5) * t39 + (-pkin(3) + t108) * t38 + t93 + t102;
t14 = t112 * t38 + t75 * t39 + t66;
t7 = t49 * t39 * qJD(4) - t112 * t36 + t75 * t35 + t63;
t6 = t72 * rSges(5,2) + t36 * pkin(3) + t35 * pkin(5) + t64 + t70;
t1 = t35 * t24 + t39 * (t73 * rSges(5,1) + t36 * rSges(5,3)) - t36 * t25 + t38 * t70 + (t38 * t72 + t39 * t74) * rSges(5,2);
t2 = [0.2e1 * m(3) * (t28 * t32 + t29 * t31) + 0.2e1 * m(4) * (t16 * t27 + t17 * t26) + (t14 * t7 + t15 * t6) * t117 + (-t78 + t47) * t95 + (-t48 - t77) * t94; m(3) * (t110 * t29 - t111 * t28 + (t110 * t32 + t111 * t31) * qJD(1)) + m(4) * (t110 * t17 - t111 * t16 + (t110 * t27 + t111 * t26) * qJD(1)) + m(5) * (t110 * t7 - t111 * t6 + (t110 * t15 + t111 * t14) * qJD(1)); 0; 0; 0; 0; m(5) * ((-t14 * t38 - t15 * t39) * t44 + (t14 * t36 - t15 * t35 - t38 * t7 - t39 * t6) * t49) + (t39 * t35 - t38 * t36) * (-Icges(5,5) * t61 - Icges(5,6) * t62) + (-t21 * t62 - t61 * t23) * t35 / 0.2e1 + (-t20 * t62 - t61 * t22) * t36 / 0.2e1 - (t82 * qJD(4) - (t73 * Icges(5,4) + t74 * Icges(5,2) + Icges(5,6) * t36) * t62 - t61 * (t73 * Icges(5,1) + t74 * Icges(5,4) + Icges(5,5) * t36) - t38 * t41 + t47 * t88 + t48 * t90) * t38 / 0.2e1 + (t81 * qJD(4) - (t71 * Icges(5,4) + t72 * Icges(5,2) + Icges(5,6) * t35) * t62 - t61 * (t71 * Icges(5,1) + t72 * Icges(5,4) + Icges(5,5) * t35) + t39 * t41 + t47 * t89 + t48 * t91) * t39 / 0.2e1; m(5) * ((-t110 * t38 + t111 * t39) * t44 + (t110 * t36 + t111 * t35 + (-t110 * t39 - t111 * t38) * qJD(1)) * t49); -m(5) * t1; ((t1 * t24 + t35 * t45) * t117 + (-t121 + t122) * t36 + (0.3e1 * t35 * t19 + t123) * t39) * t39 + ((-t82 - t19) * t36 * t39 + (0.3e1 * t36 * t18 + t123) * t38 + (t120 - t122 + (-t18 + t81) * t39) * t35 + (t25 * t1 - t45 * t36) * t117) * t38;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
