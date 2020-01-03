% Calculate time derivative of joint inertia matrix for
% S5PPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR5_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR5_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:27
% DurationCPUTime: 1.06s
% Computational Cost: add. (1138->143), mult. (2849->235), div. (0->0), fcn. (2860->6), ass. (0->77)
t108 = sin(qJ(3));
t109 = cos(qJ(3));
t94 = sin(pkin(7));
t95 = cos(pkin(7));
t48 = t95 * t108 - t94 * t109;
t40 = t48 * qJD(3);
t64 = sin(qJ(5));
t47 = -t108 * t94 - t109 * t95;
t65 = cos(qJ(5));
t92 = qJD(5) * t65;
t87 = t47 * t92;
t74 = -t40 * t64 - t87;
t81 = rSges(6,1) * t64 + rSges(6,2) * t65;
t41 = t47 * qJD(3);
t93 = qJD(5) * t64;
t88 = t48 * t93;
t71 = t41 * t65 + t88;
t86 = t48 * t92;
t72 = -t41 * t64 + t86;
t103 = t40 * t65;
t89 = t47 * t93;
t73 = t89 - t103;
t118 = -t47 * (Icges(6,5) * t72 - Icges(6,6) * t71 - Icges(6,3) * t40) - t48 * (Icges(6,5) * t74 + Icges(6,6) * t73 + Icges(6,3) * t41);
t98 = Icges(6,4) * t64;
t76 = Icges(6,2) * t65 + t98;
t23 = -Icges(6,6) * t47 + t48 * t76;
t97 = Icges(6,4) * t65;
t77 = Icges(6,1) * t64 + t97;
t25 = -Icges(6,5) * t47 + t48 * t77;
t80 = t23 * t65 + t25 * t64;
t114 = t80 * t47;
t24 = -Icges(6,6) * t48 - t47 * t76;
t26 = -Icges(6,5) * t48 - t47 * t77;
t79 = t24 * t65 + t26 * t64;
t115 = t79 * t48;
t75 = Icges(6,5) * t64 + Icges(6,6) * t65;
t21 = -Icges(6,3) * t47 + t48 * t75;
t22 = -Icges(6,3) * t48 - t47 * t75;
t117 = -t48 * t21 - t47 * t22 - t114 + t115;
t70 = t48 * t81;
t112 = 2 * m(6);
t54 = t81 * qJD(5);
t82 = rSges(6,1) * t65 - rSges(6,2) * t64;
t113 = t82 * t54 * t112;
t111 = 0.2e1 * t82;
t110 = -rSges(5,2) + pkin(3);
t105 = t40 * t47;
t99 = qJ(4) * t41 - qJD(4) * t48;
t96 = -rSges(5,3) - qJ(4);
t91 = rSges(6,3) + pkin(3) + pkin(6);
t90 = m(5) / 0.2e1 + m(6) / 0.2e1;
t28 = -t48 * rSges(6,3) - t47 * t81;
t85 = rSges(6,1) * t74 - rSges(6,2) * t103 + t41 * rSges(6,3);
t78 = -t41 * t48 + t105;
t69 = t40 * t95 - t41 * t94;
t45 = t48 * pkin(3);
t14 = -pkin(6) * t48 - qJ(4) * t47 + t28 - t45;
t43 = t48 * qJ(4);
t15 = t47 * t91 - t43 - t70;
t6 = -rSges(6,1) * t72 + rSges(6,2) * t71 + t40 * t91 + t99;
t39 = t41 * pkin(3);
t7 = t41 * pkin(6) - t40 * qJ(4) + t39 + (rSges(6,2) * t93 - qJD(4)) * t47 + t85;
t68 = -t14 * t41 - t15 * t40 - t47 * t6 + t48 * t7;
t58 = -Icges(6,1) * t65 + t98;
t57 = Icges(6,2) * t64 - t97;
t55 = t82 ^ 2;
t51 = t75 * qJD(5);
t46 = t47 ^ 2;
t31 = rSges(4,1) * t41 + rSges(4,2) * t40;
t30 = rSges(4,1) * t40 - rSges(4,2) * t41;
t27 = -t47 * rSges(6,3) + t70;
t19 = -rSges(5,3) * t48 + t110 * t47 - t43;
t18 = rSges(5,2) * t48 + t47 * t96 - t45;
t17 = -rSges(5,2) * t41 - t47 * qJD(4) + t40 * t96 + t39;
t16 = rSges(5,3) * t41 + t110 * t40 + t99;
t1 = -t40 * t28 - t47 * (rSges(6,2) * t89 + t85) + (qJD(5) * t48 * t82 - t40 * rSges(6,3)) * t48 - (t27 + t70) * t41;
t2 = [0; 0; 0; 0; m(4) * (-t30 * t95 + t31 * t94) + m(5) * (-t16 * t95 + t17 * t94) + m(6) * (-t6 * t95 + t7 * t94); t76 * t93 - t77 * t92 + (t57 * t65 + t58 * t64) * qJD(5) + 0.2e1 * m(4) * ((-rSges(4,1) * t48 + rSges(4,2) * t47) * t31 + (rSges(4,1) * t47 + rSges(4,2) * t48) * t30) + 0.2e1 * m(5) * (t16 * t19 + t17 * t18) + (t14 * t7 + t15 * t6) * t112; 0; 0.2e1 * t90 * t69; m(5) * (-t16 * t47 + t17 * t48 - t18 * t41 - t19 * t40) + m(6) * t68; 0.4e1 * t90 * t78; m(6) * t1; m(6) * (t69 * t82 + (-t47 * t95 - t48 * t94) * t54); (t80 * qJD(5) + t64 * (Icges(6,4) * t72 - Icges(6,2) * t71 - Icges(6,6) * t40) - (Icges(6,1) * t72 - Icges(6,4) * t71 - Icges(6,5) * t40) * t65) * t47 / 0.2e1 + (t64 * t23 - t25 * t65) * t40 / 0.2e1 + (t79 * qJD(5) + t64 * (Icges(6,4) * t74 + Icges(6,2) * t73 + Icges(6,6) * t41) - (Icges(6,1) * t74 + Icges(6,4) * t73 + Icges(6,5) * t41) * t65) * t48 / 0.2e1 - (t64 * t24 - t26 * t65) * t41 / 0.2e1 - t47 * (t47 * t51 + t57 * t88 - t58 * t86) / 0.2e1 - t48 * (t48 * t51 - t57 * t89 + t58 * t87) / 0.2e1 + m(6) * ((-t14 * t48 + t15 * t47) * t54 + t68 * t82) - t78 * (-Icges(6,5) * t65 + Icges(6,6) * t64); m(6) * (t105 * t111 - t46 * t54 + (-t111 * t41 - t48 * t54) * t48); -t46 * t113 + ((t1 * t27 - t41 * t55) * t112 + (-t115 + t117) * t40 + (0.3e1 * t41 * t22 - t113 + t118) * t48) * t48 + ((t80 - t22) * t40 * t48 + (-0.3e1 * t40 * t21 + t118) * t47 + (-t114 - t117 + (t21 + t79) * t48) * t41 + (-t1 * t28 + t40 * t55) * t112) * t47;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
