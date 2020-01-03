% Calculate time derivative of joint inertia matrix for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:09
% DurationCPUTime: 0.87s
% Computational Cost: add. (1776->158), mult. (1636->231), div. (0->0), fcn. (1190->8), ass. (0->95)
t60 = pkin(8) + qJ(5);
t57 = cos(t60);
t115 = rSges(6,2) * t57;
t55 = sin(t60);
t85 = rSges(6,1) * t55 + t115;
t61 = qJ(1) + pkin(7);
t58 = cos(t61);
t132 = t58 * t85;
t56 = sin(t61);
t106 = Icges(6,4) * t55;
t77 = Icges(6,2) * t57 + t106;
t26 = Icges(6,6) * t58 + t77 * t56;
t105 = Icges(6,4) * t57;
t79 = Icges(6,1) * t55 + t105;
t28 = Icges(6,5) * t58 + t79 * t56;
t82 = t26 * t57 + t28 * t55;
t131 = t82 * t58;
t50 = t58 * qJ(3);
t62 = sin(pkin(8));
t123 = pkin(4) * t62;
t69 = t85 + t123;
t120 = sin(qJ(1)) * pkin(1);
t64 = -pkin(6) - qJ(4);
t94 = -rSges(6,3) - pkin(2) + t64;
t71 = t94 * t56 - t120;
t10 = t69 * t58 + t50 + t71;
t117 = rSges(6,1) * t56;
t51 = t58 * rSges(6,3);
t30 = t56 * t115 + t55 * t117 + t51;
t59 = cos(qJ(1)) * pkin(1);
t90 = t58 * pkin(2) + t56 * qJ(3) + t59;
t11 = t56 * t123 - t58 * t64 + t30 + t90;
t130 = -t10 * t56 + t11 * t58;
t116 = rSges(6,2) * t55;
t99 = qJD(1) * t58;
t107 = qJ(3) * t99 + qJD(3) * t56;
t91 = qJD(4) * t58 + t107;
t96 = qJD(5) * t57;
t92 = t96 * t117 + t85 * t99;
t97 = qJD(5) * t56;
t6 = -t97 * t116 + (t58 * t123 + t71) * qJD(1) + t91 + t92;
t39 = rSges(6,1) * t57 - t116;
t48 = qJD(3) * t58;
t88 = -qJD(4) * t56 + t48;
t95 = qJD(5) * t58;
t7 = t39 * t95 + (-t59 + t94 * t58 + (-qJ(3) - t69) * t56) * qJD(1) + t88;
t129 = t56 * t7 - t58 * t6;
t75 = Icges(6,5) * t55 + Icges(6,6) * t57;
t128 = -Icges(6,3) * t56 + t75 * t58;
t127 = -Icges(6,6) * t56 + t77 * t58;
t126 = -Icges(6,5) * t56 + t79 * t58;
t125 = 2 * m(6);
t53 = t56 ^ 2;
t54 = t58 ^ 2;
t124 = m(6) * t39;
t119 = rSges(4,2) - pkin(2);
t112 = t55 * t26;
t111 = t55 * t127;
t110 = t56 * rSges(6,3);
t109 = t57 * t28;
t108 = t57 * t126;
t101 = rSges(5,3) + qJ(4);
t24 = Icges(6,3) * t58 + t75 * t56;
t100 = qJD(1) * t24;
t98 = qJD(5) * t55;
t93 = -pkin(2) - t101;
t35 = t85 * qJD(5);
t89 = (t53 + t54) * t35;
t86 = rSges(5,1) * t62 + rSges(5,2) * cos(pkin(8));
t81 = -t126 * t55 - t127 * t57;
t80 = Icges(6,1) * t57 - t106;
t78 = -Icges(6,2) * t55 + t105;
t76 = Icges(6,5) * t57 - Icges(6,6) * t55;
t74 = t81 * t56;
t73 = qJD(5) * t80;
t72 = qJD(5) * t78;
t68 = t58 * rSges(4,3) + t119 * t56 - t120;
t67 = t93 * t56 + t86 * t58 - t120;
t31 = t110 - t132;
t23 = -rSges(4,2) * t58 + rSges(4,3) * t56 + t90;
t22 = t50 + t68;
t21 = t48 + (-t59 + t119 * t58 + (-rSges(4,3) - qJ(3)) * t56) * qJD(1);
t20 = t68 * qJD(1) + t107;
t19 = t101 * t58 + t86 * t56 + t90;
t18 = t50 + t67;
t13 = t128 * qJD(1) + t76 * t97;
t12 = -t76 * t95 + t100;
t9 = (-t59 + t93 * t58 + (-qJ(3) - t86) * t56) * qJD(1) + t88;
t8 = t67 * qJD(1) + t91;
t5 = -t128 * t56 - t81 * t58;
t4 = t56 * t24 - t131;
t3 = -t128 * t58 + t74;
t2 = t58 * t24 + t82 * t56;
t1 = -t56 * t92 + (t53 * t116 - t39 * t54) * qJD(5) + ((-t30 + t51) * t58 + (t110 - t31 + t132) * t56) * qJD(1);
t14 = [0.2e1 * m(4) * (t20 * t23 + t21 * t22) + 0.2e1 * m(5) * (t18 * t9 + t19 * t8) - t55 * t73 - t79 * t96 - t57 * t72 + t77 * t98 + (t10 * t7 + t11 * t6) * t125; 0; 0; m(4) * (-t58 * t20 + t56 * t21 + (t22 * t58 + t23 * t56) * qJD(1)) + m(5) * (t56 * t9 - t58 * t8 + (t18 * t58 + t19 * t56) * qJD(1)) + m(6) * ((t10 * t58 + t11 * t56) * qJD(1) + t129); 0; 0; m(5) * (t56 * t8 + t58 * t9 + (-t18 * t56 + t19 * t58) * qJD(1)) + m(6) * (t130 * qJD(1) + t56 * t6 + t58 * t7); 0; 0; 0; (-t82 * qJD(5) - t55 * (t127 * qJD(1) + t56 * t72) + t57 * (t126 * qJD(1) + t56 * t73)) * t58 / 0.2e1 + (-t81 * qJD(5) - t55 * (t26 * qJD(1) - t78 * t95) + t57 * (t28 * qJD(1) - t80 * t95)) * t56 / 0.2e1 + m(6) * (t129 * t39 + t130 * t35) - (t54 / 0.2e1 + t53 / 0.2e1) * t75 * qJD(5) + ((t112 / 0.2e1 - t109 / 0.2e1 + t11 * t124) * t56 + (t111 / 0.2e1 - t108 / 0.2e1 + t10 * t124) * t58) * qJD(1); m(6) * t1; -m(6) * t89; 0; ((-t30 * t56 + t31 * t58) * t1 - t39 * t89) * t125 - qJD(1) * t56 * (t2 * t58 + t3 * t56) + t58 * ((t58 * t13 + (t3 + t131) * qJD(1)) * t58 + (-t2 * qJD(1) + (-t126 * t96 + t127 * t98) * t56 + (t12 + (t109 - t112) * qJD(5) + (-t24 + t81) * qJD(1)) * t58) * t56) + (t4 * t58 + t5 * t56) * t99 + t56 * ((t56 * t12 + (-t4 + t74) * qJD(1)) * t56 + (t5 * qJD(1) + (t26 * t98 - t28 * t96 + t100) * t58 + (t13 + (t108 - t111) * qJD(5) + t82 * qJD(1)) * t56) * t58);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
