% Calculate time derivative of joint inertia matrix for
% S4RPPR4
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:50
% DurationCPUTime: 0.71s
% Computational Cost: add. (1057->124), mult. (1316->195), div. (0->0), fcn. (968->6), ass. (0->76)
t54 = sin(qJ(4));
t56 = cos(qJ(4));
t75 = rSges(5,1) * t54 + rSges(5,2) * t56;
t53 = qJ(1) + pkin(6);
t51 = cos(t53);
t64 = t51 * t75;
t50 = sin(t53);
t90 = Icges(5,4) * t54;
t67 = Icges(5,2) * t56 + t90;
t22 = Icges(5,6) * t51 + t67 * t50;
t89 = Icges(5,4) * t56;
t69 = Icges(5,1) * t54 + t89;
t24 = Icges(5,5) * t51 + t69 * t50;
t72 = t22 * t56 + t24 * t54;
t110 = t72 * t51;
t101 = sin(qJ(1)) * pkin(1);
t80 = -rSges(5,3) - pkin(2) - pkin(5);
t62 = t80 * t50 - t101;
t81 = qJD(4) * t56;
t84 = qJD(1) * t51;
t79 = t50 * rSges(5,1) * t81 + t75 * t84;
t82 = qJD(4) * t54;
t91 = qJ(3) * t84 + qJD(3) * t50;
t6 = -t50 * rSges(5,2) * t82 + t62 * qJD(1) + t79 + t91;
t96 = t54 * rSges(5,2);
t38 = t56 * rSges(5,1) - t96;
t43 = qJD(3) * t51;
t52 = cos(qJ(1)) * pkin(1);
t83 = qJD(4) * t51;
t7 = t43 + t38 * t83 + (-t52 + t80 * t51 + (-qJ(3) - t75) * t50) * qJD(1);
t109 = t50 * t7 - t51 * t6;
t65 = Icges(5,5) * t54 + Icges(5,6) * t56;
t108 = -Icges(5,3) * t50 + t65 * t51;
t107 = -Icges(5,6) * t50 + t67 * t51;
t106 = -Icges(5,5) * t50 + t69 * t51;
t105 = 2 * m(5);
t48 = t50 ^ 2;
t49 = t51 ^ 2;
t104 = m(5) * t38;
t100 = rSges(4,2) - pkin(2);
t97 = t50 * rSges(5,3);
t46 = t51 * rSges(5,3);
t95 = t54 * t22;
t94 = t54 * t107;
t93 = t56 * t24;
t92 = t56 * t106;
t20 = Icges(5,3) * t51 + t65 * t50;
t85 = qJD(1) * t20;
t26 = t75 * t50 + t46;
t78 = t51 * pkin(2) + t50 * qJ(3) + t52;
t31 = t75 * qJD(4);
t77 = (t48 + t49) * t31;
t71 = -t106 * t54 - t107 * t56;
t70 = Icges(5,1) * t56 - t90;
t68 = -Icges(5,2) * t54 + t89;
t66 = Icges(5,5) * t56 - Icges(5,6) * t54;
t63 = t71 * t50;
t61 = qJD(4) * t70;
t60 = qJD(4) * t68;
t58 = t51 * rSges(4,3) + t100 * t50 - t101;
t45 = t51 * qJ(3);
t27 = t97 - t64;
t19 = -t51 * rSges(4,2) + t50 * rSges(4,3) + t78;
t18 = t45 + t58;
t17 = t43 + (-t52 + t100 * t51 + (-rSges(4,3) - qJ(3)) * t50) * qJD(1);
t16 = t58 * qJD(1) + t91;
t15 = t51 * pkin(5) + t26 + t78;
t14 = t45 + t64 + t62;
t9 = t66 * t50 * qJD(4) + t108 * qJD(1);
t8 = -t66 * t83 + t85;
t5 = -t108 * t50 - t71 * t51;
t4 = t50 * t20 - t110;
t3 = -t108 * t51 + t63;
t2 = t51 * t20 + t72 * t50;
t1 = -t50 * t79 + (-t38 * t49 + t48 * t96) * qJD(4) + ((-t26 + t46) * t51 + (-t27 + t64 + t97) * t50) * qJD(1);
t10 = [0.2e1 * m(4) * (t19 * t16 + t18 * t17) - t54 * t61 - t69 * t81 - t56 * t60 + t67 * t82 + (t14 * t7 + t15 * t6) * t105; 0; 0; m(4) * (-t51 * t16 + t50 * t17 + (t18 * t51 + t19 * t50) * qJD(1)) + m(5) * ((t14 * t51 + t15 * t50) * qJD(1) + t109); 0; 0; (-t72 * qJD(4) - t54 * (t107 * qJD(1) + t50 * t60) + t56 * (t106 * qJD(1) + t50 * t61)) * t51 / 0.2e1 + (-t71 * qJD(4) - t54 * (t22 * qJD(1) - t68 * t83) + t56 * (t24 * qJD(1) - t70 * t83)) * t50 / 0.2e1 + m(5) * (t109 * t38 - (t14 * t50 - t15 * t51) * t31) - (t49 / 0.2e1 + t48 / 0.2e1) * t65 * qJD(4) + ((t95 / 0.2e1 - t93 / 0.2e1 + t15 * t104) * t50 + (t94 / 0.2e1 - t92 / 0.2e1 + t14 * t104) * t51) * qJD(1); m(5) * t1; -m(5) * t77; ((-t50 * t26 + t51 * t27) * t1 - t38 * t77) * t105 - qJD(1) * t50 * (t2 * t51 + t3 * t50) + t51 * ((t51 * t9 + (t3 + t110) * qJD(1)) * t51 + (-t2 * qJD(1) + (-t106 * t81 + t107 * t82) * t50 + (t8 + (t93 - t95) * qJD(4) + (-t20 + t71) * qJD(1)) * t51) * t50) + (t4 * t51 + t5 * t50) * t84 + t50 * ((t50 * t8 + (-t4 + t63) * qJD(1)) * t50 + (t5 * qJD(1) + (t22 * t82 - t24 * t81 + t85) * t51 + (t9 + (t92 - t94) * qJD(4) + t72 * qJD(1)) * t50) * t51);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
