% Calculate time derivative of joint inertia matrix for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR4_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR4_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:17
% EndTime: 2019-12-31 17:32:19
% DurationCPUTime: 1.35s
% Computational Cost: add. (1540->137), mult. (2808->232), div. (0->0), fcn. (2826->8), ass. (0->76)
t105 = sin(qJ(3));
t106 = cos(qJ(3));
t92 = sin(pkin(7));
t93 = cos(pkin(7));
t48 = -t92 * t105 - t93 * t106;
t49 = t93 * t105 - t92 * t106;
t58 = pkin(8) + qJ(5);
t56 = sin(t58);
t57 = cos(t58);
t71 = -Icges(6,5) * t57 + Icges(6,6) * t56;
t21 = -Icges(6,3) * t48 + t71 * t49;
t114 = t49 * t21;
t22 = Icges(6,3) * t49 + t71 * t48;
t115 = t48 * t22;
t116 = -t114 + t115;
t113 = t48 ^ 2;
t112 = t49 ^ 2;
t110 = t48 * t49;
t40 = t48 * qJD(3);
t91 = qJD(5) * t56;
t86 = t49 * t91;
t65 = t40 * t57 + t86;
t39 = t49 * qJD(3);
t87 = t48 * t91;
t67 = -t39 * t57 + t87;
t107 = 2 * m(6);
t104 = rSges(6,1) * t57;
t103 = rSges(6,2) * t56;
t61 = -pkin(6) - qJ(4);
t98 = -rSges(6,3) + t61;
t97 = -t48 * rSges(6,3) - t49 * t104;
t96 = Icges(6,4) * t56;
t95 = Icges(6,4) * t57;
t94 = -rSges(5,3) - qJ(4);
t90 = qJD(5) * t57;
t89 = t48 * qJD(4);
t88 = m(5) / 0.2e1 + m(6) / 0.2e1;
t85 = t48 * t90;
t84 = t49 * t90;
t79 = t103 - t104;
t47 = -rSges(6,1) * t56 - rSges(6,2) * t57;
t76 = t39 * t48 - t40 * t49;
t73 = -Icges(6,1) * t57 + t96;
t72 = Icges(6,2) * t56 - t95;
t60 = cos(pkin(8));
t70 = rSges(5,1) * t60 - rSges(5,2) * sin(pkin(8)) + pkin(3);
t55 = pkin(4) * t60 + pkin(3);
t69 = t55 - t79;
t68 = t39 * t56 + t85;
t66 = -t40 * t56 + t84;
t62 = t65 * rSges(6,1) - t39 * rSges(6,3);
t46 = -Icges(6,1) * t56 - t95;
t45 = -Icges(6,2) * t57 - t96;
t41 = t49 * qJD(4);
t38 = t79 * qJD(5);
t35 = t71 * qJD(5);
t31 = rSges(4,1) * t40 + rSges(4,2) * t39;
t30 = rSges(4,1) * t39 - rSges(4,2) * t40;
t28 = t49 * rSges(6,3) + t48 * t79;
t27 = t49 * t103 + t97;
t26 = Icges(6,5) * t49 + t73 * t48;
t25 = -Icges(6,5) * t48 + t73 * t49;
t24 = Icges(6,6) * t49 + t72 * t48;
t23 = -Icges(6,6) * t48 + t72 * t49;
t19 = t70 * t48 + t94 * t49;
t18 = t94 * t48 - t70 * t49;
t17 = t69 * t48 + t98 * t49;
t16 = t48 * t61 + (-t55 + t103) * t49 + t97;
t11 = t65 * Icges(6,5) + t66 * Icges(6,6) - Icges(6,3) * t39;
t10 = t67 * Icges(6,5) + t68 * Icges(6,6) - Icges(6,3) * t40;
t9 = t94 * t39 + t70 * t40 - t89;
t8 = t70 * t39 - t94 * t40 - t41;
t7 = t66 * rSges(6,2) + t39 * t61 + t40 * t55 + t62 - t89;
t6 = t47 * t48 * qJD(5) + t69 * t39 - t98 * t40 - t41;
t1 = -t40 * t27 + t49 * t62 + t39 * t28 + t48 * (t67 * rSges(6,1) - t40 * rSges(6,3)) + (t48 * t68 + t49 * t66) * rSges(6,2);
t2 = [0; 0; 0; 0; m(4) * (-t30 * t93 + t31 * t92) + m(5) * (-t8 * t93 + t9 * t92) + m(6) * (-t6 * t93 + t7 * t92); -t72 * t90 - t73 * t91 + (t56 * t45 - t57 * t46) * qJD(5) + (t16 * t7 + t17 * t6) * t107 + 0.2e1 * m(5) * (t18 * t9 + t19 * t8) + 0.2e1 * m(4) * ((-rSges(4,1) * t49 + rSges(4,2) * t48) * t31 + (rSges(4,1) * t48 + rSges(4,2) * t49) * t30); 0; 0.2e1 * t88 * (t39 * t93 - t40 * t92); m(6) * (-t16 * t40 - t17 * t39 - t48 * t6 + t49 * t7) + m(5) * (-t18 * t40 - t19 * t39 - t48 * t8 + t49 * t9); 0.4e1 * t88 * t76; m(6) * t1; m(6) * ((-t92 * t39 - t93 * t40) * t47 + (-t92 * t48 + t93 * t49) * t38); -((t24 * t56 - t26 * t57) * qJD(5) - (t67 * Icges(6,4) + t68 * Icges(6,2) - Icges(6,6) * t40) * t57 - (t67 * Icges(6,1) + t68 * Icges(6,4) - Icges(6,5) * t40) * t56) * t49 / 0.2e1 + (-t24 * t57 - t26 * t56) * t40 / 0.2e1 + ((t23 * t56 - t25 * t57) * qJD(5) - (t65 * Icges(6,4) + t66 * Icges(6,2) - Icges(6,6) * t39) * t57 - (t65 * Icges(6,1) + t66 * Icges(6,4) - Icges(6,5) * t39) * t56) * t48 / 0.2e1 + (-t23 * t57 - t25 * t56) * t39 / 0.2e1 - t48 * (t48 * t35 - t45 * t84 - t46 * t86) / 0.2e1 + m(6) * ((-t16 * t48 - t17 * t49) * t38 + (-t16 * t39 + t17 * t40 - t48 * t7 - t49 * t6) * t47) + t49 * (-t49 * t35 - t45 * t85 - t46 * t87) / 0.2e1 - t76 * (-Icges(6,5) * t56 - Icges(6,6) * t57); 0; ((t27 * t49 + t28 * t48) * t1 + ((t112 + t113) * t38 + t76 * t47) * t47) * t107 + t49 * (t112 * t10 - t11 * t110) - t48 * (-t10 * t110 + t113 * t11) + (-0.3e1 * t112 * t22 + (t114 - t116) * t48) * t40 + (-0.3e1 * t113 * t21 + (t115 + t116) * t49) * t39;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
