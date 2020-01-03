% Calculate time derivative of joint inertia matrix for
% S4PRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:48
% EndTime: 2019-12-31 16:20:50
% DurationCPUTime: 0.80s
% Computational Cost: add. (1416->115), mult. (1266->179), div. (0->0), fcn. (952->6), ass. (0->76)
t47 = pkin(7) + qJ(4);
t43 = sin(t47);
t100 = qJD(2) * t43;
t45 = cos(t47);
t99 = qJD(2) * t45;
t48 = pkin(6) + qJ(2);
t44 = sin(t48);
t35 = rSges(5,1) * t43 + rSges(5,2) * t45;
t56 = t35 * qJD(4);
t98 = t44 * t56;
t46 = cos(t48);
t79 = Icges(5,4) * t45;
t61 = -Icges(5,2) * t43 + t79;
t23 = Icges(5,6) * t44 + t61 * t46;
t80 = Icges(5,4) * t43;
t63 = Icges(5,1) * t45 - t80;
t25 = Icges(5,5) * t44 + t63 * t46;
t64 = t23 * t43 - t25 * t45;
t97 = t64 * t44;
t22 = -Icges(5,6) * t46 + t61 * t44;
t24 = -Icges(5,5) * t46 + t63 * t44;
t65 = t22 * t43 - t24 * t45;
t96 = t65 * t46;
t50 = cos(pkin(7));
t42 = pkin(3) * t50 + pkin(2);
t89 = rSges(5,2) * t43;
t90 = rSges(5,1) * t45;
t68 = -t89 + t90;
t57 = -t42 - t68;
t51 = -pkin(5) - qJ(3);
t82 = rSges(5,3) - t51;
t16 = t57 * t44 + t82 * t46;
t41 = t44 * rSges(5,3);
t81 = t46 * t90 + t41;
t17 = -t44 * t51 + (t42 - t89) * t46 + t81;
t95 = t16 * t46 + t17 * t44;
t59 = Icges(5,5) * t45 - Icges(5,6) * t43;
t20 = -Icges(5,3) * t46 + t59 * t44;
t58 = rSges(4,1) * t50 - rSges(4,2) * sin(pkin(7)) + pkin(2);
t75 = rSges(4,3) + qJ(3);
t19 = t75 * t44 + t58 * t46;
t94 = 2 * m(5);
t93 = t44 ^ 2;
t92 = t46 ^ 2;
t91 = m(5) * t35;
t86 = t43 * t24;
t85 = t43 * t25;
t84 = t45 * t22;
t83 = t45 * t23;
t21 = Icges(5,3) * t44 + t59 * t46;
t74 = qJD(2) * t21;
t73 = qJD(2) * t44;
t72 = qJD(2) * t46;
t71 = qJD(4) * t43;
t70 = qJD(4) * t45;
t69 = t46 * t89;
t53 = qJD(4) * (-Icges(5,5) * t43 - Icges(5,6) * t45);
t52 = rSges(5,3) * t72 - t46 * t56 + t73 * t89;
t18 = -t58 * t44 + t75 * t46;
t40 = qJD(3) * t46;
t39 = qJD(3) * t44;
t31 = t68 * qJD(4);
t27 = -t69 + t81;
t26 = -t46 * rSges(5,3) + t68 * t44;
t15 = -qJD(2) * t19 + t40;
t14 = t18 * qJD(2) + t39;
t9 = t44 * t53 + t74;
t8 = -qJD(2) * t20 + t46 * t53;
t7 = t40 + t98 + (-t82 * t44 + t57 * t46) * qJD(2);
t6 = t39 + (-t46 * t51 + (-t42 - t90) * t44) * qJD(2) + t52;
t5 = t44 * t21 - t64 * t46;
t4 = t44 * t20 - t96;
t3 = -t46 * t21 - t97;
t2 = -t46 * t20 - t65 * t44;
t1 = (qJD(2) * t26 + t52) * t46 + (-t98 + (-t27 - t69 + t41) * qJD(2)) * t44;
t10 = [0; 0; 0.2e1 * m(4) * (t14 * t19 + t15 * t18) + (t16 * t7 + t17 * t6) * t94 + (-Icges(5,2) * t45 + t63 - t80) * t71 + (Icges(5,1) * t43 + t61 + t79) * t70; 0; m(4) * (-t46 * t14 + t44 * t15 + (t18 * t46 + t19 * t44) * qJD(2)) + m(5) * (qJD(2) * t95 + t44 * t7 - t46 * t6); 0; m(5) * t1; (-t64 * qJD(4) - t24 * t100 - t22 * t99) * t44 / 0.2e1 - (-t65 * qJD(4) + t25 * t100 + t23 * t99) * t46 / 0.2e1 + m(5) * ((-t44 * t6 - t46 * t7) * t35 - t95 * t31) + (t93 / 0.2e1 + t92 / 0.2e1) * t59 * qJD(4) + ((t83 / 0.2e1 + t85 / 0.2e1 - t17 * t91) * t46 + (t84 / 0.2e1 + t86 / 0.2e1 + t16 * t91) * t44) * qJD(2); 0; ((t26 * t44 + t27 * t46) * t1 + (t92 + t93) * t35 * t31) * t94 + (-t4 * t46 + t5 * t44) * t72 + t44 * ((t44 * t8 + (t4 + t97) * qJD(2)) * t44 + (t5 * qJD(2) + (t22 * t70 + t24 * t71) * t46 + (-t9 + (-t83 - t85) * qJD(4) + (t21 - t65) * qJD(2)) * t44) * t46) + (-t2 * t46 + t3 * t44) * t73 - t46 * ((t46 * t9 + (t3 + t96) * qJD(2)) * t46 + (t2 * qJD(2) + (-t23 * t70 - t25 * t71 + t74) * t44 + (-t8 + (t84 + t86) * qJD(4) - t64 * qJD(2)) * t46) * t44);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
