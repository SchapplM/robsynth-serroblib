% Calculate time derivative of joint inertia matrix for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:22
% EndTime: 2019-12-31 16:17:24
% DurationCPUTime: 0.88s
% Computational Cost: add. (909->112), mult. (2292->192), div. (0->0), fcn. (2340->6), ass. (0->56)
t76 = sin(pkin(6));
t77 = cos(pkin(6));
t88 = sin(qJ(3));
t89 = cos(qJ(3));
t32 = -t76 * t88 - t77 * t89;
t30 = t32 * qJD(3);
t33 = -t76 * t89 + t77 * t88;
t49 = cos(qJ(4));
t48 = sin(qJ(4));
t75 = qJD(4) * t48;
t53 = t30 * t49 + t33 * t75;
t94 = t33 ^ 2;
t58 = -Icges(5,5) * t49 + Icges(5,6) * t48;
t16 = -Icges(5,3) * t32 + t33 * t58;
t17 = Icges(5,3) * t33 + t32 * t58;
t78 = Icges(5,4) * t49;
t59 = Icges(5,2) * t48 - t78;
t19 = Icges(5,6) * t33 + t32 * t59;
t79 = Icges(5,4) * t48;
t60 = -Icges(5,1) * t49 + t79;
t21 = Icges(5,5) * t33 + t32 * t60;
t63 = t19 * t48 - t21 * t49;
t51 = t63 * t33;
t18 = -Icges(5,6) * t32 + t33 * t59;
t20 = -Icges(5,5) * t32 + t33 * t60;
t64 = t18 * t48 - t20 * t49;
t52 = t64 * t32;
t93 = -t33 * t16 + t32 * t17 - t51 - t52;
t29 = t33 * qJD(3);
t55 = -t29 * t49 + t32 * t75;
t86 = rSges(5,2) * t48;
t87 = rSges(5,1) * t49;
t65 = t86 - t87;
t38 = t65 * qJD(4);
t43 = -t48 * rSges(5,1) - rSges(5,2) * t49;
t74 = qJD(4) * t49;
t56 = t29 * t48 + t32 * t74;
t91 = 2 * m(5);
t92 = t33 * (Icges(5,5) * t55 + Icges(5,6) * t56 - Icges(5,3) * t30) + t38 * t43 * t91;
t90 = -rSges(5,3) - pkin(5);
t80 = -t32 * rSges(5,3) - t33 * t87;
t57 = pkin(3) - t65;
t54 = -t30 * t48 + t33 * t74;
t50 = t53 * rSges(5,1) - t29 * rSges(5,3);
t39 = t43 ^ 2;
t25 = rSges(4,1) * t30 + rSges(4,2) * t29;
t24 = rSges(4,1) * t29 - rSges(4,2) * t30;
t23 = t33 * rSges(5,3) + t32 * t65;
t22 = t33 * t86 + t80;
t15 = t32 * t57 + t33 * t90;
t14 = -pkin(5) * t32 + (-pkin(3) + t86) * t33 + t80;
t9 = Icges(5,5) * t53 + Icges(5,6) * t54 - Icges(5,3) * t29;
t7 = rSges(5,2) * t54 + t30 * pkin(3) - t29 * pkin(5) + t50;
t6 = qJD(4) * t32 * t43 + t29 * t57 - t30 * t90;
t1 = -t30 * t22 + t33 * t50 + t29 * t23 + t32 * (t55 * rSges(5,1) - t30 * rSges(5,3)) + (t32 * t56 + t33 * t54) * rSges(5,2);
t2 = [0; 0; 0; 0; m(4) * (-t24 * t77 + t25 * t76) + m(5) * (-t6 * t77 + t7 * t76); -t59 * t74 - t60 * t75 + ((-Icges(5,2) * t49 - t79) * t48 - (-Icges(5,1) * t48 - t78) * t49) * qJD(4) + 0.2e1 * m(4) * ((-rSges(4,1) * t33 + rSges(4,2) * t32) * t25 + (rSges(4,1) * t32 + rSges(4,2) * t33) * t24) + (t14 * t7 + t15 * t6) * t91; m(5) * t1; m(5) * ((-t29 * t76 - t30 * t77) * t43 + (-t32 * t76 + t33 * t77) * t38); m(5) * ((-t14 * t32 - t15 * t33) * t38 + (-t14 * t29 + t15 * t30 - t32 * t7 - t33 * t6) * t43) - (qJD(4) * t63 - t49 * (Icges(5,4) * t55 + Icges(5,2) * t56 - Icges(5,6) * t30) - t48 * (Icges(5,1) * t55 + Icges(5,4) * t56 - Icges(5,5) * t30)) * t33 / 0.2e1 + (-t49 * t19 - t48 * t21) * t30 / 0.2e1 + (qJD(4) * t64 - t49 * (Icges(5,4) * t53 + Icges(5,2) * t54 - Icges(5,6) * t29) - t48 * (Icges(5,1) * t53 + Icges(5,4) * t54 - Icges(5,5) * t29)) * t32 / 0.2e1 + (-t49 * t18 - t48 * t20) * t29 / 0.2e1 + (-t29 * t32 + t30 * t33) * (-Icges(5,5) * t48 - Icges(5,6) * t49) + (-t94 / 0.2e1 - t32 ^ 2 / 0.2e1) * t58 * qJD(4); ((t1 * t22 - t30 * t39) * t91 + (t51 + t93) * t29 + (-0.3e1 * t30 * t17 + t92) * t33) * t33 + (t1 * t23 * t91 - t94 * t9 + (-t32 * t9 + t92) * t32 + (-0.3e1 * t32 * t16 + t39 * t91 + (t64 + t17) * t33) * t29 + (-t52 - t93 + (t16 - t63) * t33) * t30) * t32;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
