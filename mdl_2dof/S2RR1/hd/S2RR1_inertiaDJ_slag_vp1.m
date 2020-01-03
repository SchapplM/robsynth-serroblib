% Calculate time derivative of joint inertia matrix for
% S2RR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_inertiaDJ_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR1_inertiaDJ_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR1_inertiaDJ_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:07
% EndTime: 2020-01-03 11:19:08
% DurationCPUTime: 0.84s
% Computational Cost: add. (370->91), mult. (1022->152), div. (0->0), fcn. (786->4), ass. (0->61)
t33 = sin(qJ(2));
t78 = qJD(1) * t33;
t35 = cos(qJ(2));
t77 = qJD(1) * t35;
t55 = qJD(2) * t33;
t54 = qJD(2) * t35;
t34 = sin(qJ(1));
t67 = rSges(3,2) * t33;
t68 = rSges(3,1) * t35;
t51 = -t67 + t68;
t76 = t34 * t51;
t69 = rSges(3,3) + pkin(1);
t75 = t69 * t34;
t50 = rSges(3,1) * t33 + rSges(3,2) * t35;
t74 = t50 * qJD(2);
t36 = cos(qJ(1));
t44 = Icges(3,5) * t35 - Icges(3,6) * t33;
t15 = Icges(3,3) * t36 + t44 * t34;
t61 = Icges(3,4) * t35;
t45 = -Icges(3,2) * t33 + t61;
t17 = Icges(3,6) * t36 + t45 * t34;
t62 = Icges(3,4) * t33;
t46 = Icges(3,1) * t35 - t62;
t19 = Icges(3,5) * t36 + t46 * t34;
t73 = 2 * m(3);
t72 = t34 ^ 2;
t71 = t36 ^ 2;
t70 = m(3) * t50;
t66 = t17 * t35;
t18 = -Icges(3,6) * t34 + t45 * t36;
t65 = t18 * t35;
t64 = t19 * t33;
t20 = -Icges(3,5) * t34 + t46 * t36;
t63 = t20 * t33;
t16 = -Icges(3,3) * t34 + t44 * t36;
t57 = qJD(1) * t16;
t56 = qJD(1) * t34;
t53 = t36 * t67;
t48 = t17 * t33 - t19 * t35;
t47 = t18 * t33 - t20 * t35;
t43 = t34 * rSges(3,3) + t53;
t42 = t48 * t36;
t41 = t47 * t34;
t38 = qJD(2) * (-Icges(3,5) * t33 - Icges(3,6) * t35);
t37 = t34 * t74;
t32 = t36 * t68;
t31 = t56 * t67;
t26 = t51 * qJD(2);
t22 = t32 - t43;
t21 = rSges(3,3) * t36 + t76;
t14 = -t69 * t36 - t76;
t13 = t32 - t53 - t75;
t8 = t34 * t38 + t57;
t7 = -qJD(1) * t15 + t36 * t38;
t6 = -t56 * t68 + t31 + (-t69 * qJD(1) - t74) * t36;
t5 = t37 + (-t51 * t36 + t75) * qJD(1);
t4 = -t34 * t16 - t47 * t36;
t3 = -t34 * t15 - t42;
t2 = t16 * t36 - t41;
t1 = t15 * t36 - t48 * t34;
t9 = [(t13 * t6 + t14 * t5) * t73 + (-Icges(3,2) * t35 + t46 - t62) * t55 + (Icges(3,1) * t33 + t45 + t61) * t54; m(3) * (-(t34 * t6 + t36 * t5) * t50 - (t13 * t34 + t14 * t36) * t26) + (t48 * qJD(2) - t18 * t77 - t20 * t78) * t36 / 0.2e1 - (t47 * qJD(2) + t17 * t77 + t19 * t78) * t34 / 0.2e1 - (t71 / 0.2e1 + t72 / 0.2e1) * t44 * qJD(2) + ((-t13 * t70 + t65 / 0.2e1 + t63 / 0.2e1) * t36 + (t14 * t70 + t66 / 0.2e1 + t64 / 0.2e1) * t34) * qJD(1); ((-t34 * t21 - t22 * t36) * ((-qJD(1) * t21 - t31 + (rSges(3,1) * t55 + rSges(3,2) * t54 + rSges(3,3) * qJD(1)) * t36) * t36 + (t37 + (t22 + t43) * qJD(1)) * t34) + (t71 + t72) * t50 * t26) * t73 - (t1 * t36 - t2 * t34) * t56 + t36 * ((t36 * t8 + (-t2 - t42) * qJD(1)) * t36 + (-t1 * qJD(1) + (t18 * t54 + t20 * t55 + t57) * t34 + (-t7 + (-t64 - t66) * qJD(2) + t47 * qJD(1)) * t36) * t34) - qJD(1) * t36 * (t3 * t36 - t4 * t34) - t34 * ((t34 * t7 + (-t3 - t41) * qJD(1)) * t34 + (-t4 * qJD(1) + (-t17 * t54 - t19 * t55) * t36 + (-t8 + (t63 + t65) * qJD(2) + (t16 + t48) * qJD(1)) * t34) * t36);];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t9(1), t9(2); t9(2), t9(3);];
Mq = res;
