% Calculate time derivative of joint inertia matrix for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:20
% EndTime: 2019-03-08 18:26:20
% DurationCPUTime: 0.29s
% Computational Cost: add. (531->104), mult. (1484->150), div. (0->0), fcn. (1354->6), ass. (0->60)
t38 = sin(pkin(6));
t40 = sin(qJ(1));
t41 = cos(qJ(1));
t59 = cos(pkin(6));
t60 = cos(pkin(4));
t47 = t60 * t59;
t43 = t41 * t47;
t24 = t40 * t38 - t43;
t39 = sin(pkin(4));
t64 = qJ(2) * t39;
t36 = t41 * t64;
t75 = -t24 * qJ(3) + t36;
t50 = t38 * t60;
t27 = -t40 * t50 + t41 * t59;
t74 = 0.2e1 * t39;
t73 = m(4) / 0.2e1;
t72 = m(5) / 0.2e1;
t71 = t40 * pkin(1);
t37 = t41 * pkin(1);
t70 = rSges(5,1) + pkin(3);
t69 = rSges(4,2) - pkin(2);
t68 = t39 * t40;
t67 = t39 * t41;
t56 = qJD(2) * t39;
t66 = qJD(1) * t36 + t40 * t56;
t65 = t40 * t64 + t37;
t63 = rSges(5,2) + qJ(3);
t62 = rSges(4,3) + qJ(3);
t61 = rSges(5,3) + qJ(4);
t58 = qJD(1) * t40;
t57 = qJD(1) * t41;
t55 = t73 + t72;
t54 = -pkin(2) - t61;
t53 = t27 * pkin(2) + t65;
t26 = t41 * t38 + t40 * t47;
t19 = qJD(1) * t26;
t34 = t41 * t56;
t52 = -t19 * qJ(3) - t24 * qJD(3) + t34;
t51 = t39 * t70;
t25 = t40 * t59 + t41 * t50;
t18 = qJD(1) * t25;
t46 = -t18 * pkin(2) + t26 * qJD(3) + t66;
t45 = rSges(4,1) * t67 - t71;
t44 = rSges(3,3) * t67 - t71;
t42 = t41 * t51 - t71;
t20 = t27 * qJD(1);
t17 = -qJD(1) * t43 + t38 * t58;
t14 = t27 * rSges(3,1) - t26 * rSges(3,2) + rSges(3,3) * t68 + t65;
t13 = -t25 * rSges(3,1) + t24 * rSges(3,2) + t36 + t44;
t11 = -t20 * rSges(3,1) + t19 * rSges(3,2) + t34 + (-t37 + (-rSges(3,3) - qJ(2)) * t68) * qJD(1);
t10 = -t18 * rSges(3,1) + t17 * rSges(3,2) + t44 * qJD(1) + t66;
t9 = rSges(4,1) * t68 - t27 * rSges(4,2) + t62 * t26 + t53;
t8 = -t24 * rSges(4,3) + t69 * t25 + t45 + t75;
t6 = t63 * t26 + t61 * t27 + t40 * t51 + t53;
t5 = -t24 * rSges(5,2) + t54 * t25 + t42 + t75;
t4 = -t19 * rSges(4,3) + t69 * t20 + (-t37 + (-rSges(4,1) - qJ(2)) * t68) * qJD(1) + t52;
t3 = t18 * rSges(4,2) + t45 * qJD(1) - t62 * t17 + t46;
t2 = -t19 * rSges(5,2) - t25 * qJD(4) + t54 * t20 + (-t37 + (-qJ(2) - t70) * t68) * qJD(1) + t52;
t1 = t42 * qJD(1) + t27 * qJD(4) - t63 * t17 - t61 * t18 + t46;
t7 = [0.2e1 * m(3) * (t14 * t10 + t13 * t11) + 0.2e1 * m(4) * (t9 * t3 + t8 * t4) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); (m(3) * (-t10 * t41 + t11 * t40 + t13 * t57 + t14 * t58) / 0.2e1 + (-t3 * t41 + t4 * t40 + t8 * t57 + t9 * t58) * t73 + (-t1 * t41 + t2 * t40 + t5 * t57 + t6 * t58) * t72) * t74; 0; m(4) * (-t17 * t8 + t19 * t9 + t24 * t3 + t26 * t4) + m(5) * (t24 * t1 - t17 * t5 + t19 * t6 + t26 * t2); t55 * (-t17 * t40 - t19 * t41 + (t24 * t40 + t26 * t41) * qJD(1)) * t74; 0.4e1 * t55 * (-t26 * t17 + t24 * t19); m(5) * (t25 * t1 - t18 * t5 + t27 * t2 + t20 * t6); m(5) * (-t18 * t40 - t20 * t41 + (t25 * t40 + t27 * t41) * qJD(1)) * t39; m(5) * (-t27 * t17 - t18 * t26 + t25 * t19 + t20 * t24); 0.2e1 * m(5) * (-t27 * t18 + t25 * t20);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
