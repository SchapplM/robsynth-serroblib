% Calculate time derivative of joint inertia matrix for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-01-31 13:16:39
% EndTime: 2019-01-31 13:16:40
% DurationCPUTime: 0.18s
% Computational Cost: add. (730->59), mult. (452->82), div. (0->0), fcn. (236->8), ass. (0->51)
t48 = qJ(1) + qJ(2);
t43 = pkin(7) + t48;
t41 = qJ(4) + t43;
t35 = sin(t41);
t36 = cos(t41);
t24 = t36 * rSges(5,1) - t35 * rSges(5,2);
t44 = sin(t48);
t45 = cos(t48);
t29 = t45 * rSges(3,1) - t44 * rSges(3,2);
t38 = sin(t43);
t39 = cos(t43);
t40 = pkin(2) * t45;
t62 = -t39 * rSges(4,1) + t38 * rSges(4,2) - t40;
t60 = -pkin(3) * t39 - t40;
t59 = pkin(2) * t44;
t49 = sin(qJ(1));
t58 = t49 * pkin(1);
t54 = pkin(1) * qJD(1);
t47 = qJD(1) + qJD(2);
t53 = t49 * t54;
t50 = cos(qJ(1));
t52 = t50 * t54;
t22 = t29 * t47;
t42 = qJD(4) + t47;
t14 = t24 * t42;
t51 = -pkin(3) * t38 - t59;
t28 = -t44 * rSges(3,1) - t45 * rSges(3,2);
t23 = -t35 * rSges(5,1) - t36 * rSges(5,2);
t12 = t24 - t60;
t21 = t28 * t47;
t13 = t23 * t42;
t19 = -t38 * rSges(4,1) - t39 * rSges(4,2) - t59;
t10 = t62 * t47;
t9 = t19 * t47;
t11 = t23 + t51;
t4 = t60 * t47 - t14;
t3 = t51 * t47 + t13;
t46 = t50 * pkin(1);
t26 = t29 + t46;
t25 = t28 - t58;
t18 = -t22 - t52;
t17 = t21 - t53;
t16 = -t62 + t46;
t15 = t19 - t58;
t8 = t12 + t46;
t7 = t11 - t58;
t6 = t10 - t52;
t5 = t9 - t53;
t2 = t4 - t52;
t1 = t3 - t53;
t20 = [0.2e1 * m(3) * (t26 * t17 + t25 * t18) + 0.2e1 * m(4) * (t15 * t6 + t16 * t5) + 0.2e1 * m(5) * (t8 * t1 + t7 * t2); m(3) * (t29 * t17 + t28 * t18 + t21 * t26 - t22 * t25) + m(4) * (t10 * t15 + t9 * t16 + t19 * t6 - t5 * t62) + m(5) * (t12 * t1 + t11 * t2 + t3 * t8 + t4 * t7); 0.2e1 * m(3) * (t29 * t21 - t28 * t22) + 0.2e1 * m(4) * (t19 * t10 - t62 * t9) + 0.2e1 * m(5) * (t11 * t4 + t12 * t3); 0; 0; 0; m(5) * (t24 * t1 + t13 * t8 - t14 * t7 + t23 * t2); m(5) * (-t14 * t11 + t13 * t12 + t23 * t4 + t24 * t3); 0; 0.2e1 * m(5) * (t24 * t13 - t23 * t14);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t20(1) t20(2) t20(4) t20(7); t20(2) t20(3) t20(5) t20(8); t20(4) t20(5) t20(6) t20(9); t20(7) t20(8) t20(9) t20(10);];
Mq  = res;
