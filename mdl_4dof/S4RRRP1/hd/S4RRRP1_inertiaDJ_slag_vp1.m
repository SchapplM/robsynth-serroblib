% Calculate time derivative of joint inertia matrix for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:53
% EndTime: 2019-03-08 18:35:53
% DurationCPUTime: 0.20s
% Computational Cost: add. (908->72), mult. (574->101), div. (0->0), fcn. (300->6), ass. (0->57)
t64 = rSges(5,1) + pkin(3);
t49 = qJ(1) + qJ(2);
t44 = sin(t49);
t63 = pkin(2) * t44;
t50 = sin(qJ(1));
t62 = t50 * pkin(1);
t46 = qJ(3) + t49;
t41 = sin(t46);
t48 = qJD(1) + qJD(2);
t43 = qJD(3) + t48;
t60 = t41 * t43;
t42 = cos(t46);
t59 = t42 * t43;
t58 = t44 * t48;
t45 = cos(t49);
t57 = t45 * t48;
t56 = pkin(1) * qJD(1);
t55 = pkin(2) * t58;
t54 = pkin(2) * t57;
t53 = t50 * t56;
t51 = cos(qJ(1));
t52 = t51 * t56;
t32 = t45 * rSges(3,1) - t44 * rSges(3,2);
t30 = t42 * rSges(4,1) - t41 * rSges(4,2);
t26 = -rSges(3,1) * t57 + rSges(3,2) * t58;
t20 = -rSges(4,1) * t59 + rSges(4,2) * t60;
t40 = pkin(2) * t45;
t24 = t30 + t40;
t22 = -t41 * rSges(5,2) + t64 * t42;
t31 = -t44 * rSges(3,1) - t45 * rSges(3,2);
t29 = -t41 * rSges(4,1) - t42 * rSges(4,2);
t14 = t22 + t40;
t8 = rSges(5,2) * t60 - t59 * t64;
t21 = -t42 * rSges(5,2) - t41 * t64;
t25 = t31 * t48;
t19 = t29 * t43;
t23 = t29 - t63;
t10 = t20 - t54;
t7 = t21 * t43;
t13 = t21 - t63;
t4 = t8 - t54;
t9 = t19 - t55;
t3 = t7 - t55;
t47 = t51 * pkin(1);
t28 = t32 + t47;
t27 = t31 - t62;
t18 = t26 - t52;
t17 = t25 - t53;
t16 = t24 + t47;
t15 = t23 - t62;
t12 = t14 + t47;
t11 = t13 - t62;
t6 = t10 - t52;
t5 = t9 - t53;
t2 = t4 - t52;
t1 = t3 - t53;
t33 = [0.2e1 * m(3) * (t28 * t17 + t27 * t18) + 0.2e1 * m(4) * (t15 * t6 + t16 * t5) + 0.2e1 * m(5) * (t12 * t1 + t11 * t2); m(3) * (t32 * t17 + t31 * t18 + t25 * t28 + t26 * t27) + m(4) * (t10 * t15 + t9 * t16 + t23 * t6 + t24 * t5) + m(5) * (t14 * t1 + t4 * t11 + t3 * t12 + t13 * t2); 0.2e1 * m(3) * (t32 * t25 + t31 * t26) + 0.2e1 * m(4) * (t23 * t10 + t24 * t9) + 0.2e1 * m(5) * (t13 * t4 + t14 * t3); m(4) * (t20 * t15 + t19 * t16 + t29 * t6 + t30 * t5) + m(5) * (t22 * t1 + t8 * t11 + t7 * t12 + t21 * t2); m(4) * (t29 * t10 + t19 * t24 + t20 * t23 + t30 * t9) + m(5) * (t8 * t13 + t7 * t14 + t21 * t4 + t22 * t3); 0.2e1 * m(4) * (t30 * t19 + t29 * t20) + 0.2e1 * m(5) * (t21 * t8 + t22 * t7); 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t33(1) t33(2) t33(4) t33(7); t33(2) t33(3) t33(5) t33(8); t33(4) t33(5) t33(6) t33(9); t33(7) t33(8) t33(9) t33(10);];
Mq  = res;
