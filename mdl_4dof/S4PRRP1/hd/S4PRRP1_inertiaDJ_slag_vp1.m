% Calculate time derivative of joint inertia matrix for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:49
% EndTime: 2019-03-08 18:22:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (537->41), mult. (320->58), div. (0->0), fcn. (176->4), ass. (0->34)
t36 = qJD(2) + qJD(3);
t35 = pkin(6) + qJ(2);
t34 = qJ(3) + t35;
t30 = sin(t34);
t31 = cos(t34);
t47 = rSges(5,1) + pkin(3);
t8 = t47 * t31 + (rSges(5,3) + qJ(4)) * t30;
t49 = t36 * t8;
t46 = t31 * qJ(4) - t30 * t47;
t32 = sin(t35);
t44 = pkin(2) * t32;
t42 = rSges(4,2) * t30;
t41 = t31 * t36;
t40 = pkin(2) * qJD(2);
t39 = t32 * t40;
t33 = cos(t35);
t38 = t33 * t40;
t16 = t31 * rSges(4,1) - t42;
t12 = -rSges(4,1) * t41 + t36 * t42;
t15 = -rSges(4,1) * t30 - rSges(4,2) * t31;
t7 = t31 * rSges(5,3) + t46;
t11 = t15 * t36;
t3 = rSges(5,3) * t41 + qJD(4) * t30 + t46 * t36;
t4 = qJD(4) * t31 - t49;
t29 = pkin(2) * t33;
t14 = t16 + t29;
t13 = t15 - t44;
t10 = t12 - t38;
t9 = t11 - t39;
t6 = t29 + t8;
t5 = t7 - t44;
t2 = t4 - t38;
t1 = t3 - t39;
t17 = [0; 0; 0.2e1 * m(4) * (t10 * t13 + t14 * t9) + 0.2e1 * m(5) * (t1 * t6 + t2 * t5); 0; m(4) * (t10 * t15 + t11 * t14 + t12 * t13 + t16 * t9) + m(5) * (t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5); 0.2e1 * m(4) * (t11 * t16 + t12 * t15) + 0.2e1 * m(5) * (t3 * t8 + t4 * t7); 0; m(5) * ((t36 * t5 - t1) * t31 + (t36 * t6 + t2) * t30); m(5) * ((t36 * t7 - t3) * t31 + (t4 + t49) * t30); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1) t17(2) t17(4) t17(7); t17(2) t17(3) t17(5) t17(8); t17(4) t17(5) t17(6) t17(9); t17(7) t17(8) t17(9) t17(10);];
Mq  = res;
