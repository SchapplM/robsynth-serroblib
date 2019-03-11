% Calculate time derivative of joint inertia matrix for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:37
% EndTime: 2019-03-08 18:29:37
% DurationCPUTime: 0.16s
% Computational Cost: add. (567->43), mult. (374->59), div. (0->0), fcn. (208->6), ass. (0->31)
t37 = qJ(1) + pkin(6);
t34 = qJ(3) + t37;
t30 = sin(t34);
t31 = cos(t34);
t53 = rSges(5,1) + pkin(3);
t10 = t53 * t31 + (rSges(5,3) + qJ(4)) * t30;
t36 = qJD(1) + qJD(3);
t55 = t10 * t36;
t52 = t31 * qJ(4) - t30 * t53;
t46 = pkin(2) * cos(t37) + cos(qJ(1)) * pkin(1);
t48 = t30 * rSges(4,2);
t47 = t31 * t36;
t16 = t31 * rSges(4,1) - t48;
t14 = -rSges(4,1) * t47 + t36 * t48;
t44 = -pkin(2) * sin(t37) - sin(qJ(1)) * pkin(1);
t15 = -t30 * rSges(4,1) - t31 * rSges(4,2);
t9 = t31 * rSges(5,3) + t52;
t13 = t15 * t36;
t42 = t44 * qJD(1);
t41 = t46 * qJD(1);
t3 = rSges(5,3) * t47 + qJD(4) * t30 + t52 * t36;
t4 = qJD(4) * t31 - t55;
t12 = t16 + t46;
t11 = t15 + t44;
t8 = t14 - t41;
t7 = t13 + t42;
t6 = t10 + t46;
t5 = t44 + t9;
t2 = -t41 + t4;
t1 = t42 + t3;
t17 = [0.2e1 * m(4) * (t11 * t8 + t12 * t7) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); 0; 0; m(4) * (t14 * t11 + t13 * t12 + t15 * t8 + t16 * t7) + m(5) * (t10 * t1 + t9 * t2 + t3 * t6 + t4 * t5); 0; 0.2e1 * m(4) * (t16 * t13 + t15 * t14) + 0.2e1 * m(5) * (t10 * t3 + t9 * t4); m(5) * ((t36 * t5 - t1) * t31 + (t36 * t6 + t2) * t30); 0; m(5) * ((t36 * t9 - t3) * t31 + (t4 + t55) * t30); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1) t17(2) t17(4) t17(7); t17(2) t17(3) t17(5) t17(8); t17(4) t17(5) t17(6) t17(9); t17(7) t17(8) t17(9) t17(10);];
Mq  = res;
