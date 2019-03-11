% Calculate time derivative of joint inertia matrix for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:31:39
% EndTime: 2019-03-08 18:31:40
% DurationCPUTime: 0.16s
% Computational Cost: add. (608->44), mult. (364->65), div. (0->0), fcn. (188->8), ass. (0->36)
t38 = qJ(1) + pkin(7);
t35 = qJ(3) + t38;
t31 = qJ(4) + t35;
t26 = sin(t31);
t27 = cos(t31);
t18 = t27 * rSges(5,1) - t26 * rSges(5,2);
t46 = pkin(2) * cos(t38) + cos(qJ(1)) * pkin(1);
t29 = sin(t35);
t37 = qJD(1) + qJD(3);
t48 = t29 * t37;
t30 = cos(t35);
t47 = t30 * t37;
t20 = t30 * rSges(4,1) - t29 * rSges(4,2);
t16 = -rSges(4,1) * t47 + rSges(4,2) * t48;
t34 = qJD(4) + t37;
t10 = t18 * t34;
t14 = pkin(3) * t30 + t18;
t45 = -pkin(2) * sin(t38) - sin(qJ(1)) * pkin(1);
t19 = -t29 * rSges(4,1) - t30 * rSges(4,2);
t17 = -t26 * rSges(5,1) - t27 * rSges(5,2);
t15 = t19 * t37;
t9 = t17 * t34;
t43 = t45 * qJD(1);
t42 = t46 * qJD(1);
t13 = -pkin(3) * t29 + t17;
t4 = -pkin(3) * t47 - t10;
t3 = -pkin(3) * t48 + t9;
t12 = t20 + t46;
t11 = t19 + t45;
t8 = t16 - t42;
t7 = t15 + t43;
t6 = t14 + t46;
t5 = t13 + t45;
t2 = -t42 + t4;
t1 = t43 + t3;
t21 = [0.2e1 * m(4) * (t11 * t8 + t12 * t7) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); 0; 0; m(4) * (t16 * t11 + t15 * t12 + t19 * t8 + t20 * t7) + m(5) * (t14 * t1 + t13 * t2 + t3 * t6 + t4 * t5); 0; 0.2e1 * m(4) * (t20 * t15 + t19 * t16) + 0.2e1 * m(5) * (t13 * t4 + t14 * t3); m(5) * (t18 * t1 - t10 * t5 + t17 * t2 + t9 * t6); 0; m(5) * (-t10 * t13 + t9 * t14 + t17 * t4 + t18 * t3); 0.2e1 * m(5) * (-t17 * t10 + t18 * t9);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t21(1) t21(2) t21(4) t21(7); t21(2) t21(3) t21(5) t21(8); t21(4) t21(5) t21(6) t21(9); t21(7) t21(8) t21(9) t21(10);];
Mq  = res;
