% Calculate time derivative of joint inertia matrix for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:23
% EndTime: 2018-11-14 13:44:23
% DurationCPUTime: 0.14s
% Computational Cost: add. (578->42), mult. (310->64), div. (0->0), fcn. (156->6), ass. (0->39)
t36 = pkin(7) + qJ(2);
t35 = qJ(3) + t36;
t31 = qJ(4) + t35;
t26 = sin(t31);
t27 = cos(t31);
t18 = t27 * rSges(5,1) - t26 * rSges(5,2);
t32 = sin(t36);
t44 = pkin(2) * t32;
t29 = sin(t35);
t37 = qJD(2) + qJD(3);
t42 = t29 * t37;
t30 = cos(t35);
t41 = t30 * t37;
t40 = pkin(2) * qJD(2);
t39 = t32 * t40;
t33 = cos(t36);
t38 = t33 * t40;
t20 = t30 * rSges(4,1) - t29 * rSges(4,2);
t14 = -rSges(4,1) * t41 + rSges(4,2) * t42;
t34 = qJD(4) + t37;
t10 = t18 * t34;
t12 = pkin(3) * t30 + t18;
t19 = -t29 * rSges(4,1) - t30 * rSges(4,2);
t17 = -t26 * rSges(5,1) - t27 * rSges(5,2);
t13 = t19 * t37;
t9 = t17 * t34;
t11 = -pkin(3) * t29 + t17;
t4 = -pkin(3) * t41 - t10;
t3 = -pkin(3) * t42 + t9;
t28 = pkin(2) * t33;
t16 = t20 + t28;
t15 = t19 - t44;
t8 = t14 - t38;
t7 = t13 - t39;
t6 = t12 + t28;
t5 = t11 - t44;
t2 = t4 - t38;
t1 = t3 - t39;
t21 = [0; 0; 0.2e1 * m(4) * (t15 * t8 + t16 * t7) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); 0; m(4) * (t13 * t16 + t14 * t15 + t19 * t8 + t20 * t7) + m(5) * (t12 * t1 + t11 * t2 + t3 * t6 + t4 * t5); 0.2e1 * m(4) * (t20 * t13 + t19 * t14) + 0.2e1 * m(5) * (t11 * t4 + t12 * t3); 0; m(5) * (t18 * t1 - t10 * t5 + t17 * t2 + t9 * t6); m(5) * (-t10 * t11 + t9 * t12 + t17 * t4 + t18 * t3); 0.2e1 * m(5) * (-t17 * t10 + t18 * t9);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t21(1) t21(2) t21(4) t21(7); t21(2) t21(3) t21(5) t21(8); t21(4) t21(5) t21(6) t21(9); t21(7) t21(8) t21(9) t21(10);];
Mq  = res;
