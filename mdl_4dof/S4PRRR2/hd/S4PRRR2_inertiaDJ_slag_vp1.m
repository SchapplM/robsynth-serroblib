% Calculate time derivative of joint inertia matrix for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_inertiaDJ_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.11s
% Computational Cost: add. (422->44), mult. (310->69), div. (0->0), fcn. (156->6), ass. (0->40)
t39 = sin(qJ(2));
t48 = t39 * pkin(1);
t40 = cos(qJ(2));
t47 = t40 * pkin(1);
t38 = qJ(2) + qJ(3);
t36 = qJ(4) + t38;
t30 = sin(t36);
t37 = qJD(2) + qJD(3);
t33 = qJD(4) + t37;
t46 = t30 * t33;
t31 = cos(t36);
t45 = t31 * t33;
t34 = sin(t38);
t44 = t34 * t37;
t35 = cos(t38);
t43 = t35 * t37;
t9 = rSges(5,1) * t46 + rSges(5,2) * t45;
t13 = rSges(4,1) * t44 + rSges(4,2) * t43;
t42 = pkin(1) * qJD(2);
t3 = pkin(2) * t44 + t9;
t41 = t40 * t42;
t20 = -t35 * rSges(4,1) + t34 * rSges(4,2);
t18 = -t31 * rSges(5,1) + t30 * rSges(5,2);
t14 = -rSges(4,1) * t43 + rSges(4,2) * t44;
t10 = -rSges(5,1) * t45 + rSges(5,2) * t46;
t19 = -t34 * rSges(4,1) - t35 * rSges(4,2);
t17 = -t30 * rSges(5,1) - t31 * rSges(5,2);
t12 = -pkin(2) * t35 + t18;
t11 = -pkin(2) * t34 + t17;
t4 = -pkin(2) * t43 + t10;
t32 = t39 * t42;
t16 = t20 - t47;
t15 = t19 - t48;
t8 = t14 - t41;
t7 = t32 + t13;
t6 = t12 - t47;
t5 = t11 - t48;
t2 = t4 - t41;
t1 = t32 + t3;
t21 = [0; 0; 0.2e1 * m(4) * (t15 * t8 + t16 * t7) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); 0; m(4) * (t13 * t16 + t14 * t15 + t19 * t8 + t20 * t7) + m(5) * (t12 * t1 + t11 * t2 + t3 * t6 + t4 * t5); 0.2e1 * m(4) * (t20 * t13 + t19 * t14) + 0.2e1 * m(5) * (t11 * t4 + t12 * t3); 0; m(5) * (t18 * t1 + t10 * t5 + t17 * t2 + t9 * t6); m(5) * (t10 * t11 + t9 * t12 + t17 * t4 + t18 * t3); 0.2e1 * m(5) * (t17 * t10 + t18 * t9);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t21(1), t21(2), t21(4), t21(7); t21(2), t21(3), t21(5), t21(8); t21(4), t21(5), t21(6), t21(9); t21(7), t21(8), t21(9), t21(10);];
Mq  = res;
