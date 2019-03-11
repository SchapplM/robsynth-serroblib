% Calculate time derivative of joint inertia matrix for
% S4PPRR1
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
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:38
% EndTime: 2019-03-08 18:15:38
% DurationCPUTime: 0.15s
% Computational Cost: add. (245->31), mult. (322->61), div. (0->0), fcn. (258->6), ass. (0->30)
t30 = qJ(3) + qJ(4);
t27 = sin(t30);
t28 = cos(t30);
t31 = sin(pkin(6));
t32 = cos(pkin(6));
t40 = t32 * t27 - t31 * t28;
t34 = cos(qJ(3));
t39 = t34 * pkin(3);
t33 = sin(qJ(3));
t37 = t31 * t33;
t35 = t32 * t33;
t29 = qJD(3) + qJD(4);
t13 = t40 * t29;
t17 = -t31 * t27 - t32 * t28;
t14 = t17 * t29;
t6 = t14 * rSges(5,1) + t13 * rSges(5,2);
t7 = -rSges(5,1) * t40 + t17 * rSges(5,2);
t5 = t13 * rSges(5,1) - t14 * rSges(5,2);
t8 = t17 * rSges(5,1) + rSges(5,2) * t40;
t22 = t31 * t34 - t35;
t21 = -t32 * t34 - t37;
t19 = t22 * qJD(3);
t20 = t21 * qJD(3);
t10 = t20 * rSges(4,1) - t19 * rSges(4,2);
t9 = -t19 * rSges(4,1) - t20 * rSges(4,2);
t4 = -pkin(3) * t37 - t39 * t32 + t8;
t3 = -pkin(3) * t35 + t39 * t31 + t7;
t2 = -pkin(3) * t19 + t5;
t1 = pkin(3) * t20 + t6;
t11 = [0; 0; 0; 0; m(4) * (t10 * t31 - t9 * t32) + m(5) * (t1 * t31 - t2 * t32); 0.2e1 * m(4) * ((t22 * rSges(4,1) + t21 * rSges(4,2)) * t10 + (t21 * rSges(4,1) - t22 * rSges(4,2)) * t9) + 0.2e1 * m(5) * (t3 * t1 + t4 * t2); 0; m(5) * (t6 * t31 - t5 * t32); m(5) * (t7 * t1 + t8 * t2 + t6 * t3 + t5 * t4); 0.2e1 * m(5) * (t8 * t5 + t7 * t6);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
