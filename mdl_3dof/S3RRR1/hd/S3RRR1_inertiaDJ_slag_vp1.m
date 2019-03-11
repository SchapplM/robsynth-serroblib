% Calculate time derivative of joint inertia matrix for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [3x3]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3RRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_inertiaDJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_inertiaDJ_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRR1_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRR1_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:07
% EndTime: 2019-03-08 18:08:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (422->41), mult. (310->64), div. (0->0), fcn. (156->6), ass. (0->38)
t34 = qJ(1) + qJ(2);
t31 = qJ(3) + t34;
t26 = sin(t31);
t27 = cos(t31);
t18 = t27 * rSges(4,1) - t26 * rSges(4,2);
t35 = sin(qJ(1));
t43 = t35 * pkin(1);
t29 = sin(t34);
t33 = qJD(1) + qJD(2);
t41 = t29 * t33;
t30 = cos(t34);
t40 = t30 * t33;
t39 = pkin(1) * qJD(1);
t38 = t35 * t39;
t36 = cos(qJ(1));
t37 = t36 * t39;
t20 = t30 * rSges(3,1) - t29 * rSges(3,2);
t14 = -rSges(3,1) * t40 + rSges(3,2) * t41;
t28 = qJD(3) + t33;
t10 = t18 * t28;
t12 = pkin(2) * t30 + t18;
t19 = -t29 * rSges(3,1) - t30 * rSges(3,2);
t17 = -t26 * rSges(4,1) - t27 * rSges(4,2);
t13 = t19 * t33;
t9 = t17 * t28;
t11 = -pkin(2) * t29 + t17;
t4 = -pkin(2) * t40 - t10;
t3 = -pkin(2) * t41 + t9;
t32 = t36 * pkin(1);
t16 = t20 + t32;
t15 = t19 - t43;
t8 = t14 - t37;
t7 = t13 - t38;
t6 = t12 + t32;
t5 = t11 - t43;
t2 = t4 - t37;
t1 = t3 - t38;
t21 = [0.2e1 * m(3) * (t15 * t8 + t16 * t7) + 0.2e1 * m(4) * (t6 * t1 + t5 * t2); m(3) * (t13 * t16 + t14 * t15 + t19 * t8 + t20 * t7) + m(4) * (t12 * t1 + t11 * t2 + t3 * t6 + t4 * t5); 0.2e1 * m(3) * (t20 * t13 + t19 * t14) + 0.2e1 * m(4) * (t11 * t4 + t12 * t3); m(4) * (t18 * t1 - t10 * t5 + t17 * t2 + t9 * t6); m(4) * (-t10 * t11 + t9 * t12 + t17 * t4 + t18 * t3); 0.2e1 * m(4) * (-t17 * t10 + t18 * t9);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t21(1) t21(2) t21(4); t21(2) t21(3) t21(5); t21(4) t21(5) t21(6);];
Mq  = res;
