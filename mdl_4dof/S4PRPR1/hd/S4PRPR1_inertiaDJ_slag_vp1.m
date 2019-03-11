% Calculate time derivative of joint inertia matrix for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:52
% EndTime: 2019-03-08 18:20:52
% DurationCPUTime: 0.15s
% Computational Cost: add. (455->47), mult. (482->71), div. (0->0), fcn. (386->4), ass. (0->33)
t29 = pkin(6) + qJ(2);
t28 = cos(t29);
t30 = sin(qJ(4));
t27 = sin(t29);
t31 = cos(qJ(4));
t39 = t27 * t31;
t16 = -t28 * t30 + t39;
t41 = -pkin(2) - pkin(3);
t40 = -rSges(4,1) - pkin(2);
t35 = qJD(2) * t28;
t37 = qJ(3) * t35 + qJD(3) * t27;
t24 = t27 * qJ(3);
t36 = t28 * pkin(2) + t24;
t34 = t41 * t27;
t10 = -qJD(2) * t39 + t16 * qJD(4) + t30 * t35;
t33 = t27 * t30 + t28 * t31;
t9 = (qJD(2) - qJD(4)) * t33;
t3 = t9 * rSges(5,1) - t10 * rSges(5,2);
t4 = -t10 * rSges(5,1) - t9 * rSges(5,2);
t11 = t16 * rSges(5,1) - rSges(5,2) * t33;
t12 = -rSges(5,1) * t33 - t16 * rSges(5,2);
t32 = t28 * rSges(4,3) + t40 * t27;
t25 = t28 * qJ(3);
t23 = qJD(3) * t28;
t14 = t28 * rSges(4,1) + t27 * rSges(4,3) + t36;
t13 = t25 + t32;
t8 = t23 + (t40 * t28 + (-rSges(4,3) - qJ(3)) * t27) * qJD(2);
t7 = t32 * qJD(2) + t37;
t6 = t28 * pkin(3) - t12 + t36;
t5 = t25 + t34 - t11;
t2 = t23 + (t41 * t28 - t24) * qJD(2) - t3;
t1 = qJD(2) * t34 + t37 - t4;
t15 = [0; 0; 0.2e1 * m(4) * (t13 * t8 + t14 * t7) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); 0; m(4) * (t27 * t8 - t28 * t7 + (t13 * t28 + t14 * t27) * qJD(2)) + m(5) * (-t28 * t1 + t27 * t2 + (t27 * t6 + t28 * t5) * qJD(2)); 0; 0; m(5) * (t12 * t1 + t11 * t2 + t3 * t5 + t4 * t6); m(5) * (t3 * t27 - t4 * t28 + (t11 * t28 + t12 * t27) * qJD(2)); 0.2e1 * m(5) * (t11 * t3 + t12 * t4);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t15(1) t15(2) t15(4) t15(7); t15(2) t15(3) t15(5) t15(8); t15(4) t15(5) t15(6) t15(9); t15(7) t15(8) t15(9) t15(10);];
Mq  = res;
