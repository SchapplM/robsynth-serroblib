% Calculate time derivative of joint inertia matrix for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2018-11-14 14:02
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:02:21
% EndTime: 2018-11-14 14:02:21
% DurationCPUTime: 0.12s
% Computational Cost: add. (173->21), mult. (177->33), div. (0->0), fcn. (79->6), ass. (0->22)
t20 = qJ(2) + pkin(6);
t17 = qJ(4) + t20;
t13 = sin(t17);
t14 = cos(t17);
t9 = t14 * rSges(5,1) - t13 * rSges(5,2);
t30 = 2 * m(5);
t16 = cos(t20);
t22 = cos(qJ(2));
t29 = -t22 * pkin(2) - pkin(3) * t16;
t21 = sin(qJ(2));
t27 = t21 * pkin(2);
t19 = qJD(2) + qJD(4);
t6 = t9 * t19;
t15 = sin(t20);
t25 = -pkin(3) * t15 - t27;
t8 = -t13 * rSges(5,1) - t14 * rSges(5,2);
t5 = t8 * t19;
t4 = t9 - t29;
t3 = t25 + t8;
t2 = t29 * qJD(2) - t6;
t1 = t25 * qJD(2) + t5;
t7 = [0; m(5) * t1 + (m(3) * (-t21 * rSges(3,1) - t22 * rSges(3,2)) + m(4) * (-t15 * rSges(4,1) - t16 * rSges(4,2) - t27)) * qJD(2); (t4 * t1 + t3 * t2) * t30; 0; 0; 0; m(5) * t5; m(5) * (t9 * t1 + t8 * t2 - t6 * t3 + t5 * t4); 0; (t9 * t5 - t8 * t6) * t30;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
