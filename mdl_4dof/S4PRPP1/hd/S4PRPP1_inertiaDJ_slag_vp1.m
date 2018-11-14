% Calculate time derivative of joint inertia matrix for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2018-11-14 13:41
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:59
% EndTime: 2018-11-14 13:40:59
% DurationCPUTime: 0.10s
% Computational Cost: add. (235->37), mult. (220->50), div. (0->0), fcn. (128->2), ass. (0->21)
t24 = rSges(4,2) - pkin(2);
t17 = pkin(5) + qJ(2);
t16 = cos(t17);
t13 = t16 * qJ(3);
t15 = sin(t17);
t23 = qJD(2) * t13 + qJD(3) * t15;
t22 = t16 * pkin(2) + t15 * qJ(3);
t21 = rSges(5,3) + qJ(4);
t20 = -pkin(2) - t21;
t19 = t16 * rSges(4,3) + t24 * t15;
t18 = t16 * rSges(5,2) + t20 * t15;
t11 = qJD(3) * t16;
t8 = -t16 * rSges(4,2) + t15 * rSges(4,3) + t22;
t7 = t13 + t19;
t6 = t15 * rSges(5,2) + t21 * t16 + t22;
t5 = t13 + t18;
t4 = t11 + (t24 * t16 + (-rSges(4,3) - qJ(3)) * t15) * qJD(2);
t3 = t19 * qJD(2) + t23;
t2 = -qJD(4) * t15 + t11 + ((-rSges(5,2) - qJ(3)) * t15 + t20 * t16) * qJD(2);
t1 = t18 * qJD(2) + qJD(4) * t16 + t23;
t9 = [0; 0; 0.2e1 * m(4) * (t8 * t3 + t7 * t4) + 0.2e1 * m(5) * (t6 * t1 + t5 * t2); 0; m(4) * (t15 * t4 - t16 * t3 + (t15 * t8 + t16 * t7) * qJD(2)) + m(5) * (-t16 * t1 + t15 * t2 + (t15 * t6 + t16 * t5) * qJD(2)); 0; 0; m(5) * (t15 * t1 + t16 * t2 + (-t15 * t5 + t16 * t6) * qJD(2)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t9(1) t9(2) t9(4) t9(7); t9(2) t9(3) t9(5) t9(8); t9(4) t9(5) t9(6) t9(9); t9(7) t9(8) t9(9) t9(10);];
Mq  = res;
