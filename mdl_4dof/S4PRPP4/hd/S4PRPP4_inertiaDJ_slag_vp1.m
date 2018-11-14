% Calculate time derivative of joint inertia matrix for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2018-11-14 14:09
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRPP4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:23
% EndTime: 2018-11-14 14:09:23
% DurationCPUTime: 0.08s
% Computational Cost: add. (122->18), mult. (164->30), div. (0->0), fcn. (80->4), ass. (0->13)
t12 = cos(qJ(2));
t15 = rSges(5,3) + qJ(4);
t16 = rSges(5,1) + pkin(3);
t10 = qJ(2) + pkin(5);
t7 = sin(t10);
t8 = cos(t10);
t4 = t12 * pkin(2) + t15 * t7 + t16 * t8;
t11 = sin(qJ(2));
t17 = t11 * pkin(2);
t3 = t15 * t8 - t16 * t7 - t17;
t2 = -qJD(2) * t4 + qJD(4) * t8;
t1 = qJD(2) * t3 + qJD(4) * t7;
t5 = [0; m(5) * t1 + (m(3) * (-t11 * rSges(3,1) - t12 * rSges(3,2)) + m(4) * (-t7 * rSges(4,1) - t8 * rSges(4,2) - t17)) * qJD(2); 0.2e1 * m(5) * (t4 * t1 + t3 * t2); 0; 0; 0; m(5) * qJD(2) * t7; m(5) * (-t8 * t1 + t7 * t2 + (t3 * t8 + t4 * t7) * qJD(2)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t5(1) t5(2) t5(4) t5(7); t5(2) t5(3) t5(5) t5(8); t5(4) t5(5) t5(6) t5(9); t5(7) t5(8) t5(9) t5(10);];
Mq  = res;
