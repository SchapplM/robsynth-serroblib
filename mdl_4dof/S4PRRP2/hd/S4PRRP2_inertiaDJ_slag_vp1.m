% Calculate time derivative of joint inertia matrix for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP2_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP2_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:03:19
% EndTime: 2018-11-14 14:03:19
% DurationCPUTime: 0.12s
% Computational Cost: add. (235->34), mult. (258->55), div. (0->0), fcn. (118->4), ass. (0->31)
t38 = rSges(5,1) + pkin(3);
t28 = sin(qJ(2));
t37 = t28 * pkin(2);
t27 = qJ(2) + qJ(3);
t23 = sin(t27);
t26 = qJD(2) + qJD(3);
t35 = t23 * t26;
t24 = cos(t27);
t34 = t24 * t26;
t33 = pkin(2) * qJD(2);
t32 = t28 * t33;
t29 = cos(qJ(2));
t31 = t29 * t33;
t16 = t24 * rSges(4,1) - t23 * rSges(4,2);
t10 = -rSges(4,1) * t34 + rSges(4,2) * t35;
t12 = -t23 * rSges(5,2) + t38 * t24;
t15 = -t23 * rSges(4,1) - t24 * rSges(4,2);
t4 = rSges(5,2) * t35 - t34 * t38;
t11 = -t24 * rSges(5,2) - t23 * t38;
t9 = t15 * t26;
t3 = t11 * t26;
t25 = t29 * pkin(2);
t14 = t16 + t25;
t13 = t15 - t37;
t8 = t12 + t25;
t7 = t11 - t37;
t6 = t10 - t31;
t5 = t9 - t32;
t2 = t4 - t31;
t1 = t3 - t32;
t17 = [0; m(3) * (-t28 * rSges(3,1) - t29 * rSges(3,2)) * qJD(2) + m(4) * t5 + m(5) * t1; 0.2e1 * m(4) * (t13 * t6 + t14 * t5) + 0.2e1 * m(5) * (t8 * t1 + t7 * t2); m(4) * t9 + m(5) * t3; m(4) * (t10 * t13 + t9 * t14 + t15 * t6 + t16 * t5) + m(5) * (t12 * t1 + t11 * t2 + t3 * t8 + t4 * t7); 0.2e1 * m(4) * (t15 * t10 + t16 * t9) + 0.2e1 * m(5) * (t11 * t4 + t12 * t3); 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t17(1) t17(2) t17(4) t17(7); t17(2) t17(3) t17(5) t17(8); t17(4) t17(5) t17(6) t17(9); t17(7) t17(8) t17(9) t17(10);];
Mq  = res;
