% Calculate time derivative of joint inertia matrix for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:08:16
% EndTime: 2018-11-14 14:08:16
% DurationCPUTime: 0.09s
% Computational Cost: add. (151->16), mult. (127->29), div. (0->0), fcn. (55->4), ass. (0->18)
t17 = pkin(6) + qJ(3);
t16 = qJ(4) + t17;
t12 = sin(t16);
t13 = cos(t16);
t9 = t13 * rSges(5,1) - t12 * rSges(5,2);
t22 = 2 * m(5);
t20 = pkin(3) * qJD(3);
t18 = qJD(3) + qJD(4);
t4 = t9 * t18;
t8 = -t12 * rSges(5,1) - t13 * rSges(5,2);
t3 = t8 * t18;
t15 = cos(t17);
t14 = sin(t17);
t6 = pkin(3) * t15 + t9;
t5 = -pkin(3) * t14 + t8;
t2 = -t15 * t20 - t4;
t1 = -t14 * t20 + t3;
t7 = [0; 0; 0; m(4) * (-t14 * rSges(4,1) - t15 * rSges(4,2)) * qJD(3) + m(5) * t1; 0; (t6 * t1 + t5 * t2) * t22; m(5) * t3; 0; m(5) * (t9 * t1 + t8 * t2 + t3 * t6 - t4 * t5); (t9 * t3 - t8 * t4) * t22;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t7(1) t7(2) t7(4) t7(7); t7(2) t7(3) t7(5) t7(8); t7(4) t7(5) t7(6) t7(9); t7(7) t7(8) t7(9) t7(10);];
Mq  = res;
