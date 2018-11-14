% Calculate time derivative of joint inertia matrix for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2018-11-14 13:39
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S4PPRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP1_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP1_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:39:04
% EndTime: 2018-11-14 13:39:04
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->25), mult. (316->50), div. (0->0), fcn. (278->4), ass. (0->17)
t15 = cos(pkin(5));
t18 = sin(pkin(5));
t20 = sin(qJ(3));
t21 = cos(qJ(3));
t10 = t15 * t20 - t18 * t21;
t22 = rSges(5,1) + pkin(3);
t19 = -rSges(5,3) - qJ(4);
t9 = -t15 * t21 - t18 * t20;
t8 = t9 * qJD(3);
t7 = t10 * qJD(3);
t6 = t8 * rSges(4,1) + t7 * rSges(4,2);
t5 = t7 * rSges(4,1) - t8 * rSges(4,2);
t4 = t19 * t10 + t22 * t9;
t3 = -t22 * t10 + t19 * t9;
t2 = -t9 * qJD(4) + t19 * t7 + t22 * t8;
t1 = -t10 * qJD(4) - t19 * t8 + t22 * t7;
t11 = [0; 0; 0; 0; m(4) * (-t5 * t15 + t6 * t18) + m(5) * (-t1 * t15 + t2 * t18); 0.2e1 * m(4) * ((-t10 * rSges(4,1) + t9 * rSges(4,2)) * t6 + (t9 * rSges(4,1) + t10 * rSges(4,2)) * t5) + 0.2e1 * m(5) * (t4 * t1 + t3 * t2); 0; m(5) * (t7 * t15 - t8 * t18); m(5) * (-t9 * t1 + t10 * t2 - t8 * t3 - t7 * t4); 0.2e1 * m(5) * (-t10 * t8 + t9 * t7);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t11(1) t11(2) t11(4) t11(7); t11(2) t11(3) t11(5) t11(8); t11(4) t11(5) t11(6) t11(9); t11(7) t11(8) t11(9) t11(10);];
Mq  = res;
