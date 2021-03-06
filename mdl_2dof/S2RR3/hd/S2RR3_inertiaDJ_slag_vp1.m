% Calculate time derivative of joint inertia matrix for
% S2RR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,d1,d2]';
% m [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [2x2]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:14
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S2RR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR3_inertiaDJ_slag_vp1: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR3_inertiaDJ_slag_vp1: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S2RR3_inertiaDJ_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR3_inertiaDJ_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'S2RR3_inertiaDJ_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'S2RR3_inertiaDJ_slag_vp1: Icges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-06-19 09:14:25
% EndTime: 2020-06-19 09:14:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (83->13), mult. (100->23), div. (0->0), fcn. (48->4), ass. (0->17)
t14 = qJ(1) + qJ(2);
t11 = sin(t14);
t12 = cos(t14);
t8 = t12 * rSges(3,1) - t11 * rSges(3,2);
t19 = 2 * m(3);
t17 = pkin(1) * qJD(1);
t13 = qJD(1) + qJD(2);
t4 = t8 * t13;
t7 = -rSges(3,1) * t11 - rSges(3,2) * t12;
t3 = t7 * t13;
t16 = cos(qJ(1));
t15 = sin(qJ(1));
t6 = pkin(1) * t16 + t8;
t5 = -pkin(1) * t15 + t7;
t2 = -t16 * t17 - t4;
t1 = -t15 * t17 + t3;
t9 = [(t1 * t6 + t2 * t5) * t19; m(3) * (t1 * t8 + t2 * t7 + t3 * t6 - t4 * t5); (t3 * t8 - t4 * t7) * t19;];
%% Postprocessing: Reshape Output
% From vec2symmat_2_matlab.m
res = [t9(1), t9(2); t9(2), t9(3);];
Mq = res;
