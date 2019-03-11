% Calculate time derivative of joint inertia matrix for
% S3PRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2019-03-08 18:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S3PRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR1_inertiaDJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR1_inertiaDJ_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR1_inertiaDJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR1_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRR1_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRR1_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:03:50
% EndTime: 2019-03-08 18:03:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (96->15), mult. (127->29), div. (0->0), fcn. (55->4), ass. (0->17)
t15 = qJ(2) + qJ(3);
t12 = sin(t15);
t13 = cos(t15);
t8 = t13 * rSges(4,1) - t12 * rSges(4,2);
t21 = 2 * m(4);
t19 = pkin(2) * qJD(2);
t14 = qJD(2) + qJD(3);
t4 = t8 * t14;
t7 = -t12 * rSges(4,1) - t13 * rSges(4,2);
t3 = t7 * t14;
t17 = cos(qJ(2));
t16 = sin(qJ(2));
t6 = t17 * pkin(2) + t8;
t5 = -t16 * pkin(2) + t7;
t2 = -t17 * t19 - t4;
t1 = -t16 * t19 + t3;
t9 = [0; m(3) * (-t16 * rSges(3,1) - t17 * rSges(3,2)) * qJD(2) + m(4) * t1; (t6 * t1 + t5 * t2) * t21; m(4) * t3; m(4) * (t8 * t1 + t7 * t2 + t3 * t6 - t4 * t5); (t8 * t3 - t7 * t4) * t21;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t9(1) t9(2) t9(4); t9(2) t9(3) t9(5); t9(4) t9(5) t9(6);];
Mq  = res;
