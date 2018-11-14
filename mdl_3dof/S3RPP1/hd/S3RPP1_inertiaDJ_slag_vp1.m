% Calculate time derivative of joint inertia matrix for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3RPP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_inertiaDJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_inertiaDJ_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_inertiaDJ_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPP1_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:33
% EndTime: 2018-11-14 10:13:33
% DurationCPUTime: 0.10s
% Computational Cost: add. (107->36), mult. (220->50), div. (0->0), fcn. (128->2), ass. (0->20)
t23 = rSges(3,2) - pkin(1);
t16 = cos(qJ(1));
t13 = t16 * qJ(2);
t15 = sin(qJ(1));
t22 = qJD(1) * t13 + qJD(2) * t15;
t21 = t16 * pkin(1) + t15 * qJ(2);
t20 = rSges(4,3) + qJ(3);
t19 = -pkin(1) - t20;
t18 = t16 * rSges(3,3) + t23 * t15;
t17 = t16 * rSges(4,2) + t19 * t15;
t11 = qJD(2) * t16;
t8 = -t16 * rSges(3,2) + t15 * rSges(3,3) + t21;
t7 = t13 + t18;
t6 = t15 * rSges(4,2) + t20 * t16 + t21;
t5 = t13 + t17;
t4 = t11 + (t23 * t16 + (-rSges(3,3) - qJ(2)) * t15) * qJD(1);
t3 = t18 * qJD(1) + t22;
t2 = -qJD(3) * t15 + t11 + ((-rSges(4,2) - qJ(2)) * t15 + t19 * t16) * qJD(1);
t1 = t17 * qJD(1) + qJD(3) * t16 + t22;
t9 = [0.2e1 * m(3) * (t8 * t3 + t7 * t4) + 0.2e1 * m(4) * (t6 * t1 + t5 * t2); m(3) * (t15 * t4 - t16 * t3 + (t15 * t8 + t16 * t7) * qJD(1)) + m(4) * (-t16 * t1 + t15 * t2 + (t15 * t6 + t16 * t5) * qJD(1)); 0; m(4) * (t15 * t1 + t16 * t2 + (-t15 * t5 + t16 * t6) * qJD(1)); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t9(1) t9(2) t9(4); t9(2) t9(3) t9(5); t9(4) t9(5) t9(6);];
Mq  = res;
