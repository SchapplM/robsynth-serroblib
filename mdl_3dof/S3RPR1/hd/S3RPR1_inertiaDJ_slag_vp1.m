% Calculate time derivative of joint inertia matrix for
% S3RPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
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

function Mq = S3RPR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPR1_inertiaDJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPR1_inertiaDJ_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RPR1_inertiaDJ_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPR1_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPR1_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPR1_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:14:23
% EndTime: 2018-11-14 10:14:23
% DurationCPUTime: 0.16s
% Computational Cost: add. (213->46), mult. (482->71), div. (0->0), fcn. (386->4), ass. (0->32)
t27 = sin(qJ(3));
t30 = cos(qJ(1));
t28 = sin(qJ(1));
t29 = cos(qJ(3));
t38 = t28 * t29;
t16 = -t30 * t27 + t38;
t40 = -pkin(1) - pkin(2);
t39 = -rSges(3,1) - pkin(1);
t34 = qJD(1) * t30;
t36 = qJ(2) * t34 + qJD(2) * t28;
t24 = t28 * qJ(2);
t35 = t30 * pkin(1) + t24;
t33 = t40 * t28;
t32 = t28 * t27 + t30 * t29;
t7 = (qJD(1) - qJD(3)) * t32;
t8 = -qJD(1) * t38 + t16 * qJD(3) + t27 * t34;
t4 = -t8 * rSges(4,1) - t7 * rSges(4,2);
t3 = t7 * rSges(4,1) - t8 * rSges(4,2);
t9 = t16 * rSges(4,1) - rSges(4,2) * t32;
t10 = -rSges(4,1) * t32 - t16 * rSges(4,2);
t31 = t30 * rSges(3,3) + t39 * t28;
t25 = t30 * qJ(2);
t23 = qJD(2) * t30;
t14 = t30 * rSges(3,1) + t28 * rSges(3,3) + t35;
t13 = t25 + t31;
t12 = t23 + (t39 * t30 + (-rSges(3,3) - qJ(2)) * t28) * qJD(1);
t11 = t31 * qJD(1) + t36;
t6 = t30 * pkin(2) - t10 + t35;
t5 = t25 + t33 - t9;
t2 = t23 + (t40 * t30 - t24) * qJD(1) - t3;
t1 = qJD(1) * t33 + t36 - t4;
t15 = [0.2e1 * m(3) * (t14 * t11 + t13 * t12) + 0.2e1 * m(4) * (t6 * t1 + t5 * t2); m(3) * (-t30 * t11 + t28 * t12 + (t13 * t30 + t14 * t28) * qJD(1)) + m(4) * (-t30 * t1 + t28 * t2 + (t28 * t6 + t30 * t5) * qJD(1)); 0; m(4) * (t10 * t1 + t9 * t2 + t3 * t5 + t4 * t6); m(4) * (t3 * t28 - t4 * t30 + (t10 * t28 + t30 * t9) * qJD(1)); 0.2e1 * m(4) * (t10 * t4 + t9 * t3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t15(1) t15(2) t15(4); t15(2) t15(3) t15(5); t15(4) t15(5) t15(6);];
Mq  = res;
