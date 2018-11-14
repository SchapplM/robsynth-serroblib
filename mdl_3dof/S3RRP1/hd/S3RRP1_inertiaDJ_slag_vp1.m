% Calculate time derivative of joint inertia matrix for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
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
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3RRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_inertiaDJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_inertiaDJ_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_inertiaDJ_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RRP1_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RRP1_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:07
% EndTime: 2018-11-14 10:15:08
% DurationCPUTime: 0.18s
% Computational Cost: add. (361->40), mult. (320->58), div. (0->0), fcn. (176->4), ass. (0->33)
t32 = qJD(1) + qJD(2);
t33 = qJ(1) + qJ(2);
t29 = sin(t33);
t30 = cos(t33);
t46 = rSges(4,1) + pkin(2);
t8 = t46 * t30 + (rSges(4,3) + qJ(3)) * t29;
t48 = t32 * t8;
t45 = t30 * qJ(3) - t29 * t46;
t34 = sin(qJ(1));
t43 = t34 * pkin(1);
t41 = t29 * rSges(3,2);
t40 = t30 * t32;
t39 = pkin(1) * qJD(1);
t38 = t34 * t39;
t35 = cos(qJ(1));
t37 = t35 * t39;
t16 = t30 * rSges(3,1) - t41;
t12 = -rSges(3,1) * t40 + t32 * t41;
t15 = -t29 * rSges(3,1) - t30 * rSges(3,2);
t7 = t30 * rSges(4,3) + t45;
t11 = t15 * t32;
t3 = rSges(4,3) * t40 + qJD(3) * t29 + t32 * t45;
t4 = qJD(3) * t30 - t48;
t31 = t35 * pkin(1);
t14 = t16 + t31;
t13 = t15 - t43;
t10 = t12 - t37;
t9 = t11 - t38;
t6 = t31 + t8;
t5 = t7 - t43;
t2 = t4 - t37;
t1 = t3 - t38;
t17 = [0.2e1 * m(3) * (t13 * t10 + t14 * t9) + 0.2e1 * m(4) * (t6 * t1 + t5 * t2); m(3) * (t15 * t10 + t11 * t14 + t12 * t13 + t16 * t9) + m(4) * (t8 * t1 + t7 * t2 + t3 * t6 + t4 * t5); 0.2e1 * m(3) * (t16 * t11 + t15 * t12) + 0.2e1 * m(4) * (t8 * t3 + t7 * t4); m(4) * ((t32 * t5 - t1) * t30 + (t32 * t6 + t2) * t29); m(4) * ((t32 * t7 - t3) * t30 + (t4 + t48) * t29); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t17(1) t17(2) t17(4); t17(2) t17(3) t17(5); t17(4) t17(5) t17(6);];
Mq  = res;
