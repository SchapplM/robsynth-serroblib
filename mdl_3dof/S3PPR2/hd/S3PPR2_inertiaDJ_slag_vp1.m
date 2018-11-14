% Calculate time derivative of joint inertia matrix for
% S3PPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
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
% Datum: 2018-11-14 10:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3PPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PPR2_inertiaDJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PPR2_inertiaDJ_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PPR2_inertiaDJ_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PPR2_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PPR2_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PPR2_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:11:05
% EndTime: 2018-11-14 10:11:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (7->4), mult. (34->12), div. (0->0), fcn. (12->2), ass. (0->7)
t3 = sin(qJ(3));
t4 = cos(qJ(3));
t6 = -t4 * rSges(4,1) + t3 * rSges(4,2);
t5 = t3 * rSges(4,1) + t4 * rSges(4,2);
t2 = t6 * qJD(3);
t1 = t5 * qJD(3);
t7 = [0; 0; 0; m(4) * t2; -m(4) * t1; 0.2e1 * m(4) * (t6 * t1 - t5 * t2);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t7(1) t7(2) t7(4); t7(2) t7(3) t7(5); t7(4) t7(5) t7(6);];
Mq  = res;
