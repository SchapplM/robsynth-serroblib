% Calculate time derivative of joint inertia matrix for
% S3PRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
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
% Datum: 2018-11-14 10:04
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S3PRP1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP1_inertiaDJ_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP1_inertiaDJ_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP1_inertiaDJ_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP1_inertiaDJ_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3PRP1_inertiaDJ_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3PRP1_inertiaDJ_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:04:14
% EndTime: 2018-11-14 10:04:15
% DurationCPUTime: 0.08s
% Computational Cost: add. (44->12), mult. (114->25), div. (0->0), fcn. (56->2), ass. (0->9)
t10 = rSges(4,1) + pkin(2);
t6 = sin(qJ(2));
t7 = cos(qJ(2));
t9 = rSges(4,3) + qJ(3);
t4 = t10 * t7 + t9 * t6;
t3 = -t10 * t6 + t9 * t7;
t2 = -qJD(2) * t4 + qJD(3) * t7;
t1 = qJD(2) * t3 + qJD(3) * t6;
t5 = [0; m(3) * (-t6 * rSges(3,1) - t7 * rSges(3,2)) * qJD(2) + m(4) * t1; 0.2e1 * m(4) * (t4 * t1 + t3 * t2); m(4) * qJD(2) * t6; m(4) * (-t7 * t1 + t6 * t2 + (t3 * t7 + t4 * t6) * qJD(2)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t5(1) t5(2) t5(4); t5(2) t5(3) t5(5); t5(4) t5(5) t5(6);];
Mq  = res;
