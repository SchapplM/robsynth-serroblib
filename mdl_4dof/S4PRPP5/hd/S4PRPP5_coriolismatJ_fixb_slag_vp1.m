% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:10
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4PRPP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP5_coriolismatJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:10:12
% EndTime: 2018-11-14 14:10:12
% DurationCPUTime: 0.07s
% Computational Cost: add. (120->19), mult. (221->30), div. (0->0), fcn. (147->2), ass. (0->11)
t14 = sin(qJ(2));
t15 = cos(qJ(2));
t16 = rSges(5,1) + pkin(2) + pkin(3);
t17 = rSges(4,1) + pkin(2);
t13 = t15 * qJ(3);
t5 = t15 * rSges(5,2) - t16 * t14 + t13;
t8 = t15 * rSges(4,3) - t17 * t14 + t13;
t1 = m(4) * ((t17 * t15 + (rSges(4,3) + qJ(3)) * t14) * t14 + t8 * t15) + m(5) * (((rSges(5,2) + qJ(3)) * t14 + t16 * t15) * t14 + t5 * t15);
t18 = t1 * qJD(2);
t12 = (m(4) + m(5)) * t14;
t2 = [0, t12 * qJD(3) + 0.2e1 * (m(3) * (-t14 * rSges(3,1) - t15 * rSges(3,2)) / 0.2e1 + m(4) * t8 / 0.2e1 + m(5) * t5 / 0.2e1) * qJD(2), t12 * qJD(2), 0; 0, t1 * qJD(3), t18, 0; 0, -t18, 0, 0; 0, 0, 0, 0;];
Cq  = t2;
