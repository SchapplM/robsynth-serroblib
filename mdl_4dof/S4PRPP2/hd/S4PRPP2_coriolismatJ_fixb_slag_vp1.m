% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRPP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP2_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:18:49
% EndTime: 2019-03-08 18:18:49
% DurationCPUTime: 0.06s
% Computational Cost: add. (153->15), mult. (157->26), div. (0->0), fcn. (100->4), ass. (0->14)
t11 = qJ(2) + pkin(5);
t10 = cos(t11);
t13 = cos(qJ(2));
t14 = rSges(5,3) + qJ(4);
t17 = rSges(5,1) + pkin(3);
t12 = sin(qJ(2));
t18 = t12 * pkin(2);
t9 = sin(t11);
t3 = t14 * t10 - t17 * t9 - t18;
t1 = t3 * t10 + (t13 * pkin(2) + t17 * t10 + t14 * t9) * t9;
t16 = m(5) * qJD(2);
t19 = t1 * t16;
t15 = m(5) * qJD(4);
t2 = [0, t9 * t15 + 0.2e1 * (m(3) * (-t12 * rSges(3,1) - t13 * rSges(3,2)) / 0.2e1 + m(4) * (-t9 * rSges(4,1) - t10 * rSges(4,2) - t18) / 0.2e1 + m(5) * t3 / 0.2e1) * qJD(2), 0, t9 * t16; 0, t1 * t15, 0, t19; 0, 0, 0, 0; 0, -t19, 0, 0;];
Cq  = t2;
