% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:14:29
% EndTime: 2019-03-08 18:14:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (63->24), mult. (172->39), div. (0->0), fcn. (70->2), ass. (0->19)
t15 = cos(qJ(3));
t14 = sin(qJ(3));
t18 = -t15 * rSges(5,1) + rSges(5,2) * t14;
t19 = pkin(3) * t15 - t18;
t16 = qJD(3) ^ 2;
t25 = t16 * pkin(3);
t23 = rSges(5,2) * t15;
t22 = m(4) * qJD(3);
t10 = rSges(4,1) * t14 + rSges(4,2) * t15;
t21 = qJD(3) * t10;
t12 = rSges(4,1) * t15 - rSges(4,2) * t14;
t8 = qJD(3) * t12;
t9 = rSges(5,1) * t14 + t23;
t17 = -t23 + (-rSges(5,1) - pkin(3)) * t14;
t6 = qJD(1) - t21;
t5 = qJD(2) + t8;
t2 = -t15 * t25 + t18 * t16;
t1 = -t14 * t25 - t9 * t16;
t3 = [m(5) * t2 - t8 * t22; m(5) * t1 - t21 * t22; -(-t10 * t5 - t12 * t6) * t22 + m(4) * (-t5 * t21 - t6 * t8 + (t10 * t8 - t12 * t21) * qJD(3)) + (t1 * t19 + t2 * t17 + (pkin(3) * t14 + t17 + t9) * (t19 * qJD(3) + qJD(2)) * qJD(3)) * m(5); 0;];
tauc  = t3(:);
