% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:14
% EndTime: 2019-03-08 18:25:15
% DurationCPUTime: 0.27s
% Computational Cost: add. (762->46), mult. (479->61), div. (0->0), fcn. (212->6), ass. (0->41)
t37 = pkin(7) + qJ(2);
t33 = cos(t37);
t50 = pkin(2) * qJD(2);
t46 = t33 * t50;
t35 = qJ(3) + t37;
t30 = cos(t35);
t38 = qJD(2) + qJD(3);
t53 = t30 * t38;
t31 = qJ(4) + t35;
t26 = sin(t31);
t27 = cos(t31);
t15 = t27 * rSges(5,1) - t26 * rSges(5,2);
t34 = qJD(4) + t38;
t62 = t34 * t15;
t4 = pkin(3) * t53 + t46 + t62;
t32 = sin(t37);
t47 = t32 * t50;
t29 = sin(t35);
t17 = t29 * rSges(4,1) + t30 * rSges(4,2);
t51 = t38 * t17;
t7 = -t47 - t51;
t59 = pkin(2) * t32;
t58 = pkin(2) * qJD(2) ^ 2;
t57 = pkin(3) * t38 ^ 2;
t54 = t29 * t38;
t14 = t26 * rSges(5,1) + t27 * rSges(5,2);
t11 = t34 * t14;
t49 = t32 * t58;
t48 = t33 * t58;
t18 = t30 * rSges(4,1) - t29 * rSges(4,2);
t13 = rSges(4,1) * t53 - rSges(4,2) * t54;
t45 = pkin(3) * t30 + t15;
t40 = -pkin(3) * t29 - t14;
t28 = pkin(2) * t33;
t8 = t38 * t18 + t46;
t6 = -t38 * t13 - t48;
t5 = -t38 * t51 - t49;
t3 = -pkin(3) * t54 - t11 - t47;
t2 = -t30 * t57 - t34 * t62 - t48;
t1 = -t11 * t34 - t29 * t57 - t49;
t9 = [0; m(4) * (t6 * (-t17 - t59) + t5 * (t18 + t28) + (-t13 - t46 + t8) * t7) + (t2 * (t40 - t59) + t1 * (t28 + t45)) * m(5); (t1 * t45 + t2 * t40) * m(5) + (-t8 * t51 - t7 * t13 - t6 * t17 + t5 * t18 - (-t17 * t8 - t18 * t7) * t38) * m(4); (t1 * t15 - t2 * t14 - (-t14 * t4 - t15 * t3) * t34 - t3 * t62 - t4 * t11) * m(5);];
tauc  = t9(:);
