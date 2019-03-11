% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:59
% EndTime: 2019-03-08 18:23:59
% DurationCPUTime: 0.37s
% Computational Cost: add. (382->61), mult. (438->87), div. (0->0), fcn. (189->4), ass. (0->48)
t57 = rSges(5,1) + pkin(3);
t33 = sin(qJ(2));
t46 = pkin(2) * qJD(2);
t42 = t33 * t46;
t32 = qJ(2) + qJ(3);
t27 = sin(t32);
t28 = cos(t32);
t14 = t27 * rSges(4,1) + t28 * rSges(4,2);
t31 = qJD(2) + qJD(3);
t47 = t31 * t14;
t8 = -t42 - t47;
t55 = -t27 * rSges(5,2) + t57 * t28;
t56 = t55 * t31;
t30 = t31 ^ 2;
t54 = pkin(3) * t27;
t53 = t30 * pkin(3);
t52 = t33 * pkin(2);
t34 = cos(qJ(2));
t29 = t34 * pkin(2);
t50 = t27 * t31;
t49 = t28 * rSges(5,2);
t48 = t28 * t31;
t22 = t34 * rSges(3,1) - t33 * rSges(3,2);
t45 = qJD(2) * t22;
t35 = qJD(2) ^ 2;
t44 = t35 * t52;
t43 = t35 * t29;
t41 = t57 * t27;
t16 = t28 * rSges(4,1) - t27 * rSges(4,2);
t26 = t34 * t46;
t39 = -t31 * t16 - t26;
t10 = rSges(4,1) * t48 - rSges(4,2) * t50;
t21 = t33 * rSges(3,1) + t34 * rSges(3,2);
t13 = t27 * rSges(5,1) + t49;
t37 = -t41 - t49;
t3 = qJD(1) + t26 + t56;
t4 = -t42 + (-t13 - t54) * t31;
t36 = (-t3 * t41 + (-t3 * rSges(5,2) - t4 * t57) * t28) * t31;
t19 = rSges(5,2) * t50;
t17 = t21 * qJD(2);
t12 = qJD(1) + t45;
t11 = t31 * t13;
t7 = qJD(1) - t39;
t6 = -t31 * t10 - t43;
t5 = -t31 * t47 - t44;
t2 = -t43 - t28 * t53 - t31 * (rSges(5,1) * t48 - t19);
t1 = -t13 * t30 - t27 * t53 - t44;
t9 = [-m(3) * qJD(2) * t17 + m(4) * t5 + m(5) * t1; (-t3 * (-pkin(3) * t50 - t11) + t2 * (t37 - t52) + t1 * (t29 + t55) + t36 + (t19 + t56) * t4) * m(5) + (t6 * (-t14 - t52) + t5 * (t16 + t29) + (-t39 - t10 - t26) * t8) * m(4) + (-(-qJD(2) * t12 + t22 * t35) * t21 - t12 * t17 + (-t17 * t22 + 0.2e1 * t21 * t45) * qJD(2)) * m(3); (t1 * t55 + t4 * t19 + t2 * t37 + t36 + t3 * t11 - (-t3 * t54 - t4 * t55) * t31) * m(5) + (-t8 * t10 - t6 * t14 + t5 * t16 - t7 * t47 - (-t7 * t14 - t8 * t16) * t31) * m(4); 0;];
tauc  = t9(:);
