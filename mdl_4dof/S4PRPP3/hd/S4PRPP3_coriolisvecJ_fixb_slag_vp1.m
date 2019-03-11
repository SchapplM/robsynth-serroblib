% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPP3
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:53
% EndTime: 2019-03-08 18:19:53
% DurationCPUTime: 0.34s
% Computational Cost: add. (230->74), mult. (508->94), div. (0->0), fcn. (253->2), ass. (0->45)
t44 = sin(qJ(2));
t45 = cos(qJ(2));
t26 = t45 * rSges(4,1) + t44 * rSges(4,3);
t34 = qJD(3) * t45;
t55 = t45 * pkin(2) + t44 * qJ(3);
t67 = qJD(2) * t55 - t34;
t47 = -qJD(2) * t26 - t67;
t69 = -rSges(5,1) - pkin(3);
t60 = -rSges(4,1) - pkin(2);
t66 = t60 * t44;
t25 = t45 * rSges(5,1) + t44 * rSges(5,2);
t48 = t45 * pkin(3) + t25;
t36 = t45 * qJ(3);
t20 = t44 * pkin(2) - t36;
t33 = qJD(3) * t44;
t52 = qJD(2) * t45;
t56 = qJ(3) * t52 + t33;
t65 = qJD(2) * t20 - t33 + t56;
t64 = -t45 / 0.2e1;
t51 = qJD(2) * qJD(3);
t29 = t45 * t51;
t46 = qJD(2) ^ 2;
t61 = t46 * pkin(3);
t2 = -t45 * t61 + t29 + (-qJD(2) * t25 - t67) * qJD(2);
t63 = t2 * t44;
t53 = qJD(2) * t44;
t59 = qJD(2) * (-pkin(2) * t53 + t56) + t44 * t51;
t27 = t45 * rSges(3,1) - t44 * rSges(3,2);
t54 = qJD(2) * t27;
t50 = -pkin(2) + t69;
t40 = t45 * rSges(5,2);
t49 = t69 * t44 + t40;
t23 = t44 * rSges(3,1) + t45 * rSges(3,2);
t39 = t45 * rSges(4,3);
t32 = rSges(5,2) * t52;
t31 = rSges(4,3) * t52;
t22 = t44 * rSges(4,1) - t39;
t16 = qJD(2) * t23;
t15 = qJD(1) + t54;
t9 = t33 + (-t20 + t49) * qJD(2);
t7 = t48 * qJD(2) + qJD(1) + t67;
t4 = t47 * qJD(2) + t29;
t3 = qJD(2) * (-rSges(4,1) * t53 + t31) + t59;
t1 = -t44 * t61 + qJD(2) * (-rSges(5,1) * t53 + t32) + t59;
t5 = [-m(3) * qJD(2) * t16 + m(4) * t3 + m(5) * t1; (t2 * (t36 + t40) + t1 * (t55 + t48) + t50 * t63 + (t67 + t34) * t9 + (t32 + t65) * t7 + ((t50 * t44 - t49) * t7 + (t48 + t50 * t45 + (-rSges(5,2) - qJ(3)) * t44) * t9) * qJD(2)) * m(5) + (t4 * (t36 + t39 + t66) + t3 * (t55 + t26) + (t31 + (t22 + t66) * qJD(2) + t65) * (qJD(1) - t47) + (-t47 + t34 + (t60 * t45 + (-rSges(4,3) - qJ(3)) * t44) * qJD(2)) * (t33 + (-t20 - t22) * qJD(2))) * m(4) + (-(-qJD(2) * t15 + t27 * t46) * t23 - t15 * t16 + (-t16 * t27 + 0.2e1 * t23 * t54) * qJD(2)) * m(3); 0.2e1 * (t1 * t64 + t63 / 0.2e1) * m(5) + 0.2e1 * (t3 * t64 + t4 * t44 / 0.2e1) * m(4); 0;];
tauc  = t5(:);
