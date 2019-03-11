% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
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
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:54
% EndTime: 2019-03-08 18:17:55
% DurationCPUTime: 0.31s
% Computational Cost: add. (562->72), mult. (575->89), div. (0->0), fcn. (300->2), ass. (0->45)
t69 = -rSges(5,3) - qJ(4);
t41 = pkin(5) + qJ(2);
t39 = sin(t41);
t40 = cos(t41);
t21 = t40 * pkin(2) + t39 * qJ(3);
t66 = -t39 * rSges(5,2) - t40 * rSges(5,3);
t68 = t40 * qJ(4) + t21 - t66;
t65 = t40 * rSges(4,2) - t39 * rSges(4,3);
t67 = t21 - t65;
t64 = t39 / 0.2e1;
t63 = -t40 / 0.2e1;
t62 = rSges(4,2) - pkin(2);
t60 = t40 * rSges(5,2);
t59 = t40 * rSges(4,3);
t49 = qJD(2) * qJD(3);
t52 = qJD(2) * t39;
t51 = qJD(2) * t40;
t26 = qJ(3) * t51;
t31 = qJD(3) * t39;
t57 = t26 + t31;
t58 = qJD(2) * (-pkin(2) * t52 + t57) + t39 * t49;
t56 = rSges(4,2) * t52 + rSges(4,3) * t51;
t30 = qJD(4) * t40;
t55 = t30 + t31;
t34 = t40 * qJ(3);
t18 = t39 * pkin(2) - t34;
t54 = -qJD(2) * t18 + t31;
t53 = qJD(2) ^ 2 * qJ(4);
t50 = qJD(4) * t39;
t48 = -pkin(2) + t69;
t46 = t69 * t39 + t60;
t32 = qJD(3) * t40;
t29 = rSges(5,2) * t51;
t25 = t40 * t49;
t19 = t39 * rSges(4,2) + t59;
t15 = qJD(2) * t21 - t32;
t13 = t67 * qJD(2) - t32;
t12 = t31 + (-t18 + t19) * qJD(2);
t9 = t25 + (t65 * qJD(2) - t15) * qJD(2);
t8 = qJD(2) * t56 + t58;
t7 = t68 * qJD(2) - t32 + t50;
t6 = (-t18 + t46) * qJD(2) + t55;
t2 = -t40 * t53 + t25 + (t66 * qJD(2) - t15 - 0.2e1 * t50) * qJD(2);
t1 = -t39 * t53 + (-rSges(5,3) * t52 + t29 + 0.2e1 * t30) * qJD(2) + t58;
t3 = [0; (-(t46 * qJD(2) + t30 + t54 - t6) * t7 + t2 * (t34 + t60) + t6 * t32 + t1 * t68 + t7 * (t26 + t29 + t55) + (-t6 * qJD(4) + t2 * t48) * t39 + (t6 * t48 * t40 + (t6 * (-rSges(5,2) - qJ(3)) + t7 * t48) * t39) * qJD(2)) * m(5) + (-(qJD(2) * t19 - t12 + t54) * t13 + t9 * (t62 * t39 + t34 + t59) + t12 * t32 + t8 * t67 + t13 * (t56 + t57) + (t12 * t62 * t40 + (t12 * (-rSges(4,3) - qJ(3)) - t13 * pkin(2)) * t39) * qJD(2)) * m(4); 0.2e1 * (t1 * t63 + t2 * t64) * m(5) + 0.2e1 * (t8 * t63 + t9 * t64) * m(4); m(5) * (t1 * t39 + t2 * t40);];
tauc  = t3(:);
