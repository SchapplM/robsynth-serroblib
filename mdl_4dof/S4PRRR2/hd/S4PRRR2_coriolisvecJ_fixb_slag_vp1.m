% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:22
% EndTime: 2019-07-18 13:27:22
% DurationCPUTime: 0.30s
% Computational Cost: add. (556->57), mult. (479->71), div. (0->0), fcn. (212->6), ass. (0->48)
t37 = qJD(2) + qJD(3);
t38 = qJ(2) + qJ(3);
t33 = sin(t38);
t34 = cos(t38);
t45 = t33 * rSges(4,1) + t34 * rSges(4,2);
t11 = t45 * t37;
t39 = sin(qJ(2));
t55 = pkin(1) * qJD(2);
t51 = t39 * t55;
t9 = t11 + t51;
t62 = pkin(2) * t37 ^ 2;
t61 = t39 * pkin(1);
t40 = cos(qJ(2));
t60 = t40 * pkin(1);
t35 = qJ(4) + t38;
t29 = sin(t35);
t32 = qJD(4) + t37;
t59 = t29 * t32;
t30 = cos(t35);
t58 = t30 * t32;
t57 = t33 * t37;
t56 = t34 * t37;
t7 = rSges(5,1) * t59 + rSges(5,2) * t58;
t54 = pkin(2) * t57;
t23 = pkin(2) * t56;
t41 = qJD(2) ^ 2;
t53 = t41 * t60;
t52 = t54 + t7;
t50 = t40 * t55;
t16 = t34 * rSges(4,1) - t33 * rSges(4,2);
t14 = t30 * rSges(5,1) - t29 * rSges(5,2);
t13 = -t29 * rSges(5,1) - t30 * rSges(5,2);
t49 = -t32 * t13 + t54;
t48 = -t32 * t14 - t23;
t12 = -rSges(4,1) * t56 + rSges(4,2) * t57;
t8 = -rSges(5,1) * t58 + rSges(5,2) * t59;
t44 = -pkin(2) * t34 - t14;
t10 = -t37 * t16 - t50;
t4 = t48 - t50;
t43 = -pkin(2) * t33 + t13;
t42 = t8 - t23;
t31 = t41 * t61;
t6 = t37 * t12 - t53;
t5 = t37 * t11 + t31;
t3 = t49 + t51;
t2 = t32 * t8 - t34 * t62 - t53;
t1 = t32 * t7 + t33 * t62 + t31;
t15 = [0; m(4) * (t5 * (-t16 - t60) + t6 * (-t45 - t61) + (t10 - t12 + t50) * t9) + m(5) * (t1 * (t44 - t60) + t4 * (t51 + t52) + t2 * (t43 - t61) - t3 * (t42 - t50)); (t1 * t44 + t2 * t43 + (-t49 + t52) * t4 + (-t42 + t48) * t3) * m(5) + (-(t10 * t45 + t16 * t9) * t37 + t10 * t11 - t9 * t12 - t6 * t45 - t5 * t16) * m(4); (-t1 * t14 + t2 * t13 - t3 * t8 + t4 * t7 - (-t4 * t13 + t3 * t14) * t32) * m(5);];
tauc  = t15(:);
