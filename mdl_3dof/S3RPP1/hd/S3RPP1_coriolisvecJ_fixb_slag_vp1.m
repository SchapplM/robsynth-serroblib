% Calculate vector of centrifugal and coriolis load on the joints for
% S3RPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% tauc [3x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:14
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S3RPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'S3RPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:13:33
% EndTime: 2018-11-14 10:13:34
% DurationCPUTime: 0.33s
% Computational Cost: add. (262->71), mult. (575->89), div. (0->0), fcn. (300->2), ass. (0->44)
t68 = -rSges(4,3) - qJ(3);
t39 = sin(qJ(1));
t40 = cos(qJ(1));
t21 = t40 * pkin(1) + t39 * qJ(2);
t65 = -t39 * rSges(4,2) - t40 * rSges(4,3);
t67 = t40 * qJ(3) + t21 - t65;
t64 = t40 * rSges(3,2) - t39 * rSges(3,3);
t66 = t21 - t64;
t63 = t39 / 0.2e1;
t62 = -t40 / 0.2e1;
t61 = rSges(3,2) - pkin(1);
t59 = t40 * rSges(4,2);
t58 = t40 * rSges(3,3);
t48 = qJD(1) * qJD(2);
t51 = qJD(1) * t39;
t50 = qJD(1) * t40;
t26 = qJ(2) * t50;
t31 = qJD(2) * t39;
t56 = t26 + t31;
t57 = qJD(1) * (-pkin(1) * t51 + t56) + t39 * t48;
t55 = rSges(3,2) * t51 + rSges(3,3) * t50;
t30 = qJD(3) * t40;
t54 = t30 + t31;
t34 = t40 * qJ(2);
t18 = t39 * pkin(1) - t34;
t53 = -qJD(1) * t18 + t31;
t52 = qJD(1) ^ 2 * qJ(3);
t49 = qJD(3) * t39;
t47 = -pkin(1) + t68;
t45 = t68 * t39 + t59;
t32 = qJD(2) * t40;
t29 = rSges(4,2) * t50;
t25 = t40 * t48;
t19 = t39 * rSges(3,2) + t58;
t15 = qJD(1) * t21 - t32;
t13 = t66 * qJD(1) - t32;
t12 = t31 + (-t18 + t19) * qJD(1);
t9 = t67 * qJD(1) - t32 + t49;
t8 = (-t18 + t45) * qJD(1) + t54;
t7 = t25 + (t64 * qJD(1) - t15) * qJD(1);
t6 = qJD(1) * t55 + t57;
t2 = -t40 * t52 + t25 + (t65 * qJD(1) - t15 - 0.2e1 * t49) * qJD(1);
t1 = -t39 * t52 + (-rSges(4,3) * t51 + t29 + 0.2e1 * t30) * qJD(1) + t57;
t3 = [(-(t45 * qJD(1) + t30 + t53 - t8) * t9 + t2 * (t34 + t59) + t8 * t32 + t1 * t67 + t9 * (t26 + t29 + t54) + (-t8 * qJD(3) + t2 * t47) * t39 + (t8 * t47 * t40 + (t8 * (-rSges(4,2) - qJ(2)) + t9 * t47) * t39) * qJD(1)) * m(4) + (-(qJD(1) * t19 - t12 + t53) * t13 + t7 * (t61 * t39 + t34 + t58) + t12 * t32 + t6 * t66 + t13 * (t55 + t56) + (t12 * t61 * t40 + (t12 * (-rSges(3,3) - qJ(2)) - t13 * pkin(1)) * t39) * qJD(1)) * m(3); 0.2e1 * (t1 * t62 + t2 * t63) * m(4) + 0.2e1 * (t6 * t62 + t7 * t63) * m(3); m(4) * (t1 * t39 + t2 * t40);];
tauc  = t3(:);
