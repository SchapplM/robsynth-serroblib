% Calculate vector of centrifugal and coriolis load on the joints for
% S4PRRP3
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:18
% EndTime: 2018-11-14 14:12:19
% DurationCPUTime: 0.35s
% Computational Cost: add. (387->62), mult. (438->81), div. (0->0), fcn. (189->4), ass. (0->46)
t34 = sin(qJ(2));
t49 = pkin(2) * qJD(2);
t44 = t34 * t49;
t33 = qJ(2) + qJ(3);
t29 = sin(t33);
t30 = cos(t33);
t14 = rSges(4,1) * t29 + rSges(4,2) * t30;
t32 = qJD(2) + qJD(3);
t54 = t14 * t32;
t8 = -t44 - t54;
t56 = pkin(2) * t34;
t35 = cos(qJ(2));
t31 = t35 * pkin(2);
t55 = -rSges(5,1) - pkin(3);
t53 = t29 * t32;
t52 = t30 * rSges(5,2);
t51 = t30 * t32;
t42 = t30 * rSges(5,1) - rSges(5,2) * t29;
t50 = -pkin(3) * t51 - t32 * t42;
t24 = rSges(3,1) * t35 - rSges(3,2) * t34;
t48 = qJD(2) * t24;
t36 = qJD(2) ^ 2;
t47 = t36 * t56;
t46 = t36 * t31;
t28 = t35 * t49;
t45 = -t28 + t50;
t43 = t55 * t29;
t16 = t30 * rSges(4,1) - rSges(4,2) * t29;
t41 = -(rSges(5,1) * t29 + t52) * t32 - pkin(3) * t53;
t40 = -t16 * t32 - t28;
t10 = rSges(4,1) * t51 - rSges(4,2) * t53;
t39 = pkin(3) * t30 + t42;
t23 = rSges(3,1) * t34 + rSges(3,2) * t35;
t38 = t43 - t52;
t4 = t41 - t44;
t3 = qJD(1) - t45;
t37 = (t3 * t43 + (-t3 * rSges(5,2) + t4 * t55) * t30) * t32;
t19 = rSges(5,2) * t53;
t17 = t23 * qJD(2);
t12 = qJD(1) + t48;
t7 = qJD(1) - t40;
t6 = -t10 * t32 - t46;
t5 = -t32 * t54 - t47;
t2 = -t46 + (t55 * t51 + t19) * t32;
t1 = t38 * t32 ^ 2 - t47;
t9 = [-m(3) * qJD(2) * t17 + m(4) * t5 + m(5) * t1; (t1 * (t31 + t39) - t3 * t44 + t2 * (t38 - t56) + t37 + (t19 - t28 - t3 - t45) * t4) * m(5) + (t5 * (t16 + t31) + t6 * (-t14 - t56) + (-t10 - t28 - t40) * t8) * m(4) + (-t12 * t17 + (-t17 * t24 + 0.2e1 * t23 * t48) * qJD(2) - (-qJD(2) * t12 + t24 * t36) * t23) * m(3); (t1 * t39 + t2 * t38 - t3 * t41 + t37 + (t19 - t50) * t4) * m(5) + (-(-t7 * t14 - t8 * t16) * t32 - t8 * t10 - t6 * t14 + t5 * t16 - t7 * t54) * m(4); 0;];
tauc  = t9(:);
