% Calculate vector of centrifugal and coriolis load on the joints for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPRR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:29
% EndTime: 2018-11-14 13:50:29
% DurationCPUTime: 0.48s
% Computational Cost: add. (838->50), mult. (649->67), div. (0->0), fcn. (290->8), ass. (0->41)
t43 = qJ(1) + pkin(7);
t78 = pkin(2) * cos(t43) + cos(qJ(1)) * pkin(1);
t51 = t78 * qJD(1);
t39 = qJ(3) + t43;
t34 = cos(t39);
t42 = qJD(1) + qJD(3);
t66 = t34 * t42;
t35 = qJ(4) + t39;
t29 = sin(t35);
t30 = cos(t35);
t17 = t30 * rSges(5,1) - rSges(5,2) * t29;
t38 = qJD(4) + t42;
t82 = t17 * t38;
t4 = pkin(3) * t66 + t51 + t82;
t58 = -pkin(2) * sin(t43) - sin(qJ(1)) * pkin(1);
t52 = t58 * qJD(1);
t46 = qJD(1) ^ 2;
t72 = pkin(3) * t42 ^ 2;
t33 = sin(t39);
t18 = rSges(4,1) * t33 + t34 * rSges(4,2);
t27 = t34 * rSges(4,1);
t19 = -rSges(4,2) * t33 + t27;
t8 = t42 * t19 + t51;
t71 = t18 * t8;
t67 = t33 * t42;
t16 = rSges(5,1) * t29 + rSges(5,2) * t30;
t11 = t38 * t16;
t64 = t42 * t18;
t59 = pkin(3) * t34 + t17;
t54 = t58 * t46;
t53 = t78 * t46;
t48 = -pkin(3) * t33 - t16;
t7 = t52 - t64;
t23 = rSges(4,2) * t67;
t13 = rSges(4,1) * t66 - t23;
t6 = -t42 * t13 - t53;
t5 = -t42 * t64 + t54;
t3 = -pkin(3) * t67 - t11 + t52;
t2 = -t34 * t72 - t38 * t82 - t53;
t1 = -t11 * t38 - t33 * t72 + t54;
t9 = [(t2 * (t48 + t58) + t1 * (t59 + t78)) * m(5) + (t6 * (-t18 + t58) + t7 * t23 + t5 * (t19 + t78) + (-t7 * t27 - t71) * t42 + (t8 * t58 - t7 * t78) * qJD(1)) * m(4); 0; (t1 * t59 + t2 * t48) * m(5) + (-t64 * t8 - t13 * t7 - t18 * t6 + t19 * t5 - (-t19 * t7 - t71) * t42) * m(4); (t1 * t17 - t2 * t16 - (-t16 * t4 - t17 * t3) * t38 - t3 * t82 - t4 * t11) * m(5);];
tauc  = t9(:);
