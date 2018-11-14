% Calculate vector of centrifugal and coriolis load on the joints for
% S4PRPP4
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
% tauc [4x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:09
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PRPP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:18
% EndTime: 2018-11-14 14:09:19
% DurationCPUTime: 0.40s
% Computational Cost: add. (299->40), mult. (398->58), div. (0->0), fcn. (188->4), ass. (0->31)
t36 = qJ(2) + pkin(5);
t34 = cos(t36);
t57 = rSges(5,3) + qJ(4);
t62 = t57 * t34;
t38 = cos(qJ(2));
t35 = t38 * pkin(2);
t33 = sin(t36);
t59 = t57 * t33;
t60 = rSges(5,1) + pkin(3);
t61 = t60 * t34 + t35 + t59;
t39 = qJD(2) ^ 2;
t15 = t34 * rSges(4,1) - t33 * rSges(4,2);
t37 = sin(qJ(2));
t55 = t37 * pkin(2);
t56 = -t55 + t62;
t20 = t38 * rSges(3,1) - t37 * rSges(3,2);
t50 = qJD(2) * t20;
t47 = t39 * t55;
t46 = t39 * t35;
t23 = qJD(4) * t33;
t45 = qJD(2) * t62 + t23;
t19 = t37 * rSges(3,1) + t38 * rSges(3,2);
t12 = t33 * rSges(4,1) + t34 * rSges(4,2);
t41 = -qJD(2) * t60 + qJD(4);
t24 = qJD(4) * t34;
t17 = t19 * qJD(2);
t16 = qJD(1) + t50;
t7 = -t12 * t39 - t47;
t2 = -t46 + (-qJD(2) * t59 + t41 * t34 + t24) * qJD(2);
t1 = -t47 + (t41 * t33 + t45) * qJD(2);
t3 = [-m(3) * qJD(2) * t17 + m(4) * t7 + m(5) * t1; (t7 * (t15 + t35) + (-t15 * t39 - t46) * (-t12 - t55)) * m(4) + (-t16 * t17 + (-t17 * t20 + 0.2e1 * t19 * t50) * qJD(2) - (-qJD(2) * t16 + t20 * t39) * t19) * m(3) + (t1 * t61 + t2 * (-t60 * t33 + t56) + (t45 - t23 + (-t56 - t55) * qJD(2)) * (t61 * qJD(2) + qJD(1) - t24)) * m(5); 0; 0.2e1 * (-t1 * t34 / 0.2e1 + t2 * t33 / 0.2e1) * m(5);];
tauc  = t3(:);
