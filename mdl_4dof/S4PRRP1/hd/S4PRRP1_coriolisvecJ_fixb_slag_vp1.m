% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:49
% EndTime: 2019-03-08 18:22:49
% DurationCPUTime: 0.35s
% Computational Cost: add. (830->56), mult. (543->70), div. (0->0), fcn. (268->4), ass. (0->45)
t42 = pkin(6) + qJ(2);
t41 = qJ(3) + t42;
t38 = cos(t41);
t30 = t38 * qJ(4);
t37 = sin(t41);
t69 = rSges(5,1) + pkin(3);
t74 = -t38 * rSges(5,3) + t37 * t69 - t30;
t27 = qJD(4) * t37;
t43 = qJD(2) + qJD(3);
t67 = t74 * t43;
t73 = t27 - t67;
t72 = t69 * t38;
t39 = sin(t42);
t59 = pkin(2) * qJD(2);
t54 = t39 * t59;
t18 = t37 * rSges(4,1) + t38 * rSges(4,2);
t60 = t43 * t18;
t9 = -t54 - t60;
t68 = rSges(5,3) + qJ(4);
t61 = t38 * t43;
t66 = rSges(5,3) * t61 + t43 * t30;
t65 = pkin(2) * t39;
t64 = pkin(2) * qJD(2) ^ 2;
t62 = t37 * t43;
t52 = t68 * t37 + t72;
t57 = t39 * t64;
t40 = cos(t42);
t56 = t40 * t64;
t55 = t27 + t66;
t53 = t40 * t59;
t21 = t38 * rSges(4,1) - t37 * rSges(4,2);
t12 = rSges(4,1) * t61 - rSges(4,2) * t62;
t28 = qJD(4) * t38;
t50 = t28 - t53;
t48 = -t43 * t69 + qJD(4);
t5 = -t54 + t73;
t6 = t52 * t43 - t50;
t45 = (-t5 * t72 + (-t5 * t68 - t6 * t69) * t37) * t43;
t36 = pkin(2) * t40;
t10 = t43 * t21 + t53;
t8 = -t43 * t12 - t56;
t7 = -t43 * t60 - t57;
t2 = -t56 + (t48 * t38 - t62 * t68 + t28) * t43;
t1 = -t57 + (t48 * t37 + t55) * t43;
t3 = [0; m(4) * (t8 * (-t18 - t65) + t7 * (t21 + t36) + (-t12 - t53 + t10) * t9) + (t2 * (-t74 - t65) + t5 * t50 + t1 * (t36 + t52) + t45 + (t5 + t66 + t67) * t6) * m(5); (-t2 * t74 + t45 + (t55 - t73) * t6 + (t5 * t43 + t1) * t52) * m(5) + (-t10 * t60 - t9 * t12 - t8 * t18 + t7 * t21 - (-t10 * t18 - t21 * t9) * t43) * m(4); 0.2e1 * (-t1 * t38 / 0.2e1 + t2 * t37 / 0.2e1) * m(5);];
tauc  = t3(:);
