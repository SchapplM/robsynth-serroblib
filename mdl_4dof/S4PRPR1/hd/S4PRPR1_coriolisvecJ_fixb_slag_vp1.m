% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:52
% EndTime: 2019-03-08 18:20:53
% DurationCPUTime: 0.35s
% Computational Cost: add. (752->77), mult. (775->101), div. (0->0), fcn. (544->4), ass. (0->51)
t50 = qJD(2) - qJD(4);
t49 = pkin(6) + qJ(2);
t47 = sin(t49);
t48 = cos(t49);
t51 = sin(qJ(4));
t67 = cos(qJ(4));
t54 = t47 * t51 + t48 * t67;
t58 = t47 * t67;
t74 = -t48 * t51 + t58;
t56 = -rSges(5,1) * t74 + rSges(5,2) * t54;
t77 = t50 * t56;
t28 = t48 * pkin(2) + t47 * qJ(3);
t62 = t48 * rSges(4,1) + t47 * rSges(4,3);
t76 = t28 + t62;
t69 = t48 * pkin(3);
t75 = t28 + t69;
t73 = t47 / 0.2e1;
t72 = -t48 / 0.2e1;
t71 = -pkin(2) - pkin(3);
t70 = t47 * pkin(3);
t68 = -rSges(4,1) - pkin(2);
t59 = qJD(2) * qJD(3);
t61 = qJD(2) * t47;
t39 = qJD(3) * t47;
t60 = qJD(2) * t48;
t64 = qJ(3) * t60 + t39;
t65 = qJD(2) * (-pkin(2) * t61 + t64) + t47 * t59;
t42 = t48 * qJ(3);
t25 = t47 * pkin(2) - t42;
t63 = -qJD(2) * t25 + t39;
t15 = t50 * t54;
t16 = -qJD(2) * t58 + t74 * qJD(4) + t51 * t60;
t4 = t16 * rSges(5,1) + t15 * rSges(5,2);
t3 = t15 * rSges(5,1) - t16 * rSges(5,2);
t55 = -rSges(5,1) * t54 - rSges(5,2) * t74;
t52 = qJD(2) ^ 2;
t44 = t48 * rSges(4,3);
t40 = qJD(3) * t48;
t38 = rSges(4,3) * t60;
t33 = t48 * t59;
t26 = t47 * rSges(4,1) - t44;
t20 = qJD(2) * t28 - t40;
t14 = t76 * qJD(2) - t40;
t13 = t39 + (-t25 - t26) * qJD(2);
t10 = t33 + (-t62 * qJD(2) - t20) * qJD(2);
t9 = qJD(2) * (-rSges(4,1) * t61 + t38) + t65;
t8 = t75 * qJD(2) - t55 * t50 - t40;
t7 = t77 + t39 + (-t25 - t70) * qJD(2);
t2 = -qJD(2) * t20 - t50 * t3 - t52 * t69 + t33;
t1 = t50 * t4 - t52 * t70 + t65;
t5 = [0; (-(-pkin(3) * t61 + t63 - t7 + t77) * t8 + t2 * (t71 * t47 + t42 + t56) + t7 * (-t3 + t40) + t1 * (-t55 + t75) + t8 * (t4 + t64) + (t7 * t71 * t48 + (-t7 * qJ(3) + t8 * t71) * t47) * qJD(2)) * m(5) + (-(-qJD(2) * t26 - t13 + t63) * t14 + t10 * (t68 * t47 + t42 + t44) + t13 * t40 + t9 * t76 + t14 * (t38 + t64) + (t13 * t68 * t48 + (t13 * (-rSges(4,3) - qJ(3)) + t14 * t68) * t47) * qJD(2)) * m(4); 0.2e1 * (t1 * t72 + t2 * t73) * m(5) + 0.2e1 * (t10 * t73 + t9 * t72) * m(4); (t1 * t55 - t2 * t56 + t7 * t3 - t8 * t4 - (-t7 * t55 - t56 * t8) * t50) * m(5);];
tauc  = t5(:);
