% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:32:58
% EndTime: 2019-03-08 18:32:59
% DurationCPUTime: 0.58s
% Computational Cost: add. (1056->94), mult. (803->108), div. (0->0), fcn. (388->6), ass. (0->69)
t54 = qJ(1) + qJ(2);
t49 = sin(t54);
t90 = pkin(2) * t49;
t48 = pkin(6) + t54;
t46 = cos(t48);
t38 = t46 * qJ(4);
t45 = sin(t48);
t87 = -rSges(5,1) - pkin(3);
t97 = t46 * rSges(5,3) + t87 * t45 + t38;
t60 = -t90 + t97;
t98 = rSges(5,3) + qJ(4);
t50 = cos(t54);
t47 = pkin(2) * t50;
t96 = -t46 * rSges(4,1) - t47;
t55 = sin(qJ(1));
t82 = pkin(1) * qJD(1);
t75 = t55 * t82;
t26 = rSges(3,1) * t49 + rSges(3,2) * t50;
t53 = qJD(1) + qJD(2);
t86 = t26 * t53;
t13 = -t75 - t86;
t95 = t87 * t46 - t47;
t94 = -rSges(4,2) * t45 - t96;
t84 = t46 * t53;
t93 = rSges(5,3) * t84 + t53 * t38;
t92 = t45 * t98 - t95;
t52 = t53 ^ 2;
t91 = pkin(1) * t55;
t89 = pkin(2) * t52;
t56 = cos(qJ(1));
t51 = t56 * pkin(1);
t85 = t45 * t53;
t83 = t49 * t53;
t44 = t50 * rSges(3,1);
t36 = qJD(4) * t46;
t74 = t56 * t82;
t69 = t36 - t74;
t6 = t53 * t92 - t69;
t80 = t6 * t90;
t57 = qJD(1) ^ 2;
t79 = t57 * t91;
t78 = t57 * t51;
t35 = qJD(4) * t45;
t77 = t97 * t53 + t35;
t76 = t35 + t93;
t22 = rSges(4,1) * t45 + rSges(4,2) * t46;
t61 = -t22 - t90;
t27 = -rSges(3,2) * t49 + t44;
t19 = -rSges(3,2) * t83 + t53 * t44;
t70 = t35 - t75;
t66 = t87 * t53 + qJD(4);
t65 = -t49 * t89 - t79;
t64 = -t50 * t89 - t78;
t62 = -pkin(2) * t83 - t75;
t10 = t53 * t94 + t74;
t9 = t61 * t53 - t75;
t59 = (t10 * t61 + t9 * t96) * t53;
t5 = t60 * t53 + t70;
t58 = (t5 * t95 - t80 + (-t5 * t98 + t6 * t87) * t45) * t53;
t30 = rSges(4,2) * t85;
t17 = t53 * t22;
t14 = t27 * t53 + t74;
t12 = -t19 * t53 - t78;
t11 = -t53 * t86 - t79;
t8 = -t53 * (rSges(4,1) * t84 - t30) + t64;
t7 = -t22 * t52 + t65;
t2 = (t66 * t46 - t85 * t98 + t36) * t53 + t64;
t1 = (t66 * t45 + t76) * t53 + t65;
t3 = [m(3) * (t12 * (-t26 - t91) + t11 * (t27 + t51) + (-t19 - t74 + t14) * t13) + (t2 * (t60 - t91) + t5 * t69 + t1 * (t51 + t92) + t58 + (t70 + t5 - t62 - t77 + t93) * t6) * m(5) + (t8 * (t61 - t91) + t9 * (t30 - t74) + t7 * (t51 + t94) + t59 + (-t75 + t17 - t62 + t9) * t10) * m(4); (t10 * t17 - (-t10 * t90 - t9 * t94) * t53 + t9 * t30 + t8 * t61 + t7 * t94 + t59) * m(4) + (-(-t13 * t27 - t14 * t26) * t53 + t11 * t27 - t12 * t26 - t13 * t19 - t14 * t86) * m(3) + (-(-t5 * t92 - t80) * t53 + t1 * t92 + t2 * t60 + t58 + (-t77 + t76) * t6) * m(5); 0; 0.2e1 * (-t1 * t46 / 0.2e1 + t2 * t45 / 0.2e1) * m(5);];
tauc  = t3(:);
