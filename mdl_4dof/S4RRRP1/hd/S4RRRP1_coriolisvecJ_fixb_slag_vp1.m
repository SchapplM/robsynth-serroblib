% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:53
% EndTime: 2019-03-08 18:35:53
% DurationCPUTime: 0.54s
% Computational Cost: add. (1122->98), mult. (833->109), div. (0->0), fcn. (374->6), ass. (0->73)
t92 = rSges(5,1) + pkin(3);
t47 = sin(qJ(1));
t74 = pkin(1) * qJD(1);
t69 = t47 * t74;
t46 = qJ(1) + qJ(2);
t40 = sin(t46);
t41 = cos(t46);
t23 = t40 * rSges(3,1) + t41 * rSges(3,2);
t45 = qJD(1) + qJD(2);
t75 = t45 * t23;
t13 = -t69 - t75;
t39 = qJD(3) + t45;
t42 = qJ(3) + t46;
t36 = sin(t42);
t37 = cos(t42);
t89 = -t36 * rSges(5,2) + t92 * t37;
t91 = t89 * t39;
t77 = t40 * t45;
t73 = pkin(2) * t77;
t57 = -t69 - t73;
t88 = t57 + t69;
t38 = t39 ^ 2;
t87 = pkin(2) * t40;
t86 = pkin(2) * t45 ^ 2;
t85 = pkin(3) * t36;
t84 = pkin(3) * t38;
t83 = t47 * pkin(1);
t48 = cos(qJ(1));
t43 = t48 * pkin(1);
t81 = t36 * t39;
t80 = t37 * rSges(5,2);
t79 = t37 * t39;
t20 = t36 * rSges(4,1) + t37 * rSges(4,2);
t16 = t39 * t20;
t22 = t37 * rSges(4,1) - t36 * rSges(4,2);
t78 = t39 * t22;
t76 = t41 * t45;
t72 = pkin(2) * t76;
t49 = qJD(1) ^ 2;
t71 = t49 * t83;
t70 = t49 * t43;
t68 = t48 * t74;
t67 = t92 * t36;
t24 = t41 * rSges(3,1) - t40 * rSges(3,2);
t18 = rSges(3,1) * t76 - rSges(3,2) * t77;
t12 = rSges(4,1) * t79 - rSges(4,2) * t81;
t35 = pkin(2) * t41;
t65 = t22 + t35;
t19 = t36 * rSges(5,1) + t80;
t62 = t35 + t89;
t61 = -t67 - t80;
t60 = -t40 * t86 - t71;
t59 = -t41 * t86 - t70;
t56 = t68 + t72;
t55 = -t20 - t87;
t15 = t39 * t19;
t54 = -pkin(3) * t81 - t15 - t73;
t53 = -t12 - t72;
t52 = t61 - t87;
t3 = (-t19 - t85) * t39 + t57;
t4 = t56 + t91;
t50 = (-t4 * t67 + (-t4 * rSges(5,2) - t3 * t92) * t37) * t39;
t25 = rSges(5,2) * t81;
t14 = t45 * t24 + t68;
t10 = -t45 * t18 - t70;
t9 = -t45 * t75 - t71;
t8 = t56 + t78;
t7 = t57 - t16;
t6 = -t39 * t12 + t59;
t5 = -t16 * t39 + t60;
t2 = -t37 * t84 - t39 * (rSges(5,1) * t79 - t25) + t59;
t1 = -t19 * t38 - t36 * t84 + t60;
t11 = [m(3) * (t10 * (-t23 - t83) + t9 * (t24 + t43) + (-t18 - t68 + t14) * t13) + (t2 * (t52 - t83) + t3 * (t25 - t56) + t1 * (t43 + t62) + t50 + (t3 - t54 + t88) * t4) * m(5) + (t6 * (t55 - t83) + t7 * (t53 - t68) + t5 * (t43 + t65) + (t7 - t88 - t73) * t8) * m(4); (t1 * t62 + t2 * t52 + t50 + (-t54 - t73) * t4 + (t25 + t91) * t3) * m(5) + (t5 * t65 + t6 * t55 + (t72 + t78 + t53) * t7) * m(4) + (-(-t13 * t24 - t14 * t23) * t45 - t10 * t23 - t13 * t18 - t14 * t75 + t9 * t24) * m(3); (t1 * t89 + t2 * t61 + t3 * t25 + t50 + t4 * t15 - (-t3 * t89 - t4 * t85) * t39) * m(5) + (-t8 * t16 - t7 * t12 - t6 * t20 + t5 * t22 - (-t20 * t8 - t22 * t7) * t39) * m(4); 0;];
tauc  = t11(:);
