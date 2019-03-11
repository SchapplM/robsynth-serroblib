% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:59
% EndTime: 2019-03-08 18:34:00
% DurationCPUTime: 0.55s
% Computational Cost: add. (1064->104), mult. (951->122), div. (0->0), fcn. (480->4), ass. (0->73)
t60 = sin(qJ(1));
t83 = pkin(1) * qJD(1);
t77 = t60 * t83;
t59 = qJ(1) + qJ(2);
t54 = sin(t59);
t55 = cos(t59);
t27 = rSges(3,1) * t54 + rSges(3,2) * t55;
t58 = qJD(1) + qJD(2);
t89 = t27 * t58;
t17 = -t77 - t89;
t28 = t55 * pkin(2) + t54 * qJ(3);
t29 = t55 * rSges(5,1) + t54 * rSges(5,2);
t52 = t55 * pkin(3);
t96 = t52 + t28 + t29;
t102 = t96 * t58;
t48 = t55 * rSges(5,2);
t25 = rSges(5,1) * t54 - t48;
t87 = t55 * t58;
t38 = rSges(5,2) * t87;
t101 = t58 * t25 + t38;
t47 = t55 * rSges(4,3);
t26 = rSges(4,1) * t54 - t47;
t37 = rSges(4,3) * t87;
t100 = t58 * t26 + t37;
t99 = -t55 * rSges(4,1) - t54 * rSges(4,3);
t44 = t55 * qJ(3);
t24 = pkin(2) * t54 - t44;
t21 = t58 * t24;
t35 = t58 * t44;
t98 = t35 + t21;
t41 = qJD(3) * t54;
t86 = t35 + t41;
t97 = t86 - t41 + t21;
t95 = -t55 / 0.2e1;
t94 = pkin(1) * t60;
t93 = pkin(3) * t54;
t42 = qJD(3) * t55;
t16 = t28 * t58 - t42;
t57 = t58 ^ 2;
t61 = cos(qJ(1));
t56 = t61 * pkin(1);
t62 = qJD(1) ^ 2;
t79 = t62 * t56;
t82 = qJD(3) * t58;
t75 = t55 * t82 - t79;
t2 = -t57 * t52 + (-t29 * t58 - t16) * t58 + t75;
t92 = t2 * t54;
t90 = -rSges(4,1) - pkin(2);
t88 = t54 * t58;
t76 = t28 - t99;
t84 = t44 + t48;
t81 = -rSges(5,1) - pkin(2) - pkin(3);
t80 = t62 * t94;
t78 = t61 * t83;
t31 = t55 * rSges(3,1) - rSges(3,2) * t54;
t20 = rSges(3,1) * t87 - rSges(3,2) * t88;
t73 = t41 - t77;
t72 = -t42 + t78;
t69 = t58 * (-pkin(2) * t88 + t86) + t54 * t82 - t80;
t66 = t90 * t54 + t44 + t47;
t11 = (-t24 - t26) * t58 + t73;
t12 = t76 * t58 + t72;
t64 = (t11 * t90 * t55 + (t11 * (-rSges(4,3) - qJ(3)) + t12 * t90) * t54) * t58;
t7 = (-t24 - t25 - t93) * t58 + t73;
t8 = t72 + t102;
t63 = (t7 * t81 * t55 + (t7 * (-rSges(5,2) - qJ(3)) + t8 * t81) * t54) * t58 + t81 * t92;
t18 = t31 * t58 + t78;
t15 = -t20 * t58 - t79;
t14 = -t58 * t89 - t80;
t4 = (t99 * t58 - t16) * t58 + t75;
t3 = t58 * (-rSges(4,1) * t88 + t37) + t69;
t1 = -t57 * t93 + t58 * (-rSges(5,1) * t88 + t38) + t69;
t5 = [m(3) * (t15 * (-t27 - t94) + t14 * (t31 + t56) + (-t20 - t78 + t18) * t17) + (t2 * (t84 - t94) - t7 * t72 + t1 * (t56 + t96) + t63 + (pkin(3) * t88 + t101 + t7 + t98) * t8) * m(5) + (t4 * (t66 - t94) - t11 * t72 + t3 * (t56 + t76) + t64 + (t11 + t98 + t100) * t12) * m(4); (t4 * t66 + t64 + (t11 * t58 + t3) * t76 + (t97 + t100) * t12) * m(4) + (-(-t17 * t31 - t18 * t27) * t58 + t14 * t31 - t15 * t27 - t17 * t20 - t18 * t89) * m(3) + (t7 * t102 + t1 * t96 + t2 * t84 + t63 + (t93 * t58 + t101 + t97) * t8) * m(5); 0.2e1 * (t1 * t95 + t92 / 0.2e1) * m(5) + 0.2e1 * (t3 * t95 + t4 * t54 / 0.2e1) * m(4); 0;];
tauc  = t5(:);
