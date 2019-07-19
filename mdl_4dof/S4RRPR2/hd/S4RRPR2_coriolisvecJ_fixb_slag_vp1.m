% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR2
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:31
% EndTime: 2019-07-18 18:16:32
% DurationCPUTime: 0.68s
% Computational Cost: add. (1535->122), mult. (1375->137), div. (0->0), fcn. (900->6), ass. (0->83)
t70 = qJD(1) + qJD(2);
t65 = -qJD(4) + t70;
t108 = cos(qJ(4));
t71 = qJ(1) + qJ(2);
t66 = sin(t71);
t67 = cos(t71);
t72 = sin(qJ(4));
t35 = -t66 * t108 + t67 * t72;
t80 = t67 * t108 + t66 * t72;
t86 = rSges(5,1) * t35 + rSges(5,2) * t80;
t120 = t65 * t86;
t38 = rSges(3,1) * t66 + rSges(3,2) * t67;
t106 = t38 * t70;
t100 = pkin(1) * qJD(1);
t73 = sin(qJ(1));
t96 = t73 * t100;
t26 = -t96 - t106;
t58 = t67 * qJ(3);
t36 = pkin(2) * t66 - t58;
t32 = t70 * t36;
t55 = qJD(3) * t66;
t101 = t55 - t32;
t105 = t66 * t70;
t18 = t65 * t80;
t19 = t65 * t35;
t6 = t19 * rSges(5,1) + t18 * rSges(5,2);
t119 = pkin(3) * t105 - t101 - t120 + t6;
t118 = -rSges(5,1) * t80 + t35 * rSges(5,2);
t60 = t67 * rSges(4,3);
t37 = rSges(4,1) * t66 - t60;
t104 = t67 * t70;
t48 = rSges(4,3) * t104;
t117 = t70 * t37 + t48;
t39 = t67 * pkin(2) + t66 * qJ(3);
t63 = t67 * pkin(3);
t93 = t39 + t63;
t116 = -t67 * rSges(4,1) - t66 * rSges(4,3);
t115 = -t118 * t65 + t93 * t70;
t114 = t66 / 0.2e1;
t113 = -t67 / 0.2e1;
t112 = -pkin(2) - pkin(3);
t111 = pkin(1) * t73;
t110 = pkin(3) * t66;
t74 = cos(qJ(1));
t68 = t74 * pkin(1);
t109 = -rSges(4,1) - pkin(2);
t92 = t39 - t116;
t46 = t70 * t58;
t102 = t46 + t55;
t99 = qJD(3) * t70;
t75 = qJD(1) ^ 2;
t98 = t75 * t111;
t97 = t75 * t68;
t95 = t74 * t100;
t41 = t67 * rSges(3,1) - rSges(3,2) * t66;
t91 = t67 * t99 - t97;
t29 = rSges(3,1) * t104 - rSges(3,2) * t105;
t90 = t55 - t96;
t56 = qJD(3) * t67;
t89 = -t56 + t95;
t88 = -t118 + t93;
t5 = rSges(5,1) * t18 - rSges(5,2) * t19;
t85 = t70 * (-pkin(2) * t105 + t102) + t66 * t99 - t98;
t84 = t46 + t90;
t83 = t109 * t66 + t58 + t60;
t82 = -t5 + t56;
t78 = t112 * t66 + t58 + t86;
t7 = t120 + (-t36 - t110) * t70 + t90;
t8 = t115 + t89;
t77 = (t7 * t112 * t67 + (-t7 * qJ(3) + t8 * t112) * t66) * t70;
t13 = (-t36 - t37) * t70 + t90;
t14 = t92 * t70 + t89;
t76 = (t13 * t109 * t67 + (t13 * (-rSges(4,3) - qJ(3)) + t14 * t109) * t66) * t70;
t69 = t70 ^ 2;
t27 = t41 * t70 + t95;
t25 = t39 * t70 - t56;
t24 = -t29 * t70 - t97;
t23 = -t106 * t70 - t98;
t10 = (t116 * t70 - t25) * t70 + t91;
t9 = t70 * (-rSges(4,1) * t105 + t48) + t85;
t2 = -t70 * t25 - t5 * t65 - t69 * t63 + t91;
t1 = -t69 * t110 + t6 * t65 + t85;
t3 = [m(3) * (t24 * (-t38 - t111) + t23 * (t41 + t68) + (-t29 - t95 + t27) * t26) + (t2 * (t78 - t111) + t7 * (t82 - t95) + t1 * (t68 + t88) + t77 + (t84 + t7 + t96 + t119) * t8) * m(5) + (t10 * (t83 - t111) - t13 * t89 + t9 * (t68 + t92) + t76 + (t84 + t13 + t32 - t90 + t117) * t14) * m(4); (t1 * t88 + t2 * t78 + t77 + (t102 + t119) * t8 + (t115 - t56 + t82) * t7) * m(5) + (t10 * t83 + t76 + (-t101 + t102 + t117) * t14 + (t13 * t70 + t9) * t92) * m(4) + (-(-t26 * t41 - t27 * t38) * t70 + t23 * t41 - t24 * t38 - t26 * t29 - t27 * t106) * m(3); 0.2e1 * (t1 * t113 + t2 * t114) * m(5) + 0.2e1 * (t10 * t114 + t9 * t113) * m(4); (t1 * t118 - t2 * t86 + t7 * t5 - t8 * t6 - (-t7 * t118 - t8 * t86) * t65) * m(5);];
tauc  = t3(:);
