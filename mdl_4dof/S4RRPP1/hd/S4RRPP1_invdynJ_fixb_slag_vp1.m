% Calculate vector of inverse dynamics joint torques for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:32:58
% EndTime: 2019-03-08 18:32:59
% DurationCPUTime: 0.77s
% Computational Cost: add. (1393->125), mult. (1017->131), div. (0->0), fcn. (518->6), ass. (0->74)
t71 = qJ(1) + qJ(2);
t65 = sin(t71);
t104 = pkin(2) * t65;
t102 = -rSges(5,1) - pkin(3);
t64 = pkin(6) + t71;
t58 = cos(t64);
t50 = t58 * qJ(4);
t57 = sin(t64);
t117 = t58 * rSges(5,3) + t102 * t57 + t50;
t122 = -t104 + t117;
t72 = sin(qJ(1));
t105 = pkin(1) * t72;
t73 = cos(qJ(1));
t67 = t73 * pkin(1);
t74 = qJD(1) ^ 2;
t121 = qJDD(1) * t67 - t105 * t74 - g(2);
t66 = cos(t71);
t59 = pkin(2) * t66;
t70 = qJD(1) + qJD(2);
t68 = t70 ^ 2;
t69 = qJDD(1) + qJDD(2);
t119 = t104 * t68 - t69 * t59 - t121;
t118 = (-qJDD(1) * t72 - t73 * t74) * pkin(1) - g(1);
t116 = rSges(5,3) + qJ(4);
t36 = rSges(3,1) * t65 + rSges(3,2) * t66;
t101 = t36 * t70;
t95 = pkin(1) * qJD(1);
t88 = t72 * t95;
t19 = -t88 - t101;
t100 = t57 * t70;
t48 = qJD(4) * t58;
t80 = t102 * t70 + qJD(4);
t97 = t66 * t68;
t115 = -pkin(2) * t97 + qJDD(4) * t57 + t122 * t69 + (-t100 * t116 + t58 * t80 + t48) * t70 + t118;
t32 = rSges(4,1) * t57 + rSges(4,2) * t58;
t41 = rSges(4,2) * t100;
t99 = t58 * t70;
t114 = -t69 * t32 - t70 * (rSges(4,1) * t99 - t41) + (-t65 * t69 - t97) * pkin(2) + t118;
t56 = t66 * rSges(3,1);
t98 = t65 * t70;
t27 = -rSges(3,2) * t98 + t56 * t70;
t113 = -t27 * t70 - t36 * t69 + t118;
t54 = t58 * rSges(4,1);
t35 = -rSges(4,2) * t57 + t54;
t112 = -t32 * t68 + t69 * t35 - t119;
t108 = rSges(5,3) * t99 + t70 * t50;
t47 = qJD(4) * t57;
t90 = t47 + t108;
t107 = t102 * t58;
t96 = t116 * t57 - t107;
t111 = qJDD(4) * t58 - t96 * t69 - (t57 * t80 + t90) * t70 + t119;
t37 = -rSges(3,2) * t65 + t56;
t110 = -t101 * t70 + t37 * t69 + t121;
t109 = t35 + t59;
t106 = t59 + t96;
t87 = t73 * t95;
t81 = t48 - t87;
t8 = t106 * t70 - t81;
t92 = t8 * t104;
t91 = t117 * t70 + t47;
t89 = (Icges(5,2) + Icges(3,3) + Icges(4,3)) * t69;
t21 = -t32 - t104;
t82 = t47 - t88;
t46 = rSges(2,1) * t73 - rSges(2,2) * t72;
t45 = rSges(2,1) * t72 + rSges(2,2) * t73;
t78 = -pkin(2) * t98 - t88;
t15 = t21 * t70 - t88;
t16 = t109 * t70 + t87;
t76 = (t15 * (-t54 - t59) + t16 * t21) * t70;
t7 = t122 * t70 + t82;
t75 = (t7 * (-t59 + t107) - t92 + (t102 * t8 - t116 * t7) * t57) * t70;
t25 = t70 * t32;
t20 = t37 * t70 + t87;
t1 = [Icges(2,3) * qJDD(1) + t89 + (t110 * (t37 + t67) + t113 * (-t36 - t105) + (-t27 - t87 + t20) * t19) * m(3) + (g(1) * t45 - g(2) * t46 + (t45 ^ 2 + t46 ^ 2) * qJDD(1)) * m(2) + (t7 * t81 + t75 + (t7 - t78 - t91 + t82 + t108) * t8 - t111 * (t67 + t106) + t115 * (t122 - t105)) * m(5) + (t15 * (t41 - t87) + t76 + t112 * (t109 + t67) + t114 * (t21 - t105) + (t15 + t25 - t78 - t88) * t16) * m(4); t89 + (t16 * t25 - (-t104 * t16 - t109 * t15) * t70 + t15 * t41 + t76 + t112 * t109 + t114 * t21) * m(4) + (-t19 * t27 - t20 * t101 + (t19 * t70 + t110) * t37 + (t20 * t70 - t113) * t36) * m(3) + (-(-t106 * t7 - t92) * t70 + t75 + (-t91 + t90) * t8 - t111 * t106 + t115 * t122) * m(5); (m(4) + m(5)) * (-g(3) + qJDD(3)); (t111 * t58 + t115 * t57) * m(5);];
tau  = t1;
