% Calculate vector of inverse dynamics joint torques for
% S4RRPP2
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:33:59
% EndTime: 2019-03-08 18:34:00
% DurationCPUTime: 0.74s
% Computational Cost: add. (1389->137), mult. (1196->150), div. (0->0), fcn. (636->4), ass. (0->83)
t78 = qJ(1) + qJ(2);
t72 = sin(t78);
t73 = cos(t78);
t38 = rSges(3,1) * t72 + rSges(3,2) * t73;
t77 = qJD(1) + qJD(2);
t105 = t38 * t77;
t79 = sin(qJ(1));
t99 = pkin(1) * qJD(1);
t95 = t79 * t99;
t25 = -t95 - t105;
t39 = t73 * pkin(2) + t72 * qJ(3);
t56 = qJD(3) * t73;
t22 = t39 * t77 - t56;
t40 = t73 * rSges(5,1) + t72 * rSges(5,2);
t66 = t73 * pkin(3);
t75 = t77 ^ 2;
t76 = qJDD(1) + qJDD(2);
t80 = cos(qJ(1));
t81 = qJD(1) ^ 2;
t86 = (-qJDD(1) * t79 - t80 * t81) * pkin(1);
t98 = qJD(3) * t77;
t85 = qJDD(3) * t72 + t73 * t98 + t86;
t108 = pkin(3) * t72;
t58 = t73 * qJ(3);
t35 = pkin(2) * t72 - t58;
t62 = t73 * rSges(5,2);
t36 = rSges(5,1) * t72 - t62;
t92 = -t35 - t36 - t108;
t1 = -t75 * t66 + (-t40 * t77 - t22) * t77 + t92 * t76 + t85;
t123 = -g(1) + t1;
t61 = t73 * rSges(4,3);
t37 = rSges(4,1) * t72 - t61;
t102 = -t35 - t37;
t41 = t73 * rSges(4,1) + t72 * rSges(4,3);
t3 = t102 * t76 + (-t41 * t77 - t22) * t77 + t85;
t122 = -g(1) + t3;
t104 = t72 * t77;
t103 = t73 * t77;
t51 = rSges(5,2) * t103;
t48 = t77 * t58;
t55 = qJD(3) * t72;
t101 = t48 + t55;
t109 = pkin(1) * t79;
t74 = t80 * pkin(1);
t93 = qJDD(1) * t74 - t81 * t109;
t84 = -qJDD(3) * t73 + t77 * (-pkin(2) * t104 + t101) + t76 * t39 + t72 * t98 + t93;
t2 = t76 * t40 + t77 * (-rSges(5,1) * t104 + t51) + (-t72 * t75 + t73 * t76) * pkin(3) + t84;
t121 = -g(2) + t2;
t50 = rSges(4,3) * t103;
t4 = t76 * t41 + t77 * (-rSges(4,1) * t104 + t50) + t84;
t120 = -g(2) + t4;
t28 = rSges(3,1) * t103 - rSges(3,2) * t104;
t119 = -t28 * t77 - t38 * t76 - g(1) + t86;
t42 = t73 * rSges(3,1) - rSges(3,2) * t72;
t118 = -t105 * t77 + t42 * t76 - g(2) + t93;
t112 = t66 + t39 + t40;
t117 = t112 * t77;
t116 = t77 * t36 + t51;
t115 = t77 * t37 + t50;
t32 = t77 * t35;
t114 = t48 + t32;
t113 = t101 - t55 + t32;
t111 = t72 / 0.2e1;
t110 = -t73 / 0.2e1;
t106 = -rSges(4,1) - pkin(2);
t24 = t39 + t41;
t97 = -rSges(5,1) - pkin(2) - pkin(3);
t96 = (Icges(4,2) + Icges(3,3) + Icges(5,3)) * t76;
t94 = t80 * t99;
t90 = t55 - t95;
t89 = -t56 + t94;
t53 = rSges(2,1) * t80 - rSges(2,2) * t79;
t52 = rSges(2,1) * t79 + rSges(2,2) * t80;
t23 = t106 * t72 + t58 + t61;
t17 = t97 * t72 + t58 + t62;
t13 = t102 * t77 + t90;
t14 = t24 * t77 + t89;
t83 = (t13 * t106 * t73 + (t13 * (-rSges(4,3) - qJ(3)) + t14 * t106) * t72) * t77;
t7 = t92 * t77 + t90;
t8 = t89 + t117;
t82 = (t7 * t97 * t73 + (t7 * (-rSges(5,2) - qJ(3)) + t8 * t97) * t72) * t77;
t26 = t42 * t77 + t94;
t5 = [Icges(2,3) * qJDD(1) + t96 + (t118 * (t42 + t74) + t119 * (-t38 - t109) + (-t28 - t94 + t26) * t25) * m(3) + (g(1) * t52 - g(2) * t53 + (t52 ^ 2 + t53 ^ 2) * qJDD(1)) * m(2) + (-t7 * t89 + t82 + (pkin(3) * t104 + t114 + t116 + t7) * t8 + t121 * (t74 + t112) + t123 * (t17 - t109)) * m(5) + (-t13 * t89 + t83 + t120 * (t74 + t24) + t122 * (t23 - t109) + (t13 + t114 + t115) * t14) * m(4); t96 + (t83 + t122 * t23 + (t113 + t115) * t14 + (t13 * t77 + t120) * t24) * m(4) + (-t25 * t28 - t26 * t105 + (t25 * t77 + t118) * t42 + (t26 * t77 - t119) * t38) * m(3) + (t7 * t117 + t82 + (t108 * t77 + t113 + t116) * t8 + t121 * t112 + t123 * t17) * m(5); (-m(4) - m(5)) * (g(1) * t72 - g(2) * t73) + 0.2e1 * (t1 * t111 + t2 * t110) * m(5) + 0.2e1 * (t4 * t110 + t3 * t111) * m(4); (g(3) + qJDD(4)) * m(5);];
tau  = t5;
