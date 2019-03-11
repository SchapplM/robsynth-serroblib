% Calculate vector of inverse dynamics joint torques for
% S4RRRP1
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:35:53
% EndTime: 2019-03-08 18:35:54
% DurationCPUTime: 0.72s
% Computational Cost: add. (1483->128), mult. (1058->137), div. (0->0), fcn. (504->6), ass. (0->79)
t67 = sin(qJ(1));
t68 = cos(qJ(1));
t69 = qJD(1) ^ 2;
t113 = (-qJDD(1) * t67 - t68 * t69) * pkin(1) - g(1);
t62 = t68 * pkin(1);
t99 = pkin(1) * t67;
t112 = qJDD(1) * t62 - t69 * t99 - g(2);
t66 = qJ(1) + qJ(2);
t59 = sin(t66);
t60 = cos(t66);
t65 = qJD(1) + qJD(2);
t63 = t65 ^ 2;
t64 = qJDD(1) + qJDD(2);
t111 = (-t59 * t64 - t60 * t63) * pkin(2) + t113;
t51 = pkin(2) * t60;
t98 = pkin(2) * t59;
t110 = t64 * t51 - t63 * t98 + t112;
t87 = pkin(1) * qJD(1);
t83 = t67 * t87;
t35 = rSges(3,1) * t59 + rSges(3,2) * t60;
t93 = t35 * t65;
t19 = -t83 - t93;
t61 = qJ(3) + t66;
t53 = cos(t61);
t58 = qJD(3) + t65;
t91 = t53 * t58;
t52 = sin(t61);
t92 = t52 * t58;
t18 = rSges(4,1) * t91 - rSges(4,2) * t92;
t32 = rSges(4,1) * t52 + rSges(4,2) * t53;
t57 = qJDD(3) + t64;
t109 = -t18 * t58 - t32 * t57 + t111;
t89 = t60 * t65;
t90 = t59 * t65;
t28 = rSges(3,1) * t89 - rSges(3,2) * t90;
t108 = -t28 * t65 - t35 * t64 + t113;
t26 = t58 * t32;
t34 = t53 * rSges(4,1) - rSges(4,2) * t52;
t107 = -t26 * t58 + t34 * t57 + t110;
t36 = t60 * rSges(3,1) - rSges(3,2) * t59;
t106 = t36 * t64 - t65 * t93 + t112;
t95 = rSges(5,2) * t53;
t31 = rSges(5,1) * t52 + t95;
t37 = rSges(5,2) * t92;
t56 = t58 ^ 2;
t105 = -t57 * t31 - t58 * (rSges(5,1) * t91 - t37) + (-t52 * t57 - t53 * t56) * pkin(3) + t111;
t33 = t53 * rSges(5,1) - rSges(5,2) * t52;
t104 = t57 * t33 - t31 * t56 + (-t52 * t56 + t53 * t57) * pkin(3) + t110;
t101 = pkin(3) * t53 + t33;
t103 = t58 * t101;
t86 = pkin(2) * t90;
t77 = -t83 - t86;
t100 = t77 + t83;
t97 = pkin(3) * t52;
t96 = -rSges(5,1) - pkin(3);
t94 = t34 * t58;
t88 = (Icges(4,3) + Icges(5,3)) * t57;
t85 = pkin(2) * t89;
t84 = Icges(3,3) * t64 + t88;
t82 = t68 * t87;
t81 = t96 * t52;
t24 = t34 + t51;
t44 = rSges(2,1) * t68 - rSges(2,2) * t67;
t43 = rSges(2,1) * t67 + rSges(2,2) * t68;
t14 = t101 + t51;
t21 = t81 - t95;
t76 = t82 + t85;
t23 = -t32 - t98;
t25 = t58 * t31;
t75 = -pkin(3) * t92 - t25 - t86;
t74 = -t18 - t85;
t13 = t21 - t98;
t7 = (-t31 - t97) * t58 + t77;
t8 = t76 + t103;
t71 = (t8 * t81 + (-t8 * rSges(5,2) + t7 * t96) * t53) * t58;
t20 = t36 * t65 + t82;
t10 = t76 + t94;
t9 = t77 - t26;
t1 = [Icges(2,3) * qJDD(1) + t84 + (t106 * (t36 + t62) + t108 * (-t35 - t99) + (-t28 - t82 + t20) * t19) * m(3) + (g(1) * t43 - g(2) * t44 + (t43 ^ 2 + t44 ^ 2) * qJDD(1)) * m(2) + (t7 * (t37 - t76) + t71 + (t7 - t75 + t100) * t8 + t104 * (t14 + t62) + t105 * (t13 - t99)) * m(5) + (t9 * (t74 - t82) + t107 * (t24 + t62) + t109 * (t23 - t99) + (t9 - t100 - t86) * t10) * m(4); t84 + (t71 + (-t75 - t86) * t8 + t104 * t14 + t105 * t13 + (t37 + t103) * t7) * m(5) + ((t85 + t94 + t74) * t9 + t107 * t24 + t109 * t23) * m(4) + (-t19 * t28 - t20 * t93 + (t19 * t65 + t106) * t36 + (t20 * t65 - t108) * t35) * m(3); t88 + (t7 * t37 + t71 + t8 * t25 - (-t101 * t7 - t8 * t97) * t58 + t104 * t101 + t105 * t21) * m(5) + (-t10 * t26 - t9 * t18 + (t58 * t9 + t107) * t34 + (t10 * t58 - t109) * t32) * m(4); (-g(3) + qJDD(4)) * m(5);];
tau  = t1;
