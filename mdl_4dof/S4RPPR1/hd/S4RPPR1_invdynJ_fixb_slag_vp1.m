% Calculate vector of inverse dynamics joint torques for
% S4RPPR1
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau = S4RPPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR1_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:44
% EndTime: 2018-11-14 13:46:45
% DurationCPUTime: 0.91s
% Computational Cost: add. (1095->129), mult. (1180->132), div. (0->0), fcn. (802->6), ass. (0->73)
t74 = cos(qJ(1));
t68 = t74 * pkin(1);
t73 = sin(qJ(1));
t104 = t73 * pkin(1);
t71 = qJ(1) + pkin(6);
t66 = sin(t71);
t67 = cos(t71);
t37 = rSges(3,1) * t66 + rSges(3,2) * t67;
t29 = -t37 - t104;
t64 = t67 * pkin(2);
t111 = t66 * qJ(3) + t64;
t57 = qJD(3) * t67;
t24 = qJD(1) * t111 - t57;
t69 = qJDD(1) - qJDD(4);
t70 = qJD(1) - qJD(4);
t102 = cos(qJ(4));
t72 = sin(qJ(4));
t78 = t67 * t102 + t66 * t72;
t16 = t70 * t78;
t88 = t66 * t102;
t32 = t67 * t72 - t88;
t92 = qJD(1) * t67;
t17 = -qJD(1) * t88 - qJD(4) * t32 + t72 * t92;
t7 = rSges(5,1) * t16 - rSges(5,2) * t17;
t75 = qJD(1) ^ 2;
t106 = pkin(2) * t66;
t59 = t67 * qJ(3);
t35 = -t59 + t106;
t79 = -t66 * pkin(3) - t104 - t35;
t105 = pkin(3) * t67;
t82 = t68 + t105;
t90 = qJD(1) * qJD(3);
t97 = qJDD(3) * t66 + t67 * t90;
t98 = t32 * rSges(5,1) + rSges(5,2) * t78;
t1 = -qJD(1) * t24 + t79 * qJDD(1) + t69 * t98 - t7 * t70 - t82 * t75 + t97;
t120 = -g(1) + t1;
t18 = -rSges(5,1) * t78 + rSges(5,2) * t32;
t85 = qJDD(1) * t68 - t75 * t104;
t93 = qJD(1) * t66;
t56 = qJD(3) * t66;
t96 = qJ(3) * t92 + t56;
t76 = -qJDD(3) * t67 + qJD(1) * (-pkin(2) * t93 + t96) + qJDD(1) * t111 + t66 * t90 + t85;
t8 = rSges(5,1) * t17 + rSges(5,2) * t16;
t2 = -t18 * t69 + t70 * t8 + (qJDD(1) * t67 - t66 * t75) * pkin(3) + t76;
t119 = -g(2) + t2;
t118 = t75 * t68;
t117 = qJD(1) * t35 - t56 + t96;
t116 = (-rSges(4,1) - pkin(2)) * t66;
t115 = -(t111 + t82) * qJD(1) + t57;
t114 = t70 * t98;
t62 = t67 * rSges(4,1);
t39 = t66 * rSges(4,3) + t62;
t89 = t68 + t111;
t22 = t89 + t39;
t40 = t67 * rSges(3,1) - rSges(3,2) * t66;
t30 = t40 + t68;
t109 = t66 / 0.2e1;
t108 = -t67 / 0.2e1;
t107 = -m(4) - m(5);
t94 = Icges(5,3) * t69;
t61 = t67 * rSges(4,3);
t87 = -rSges(4,1) * t66 - t104 + t61;
t86 = t59 - t104;
t84 = -t35 + t87;
t47 = rSges(2,1) * t74 - rSges(2,2) * t73;
t46 = rSges(2,1) * t73 + rSges(2,2) * t74;
t54 = rSges(4,3) * t92;
t14 = t84 * qJD(1) + t56;
t10 = -t70 * t18 - t115;
t9 = t79 * qJD(1) + t114 + t56;
t4 = qJDD(1) * t39 + qJD(1) * (-rSges(4,1) * t93 + t54) + t76;
t3 = -t118 + t84 * qJDD(1) + (-t39 * qJD(1) - t24) * qJD(1) + t97;
t5 = [-m(2) * (-g(1) * t46 + g(2) * t47) + t94 + ((-t37 * t75 - g(2) + t85) * t30 + (-g(1) - t118 + (-0.2e1 * t40 - t68 + t30) * t75) * t29) * m(3) + (Icges(2,3) + Icges(3,3) + Icges(4,2) + m(2) * (t46 ^ 2 + t47 ^ 2) + m(3) * (t29 ^ 2 + t40 * t30)) * qJDD(1) + ((t115 - t7) * t9 + t120 * ((-pkin(2) - pkin(3)) * t66 + t86 + t98) + t119 * (-t18 + t89 + t105) + (-t106 * qJD(1) - t114 + t117 + t8 + t9) * t10) * m(5) + ((-g(2) + t4) * t22 + (-g(1) + t3) * (t61 + t86 + t116) + (t57 + (-t62 - t68 - t64 + (-rSges(4,3) - qJ(3)) * t66) * qJD(1)) * t14 + (t14 + t54 + (-t104 - t87 + t116) * qJD(1) + t117) * (t22 * qJD(1) - t57)) * m(4); (-g(3) + qJDD(2)) * (m(3) - t107); t107 * (g(1) * t66 - g(2) * t67) + 0.2e1 * (t1 * t109 + t2 * t108) * m(5) + 0.2e1 * (t4 * t108 + t3 * t109) * m(4); -t94 + (-t10 * t8 + t9 * t7 - (-t10 * t70 + t120) * t98 + (t9 * t70 + t119) * t18) * m(5);];
tau  = t5;
