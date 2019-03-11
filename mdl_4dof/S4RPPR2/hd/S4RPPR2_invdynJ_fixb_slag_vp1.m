% Calculate vector of inverse dynamics joint torques for
% S4RPPR2
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:26
% EndTime: 2019-03-08 18:28:28
% DurationCPUTime: 0.98s
% Computational Cost: add. (1001->159), mult. (1541->174), div. (0->0), fcn. (1142->6), ass. (0->90)
t114 = qJD(1) * qJD(2);
t93 = sin(qJ(1));
t94 = cos(qJ(1));
t122 = qJDD(2) * t93 + t94 * t114;
t133 = qJD(1) ^ 2 * pkin(2);
t102 = -t133 * t94 + t122;
t135 = pkin(2) * t93;
t85 = t94 * qJ(2);
t56 = pkin(1) * t93 - t85;
t111 = -t56 - t135;
t92 = sin(pkin(6));
t129 = t92 * t94;
t118 = cos(pkin(6));
t74 = pkin(3) * t118 + pkin(2);
t107 = pkin(3) * t129 - t93 * t74;
t30 = t107 + t135;
t105 = t111 + t30;
t113 = pkin(6) + qJ(4);
t108 = cos(t113);
t101 = t93 * t108;
t80 = sin(t113);
t139 = -t80 * t94 + t101;
t97 = t108 * t94 + t93 * t80;
t127 = -rSges(5,1) * t139 + rSges(5,2) * t97;
t128 = t93 * t92;
t123 = pkin(3) * t128 + t94 * t74;
t134 = pkin(2) * t94;
t31 = t123 - t134;
t59 = t94 * pkin(1) + t93 * qJ(2);
t83 = qJD(2) * t94;
t41 = t59 * qJD(1) - t83;
t91 = qJD(1) - qJD(4);
t19 = t91 * t97;
t116 = qJD(1) * t94;
t20 = -qJD(1) * t101 + t139 * qJD(4) + t116 * t80;
t9 = rSges(5,1) * t19 - rSges(5,2) * t20;
t90 = qJDD(1) - qJDD(4);
t1 = t90 * t127 - t91 * t9 + t105 * qJDD(1) + (-t31 * qJD(1) - t41) * qJD(1) + t102;
t146 = -g(1) + t1;
t10 = rSges(5,1) * t20 + t19 * rSges(5,2);
t117 = qJD(1) * t93;
t21 = -rSges(5,1) * t97 - rSges(5,2) * t139;
t112 = t92 * t116;
t65 = pkin(3) * t112;
t82 = qJD(2) * t93;
t121 = qJ(2) * t116 + t82;
t100 = -qJDD(2) * t94 + qJD(1) * (-pkin(1) * t117 + t121) + qJDD(1) * t59 + t93 * t114;
t96 = qJDD(1) * t134 - t133 * t93 + t100;
t2 = qJDD(1) * t31 + t91 * t10 - t90 * t21 + (t65 + (pkin(2) - t74) * t117) * qJD(1) + t96;
t145 = -g(2) + t2;
t144 = t91 * t127;
t141 = m(4) + m(5);
t109 = t93 * t118;
t49 = -t109 + t129;
t99 = t118 * t94 + t128;
t125 = t49 * rSges(4,1) + rSges(4,2) * t99;
t60 = t94 * rSges(3,1) + t93 * rSges(3,3);
t110 = t59 + t134;
t27 = rSges(4,1) * t99 - t49 * rSges(4,2);
t140 = t110 + t27;
t138 = t93 / 0.2e1;
t137 = -t94 / 0.2e1;
t136 = -pkin(1) - pkin(2);
t132 = -rSges(3,1) - pkin(1);
t131 = -pkin(1) - t74;
t42 = t99 * qJD(1);
t43 = -qJD(1) * t109 + t112;
t126 = t43 * rSges(4,1) + t42 * rSges(4,2);
t87 = t94 * rSges(3,3);
t57 = rSges(3,1) * t93 - t87;
t124 = -t56 - t57;
t36 = t59 + t60;
t120 = -qJD(1) * t56 + t82;
t119 = Icges(5,3) * t90;
t106 = t111 + t125;
t61 = rSges(2,1) * t94 - rSges(2,2) * t93;
t58 = rSges(2,1) * t93 + rSges(2,2) * t94;
t103 = -rSges(4,1) * t42 + rSges(4,2) * t43;
t79 = rSges(3,3) * t116;
t29 = qJD(1) * t36 - t83;
t28 = qJD(1) * t124 + t82;
t17 = t140 * qJD(1) - t83;
t16 = qJD(1) * t106 + t82;
t12 = qJDD(1) * t60 + qJD(1) * (-rSges(3,1) * t117 + t79) + t100;
t11 = t124 * qJDD(1) + (-t60 * qJD(1) - t41) * qJD(1) + t122;
t8 = -t91 * t21 - t83 + (t110 + t31) * qJD(1);
t7 = qJD(1) * t105 + t144 + t82;
t4 = qJD(1) * t126 + qJDD(1) * t27 + t96;
t3 = (t103 - t41) * qJD(1) + t106 * qJDD(1) + t102;
t5 = [-m(2) * (-g(1) * t58 + g(2) * t61) + t119 + (-(-t7 + t144 + (t30 - t135) * qJD(1) + t120) * t8 + t7 * (t83 - t9) + t8 * (t10 + t65 + t121) + (t7 * t131 * t94 + (t7 * (-pkin(3) * t92 - qJ(2)) + t8 * t131) * t93) * qJD(1) + t145 * (-t21 + t59 + t123) + t146 * (t107 - t56 + t127)) * m(5) + (-(-t16 + (t125 - t135) * qJD(1) + t120) * t17 + t16 * (t103 + t83) + t17 * (t121 + t126) + (t16 * t136 * t94 + (-t16 * qJ(2) + t136 * t17) * t93) * qJD(1) + (-g(2) + t4) * t140 + (-g(1) + t3) * (t136 * t93 + t125 + t85)) * m(4) + (-(-qJD(1) * t57 + t120 - t28) * t29 + t28 * t83 + t29 * (t79 + t121) + (t28 * t132 * t94 + (t28 * (-rSges(3,3) - qJ(2)) + t29 * t132) * t93) * qJD(1) + (-g(2) + t12) * t36 + (-g(1) + t11) * (t132 * t93 + t85 + t87)) * m(3) + (Icges(2,3) + Icges(3,2) + Icges(4,3) + m(2) * (t58 ^ 2 + t61 ^ 2)) * qJDD(1); (-m(3) - t141) * (g(1) * t93 - g(2) * t94) + 0.2e1 * (t1 * t138 + t2 * t137) * m(5) + 0.2e1 * (t4 * t137 + t138 * t3) * m(4) + 0.2e1 * (t11 * t138 + t12 * t137) * m(3); t141 * (g(3) + qJDD(3)); -t119 + (-t8 * t10 + t7 * t9 - (-t91 * t8 + t146) * t127 + (t7 * t91 + t145) * t21) * m(5);];
tau  = t5;
