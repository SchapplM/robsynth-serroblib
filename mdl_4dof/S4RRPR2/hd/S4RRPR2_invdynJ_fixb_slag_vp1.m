% Calculate vector of inverse dynamics joint torques for
% S4RRPR2
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
% Datum: 2019-05-28 15:34
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-28 15:34:02
% EndTime: 2019-05-28 15:34:04
% DurationCPUTime: 1.08s
% Computational Cost: add. (1968->159), mult. (1704->168), div. (0->0), fcn. (1158->6), ass. (0->92)
t94 = sin(qJ(1));
t95 = cos(qJ(1));
t96 = qJD(1) ^ 2;
t102 = (-qJDD(1) * t94 - t95 * t96) * pkin(1);
t91 = qJD(1) + qJD(2);
t117 = qJD(3) * t91;
t92 = qJ(1) + qJ(2);
t86 = sin(t92);
t87 = cos(t92);
t101 = qJDD(3) * t86 + t87 * t117 + t102;
t128 = cos(qJ(4));
t93 = sin(qJ(4));
t104 = t87 * t128 + t86 * t93;
t135 = t86 * t128 - t87 * t93;
t122 = -rSges(5,1) * t135 + rSges(5,2) * t104;
t144 = qJD(4) - t91;
t51 = t87 * pkin(2) + t86 * qJ(3);
t72 = qJD(3) * t87;
t30 = t51 * t91 - t72;
t74 = t87 * qJ(3);
t48 = pkin(2) * t86 - t74;
t22 = t144 * t104;
t23 = t144 * t135;
t7 = -rSges(5,1) * t22 - rSges(5,2) * t23;
t90 = qJDD(1) + qJDD(2);
t84 = -qJDD(4) + t90;
t89 = t91 ^ 2;
t1 = t122 * t84 - t91 * t30 - t90 * t48 + t7 * t144 + (-t86 * t90 - t87 * t89) * pkin(3) + t101;
t148 = -g(1) + t1;
t138 = -rSges(5,1) * t104 - rSges(5,2) * t135;
t8 = t23 * rSges(5,1) - t22 * rSges(5,2);
t130 = pkin(1) * t94;
t88 = t95 * pkin(1);
t111 = qJDD(1) * t88 - t96 * t130;
t60 = t91 * t74;
t71 = qJD(3) * t86;
t120 = t60 + t71;
t125 = t86 * t91;
t99 = -qJDD(3) * t87 + t111 + t91 * (-pkin(2) * t125 + t120) + t90 * t51 + t86 * t117;
t2 = -t138 * t84 - t8 * t144 + (-t86 * t89 + t87 * t90) * pkin(3) + t99;
t142 = -g(2) + t2;
t146 = t144 * t122;
t131 = -pkin(2) - pkin(3);
t145 = t131 * t86 + t122 + t74;
t118 = pkin(1) * qJD(1);
t115 = t94 * t118;
t50 = rSges(3,1) * t86 + rSges(3,2) * t87;
t126 = t50 * t91;
t33 = -t115 - t126;
t76 = t87 * rSges(4,3);
t49 = rSges(4,1) * t86 - t76;
t121 = -t48 - t49;
t52 = t87 * rSges(4,1) + t86 * rSges(4,3);
t3 = t121 * t90 + (-t52 * t91 - t30) * t91 + t101;
t143 = -g(1) + t3;
t124 = t87 * t91;
t62 = rSges(4,3) * t124;
t4 = t90 * t52 + t91 * (-rSges(4,1) * t125 + t62) + t99;
t141 = -g(2) + t4;
t36 = rSges(3,1) * t124 - rSges(3,2) * t125;
t140 = -t36 * t91 - t50 * t90 - g(1) + t102;
t53 = t87 * rSges(3,1) - rSges(3,2) * t86;
t139 = -t126 * t91 + t53 * t90 - g(2) + t111;
t137 = t91 * t49 + t62;
t112 = t87 * pkin(3) + t51;
t44 = t91 * t48;
t119 = t71 - t44;
t136 = pkin(3) * t125 - t119 + t146 + t8;
t134 = t112 * t91 + t138 * t144;
t133 = t86 / 0.2e1;
t132 = -t87 / 0.2e1;
t129 = -rSges(4,1) - pkin(2);
t32 = t51 + t52;
t69 = Icges(5,3) * t84;
t116 = t69 + (Icges(4,2) + Icges(3,3)) * t90;
t114 = t95 * t118;
t110 = t71 - t115;
t109 = -t72 + t114;
t14 = -t138 + t112;
t64 = rSges(2,1) * t95 - rSges(2,2) * t94;
t63 = rSges(2,1) * t94 + rSges(2,2) * t95;
t107 = t110 + t60;
t31 = t129 * t86 + t74 + t76;
t106 = -t7 + t72;
t10 = t109 + t134;
t9 = -t146 + (-pkin(3) * t86 - t48) * t91 + t110;
t98 = (t9 * t131 * t87 + (-t9 * qJ(3) + t10 * t131) * t86) * t91;
t17 = t121 * t91 + t110;
t18 = t32 * t91 + t109;
t97 = (t17 * t129 * t87 + (t17 * (-rSges(4,3) - qJ(3)) + t18 * t129) * t86) * t91;
t34 = t53 * t91 + t114;
t5 = [Icges(2,3) * qJDD(1) + t116 + (t139 * (t53 + t88) + t140 * (-t50 - t130) + (-t36 - t114 + t34) * t33) * m(3) + (g(1) * t63 - g(2) * t64 + (t63 ^ 2 + t64 ^ 2) * qJDD(1)) * m(2) + (t9 * (t106 - t114) + t98 + t142 * (t88 + t14) + (t9 + t115 + t107 + t136) * t10 + t148 * (t145 - t130)) * m(5) + (-t17 * t109 + t97 + t141 * (t88 + t32) + t143 * (t31 - t130) + (-t110 + t17 + t44 + t107 + t137) * t18) * m(4); t116 + (t98 + (t106 + t134 - t72) * t9 + t142 * t14 + (t120 + t136) * t10 + t148 * t145) * m(5) + (t97 + t143 * t31 + (-t119 + t120 + t137) * t18 + (t17 * t91 + t141) * t32) * m(4) + (-t33 * t36 - t34 * t126 + (t33 * t91 + t139) * t53 + (t34 * t91 - t140) * t50) * m(3); (-m(4) - m(5)) * (g(1) * t86 - g(2) * t87) + 0.2e1 * (t1 * t133 + t2 * t132) * m(5) + 0.2e1 * (t4 * t132 + t3 * t133) * m(4); -t69 + (-t10 * t8 + t9 * t7 - (t10 * t144 + t148) * t122 + (-t144 * t9 + t142) * t138) * m(5);];
tau  = t5;
