% Calculate vector of inverse dynamics joint torques for
% S4PPRR3
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:17:21
% EndTime: 2019-12-31 16:17:25
% DurationCPUTime: 2.62s
% Computational Cost: add. (2713->259), mult. (6385->368), div. (0->0), fcn. (6846->6), ass. (0->140)
t107 = sin(pkin(6));
t108 = cos(pkin(6));
t184 = sin(qJ(3));
t185 = cos(qJ(3));
t117 = t107 * t184 + t108 * t185;
t109 = sin(qJ(4));
t110 = cos(qJ(4));
t163 = Icges(5,4) * t110;
t126 = -Icges(5,2) * t109 + t163;
t87 = -t107 * t185 + t108 * t184;
t48 = Icges(5,6) * t87 + t117 * t126;
t164 = Icges(5,4) * t109;
t128 = Icges(5,1) * t110 - t164;
t51 = Icges(5,5) * t87 + t117 * t128;
t131 = t109 * t48 - t110 * t51;
t124 = Icges(5,5) * t110 - Icges(5,6) * t109;
t45 = Icges(5,3) * t87 + t117 * t124;
t17 = -t117 * t131 + t87 * t45;
t44 = Icges(5,3) * t117 - t124 * t87;
t202 = t87 * t44;
t125 = Icges(5,2) * t110 + t164;
t97 = Icges(5,1) * t109 + t163;
t129 = -t109 * t125 + t110 * t97;
t123 = Icges(5,5) * t109 + Icges(5,6) * t110;
t57 = t123 * t87;
t35 = -t117 * t129 - t57;
t200 = qJD(3) * t35;
t135 = t109 * rSges(5,1) + rSges(5,2) * t110;
t196 = qJD(4) * t135;
t198 = t117 * t196;
t132 = t109 * t51 + t110 * t48;
t69 = -rSges(4,1) * t117 + rSges(4,2) * t87;
t197 = qJD(3) * t69;
t58 = t123 * t117;
t195 = -pkin(3) * t117 - t87 * pkin(5);
t80 = t87 * qJD(3);
t81 = t117 * qJD(3);
t193 = -t81 * pkin(3) - t80 * pkin(5);
t165 = t110 * t87;
t168 = t109 * t87;
t53 = -rSges(5,1) * t165 + rSges(5,2) * t168 + rSges(5,3) * t117;
t31 = -t87 * pkin(3) + pkin(5) * t117 + t53;
t192 = qJD(3) * t31 - t198;
t50 = Icges(5,5) * t117 - t128 * t87;
t191 = t117 * t50 + t51 * t87;
t188 = m(3) + m(4);
t186 = -rSges(5,3) - pkin(5);
t82 = t87 * rSges(5,3);
t182 = -t117 * t44 + t50 * t165;
t166 = t110 * t117;
t181 = t50 * t166 + t202;
t171 = t109 * rSges(5,2);
t100 = -rSges(5,1) * t110 + t171;
t54 = -t100 * t117 + t82;
t178 = -t54 + t195;
t167 = t110 * t81;
t177 = -rSges(5,1) * t167 - t80 * rSges(5,3);
t176 = -t125 + t128;
t175 = -t126 - t97;
t47 = Icges(5,6) * t117 - t126 * t87;
t170 = t109 * t47;
t169 = t109 * t117;
t162 = qJD(4) * t117;
t161 = qJD(4) * t87;
t91 = t100 * qJD(4);
t160 = qJD(4) * t91;
t34 = t129 * t87 - t58;
t158 = t34 * qJD(3);
t157 = t124 * qJD(3);
t156 = qJD(2) * t108;
t155 = qJD(4) * t109;
t154 = qJD(4) * t110;
t153 = -m(5) - t188;
t152 = qJDD(2) * t108;
t149 = t87 * t155;
t146 = -t162 / 0.2e1;
t145 = -t161 / 0.2e1;
t144 = t161 / 0.2e1;
t52 = -rSges(5,1) * t166 + rSges(5,2) * t169 - t82;
t143 = t52 + t195;
t134 = t109 * t50 + t110 * t47;
t118 = t149 - t167;
t119 = t109 * t81 + t154 * t87;
t26 = Icges(5,4) * t118 + Icges(5,2) * t119 - Icges(5,6) * t80;
t28 = Icges(5,1) * t118 + Icges(5,4) * t119 - Icges(5,5) * t80;
t113 = qJD(4) * t134 + t109 * t26 - t110 * t28;
t120 = -t110 * t80 - t117 * t155;
t121 = t109 * t80 - t117 * t154;
t25 = Icges(5,4) * t120 + Icges(5,2) * t121 + Icges(5,6) * t81;
t27 = Icges(5,1) * t120 + Icges(5,4) * t121 + Icges(5,5) * t81;
t114 = qJD(4) * t132 + t109 * t25 - t110 * t27;
t133 = -t110 * t50 + t170;
t23 = Icges(5,5) * t120 + Icges(5,6) * t121 + Icges(5,3) * t81;
t24 = Icges(5,5) * t118 + Icges(5,6) * t119 - Icges(5,3) * t80;
t142 = (-t113 * t117 + t133 * t80 + t87 * t24 + t81 * t44) * t117 + (-t114 * t117 + t131 * t80 + t87 * t23 + t81 * t45) * t87;
t141 = (t113 * t87 + t117 * t24 + t133 * t81 - t80 * t44) * t117 + (t114 * t87 + t117 * t23 + t131 * t81 - t80 * t45) * t87;
t14 = t168 * t47 - t182;
t15 = t117 * t45 - t51 * t165 + t168 * t48;
t140 = t117 * t14 + t15 * t87;
t16 = -t169 * t47 + t181;
t139 = t117 * t16 + t17 * t87;
t106 = qJD(2) * t107;
t19 = t106 + t192;
t20 = qJD(3) * t178 + t87 * t196 - t156;
t138 = t117 * t19 - t20 * t87;
t29 = rSges(5,1) * t120 + rSges(5,2) * t121 + t81 * rSges(5,3);
t30 = rSges(5,1) * t149 + rSges(5,2) * t119 + t177;
t137 = -t117 * t29 + t30 * t87;
t136 = -t117 * t54 + t53 * t87;
t130 = -t109 * t97 - t110 * t125;
t122 = pkin(3) - t100;
t116 = t117 * t47 + t48 * t87;
t115 = (t109 * t175 + t110 * t176) * qJD(3);
t89 = t126 * qJD(4);
t90 = t128 * qJD(4);
t112 = qJD(4) * t130 - t109 * t89 + t110 * t90;
t111 = t109 * t191 + t116 * t110;
t105 = qJDD(2) * t107;
t88 = t124 * qJD(4);
t70 = -rSges(4,1) * t87 - rSges(4,2) * t117;
t68 = -qJD(4) * t80 + qJDD(4) * t117;
t67 = qJD(4) * t81 + qJDD(4) * t87;
t66 = t135 * t117;
t65 = t135 * t87;
t64 = -rSges(4,1) * t81 + rSges(4,2) * t80;
t63 = -rSges(4,1) * t80 - rSges(4,2) * t81;
t56 = -t156 + t197;
t33 = -qJD(3) * t63 + qJDD(3) * t69 - t152;
t32 = qJD(3) * t64 + qJDD(3) * t70 + t105;
t18 = qJD(4) * t136 + qJD(1);
t13 = t112 * t87 - t117 * t88 + t123 * t80 + t129 * t81;
t12 = -t112 * t117 - t123 * t81 + t129 * t80 - t87 * t88;
t11 = -t87 * t160 - t152 + t67 * t135 + t178 * qJDD(3) + (t80 * pkin(3) - t81 * pkin(5) - t29) * qJD(3);
t10 = t117 * t160 - t68 * t135 + t105 + t31 * qJDD(3) + (t193 + t30) * qJD(3);
t9 = qJD(4) * t131 - t109 * t27 - t110 * t25;
t8 = qJD(4) * t133 - t109 * t28 - t110 * t26;
t7 = qJD(4) * t137 + t53 * t67 - t54 * t68 + qJDD(1);
t6 = qJD(4) * t139 - t200;
t5 = qJD(4) * t140 - t158;
t1 = [m(5) * t7 + (m(2) + t188) * qJDD(1) + (-m(2) + t153) * g(3); t153 * (g(1) * t107 - g(2) * t108) + m(4) * (t107 * t32 - t108 * t33) + m(5) * (t10 * t107 - t108 * t11) + m(3) * (t107 ^ 2 + t108 ^ 2) * qJDD(2); t5 * t144 - qJD(3) * (-qJD(4) * t129 - t109 * t90 - t110 * t89) - (-t132 + t35) * t67 / 0.2e1 - (-t134 + t34) * t68 / 0.2e1 + (Icges(4,3) - t130) * qJDD(3) + (t9 + t12 - t158 + ((t15 - t16 + t181) * t87 - t182 * t117) * qJD(4)) * t145 + (-g(2) * t143 - t18 * (t52 + t54) * t161 + (-g(1) + t10) * t31 + (-t117 * t122 + t186 * t87) * t11 + (t122 * t80 + t186 * t81 + t192 + t198) * t20 + (-qJD(3) * t143 + t171 * t81 + t177 + t193) * t19) * m(5) + (-t56 * t63 + (t64 - t197) * (qJD(3) * t70 + t106) + (qJD(3) * t56 - g(1) + t32) * t70 + (-g(2) + t33) * t69) * m(4) + (t6 + t8 + t13 + ((-t14 + (-t45 + t170) * t87 - t182) * t87 - (-t15 - (t133 - t45) * t117 + t202) * t117) * qJD(4) + t200) * t146; t81 * t6 / 0.2e1 + t87 * (-qJD(3) * t12 + qJD(4) * t142 - qJDD(3) * t35 + t16 * t68 + t17 * t67) / 0.2e1 + t67 * t139 / 0.2e1 + (-t16 * t80 + t17 * t81 + t142) * t144 - t80 * t5 / 0.2e1 + t117 * (-qJD(3) * t13 + qJD(4) * t141 - qJDD(3) * t34 + t14 * t68 + t15 * t67) / 0.2e1 + t68 * t140 / 0.2e1 + (-t14 * t80 + t15 * t81 + t141) * t162 / 0.2e1 - qJDD(3) * (-t117 * t134 - t132 * t87) / 0.2e1 - qJD(3) * (t117 * t8 - t132 * t81 + t134 * t80 + t87 * t9) / 0.2e1 + ((-t161 * t58 + t157) * t87 - (-t115 + (-t87 * t57 + t111) * qJD(4)) * t117) * t145 + (-(-t162 * t57 - t157) * t117 + (-t115 + (-t117 * t58 + t111) * qJD(4)) * t87) * t146 + qJD(3) * (-(t176 * t109 - t175 * t110) * qJD(3) + (t116 * t109 - t110 * t191) * qJD(4)) / 0.2e1 + (t7 * t136 + t18 * (t53 * t81 + t54 * t80 + t137) + t138 * t91 - (t10 * t117 - t11 * t87 - t19 * t80 - t20 * t81) * t135 - (t19 * t65 + t20 * t66) * qJD(3) - (t18 * (t117 * t66 + t65 * t87) + t138 * t100) * qJD(4) + g(1) * t66 - g(2) * t65 - g(3) * t100) * m(5);];
tau = t1;
