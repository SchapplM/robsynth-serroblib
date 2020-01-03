% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:59
% DurationCPUTime: 2.03s
% Computational Cost: add. (2770->238), mult. (3341->339), div. (0->0), fcn. (2552->4), ass. (0->135)
t113 = sin(qJ(4));
t114 = cos(qJ(4));
t135 = rSges(5,1) * t113 + rSges(5,2) * t114;
t163 = Icges(5,4) * t114;
t129 = Icges(5,1) * t113 + t163;
t86 = -Icges(5,2) * t113 + t163;
t175 = t86 + t129;
t164 = Icges(5,4) * t113;
t128 = Icges(5,2) * t114 + t164;
t88 = Icges(5,1) * t114 - t164;
t176 = -t128 + t88;
t199 = (t113 * t175 - t114 * t176) * qJD(2);
t198 = 2 * qJD(4);
t112 = pkin(6) + qJ(2);
t111 = cos(t112);
t108 = t111 * rSges(5,3);
t110 = sin(t112);
t158 = t110 * t114;
t159 = t110 * t113;
t55 = rSges(5,1) * t159 + rSges(5,2) * t158 + t108;
t74 = t111 * pkin(2) + t110 * qJ(3);
t197 = t55 + t74;
t138 = -rSges(4,2) * t111 + t110 * rSges(4,3);
t196 = t138 + t74;
t52 = -Icges(5,6) * t110 + t111 * t128;
t54 = -Icges(5,5) * t110 + t111 * t129;
t132 = t113 * t54 + t114 * t52;
t195 = t132 * t111;
t127 = Icges(5,5) * t113 + Icges(5,6) * t114;
t50 = -Icges(5,3) * t110 + t111 * t127;
t160 = qJD(2) * t50;
t24 = t113 * t52 - t114 * t54;
t51 = Icges(5,6) * t111 + t110 * t128;
t62 = t86 * t111;
t30 = qJD(2) * t51 - qJD(4) * t62;
t162 = Icges(5,5) * t111;
t64 = t88 * t111;
t32 = -qJD(4) * t64 + (t110 * t129 + t162) * qJD(2);
t194 = qJD(4) * t24 + t113 * t32 + t114 * t30 + t160;
t77 = t128 * qJD(4);
t78 = t129 * qJD(4);
t84 = Icges(5,5) * t114 - Icges(5,6) * t113;
t193 = qJD(2) * t84 + (t113 * t86 - t114 * t88) * qJD(4) + t113 * t78 + t114 * t77;
t91 = Icges(5,4) * t158;
t53 = Icges(5,1) * t159 + t162 + t91;
t133 = t113 * t51 - t114 * t53;
t49 = Icges(5,3) * t111 + t110 * t127;
t161 = qJD(2) * t49;
t153 = qJD(4) * t110;
t31 = qJD(2) * t52 + t153 * t86;
t63 = t88 * t110;
t33 = qJD(2) * t54 + qJD(4) * t63;
t192 = qJD(4) * t133 - t113 * t33 - t114 * t31 + t161;
t178 = t54 + t62;
t180 = t52 - t64;
t190 = t113 * t180 - t114 * t178;
t179 = -Icges(5,2) * t159 + t53 + t91;
t181 = t51 - t63;
t189 = t113 * t181 - t114 * t179;
t188 = t110 / 0.2e1;
t187 = -t111 / 0.2e1;
t185 = rSges(4,2) - pkin(2);
t184 = pkin(5) * t110;
t183 = pkin(5) * qJD(2) ^ 2;
t182 = -qJD(2) / 0.2e1;
t151 = qJD(2) * qJD(3);
t155 = qJD(2) * t110;
t154 = qJD(2) * t111;
t99 = qJD(3) * t110;
t174 = qJ(3) * t154 + t99;
t177 = qJD(2) * (-pkin(2) * t155 + t174) + t110 * t151;
t173 = rSges(4,2) * t155 + rSges(4,3) * t154;
t102 = t111 * qJ(3);
t71 = pkin(2) * t110 - t102;
t172 = -qJD(2) * t71 + t99;
t170 = rSges(5,1) * t114;
t168 = rSges(5,2) * t113;
t166 = rSges(4,3) * t111;
t59 = t110 * t84;
t165 = t111 * t84;
t131 = t113 * t88 + t114 * t86;
t35 = t111 * t131 - t59;
t157 = t35 * qJD(2);
t156 = t127 * qJD(2);
t152 = qJD(4) * t111;
t150 = -rSges(5,3) - pkin(2) - pkin(5);
t14 = t111 * t49 + t158 * t51 + t159 * t53;
t15 = -t111 * t50 - t158 * t52 - t159 * t54;
t149 = t135 * t154 + t153 * t170;
t148 = qJD(4) * t168;
t90 = -t168 + t170;
t68 = t90 * t153;
t147 = t90 * t152;
t106 = t110 * rSges(5,3);
t56 = t111 * t135 - t106;
t146 = t56 - t184;
t143 = -t153 / 0.2e1;
t141 = -t152 / 0.2e1;
t134 = t113 * t53 + t114 * t51;
t66 = t90 * t111;
t125 = (t110 * t15 + t111 * t14) * qJD(4);
t46 = t110 * t49;
t16 = -t111 * t134 + t46;
t17 = -t110 * t50 + t195;
t124 = (t110 * t17 + t111 * t16) * qJD(4);
t120 = -qJD(2) * t132 - qJD(4) * t165 + t161;
t119 = qJD(2) * t134 + qJD(4) * t59 + t160;
t118 = t131 * qJD(2) - qJD(4) * t127;
t100 = qJD(3) * t111;
t95 = t111 * t151;
t79 = t135 * qJD(4);
t72 = rSges(4,2) * t110 + t166;
t65 = t90 * t110;
t58 = qJD(2) * t74 - t100;
t41 = qJD(2) * t196 - t100;
t40 = t99 + (-t71 + t72) * qJD(2);
t37 = (-rSges(5,3) * qJD(2) - t148) * t110 + t149;
t36 = -qJD(4) * t66 + (t110 * t135 + t108) * qJD(2);
t34 = t110 * t131 + t165;
t27 = t95 + (-qJD(2) * t138 - t58) * qJD(2);
t26 = qJD(2) * t173 + t177;
t25 = t34 * qJD(2);
t22 = qJD(1) + (-t110 * t55 - t111 * t56) * qJD(4);
t21 = -t147 - t100 + (pkin(5) * t111 + t197) * qJD(2);
t20 = t68 + t99 + (t146 - t71) * qJD(2);
t13 = -t111 * t183 - t79 * t153 + t95 + (-t36 - t58 + t147) * qJD(2);
t12 = -t110 * t183 + t79 * t152 + (t37 + t68) * qJD(2) + t177;
t11 = -t110 * t193 + t118 * t111;
t10 = t118 * t110 + t111 * t193;
t9 = qJD(4) * t132 - t113 * t30 + t114 * t32;
t8 = -qJD(4) * t134 - t113 * t31 + t114 * t33;
t7 = (-t110 * t37 + t111 * t36 + (t110 * t56 - t111 * t55) * qJD(2)) * qJD(4);
t6 = t124 - t157;
t5 = t25 + t125;
t1 = [m(5) * t7; (-qJD(4) * t131 + t113 * t77 - t114 * t78) * qJD(2) + (t25 + ((-t16 + t46 + t15) * t110 + (t17 - t195 + (-t134 + t50) * t110 + t14) * t111) * qJD(4)) * t143 + (t6 + t157 + (t110 ^ 2 * t50 + (-t46 + t15 + (t134 + t50) * t111) * t111) * qJD(4)) * t141 + (t13 * (-t106 - t71 - t184) + t20 * t100 + t12 * t197 + t21 * (-t110 * t148 + t149 + t174) + (qJD(4) * t20 * t90 + t12 * pkin(5) + t13 * t135) * t111 + (t20 * t150 * t111 + (t20 * (-qJ(3) - t135) + t21 * t150) * t110) * qJD(2) - (qJD(2) * t146 + t172 - t20 + t68) * t21) * m(5) + (t27 * (t110 * t185 + t102 + t166) + t40 * t100 + t26 * t196 + t41 * (t173 + t174) + (t40 * t185 * t111 + (t40 * (-rSges(4,3) - qJ(3)) - t41 * pkin(2)) * t110) * qJD(2) - (qJD(2) * t72 + t172 - t40) * t41) * m(4) + (t9 + t10 + t5) * t153 / 0.2e1 + (qJD(2) * t24 + t11 + t8) * t152 / 0.2e1 + (t111 * t35 + (-t133 + t34) * t110) * qJD(4) * t182; 0.2e1 * (t12 * t187 + t13 * t188) * m(5) + 0.2e1 * (t187 * t26 + t188 * t27) * m(4); qJD(2) * (t9 * t110 + t8 * t111 + (t110 * t133 + t24 * t111) * qJD(2)) / 0.2e1 + ((t59 * t152 - t156) * t111 + (-t199 + (t190 * t110 + (-t189 - t165) * t111) * qJD(4)) * t110) * t141 + ((-t153 * t165 - t156) * t110 + (t199 + (t189 * t111 + (-t190 + t59) * t110) * qJD(4)) * t111) * t143 + ((-t113 * t176 - t114 * t175) * qJD(2) + ((t110 * t180 - t111 * t181) * t114 + (t110 * t178 - t111 * t179) * t113) * qJD(4)) * t182 + (t10 * qJD(2) + ((t119 * t110 + t111 * t192) * t111 + t110 * (t120 * t110 - t111 * t194) + (-t16 * t110 + t17 * t111) * qJD(2)) * t198) * t188 + (t11 * qJD(2) + (t110 * (t110 * t194 + t120 * t111) + t111 * (-t110 * t192 + t119 * t111) + (-t14 * t110 + t15 * t111) * qJD(2)) * t198) * t111 / 0.2e1 - (t5 + t125) * t155 / 0.2e1 + (t6 + t124) * t154 / 0.2e1 + ((t21 * t79 - t7 * t56 + t22 * (-qJD(2) * t55 + t36) + (qJD(2) * t20 - t12) * t90) * t111 + (-t20 * t79 - t7 * t55 + t22 * (qJD(2) * t56 - t37) + (qJD(2) * t21 + t13) * t90) * t110 - (t20 * t66 + t21 * t65) * qJD(2) - (t22 * (-t110 * t65 - t111 * t66) - (t110 * t20 - t111 * t21) * t135) * qJD(4)) * m(5);];
tauc = t1(:);
