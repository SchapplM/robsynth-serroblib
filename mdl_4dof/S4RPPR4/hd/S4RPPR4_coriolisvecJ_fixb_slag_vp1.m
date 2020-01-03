% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:53
% DurationCPUTime: 3.19s
% Computational Cost: add. (2848->251), mult. (3515->347), div. (0->0), fcn. (2632->6), ass. (0->142)
t117 = qJ(1) + pkin(6);
t114 = sin(t117);
t115 = cos(t117);
t152 = -rSges(4,2) * t115 + rSges(4,3) * t114;
t121 = cos(qJ(1));
t116 = t121 * pkin(1);
t216 = pkin(2) * t115 + qJ(3) * t114;
t162 = t116 + t216;
t220 = t152 + t162;
t122 = qJD(1) ^ 2;
t118 = sin(qJ(4));
t120 = cos(qJ(4));
t146 = rSges(5,1) * t118 + rSges(5,2) * t120;
t182 = Icges(5,4) * t120;
t138 = Icges(5,1) * t118 + t182;
t87 = -Icges(5,2) * t118 + t182;
t194 = t87 + t138;
t183 = Icges(5,4) * t118;
t137 = Icges(5,2) * t120 + t183;
t89 = Icges(5,1) * t120 - t183;
t195 = -t137 + t89;
t219 = (t118 * t194 - t120 * t195) * qJD(1);
t218 = 2 * qJD(4);
t52 = -Icges(5,6) * t114 + t115 * t137;
t54 = -Icges(5,5) * t114 + t115 * t138;
t141 = t118 * t54 + t120 * t52;
t217 = t141 * t115;
t153 = rSges(3,1) * t115 - rSges(3,2) * t114;
t215 = t116 + t153;
t119 = sin(qJ(1));
t203 = pkin(1) * t119;
t214 = rSges(4,3) * t115 - t203;
t136 = Icges(5,5) * t118 + Icges(5,6) * t120;
t50 = -Icges(5,3) * t114 + t115 * t136;
t179 = qJD(1) * t50;
t24 = t118 * t52 - t120 * t54;
t51 = Icges(5,6) * t115 + t114 * t137;
t64 = t87 * t115;
t32 = qJD(1) * t51 - qJD(4) * t64;
t181 = Icges(5,5) * t115;
t66 = t89 * t115;
t34 = -qJD(4) * t66 + (t114 * t138 + t181) * qJD(1);
t212 = qJD(4) * t24 + t118 * t34 + t120 * t32 + t179;
t77 = t137 * qJD(4);
t78 = t138 * qJD(4);
t85 = Icges(5,5) * t120 - Icges(5,6) * t118;
t211 = qJD(1) * t85 + qJD(4) * (t118 * t87 - t120 * t89) + t118 * t78 + t120 * t77;
t178 = t114 * t118;
t177 = t114 * t120;
t94 = Icges(5,4) * t177;
t53 = Icges(5,1) * t178 + t181 + t94;
t142 = t118 * t51 - t120 * t53;
t49 = Icges(5,3) * t115 + t114 * t136;
t180 = qJD(1) * t49;
t171 = qJD(4) * t114;
t33 = qJD(1) * t52 + t171 * t87;
t65 = t89 * t114;
t35 = qJD(1) * t54 + qJD(4) * t65;
t210 = qJD(4) * t142 - t118 * t35 - t120 * t33 + t180;
t197 = t54 + t64;
t199 = t52 - t66;
t208 = t118 * t199 - t120 * t197;
t198 = -Icges(5,2) * t178 + t53 + t94;
t200 = t51 - t65;
t207 = t118 * t200 - t120 * t198;
t206 = t114 / 0.2e1;
t205 = -t115 / 0.2e1;
t202 = pkin(2) * t114;
t201 = -qJD(1) / 0.2e1;
t169 = qJD(1) * qJD(3);
t173 = qJD(1) * t114;
t102 = qJD(3) * t114;
t172 = qJD(1) * t115;
t184 = qJ(3) * t172 + t102;
t196 = qJD(1) * (-pkin(2) * t173 + t184) + t114 * t169;
t192 = rSges(5,1) * t120;
t189 = rSges(5,2) * t118;
t61 = t114 * t85;
t111 = t115 * rSges(5,3);
t186 = t115 * t85;
t105 = t115 * qJ(3);
t72 = -t105 + t202;
t185 = -qJD(1) * t72 + t102;
t140 = t118 * t89 + t120 * t87;
t37 = t115 * t140 - t61;
t176 = t37 * qJD(1);
t175 = t136 * qJD(1);
t174 = rSges(4,2) * t173 + rSges(4,3) * t172;
t170 = qJD(4) * t115;
t168 = -rSges(5,3) - pkin(2) - pkin(5);
t14 = t115 * t49 + t177 * t51 + t178 * t53;
t15 = -t115 * t50 - t177 * t52 - t178 * t54;
t167 = t146 * t172 + t171 * t192;
t166 = t122 * t203;
t165 = t122 * t116;
t55 = rSges(5,1) * t178 + rSges(5,2) * t177 + t111;
t164 = qJD(4) * t189;
t92 = -t189 + t192;
t69 = t92 * t171;
t163 = t92 * t170;
t158 = -t171 / 0.2e1;
t156 = -t170 / 0.2e1;
t150 = -pkin(5) * t114 - t203;
t149 = -pkin(5) * t115 - t116;
t74 = rSges(3,1) * t114 + rSges(3,2) * t115;
t143 = t118 * t53 + t120 * t51;
t109 = t114 * rSges(5,3);
t56 = t115 * t146 - t109;
t135 = t150 + t56;
t68 = t92 * t115;
t133 = (t114 * t15 + t115 * t14) * qJD(4);
t46 = t114 * t49;
t16 = -t115 * t143 + t46;
t17 = -t114 * t50 + t217;
t132 = (t114 * t17 + t115 * t16) * qJD(4);
t127 = -qJD(1) * t141 - qJD(4) * t186 + t180;
t126 = qJD(1) * t143 + qJD(4) * t61 + t179;
t125 = t140 * qJD(1) - qJD(4) * t136;
t103 = qJD(3) * t115;
t98 = t115 * t169;
t79 = t146 * qJD(4);
t67 = t92 * t114;
t58 = qJD(1) * t216 - t103;
t39 = (-rSges(5,3) * qJD(1) - t164) * t114 + t167;
t38 = -qJD(4) * t68 + (t114 * t146 + t111) * qJD(1);
t36 = t114 * t140 + t186;
t27 = t36 * qJD(1);
t26 = -t165 + t98 + (-qJD(1) * t152 - t58) * qJD(1);
t25 = qJD(1) * t174 - t166 + t196;
t22 = qJD(2) + (-t114 * t55 - t115 * t56) * qJD(4);
t21 = -t163 - t103 + (-t149 + t55 + t216) * qJD(1);
t20 = t102 + t69 + (t135 - t72) * qJD(1);
t13 = -t79 * t171 + t98 + t149 * t122 + (-t38 - t58 + t163) * qJD(1);
t12 = t79 * t170 + t150 * t122 + (t39 + t69) * qJD(1) + t196;
t11 = -t114 * t211 + t115 * t125;
t10 = t114 * t125 + t115 * t211;
t9 = qJD(4) * t141 - t118 * t32 + t120 * t34;
t8 = -qJD(4) * t143 - t118 * t33 + t120 * t35;
t7 = (-t114 * t39 + t115 * t38 + (t114 * t56 - t115 * t55) * qJD(1)) * qJD(4);
t6 = t132 - t176;
t5 = t27 + t133;
t1 = [(-qJD(4) * t140 + t118 * t77 - t120 * t78) * qJD(1) + (t27 + ((-t16 + t46 + t15) * t114 + (t17 - t217 + (-t143 + t50) * t114 + t14) * t115) * qJD(4)) * t158 + m(3) * ((-t122 * t74 - t166) * t215 + (-t165 + (-0.2e1 * t153 - t116 + t215) * t122) * (-t74 - t203)) + (t6 + t176 + (t114 ^ 2 * t50 + (-t46 + t15 + (t143 + t50) * t115) * t115) * qJD(4)) * t156 + (t13 * (-t109 + t150 - t72) + t20 * t103 + t12 * (t162 + t55) + t21 * (-t114 * t164 + t167 + t184) + (qJD(4) * t20 * t92 + t12 * pkin(5) + t13 * t146) * t115 + ((-t119 * t21 - t121 * t20) * pkin(1) + t20 * t168 * t115 + (t20 * (-qJ(3) - t146) + t21 * t168) * t114) * qJD(1) - (qJD(1) * t135 + t185 - t20 + t69) * t21) * m(5) + (t26 * (t105 + (rSges(4,2) - pkin(2)) * t114 + t214) + t25 * t220 + (t174 + t184 - t185 + (-rSges(4,2) * t114 - t202 - t203 - t214) * qJD(1)) * (qJD(1) * t220 - t103)) * m(4) + (t9 + t10 + t5) * t171 / 0.2e1 + (qJD(1) * t24 + t11 + t8) * t170 / 0.2e1 + (t115 * t37 + (-t142 + t36) * t114) * qJD(4) * t201; m(5) * t7; 0.2e1 * (t12 * t205 + t13 * t206) * m(5) + 0.2e1 * (t205 * t25 + t206 * t26) * m(4); qJD(1) * (t9 * t114 + t8 * t115 + (t114 * t142 + t24 * t115) * qJD(1)) / 0.2e1 + ((t170 * t61 - t175) * t115 + (-t219 + (t208 * t114 + (-t207 - t186) * t115) * qJD(4)) * t114) * t156 + ((-t171 * t186 - t175) * t114 + (t219 + (t207 * t115 + (-t208 + t61) * t114) * qJD(4)) * t115) * t158 + ((-t118 * t195 - t120 * t194) * qJD(1) + ((t114 * t199 - t115 * t200) * t120 + (t114 * t197 - t115 * t198) * t118) * qJD(4)) * t201 + (t10 * qJD(1) + ((t114 * t126 + t115 * t210) * t115 + t114 * (t114 * t127 - t115 * t212) + (-t16 * t114 + t17 * t115) * qJD(1)) * t218) * t206 + (t11 * qJD(1) + (t114 * (t114 * t212 + t115 * t127) + t115 * (-t114 * t210 + t115 * t126) + (-t14 * t114 + t15 * t115) * qJD(1)) * t218) * t115 / 0.2e1 - (t5 + t133) * t173 / 0.2e1 + (t6 + t132) * t172 / 0.2e1 + ((t21 * t79 - t7 * t56 + t22 * (-qJD(1) * t55 + t38) + (qJD(1) * t20 - t12) * t92) * t115 + (-t20 * t79 - t7 * t55 + t22 * (qJD(1) * t56 - t39) + (qJD(1) * t21 + t13) * t92) * t114 - (t20 * t68 + t21 * t67) * qJD(1) - (t22 * (-t114 * t67 - t115 * t68) - (t114 * t20 - t115 * t21) * t146) * qJD(4)) * m(5);];
tauc = t1(:);
