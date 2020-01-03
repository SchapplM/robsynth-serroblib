% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR6
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:42
% DurationCPUTime: 3.62s
% Computational Cost: add. (2732->337), mult. (7410->455), div. (0->0), fcn. (7192->6), ass. (0->161)
t155 = sin(qJ(1));
t153 = cos(pkin(6));
t154 = sin(qJ(4));
t152 = sin(pkin(6));
t231 = cos(qJ(4));
t195 = t152 * t231;
t118 = -t153 * t154 + t195;
t156 = cos(qJ(1));
t108 = t118 * t156;
t168 = t152 * t154 + t153 * t231;
t109 = t168 * t156;
t215 = t153 * t155;
t106 = t154 * t215 - t155 * t195;
t107 = t168 * t155;
t99 = Icges(5,4) * t107;
t57 = Icges(5,2) * t106 - Icges(5,6) * t156 - t99;
t98 = Icges(5,4) * t106;
t59 = Icges(5,1) * t107 + Icges(5,5) * t156 - t98;
t228 = -t108 * t57 + t109 * t59;
t53 = Icges(5,5) * t107 - Icges(5,6) * t106 + Icges(5,3) * t156;
t16 = -t155 * t53 + t228;
t249 = t106 * t57 + t107 * t59;
t248 = t118 * t59 + t168 * t57;
t230 = pkin(5) * t156;
t120 = pkin(3) * t215 + t230;
t62 = rSges(5,1) * t107 - rSges(5,2) * t106 + rSges(5,3) * t156;
t246 = -t120 - t62;
t14 = t156 * t53 + t249;
t218 = Icges(5,4) * t109;
t58 = Icges(5,2) * t108 - Icges(5,6) * t155 + t218;
t100 = Icges(5,4) * t108;
t61 = Icges(5,1) * t109 - Icges(5,5) * t155 + t100;
t227 = t108 * t58 + t109 * t61;
t55 = Icges(5,5) * t109 + Icges(5,6) * t108 - Icges(5,3) * t155;
t182 = t155 * t55 - t227;
t245 = t182 - t14;
t89 = Icges(5,5) * t118 - Icges(5,6) * t168;
t217 = Icges(5,4) * t118;
t91 = -Icges(5,2) * t168 + t217;
t116 = Icges(5,4) * t168;
t93 = Icges(5,1) * t118 - t116;
t164 = t106 * t91 - t107 * t93 - t156 * t89;
t244 = t164 * qJD(1);
t243 = 2 * qJD(4);
t242 = -t106 * t58 + t107 * t61;
t124 = t156 * pkin(1) + t155 * qJ(2);
t214 = t153 * t156;
t216 = t152 * t156;
t170 = rSges(3,1) * t214 - rSges(3,2) * t216 + t155 * rSges(3,3);
t241 = t124 + t170;
t135 = pkin(3) * t214;
t183 = pkin(2) * t214 + qJ(3) * t216 + t124;
t240 = t135 + t183;
t196 = rSges(4,1) * t214 + t155 * rSges(4,2) + rSges(4,3) * t216;
t239 = t183 + t196;
t238 = qJD(4) * t118;
t237 = t155 * (-Icges(5,2) * t109 + t100 + t61) - t156 * (-Icges(5,2) * t107 + t59 - t98);
t235 = t155 / 0.2e1;
t234 = -t156 / 0.2e1;
t232 = -rSges(5,3) - pkin(5);
t229 = -qJD(1) / 0.2e1;
t112 = t168 * qJD(4);
t205 = qJD(1) * t155;
t69 = -t112 * t156 - t118 * t205;
t70 = t156 * t238 - t168 * t205;
t224 = t70 * rSges(5,1) + t69 * rSges(5,2);
t223 = -Icges(5,2) * t118 - t116 + t93;
t222 = -Icges(5,1) * t168 - t217 - t91;
t221 = rSges(3,2) * t152;
t220 = -rSges(4,3) - qJ(3);
t144 = qJD(2) * t156;
t111 = qJD(1) * t124 - t144;
t175 = pkin(2) * t153 + qJ(3) * t152;
t203 = qJD(3) * t152;
t194 = t155 * t203;
t204 = qJD(1) * t156;
t219 = -t175 * t204 - t111 - t194;
t213 = t109 * rSges(5,1) + t108 * rSges(5,2);
t200 = qJD(1) * qJD(2);
t143 = qJD(2) * t155;
t207 = qJ(2) * t204 + t143;
t212 = qJD(1) * (-pkin(1) * t205 + t207) + t155 * t200;
t114 = t175 * t155;
t146 = t156 * qJ(2);
t122 = pkin(1) * t155 - t146;
t211 = -t114 - t122;
t131 = t155 * t221;
t210 = rSges(3,3) * t204 + qJD(1) * t131;
t129 = t156 * t203;
t209 = t129 + t143;
t208 = t156 * rSges(3,3) + t131;
t206 = -qJD(1) * t122 + t143;
t202 = qJD(4) * t155;
t201 = qJD(4) * t156;
t15 = t156 * t55 + t242;
t199 = rSges(3,1) * t215;
t95 = rSges(5,1) * t118 - rSges(5,2) * t168;
t87 = t95 * t201;
t197 = t129 + t207;
t193 = -rSges(3,1) * t153 - pkin(1);
t189 = t202 / 0.2e1;
t188 = -t201 / 0.2e1;
t185 = t212 + (-t175 * t205 + 0.2e1 * t129) * qJD(1);
t184 = -qJD(1) * t114 + t129 + t206;
t179 = -pkin(1) + (-rSges(4,1) - pkin(2)) * t153;
t71 = qJD(1) * t108 - t112 * t155;
t72 = qJD(1) * t109 + t155 * t238;
t178 = -rSges(5,1) * t72 - rSges(5,2) * t71;
t176 = rSges(4,1) * t153 + rSges(4,3) * t152;
t174 = t155 * (Icges(5,5) * t108 - Icges(5,6) * t109) - t156 * (-Icges(5,5) * t106 - Icges(5,6) * t107);
t171 = qJD(4) * t95 + t203;
t169 = -pkin(1) - t175;
t166 = (t14 * t156 - t15 * t155) * qJD(4);
t165 = (t155 * t182 + t156 * t16) * qJD(4);
t163 = -(Icges(5,1) * t108 - t218 - t58) * t155 + (-Icges(5,1) * t106 + t57 - t99) * t156;
t161 = -pkin(3) * t153 + t169;
t150 = t156 * rSges(4,2);
t142 = rSges(4,2) * t204;
t139 = t156 * t200;
t104 = t199 - t208;
t103 = t155 * t176 - t150;
t88 = -Icges(5,5) * t168 - Icges(5,6) * t118;
t86 = -rSges(5,1) * t112 - rSges(5,2) * t238;
t85 = -Icges(5,1) * t112 - Icges(5,4) * t238;
t84 = -Icges(5,4) * t112 - Icges(5,2) * t238;
t83 = -Icges(5,5) * t112 - Icges(5,6) * t238;
t82 = qJD(1) * t241 - t144;
t81 = t143 + (-t104 - t122) * qJD(1);
t80 = rSges(5,1) * t108 - rSges(5,2) * t109;
t79 = -rSges(5,1) * t106 - rSges(5,2) * t107;
t64 = -rSges(5,3) * t155 + t213;
t51 = t139 + (-qJD(1) * t170 - t111) * qJD(1);
t50 = qJD(1) * (-qJD(1) * t199 + t210) + t212;
t49 = qJD(1) * t239 - t144 + t194;
t48 = (-t103 + t211) * qJD(1) + t209;
t38 = t139 + (-qJD(1) * t196 - t194 + t219) * qJD(1);
t37 = (-t176 * t205 + t142) * qJD(1) + t185;
t36 = -rSges(5,3) * t205 - t178;
t35 = -rSges(5,3) * t204 + t224;
t34 = Icges(5,1) * t72 + Icges(5,4) * t71 - Icges(5,5) * t205;
t33 = Icges(5,1) * t70 + Icges(5,4) * t69 - Icges(5,5) * t204;
t32 = Icges(5,4) * t72 + Icges(5,2) * t71 - Icges(5,6) * t205;
t31 = Icges(5,4) * t70 + Icges(5,2) * t69 - Icges(5,6) * t204;
t30 = Icges(5,5) * t72 + Icges(5,6) * t71 - Icges(5,3) * t205;
t29 = Icges(5,5) * t70 + Icges(5,6) * t69 - Icges(5,3) * t204;
t28 = -qJD(3) * t153 + (-t155 * t62 - t156 * t64) * qJD(4);
t27 = t118 * t61 - t168 * t58;
t25 = t108 * t91 + t109 * t93 - t155 * t89;
t23 = t25 * qJD(1);
t22 = -t144 + t171 * t155 + (-pkin(5) * t155 + t240 + t64) * qJD(1);
t21 = t87 + (t211 + t246) * qJD(1) + t209;
t13 = t86 * t201 + t139 + (-qJD(1) * t135 - t36 + (pkin(5) * qJD(1) - t171) * t155 + t219) * qJD(1);
t12 = t86 * t202 + (-t120 * qJD(1) + t35 + t87) * qJD(1) + t185;
t11 = (-t155 * t36 - t156 * t35 + (t155 * t64 - t156 * t62) * qJD(1)) * qJD(4);
t10 = -t106 * t84 + t107 * t85 + t156 * t83 - t205 * t89 + t71 * t91 + t72 * t93;
t9 = t108 * t84 + t109 * t85 - t155 * t83 - t204 * t89 + t69 * t91 + t70 * t93;
t8 = -t112 * t61 + t118 * t33 - t168 * t31 - t238 * t58;
t7 = -t112 * t59 + t118 * t34 - t168 * t32 + t238 * t57;
t6 = t23 + t165;
t5 = t166 - t244;
t1 = [(-t112 * t93 + t118 * t85 - t168 * t84 - t238 * t91) * qJD(1) + (t23 + (t228 * t156 + (t245 + t249) * t155) * qJD(4)) * t188 - (t8 + t9) * t202 / 0.2e1 + (t5 + t244 + ((t227 + t245) * t156 + t242 * t155) * qJD(4)) * t189 + (t13 * (t146 - t62 - t230) + t21 * (t144 + t178) + t12 * (t213 + t240) + t22 * (t197 + t224) + (t12 * t232 + t13 * t161 - t21 * t203) * t155 + ((t161 * t21 + t22 * t232) * t156 + (t21 * (-qJ(2) - t232) + t22 * t161) * t155) * qJD(1) - (qJD(1) * t246 + t184 - t21 + t87) * t22) * m(5) + (t38 * (t146 + t150) + t48 * t144 + t37 * t239 + t49 * (t142 + t197) + (t38 * t179 + (-t48 * qJD(3) + t220 * t38) * t152) * t155 + (t48 * (t152 * t220 + t179) * t156 + (t48 * (-rSges(4,2) - qJ(2)) + t49 * (t169 - t176)) * t155) * qJD(1) - (-qJD(1) * t103 + t184 - t48) * t49) * m(4) + (t51 * (t155 * t193 + t146 + t208) + t81 * t144 + t50 * t241 + t82 * (t207 + t210) + (t81 * (t193 + t221) * t156 + (t81 * (-rSges(3,3) - qJ(2)) + t82 * t193) * t155) * qJD(1) - (-qJD(1) * t104 + t206 - t81) * t82) * m(3) + (t7 + t10 + t6) * t201 / 0.2e1 + ((t27 + t25) * t156 + (-t164 + t248) * t155) * qJD(4) * t229; 0.2e1 * (t12 * t234 + t13 * t235) * m(5) + 0.2e1 * (t234 * t37 + t235 * t38) * m(4) + 0.2e1 * (t234 * t50 + t235 * t51) * m(3); -t11 * t153 * m(5) + 0.2e1 * (m(4) * (t155 * t37 + t156 * t38) / 0.2e1 + m(5) * (t12 * t155 + t13 * t156) / 0.2e1) * t152; qJD(1) * (-t155 * t8 + t156 * t7 + (-t248 * t155 - t156 * t27) * qJD(1)) / 0.2e1 + ((t108 * t223 + t109 * t222 - t155 * t88) * qJD(1) + (-t108 * t237 + t109 * t163 + t155 * t174) * qJD(4)) * t189 + ((-t106 * t223 + t107 * t222 + t156 * t88) * qJD(1) + (t106 * t237 + t163 * t107 - t174 * t156) * qJD(4)) * t188 + ((t118 * t222 - t168 * t223) * qJD(1) + (t118 * t163 + t168 * t237) * qJD(4)) * t229 - (qJD(1) * t9 + ((t108 * t32 + t109 * t34 - t155 * t30 - t204 * t53 - t57 * t69 + t59 * t70) * t156 - t155 * (t108 * t31 + t109 * t33 - t155 * t29 - t204 * t55 + t58 * t69 + t61 * t70) + (-t16 * t155 + t156 * t182) * qJD(1)) * t243) * t155 / 0.2e1 + (qJD(1) * t10 + (-t155 * (-t106 * t31 + t107 * t33 + t156 * t29 - t205 * t55 + t58 * t71 + t61 * t72) + t156 * (-t106 * t32 + t107 * t34 + t156 * t30 - t205 * t53 - t57 * t71 + t59 * t72) + (-t14 * t155 - t15 * t156) * qJD(1)) * t243) * t156 / 0.2e1 - (t5 + t166) * t205 / 0.2e1 - (t6 + t165) * t204 / 0.2e1 + ((t21 * t86 - t11 * t64 + t28 * (-qJD(1) * t62 - t35) + (qJD(1) * t22 + t13) * t95) * t156 + (t22 * t86 - t11 * t62 + t28 * (qJD(1) * t64 - t36) + (-qJD(1) * t21 + t12) * t95) * t155 - (-t21 * t79 + t22 * t80) * qJD(1) - (t28 * (-t155 * t79 - t156 * t80) + (t22 * t155 + t21 * t156) * (-rSges(5,1) * t168 - rSges(5,2) * t118)) * qJD(4)) * m(5);];
tauc = t1(:);
