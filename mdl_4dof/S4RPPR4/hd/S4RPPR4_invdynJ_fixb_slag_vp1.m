% Calculate vector of inverse dynamics joint torques for
% S4RPPR4
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:49
% EndTime: 2019-12-31 16:38:53
% DurationCPUTime: 3.64s
% Computational Cost: add. (3173->306), mult. (3818->398), div. (0->0), fcn. (2856->6), ass. (0->161)
t136 = cos(qJ(1));
t131 = t136 * pkin(1);
t134 = sin(qJ(1));
t216 = pkin(1) * t134;
t132 = qJ(1) + pkin(6);
t129 = sin(t132);
t130 = cos(t132);
t83 = rSges(3,1) * t129 + rSges(3,2) * t130;
t73 = -t83 - t216;
t137 = qJD(1) ^ 2;
t233 = t137 * t131;
t229 = t130 * pkin(2) + t129 * qJ(3);
t172 = t131 + t229;
t85 = -rSges(4,2) * t130 + t129 * rSges(4,3);
t52 = t85 + t172;
t133 = sin(qJ(4));
t135 = cos(qJ(4));
t158 = rSges(5,1) * t133 + rSges(5,2) * t135;
t195 = Icges(5,4) * t135;
t100 = -Icges(5,2) * t133 + t195;
t150 = Icges(5,1) * t133 + t195;
t186 = t100 + t150;
t196 = Icges(5,4) * t133;
t102 = Icges(5,1) * t135 - t196;
t149 = Icges(5,2) * t135 + t196;
t198 = t102 - t149;
t232 = (t133 * t186 - t135 * t198) * qJD(1);
t86 = t130 * rSges(3,1) - rSges(3,2) * t129;
t74 = t131 + t86;
t56 = -Icges(5,6) * t129 + t130 * t149;
t58 = -Icges(5,5) * t129 + t130 * t150;
t151 = t133 * t58 + t135 * t56;
t231 = t151 * t130;
t178 = qJD(1) * qJD(3);
t182 = qJD(1) * t129;
t116 = qJD(3) * t129;
t181 = qJD(1) * t130;
t184 = qJ(3) * t181 + t116;
t228 = t129 * t178 + qJD(1) * (-pkin(2) * t182 + t184) + qJDD(1) * t229;
t148 = Icges(5,5) * t133 + Icges(5,6) * t135;
t54 = -Icges(5,3) * t129 + t130 * t148;
t192 = qJD(1) * t54;
t26 = t133 * t56 - t135 * t58;
t55 = Icges(5,6) * t130 + t129 * t149;
t68 = t100 * t130;
t32 = qJD(1) * t55 - qJD(4) * t68;
t194 = Icges(5,5) * t130;
t70 = t102 * t130;
t34 = -qJD(4) * t70 + (t129 * t150 + t194) * qJD(1);
t227 = t26 * qJD(4) + t133 * t34 + t135 * t32 + t192;
t146 = t133 * t100 - t102 * t135;
t89 = t149 * qJD(4);
t90 = t150 * qJD(4);
t98 = Icges(5,5) * t135 - Icges(5,6) * t133;
t226 = qJD(1) * t98 + qJD(4) * t146 + t133 * t90 + t135 * t89;
t190 = t129 * t135;
t107 = Icges(5,4) * t190;
t191 = t129 * t133;
t57 = Icges(5,1) * t191 + t107 + t194;
t152 = t133 * t55 - t135 * t57;
t53 = Icges(5,3) * t130 + t129 * t148;
t193 = qJD(1) * t53;
t179 = t129 * qJD(4);
t33 = qJD(1) * t56 + t100 * t179;
t69 = t102 * t129;
t35 = qJD(1) * t58 + qJD(4) * t69;
t225 = qJD(4) * t152 - t133 * t35 - t135 * t33 + t193;
t208 = t58 + t68;
t210 = t56 - t70;
t223 = t133 * t210 - t135 * t208;
t209 = -Icges(5,2) * t191 + t107 + t57;
t211 = t55 - t69;
t222 = t133 * t211 - t135 * t209;
t177 = qJD(1) * qJD(4);
t76 = qJDD(4) * t129 + t130 * t177;
t221 = t76 / 0.2e1;
t77 = qJDD(4) * t130 - t129 * t177;
t220 = t77 / 0.2e1;
t219 = -m(4) - m(5);
t218 = -pkin(2) - pkin(5);
t217 = t129 / 0.2e1;
t215 = pkin(2) * t129;
t214 = pkin(5) * t130;
t203 = rSges(5,2) * t133;
t206 = rSges(5,1) * t135;
t105 = -t203 + t206;
t163 = -pkin(5) * t129 - t216;
t60 = -t129 * rSges(5,3) + t130 * t158;
t145 = t163 + t60;
t119 = t130 * qJ(3);
t81 = -t119 + t215;
t144 = t145 - t81;
t162 = -t131 - t214;
t185 = qJDD(3) * t129 + t130 * t178;
t125 = t130 * rSges(5,3);
t72 = t105 * t130;
t38 = -qJD(4) * t72 + (t129 * t158 + t125) * qJD(1);
t117 = qJD(3) * t130;
t62 = qJD(1) * t229 - t117;
t91 = t158 * qJD(4);
t7 = -t91 * t179 + t105 * t76 + t162 * t137 + (-t38 - t62) * qJD(1) + t144 * qJDD(1) + t185;
t213 = t7 * t129;
t128 = qJDD(1) * t131;
t173 = qJD(4) * t203;
t174 = qJD(4) * t206;
t175 = t129 * t174 + t158 * t181;
t39 = (-rSges(5,3) * qJD(1) - t173) * t129 + t175;
t59 = rSges(5,1) * t191 + rSges(5,2) * t190 + t125;
t8 = qJD(1) * t39 + qJDD(1) * t59 - t105 * t77 + t128 + t163 * t137 + (pkin(5) * qJDD(1) + qJD(4) * t91 - qJDD(3)) * t130 + t228;
t212 = t8 * t130;
t201 = rSges(4,3) * t130;
t65 = t129 * t98;
t200 = t130 * t98;
t75 = t105 * t179;
t22 = qJD(1) * t144 + t116 + t75;
t199 = t22 * t130;
t197 = -qJD(1) * t81 + t116;
t147 = t135 * t100 + t133 * t102;
t37 = t130 * t147 - t65;
t188 = t37 * qJD(1);
t187 = t148 * qJD(1);
t183 = rSges(4,2) * t182 + rSges(4,3) * t181;
t180 = qJD(4) * t130;
t176 = -rSges(5,3) + t218;
t16 = t130 * t53 + t55 * t190 + t57 * t191;
t17 = -t130 * t54 - t56 * t190 - t58 * t191;
t171 = rSges(4,2) * t129 + t201 - t216;
t170 = -t180 / 0.2e1;
t169 = t180 / 0.2e1;
t168 = -t179 / 0.2e1;
t167 = t179 / 0.2e1;
t166 = t119 - t216;
t164 = -t137 * t216 + t128;
t153 = t133 * t57 + t135 * t55;
t139 = qJD(1) * t153 + qJD(4) * t65 + t192;
t140 = -qJD(1) * t151 - qJD(4) * t200 + t193;
t161 = (t139 * t129 + t130 * t225) * t130 + t129 * (t140 * t129 - t130 * t227);
t160 = t129 * (t129 * t227 + t140 * t130) + t130 * (-t129 * t225 + t139 * t130);
t106 = rSges(2,1) * t136 - rSges(2,2) * t134;
t104 = rSges(2,1) * t134 + rSges(2,2) * t136;
t155 = t129 * t17 + t130 * t16;
t48 = t129 * t53;
t18 = -t153 * t130 + t48;
t19 = -t129 * t54 + t231;
t154 = t129 * t19 + t130 * t18;
t138 = t147 * qJD(1) - t148 * qJD(4);
t71 = t105 * t129;
t36 = t129 * t147 + t200;
t27 = t36 * qJD(1);
t24 = qJD(2) + (-t129 * t59 - t130 * t60) * qJD(4);
t23 = -t105 * t180 - t117 + (-t162 + t59 + t229) * qJD(1);
t15 = qJD(1) * t183 + qJDD(1) * t85 - qJDD(3) * t130 + t164 + t228;
t14 = -t233 + (t171 - t81) * qJDD(1) + (-qJD(1) * t85 - t62) * qJD(1) + t185;
t13 = -t129 * t226 + t138 * t130;
t12 = t138 * t129 + t130 * t226;
t11 = t151 * qJD(4) - t133 * t32 + t135 * t34;
t10 = -qJD(4) * t153 - t133 * t33 + t135 * t35;
t9 = -t59 * t76 - t60 * t77 + qJDD(2) + (-t129 * t39 + t130 * t38) * qJD(4);
t6 = qJD(4) * t154 - t188;
t5 = qJD(4) * t155 + t27;
t1 = [(t27 + ((-t18 + t48 + t17) * t129 + (t19 - t231 + (-t153 + t54) * t129 + t16) * t130) * qJD(4)) * t168 + (-t147 * qJD(4) + t133 * t89 - t135 * t90) * qJD(1) - m(2) * (-g(1) * t104 + g(2) * t106) - t76 * t37 / 0.2e1 + t26 * t221 + (-t152 + t36) * t220 + (t188 + (t129 ^ 2 * t54 + (-t48 + t17 + (t153 + t54) * t130) * t130) * qJD(4) + t6) * t170 + (t10 + t13) * t169 + ((-t137 * t83 - g(2) + t164) * t74 + (-t233 + (-0.2e1 * t86 - t131 + t74) * t137 - g(1)) * t73) * m(3) + (t11 + t12 + t5) * t167 + (-(qJD(1) * t145 + t197 - t22 + t75) * t23 + t22 * (t117 + (-t173 + t174) * t130) + t23 * (-t129 * t173 + t175 + t184) + ((-t23 * t134 - t22 * t136) * pkin(1) + t176 * t199 + (t22 * (-qJ(3) - t158) + t23 * t176) * t129) * qJD(1) + (t8 - g(2)) * (t172 + t59 + t214) + (t7 - g(1)) * (t129 * t218 + t166 + t60)) * m(5) + ((t14 - g(1)) * (t201 + (rSges(4,2) - pkin(2)) * t129 + t166) + (t15 - g(2)) * t52 + (t183 + t184 - t197 + (-t215 - t216 - t171) * qJD(1)) * (t52 * qJD(1) - t117)) * m(4) + (m(3) * (t73 ^ 2 + t86 * t74) + m(2) * (t104 ^ 2 + t106 ^ 2) - t146 + Icges(2,3) + Icges(3,3) + Icges(4,1)) * qJDD(1); m(5) * t9 + (m(3) + m(4)) * qJDD(2) + (-m(3) + t219) * g(3); t219 * (g(1) * t129 - g(2) * t130) + 0.2e1 * (t213 / 0.2e1 - t212 / 0.2e1) * m(5) + 0.2e1 * (t14 * t217 - t15 * t130 / 0.2e1) * m(4); -t5 * t182 / 0.2e1 + t130 * (t13 * qJD(1) + t160 * qJD(4) + t36 * qJDD(1) + t16 * t77 + t17 * t76) / 0.2e1 + t155 * t220 + ((-t16 * t129 + t17 * t130) * qJD(1) + t160) * t169 + t6 * t181 / 0.2e1 + (t12 * qJD(1) + qJD(4) * t161 - t37 * qJDD(1) + t18 * t77 + t19 * t76) * t217 + t154 * t221 + ((-t18 * t129 + t19 * t130) * qJD(1) + t161) * t167 + qJDD(1) * (t26 * t129 - t130 * t152) / 0.2e1 + qJD(1) * (t10 * t130 + t11 * t129 + (t129 * t152 + t26 * t130) * qJD(1)) / 0.2e1 + ((t65 * t180 - t187) * t130 + (-t232 + (t223 * t129 + (-t222 - t200) * t130) * qJD(4)) * t129) * t170 + ((-t179 * t200 - t187) * t129 + (t232 + (t222 * t130 + (-t223 + t65) * t129) * qJD(4)) * t130) * t168 - qJD(1) * ((-t198 * t133 - t186 * t135) * qJD(1) + ((t129 * t210 - t130 * t211) * t135 + (t129 * t208 - t130 * t209) * t133) * qJD(4)) / 0.2e1 + ((t23 * t91 - t9 * t60 + t24 * (-qJD(1) * t59 + t38)) * t130 + (-t22 * t91 - t9 * t59 + t24 * (qJD(1) * t60 - t39)) * t129 + (t213 - t212 + (t23 * t129 + t199) * qJD(1)) * t105 - (t22 * t72 + t23 * t71) * qJD(1) - (t24 * (-t129 * t71 - t130 * t72) - (t22 * t129 - t130 * t23) * t158) * qJD(4) - g(1) * t71 + g(2) * t72 + g(3) * t158) * m(5);];
tau = t1;
