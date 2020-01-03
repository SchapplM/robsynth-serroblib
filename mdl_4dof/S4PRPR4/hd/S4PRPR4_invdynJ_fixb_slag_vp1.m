% Calculate vector of inverse dynamics joint torques for
% S4PRPR4
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:59
% DurationCPUTime: 2.27s
% Computational Cost: add. (3074->282), mult. (3602->381), div. (0->0), fcn. (2750->4), ass. (0->149)
t124 = sin(qJ(4));
t125 = cos(qJ(4));
t143 = rSges(5,1) * t124 + rSges(5,2) * t125;
t174 = Icges(5,4) * t125;
t135 = Icges(5,1) * t124 + t174;
t96 = -Icges(5,2) * t124 + t174;
t185 = t96 + t135;
t175 = Icges(5,4) * t124;
t134 = Icges(5,2) * t125 + t175;
t98 = Icges(5,1) * t125 - t175;
t186 = -t134 + t98;
t212 = (t124 * t185 - t125 * t186) * qJD(2);
t123 = pkin(6) + qJ(2);
t122 = cos(t123);
t121 = sin(t123);
t158 = qJD(2) * qJD(3);
t163 = qJD(2) * t121;
t110 = qJD(3) * t121;
t162 = qJD(2) * t122;
t165 = qJ(3) * t162 + t110;
t82 = t122 * pkin(2) + t121 * qJ(3);
t154 = t121 * t158 + qJD(2) * (-pkin(2) * t163 + t165) + qJDD(2) * t82;
t164 = rSges(4,2) * t163 + rSges(4,3) * t162;
t83 = -rSges(4,2) * t122 + t121 * rSges(4,3);
t211 = qJD(2) * t164 + qJDD(2) * t83 - qJDD(3) * t122 - g(2) + t154;
t210 = g(1) * t121;
t56 = -Icges(5,6) * t121 + t122 * t134;
t58 = -Icges(5,5) * t121 + t122 * t135;
t138 = t124 * t58 + t125 * t56;
t209 = t138 * t122;
t181 = rSges(5,2) * t124;
t183 = rSges(5,1) * t125;
t100 = -t181 + t183;
t60 = -t121 * rSges(5,3) + t122 * t143;
t151 = -pkin(5) * t121 + t60;
t113 = t122 * qJ(3);
t79 = pkin(2) * t121 - t113;
t146 = t151 - t79;
t161 = qJD(4) * t121;
t166 = qJDD(3) * t121 + t122 * t158;
t194 = pkin(5) * qJD(2) ^ 2;
t119 = t122 * rSges(5,3);
t70 = t100 * t122;
t36 = -qJD(4) * t70 + (t121 * t143 + t119) * qJD(2);
t111 = qJD(3) * t122;
t62 = t82 * qJD(2) - t111;
t157 = qJD(2) * qJD(4);
t74 = qJDD(4) * t121 + t122 * t157;
t89 = t143 * qJD(4);
t8 = -t122 * t194 - t89 * t161 + t100 * t74 + (-t36 - t62) * qJD(2) + t146 * qJDD(2) + t166;
t152 = qJD(4) * t181;
t153 = qJD(4) * t183;
t155 = t121 * t153 + t143 * t162;
t37 = (-rSges(5,3) * qJD(2) - t152) * t121 + t155;
t169 = t121 * t125;
t170 = t121 * t124;
t59 = rSges(5,1) * t170 + rSges(5,2) * t169 + t119;
t75 = qJDD(4) * t122 - t121 * t157;
t9 = -t121 * t194 + qJD(2) * t37 + qJDD(2) * t59 - t100 * t75 + (pkin(5) * qJDD(2) + qJD(4) * t89 - qJDD(3)) * t122 + t154;
t208 = t8 * t121 - t9 * t122;
t207 = pkin(5) * t122 + t59 + t82;
t133 = Icges(5,5) * t124 + Icges(5,6) * t125;
t54 = -Icges(5,3) * t121 + t122 * t133;
t171 = qJD(2) * t54;
t26 = t124 * t56 - t125 * t58;
t55 = Icges(5,6) * t122 + t121 * t134;
t66 = t96 * t122;
t30 = qJD(2) * t55 - qJD(4) * t66;
t173 = Icges(5,5) * t122;
t68 = t98 * t122;
t32 = -qJD(4) * t68 + (t121 * t135 + t173) * qJD(2);
t206 = t26 * qJD(4) + t124 * t32 + t125 * t30 + t171;
t136 = t124 * t96 - t125 * t98;
t87 = t134 * qJD(4);
t88 = t135 * qJD(4);
t94 = Icges(5,5) * t125 - Icges(5,6) * t124;
t205 = qJD(2) * t94 + qJD(4) * t136 + t124 * t88 + t125 * t87;
t101 = Icges(5,4) * t169;
t57 = Icges(5,1) * t170 + t101 + t173;
t139 = t124 * t55 - t125 * t57;
t53 = Icges(5,3) * t122 + t121 * t133;
t172 = qJD(2) * t53;
t31 = qJD(2) * t56 + t161 * t96;
t67 = t98 * t121;
t33 = qJD(2) * t58 + qJD(4) * t67;
t204 = qJD(4) * t139 - t124 * t33 - t125 * t31 + t172;
t188 = t58 + t66;
t190 = t56 - t68;
t202 = t124 * t190 - t125 * t188;
t189 = -Icges(5,2) * t170 + t101 + t57;
t191 = t55 - t67;
t201 = t124 * t191 - t125 * t189;
t200 = t74 / 0.2e1;
t199 = t75 / 0.2e1;
t198 = -pkin(2) - pkin(5);
t197 = t121 / 0.2e1;
t196 = rSges(4,2) - pkin(2);
t179 = rSges(4,3) * t122;
t80 = rSges(4,2) * t121 + t179;
t187 = -t79 + t80;
t52 = t82 + t83;
t63 = t121 * t94;
t178 = t122 * t94;
t73 = t100 * t161;
t22 = qJD(2) * t146 + t110 + t73;
t177 = t22 * t122;
t176 = -qJD(2) * t79 + t110;
t137 = t124 * t98 + t125 * t96;
t35 = t122 * t137 - t63;
t168 = t35 * qJD(2);
t167 = t133 * qJD(2);
t160 = qJD(4) * t122;
t159 = m(2) + m(3) + m(4);
t156 = -rSges(5,3) + t198;
t14 = t122 * t53 + t55 * t169 + t57 * t170;
t15 = -t122 * t54 - t56 * t169 - t58 * t170;
t150 = -t161 / 0.2e1;
t149 = t161 / 0.2e1;
t148 = -t160 / 0.2e1;
t147 = t160 / 0.2e1;
t140 = t124 * t57 + t125 * t55;
t128 = qJD(2) * t140 + qJD(4) * t63 + t171;
t129 = -qJD(2) * t138 - qJD(4) * t178 + t172;
t145 = (t128 * t121 + t204 * t122) * t122 + t121 * (t129 * t121 - t206 * t122);
t144 = t121 * (t206 * t121 + t129 * t122) + t122 * (-t204 * t121 + t128 * t122);
t84 = rSges(3,1) * t122 - rSges(3,2) * t121;
t81 = rSges(3,1) * t121 + rSges(3,2) * t122;
t142 = t121 * t15 + t122 * t14;
t48 = t121 * t53;
t16 = -t122 * t140 + t48;
t17 = -t121 * t54 + t209;
t141 = t121 * t17 + t122 * t16;
t127 = qJD(2) * t137 - t133 * qJD(4);
t69 = t100 * t121;
t43 = qJD(2) * t52 - t111;
t42 = qJD(2) * t187 + t110;
t34 = t121 * t137 + t178;
t27 = t34 * qJD(2);
t24 = qJD(1) + (-t121 * t59 - t122 * t60) * qJD(4);
t23 = t207 * qJD(2) - t100 * t160 - t111;
t18 = t187 * qJDD(2) + (-t83 * qJD(2) - t62) * qJD(2) + t166;
t13 = -t205 * t121 + t127 * t122;
t12 = t127 * t121 + t205 * t122;
t11 = qJD(4) * t138 - t124 * t30 + t125 * t32;
t10 = -qJD(4) * t140 - t124 * t31 + t125 * t33;
t7 = -t59 * t74 - t60 * t75 + qJDD(1) + (-t121 * t37 + t122 * t36) * qJD(4);
t6 = qJD(4) * t141 - t168;
t5 = qJD(4) * t142 + t27;
t1 = [m(5) * t7 + t159 * qJDD(1) + (-m(5) - t159) * g(3); (-qJD(4) * t137 + t124 * t87 - t125 * t88) * qJD(2) + (t27 + ((-t16 + t48 + t15) * t121 + (t17 - t209 + (-t140 + t54) * t121 + t14) * t122) * qJD(4)) * t150 - m(3) * (-g(1) * t81 + g(2) * t84) - t74 * t35 / 0.2e1 + t26 * t200 + (-t139 + t34) * t199 + (t168 + (t121 ^ 2 * t54 + (-t48 + t15 + (t140 + t54) * t122) * t122) * qJD(4) + t6) * t148 + (t10 + t13) * t147 + (t11 + t12 + t5) * t149 + (-(qJD(2) * t151 + t176 - t22 + t73) * t23 + t22 * (t111 + (-t152 + t153) * t122) + t23 * (-t121 * t152 + t155 + t165) + (t156 * t177 + (t22 * (-qJ(3) - t143) + t23 * t156) * t121) * qJD(2) + (t9 - g(2)) * t207 + (t8 - g(1)) * (t121 * t198 + t113 + t60)) * m(5) + (-(qJD(2) * t80 + t176 - t42) * t43 + t42 * t111 + t43 * (t164 + t165) + (t42 * t196 * t122 + (t42 * (-rSges(4,3) - qJ(3)) - t43 * pkin(2)) * t121) * qJD(2) + t211 * t52 + (t18 - g(1)) * (t121 * t196 + t113 + t179)) * m(4) + (-t136 + m(3) * (t81 ^ 2 + t84 ^ 2) + Icges(3,3) + Icges(4,1)) * qJDD(2); (-t122 * t211 + 0.2e1 * t18 * t197 - t210) * m(4) + (g(2) * t122 + t208 - t210) * m(5); -t5 * t163 / 0.2e1 + t122 * (t13 * qJD(2) + qJD(4) * t144 + t34 * qJDD(2) + t14 * t75 + t15 * t74) / 0.2e1 + t142 * t199 + ((-t14 * t121 + t15 * t122) * qJD(2) + t144) * t147 + t6 * t162 / 0.2e1 + (t12 * qJD(2) + qJD(4) * t145 - t35 * qJDD(2) + t16 * t75 + t17 * t74) * t197 + t141 * t200 + ((-t16 * t121 + t17 * t122) * qJD(2) + t145) * t149 + qJDD(2) * (t26 * t121 - t122 * t139) / 0.2e1 + qJD(2) * (t10 * t122 + t11 * t121 + (t121 * t139 + t26 * t122) * qJD(2)) / 0.2e1 + ((t160 * t63 - t167) * t122 + (-t212 + (t202 * t121 + (-t201 - t178) * t122) * qJD(4)) * t121) * t148 + ((-t161 * t178 - t167) * t121 + (t212 + (t201 * t122 + (-t202 + t63) * t121) * qJD(4)) * t122) * t150 - qJD(2) * ((-t124 * t186 - t125 * t185) * qJD(2) + ((t121 * t190 - t122 * t191) * t125 + (t121 * t188 - t122 * t189) * t124) * qJD(4)) / 0.2e1 + ((t23 * t89 - t7 * t60 + t24 * (-qJD(2) * t59 + t36)) * t122 + (-t22 * t89 - t7 * t59 + t24 * (qJD(2) * t60 - t37)) * t121 + ((t23 * t121 + t177) * qJD(2) + t208) * t100 - (t22 * t70 + t23 * t69) * qJD(2) - (t24 * (-t121 * t69 - t122 * t70) - (t22 * t121 - t122 * t23) * t143) * qJD(4) - g(1) * t69 + g(2) * t70 + g(3) * t143) * m(5);];
tau = t1;
