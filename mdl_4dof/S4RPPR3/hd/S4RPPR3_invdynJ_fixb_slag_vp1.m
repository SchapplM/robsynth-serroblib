% Calculate vector of inverse dynamics joint torques for
% S4RPPR3
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:50
% EndTime: 2019-12-31 16:37:55
% DurationCPUTime: 4.00s
% Computational Cost: add. (4572->321), mult. (4030->421), div. (0->0), fcn. (3046->8), ass. (0->176)
t147 = cos(qJ(1));
t140 = t147 * pkin(1);
t146 = sin(qJ(1));
t234 = pkin(1) * t146;
t142 = qJ(1) + pkin(6);
t137 = sin(t142);
t139 = cos(t142);
t93 = rSges(3,1) * t137 + rSges(3,2) * t139;
t75 = -t93 - t234;
t148 = qJD(1) ^ 2;
t250 = t148 * t140;
t144 = cos(pkin(7));
t134 = pkin(3) * t144 + pkin(2);
t145 = -pkin(5) - qJ(3);
t198 = -t137 * t134 - t139 * t145;
t141 = pkin(7) + qJ(4);
t136 = sin(t141);
t206 = t136 * t137;
t106 = rSges(5,2) * t206;
t138 = cos(t141);
t204 = t137 * t138;
t59 = rSges(5,1) * t204 - t139 * rSges(5,3) - t106;
t42 = t198 - t59 - t234;
t123 = t139 * qJ(3);
t233 = pkin(2) * t137;
t92 = -t123 + t233;
t176 = t92 + t42;
t162 = t176 - t92;
t120 = qJD(3) * t137;
t191 = qJD(4) * t139;
t220 = rSges(5,2) * t138;
t91 = rSges(5,1) * t136 + t220;
t173 = -t191 * t91 + t120;
t22 = qJD(1) * t162 + t173;
t121 = qJD(3) * t139;
t245 = t139 * pkin(2) + t137 * qJ(3);
t184 = t245 + t140;
t192 = qJD(4) * t137;
t177 = t139 * t134 - t137 * t145;
t128 = t137 * rSges(5,3);
t202 = t138 * t139;
t205 = t136 * t139;
t60 = rSges(5,1) * t202 - rSges(5,2) * t205 + t128;
t247 = t177 + t60;
t229 = -t245 + t247;
t23 = -t91 * t192 - t121 + (t184 + t229) * qJD(1);
t71 = t91 * t137;
t72 = t91 * t139;
t249 = t22 * t71 - t23 * t72;
t96 = t139 * rSges(3,1) - rSges(3,2) * t137;
t76 = t140 + t96;
t143 = sin(pkin(7));
t221 = rSges(4,2) * t143;
t224 = rSges(4,1) * t144;
t62 = t137 * rSges(4,3) + (-t221 + t224) * t139;
t45 = t184 + t62;
t127 = Icges(5,4) * t138;
t160 = -Icges(5,2) * t136 + t127;
t89 = Icges(5,1) * t136 + t127;
t109 = t137 * t221;
t244 = t139 * rSges(4,3) + t109 - t234;
t212 = Icges(5,4) * t136;
t87 = Icges(5,2) * t138 + t212;
t165 = t136 * t87 - t138 * t89;
t86 = Icges(5,5) * t138 - Icges(5,6) * t136;
t243 = t165 * qJD(1) + t86 * qJD(4);
t85 = Icges(5,5) * t136 + Icges(5,6) * t138;
t156 = qJD(4) * t85;
t56 = Icges(5,6) * t137 + t139 * t160;
t218 = t136 * t56;
t90 = Icges(5,1) * t138 - t212;
t58 = Icges(5,5) * t137 + t139 * t90;
t167 = -t138 * t58 + t218;
t208 = Icges(5,3) * t139;
t242 = -t139 * t156 + (-t137 * t86 + t167 + t208) * qJD(1);
t209 = Icges(5,6) * t139;
t55 = Icges(5,4) * t204 - Icges(5,2) * t206 - t209;
t219 = t136 * t55;
t103 = Icges(5,4) * t206;
t211 = Icges(5,5) * t139;
t57 = Icges(5,1) * t204 - t103 - t211;
t168 = -t138 * t57 + t219;
t54 = Icges(5,3) * t137 + t139 * t86;
t207 = qJD(1) * t54;
t241 = qJD(1) * t168 - t137 * t156 + t207;
t53 = Icges(5,5) * t204 - Icges(5,6) * t206 - t208;
t18 = -t137 * t168 - t139 * t53;
t240 = t137 * (-t87 * t139 + t58) - t139 * (-Icges(5,2) * t204 - t103 + t57);
t189 = qJD(1) * qJD(4);
t81 = qJDD(4) * t137 + t139 * t189;
t239 = t81 / 0.2e1;
t82 = -qJDD(4) * t139 + t137 * t189;
t238 = t82 / 0.2e1;
t237 = -m(4) - m(5);
t236 = t137 / 0.2e1;
t235 = -t139 / 0.2e1;
t231 = -t137 * t53 - t57 * t202;
t230 = t137 * t54 + t58 * t202;
t226 = -t87 + t90;
t225 = t89 + t160;
t223 = rSges(5,1) * t138;
t217 = t137 * t85;
t215 = t139 * t85;
t193 = qJD(1) * t139;
t214 = rSges(5,3) * t193 + qJD(1) * t106;
t28 = -t137 * t165 - t215;
t200 = t28 * qJD(1);
t199 = t86 * qJD(1);
t197 = rSges(4,3) * t193 + qJD(1) * t109;
t116 = qJ(3) * t193;
t195 = t116 + t120;
t194 = qJD(1) * t137;
t190 = qJD(1) * qJD(3);
t188 = t137 * t224;
t186 = -pkin(2) - t224;
t185 = -t188 + t244;
t183 = -t192 / 0.2e1;
t182 = t192 / 0.2e1;
t181 = -t191 / 0.2e1;
t180 = t191 / 0.2e1;
t179 = -t53 + t218;
t46 = t58 * t204;
t178 = t139 * t54 - t46;
t175 = t185 - t92;
t174 = qJDD(1) * t140 - t148 * t234;
t26 = t136 * t58 + t138 * t56;
t157 = qJD(4) * t87;
t34 = -t139 * t157 + (-t137 * t160 + t209) * qJD(1);
t158 = qJD(4) * t89;
t36 = -t139 * t158 + (-t137 * t90 + t211) * qJD(1);
t151 = -qJD(4) * t26 - t136 * t34 + t138 * t36 + t207;
t25 = t136 * t57 + t138 * t55;
t35 = qJD(1) * t56 - t137 * t157;
t37 = qJD(1) * t58 - t137 * t158;
t152 = qJD(1) * t53 - qJD(4) * t25 - t136 * t35 + t138 * t37;
t172 = -(t241 * t137 + t152 * t139) * t139 + t137 * (t242 * t137 + t151 * t139);
t171 = t137 * (t151 * t137 - t242 * t139) - t139 * (t152 * t137 - t241 * t139);
t112 = rSges(2,1) * t147 - rSges(2,2) * t146;
t111 = rSges(2,1) * t146 + rSges(2,2) * t147;
t94 = -rSges(5,2) * t136 + t223;
t166 = t136 * t89 + t138 * t87;
t19 = -t206 * t56 - t178;
t164 = t137 * t19 - t139 * t18;
t20 = -t205 * t55 - t231;
t21 = -t205 * t56 + t230;
t163 = t137 * t21 - t139 * t20;
t159 = qJDD(3) * t137 + t139 * t190 - t250;
t155 = -t137 * t56 + t139 * t55;
t154 = (-t136 * t225 + t138 * t226) * qJD(1);
t153 = -qJDD(3) * t139 + t137 * t190 + t174 + qJD(1) * (-pkin(2) * t194 + t195) + qJDD(1) * t245;
t78 = t160 * qJD(4);
t79 = t90 * qJD(4);
t150 = qJD(1) * t85 - qJD(4) * t166 - t136 * t78 + t138 * t79;
t149 = -t240 * t136 + t155 * t138;
t84 = qJD(1) * t92;
t80 = t94 * qJD(4);
t70 = qJD(1) * t245 - t121;
t41 = qJD(1) * t45 - t121;
t40 = qJD(1) * t175 + t120;
t39 = -qJD(4) * t71 + (t139 * t94 + t128) * qJD(1);
t38 = -t191 * t220 + (-t136 * t191 - t138 * t194) * rSges(5,1) + t214;
t29 = -t139 * t165 + t217;
t27 = t29 * qJD(1);
t24 = qJD(2) + (t137 * t59 + t139 * t60) * qJD(4);
t15 = qJDD(1) * t62 + qJD(1) * (-qJD(1) * t188 + t197) + t153;
t14 = t175 * qJDD(1) + (-t62 * qJD(1) - t70) * qJD(1) + t159;
t13 = t150 * t137 - t243 * t139;
t12 = t243 * t137 + t150 * t139;
t11 = -qJD(4) * t167 + t136 * t36 + t138 * t34;
t10 = -t168 * qJD(4) + t136 * t37 + t138 * t35;
t9 = t59 * t81 - t60 * t82 + qJDD(2) + (t137 * t39 + t139 * t38) * qJD(4);
t8 = -t80 * t192 - t81 * t91 + t229 * qJDD(1) + (-t116 + t38 + (t233 + t198) * qJD(1)) * qJD(1) + t153;
t7 = -t80 * t191 + t82 * t91 + t162 * qJDD(1) + (-t70 - t39 + (t245 - t177) * qJD(1)) * qJD(1) + t159;
t6 = qJD(4) * t163 + t27;
t5 = qJD(4) * t164 + t200;
t1 = [(-t165 * qJD(4) + t136 * t79 + t138 * t78) * qJD(1) - m(2) * (-g(1) * t111 + g(2) * t112) + (t27 + ((t19 - t46 + (t54 + t219) * t139 + t231) * t139 + t230 * t137) * qJD(4)) * t180 + (t29 + t26) * t239 + (t28 + t25) * t238 + (-t200 + ((t139 * t179 + t21 - t230) * t139 + (t137 * t179 + t178 + t20) * t137) * qJD(4) + t5) * t183 + (t11 + t12) * t182 + ((-t93 * t148 - g(2) + t174) * t76 + (-t250 + (-0.2e1 * t96 - t140 + t76) * t148 - g(1)) * t75) * m(3) + (t10 + t13 + t6) * t181 + (t22 * t121 + t23 * (t120 + t214) + t249 * qJD(4) + ((-t23 * t146 - t22 * t147) * pkin(1) + (t22 * (-t134 - t94) - t23 * t145) * t139 + (t22 * (-rSges(5,3) + t145) + t23 * (-t134 - t223)) * t137) * qJD(1) - (qJD(1) * t176 + t173 - t22 - t84) * t23 + (-g(2) + t8) * (t140 + t247) + (-g(1) + t7) * t42) * m(5) + (t40 * t121 + t41 * (t195 + t197) + ((-t41 * t146 - t40 * t147) * pkin(1) + t40 * (t186 + t221) * t139 + (t40 * (-rSges(4,3) - qJ(3)) + t41 * t186) * t137) * qJD(1) - (qJD(1) * t185 + t120 - t40 - t84) * t41 + (t15 - g(2)) * t45 + (t14 - g(1)) * (t186 * t137 + t123 + t244)) * m(4) + (t166 + Icges(4,2) * t144 ^ 2 + (Icges(4,1) * t143 + 0.2e1 * Icges(4,4) * t144) * t143 + m(2) * (t111 ^ 2 + t112 ^ 2) + m(3) * (t75 ^ 2 + t96 * t76) + Icges(2,3) + Icges(3,3)) * qJDD(1); m(5) * t9 + (m(3) + m(4)) * qJDD(2) + (-m(3) + t237) * g(3); t237 * (g(1) * t137 - g(2) * t139) + 0.2e1 * (t235 * t8 + t236 * t7) * m(5) + 0.2e1 * (t14 * t236 + t15 * t235) * m(4); t6 * t193 / 0.2e1 + (qJD(1) * t12 + qJD(4) * t172 + qJDD(1) * t29 + t20 * t82 + t21 * t81) * t236 + t163 * t239 + ((t20 * t137 + t21 * t139) * qJD(1) + t172) * t182 + t5 * t194 / 0.2e1 + (qJD(1) * t13 + qJD(4) * t171 + qJDD(1) * t28 + t18 * t82 + t19 * t81) * t235 + t164 * t238 + ((t18 * t137 + t19 * t139) * qJD(1) + t171) * t181 + qJDD(1) * (t137 * t26 - t139 * t25) / 0.2e1 + qJD(1) * (-t10 * t139 + t11 * t137 + (t25 * t137 + t139 * t26) * qJD(1)) / 0.2e1 + ((-t192 * t215 + t199) * t137 + (t154 + (t137 * t217 + t149) * qJD(4)) * t139) * t183 + ((-t191 * t217 - t199) * t139 + (t154 + (t139 * t215 + t149) * qJD(4)) * t137) * t180 - qJD(1) * ((t226 * t136 + t225 * t138) * qJD(1) + (t155 * t136 + t240 * t138) * qJD(4)) / 0.2e1 + ((-t22 * t80 + t9 * t60 + t24 * (qJD(1) * t59 + t38) + (-qJD(1) * t23 - t7) * t91) * t139 + (-t23 * t80 + t9 * t59 + t24 * (-qJD(1) * t60 + t39) + (qJD(1) * t22 - t8) * t91) * t137 - t249 * qJD(1) - (t24 * (-t137 * t71 - t139 * t72) + (-t23 * t137 - t22 * t139) * t94) * qJD(4) + g(1) * t72 + g(2) * t71 - g(3) * t94) * m(5);];
tau = t1;
