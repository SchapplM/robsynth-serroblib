% Calculate vector of inverse dynamics joint torques for
% S4PRPR3
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:47
% EndTime: 2019-12-31 16:20:52
% DurationCPUTime: 2.48s
% Computational Cost: add. (4473->296), mult. (3814->396), div. (0->0), fcn. (2940->6), ass. (0->164)
t133 = pkin(6) + qJ(2);
t129 = sin(t133);
t114 = qJD(3) * t129;
t131 = cos(t133);
t174 = qJD(4) * t131;
t132 = pkin(7) + qJ(4);
t128 = sin(t132);
t130 = cos(t132);
t203 = rSges(5,2) * t130;
t89 = rSges(5,1) * t128 + t203;
t159 = -t174 * t89 + t114;
t135 = cos(pkin(7));
t127 = pkin(3) * t135 + pkin(2);
t136 = -pkin(5) - qJ(3);
t182 = -t129 * t127 - t131 * t136;
t189 = t128 * t129;
t102 = rSges(5,2) * t189;
t187 = t129 * t130;
t59 = rSges(5,1) * t187 - t131 * rSges(5,3) - t102;
t227 = -t59 + t182;
t117 = t131 * qJ(3);
t217 = pkin(2) * t129;
t90 = -t117 + t217;
t213 = t90 + t227;
t170 = -t90 + t213;
t22 = qJD(2) * t170 + t159;
t115 = qJD(3) * t131;
t175 = qJD(4) * t129;
t226 = t131 * pkin(2) + t129 * qJ(3);
t160 = t131 * t127 - t129 * t136;
t122 = t129 * rSges(5,3);
t185 = t130 * t131;
t188 = t128 * t131;
t60 = rSges(5,1) * t185 - rSges(5,2) * t188 + t122;
t228 = t160 + t60;
t212 = -t226 + t228;
t23 = -t89 * t175 - t115 + (t226 + t212) * qJD(2);
t71 = t89 * t129;
t72 = t89 * t131;
t230 = t22 * t71 - t23 * t72;
t121 = Icges(5,4) * t130;
t148 = -Icges(5,2) * t128 + t121;
t87 = Icges(5,1) * t128 + t121;
t134 = sin(pkin(7));
t204 = rSges(4,2) * t134;
t206 = rSges(4,1) * t135;
t62 = t129 * rSges(4,3) + (-t204 + t206) * t131;
t195 = Icges(5,4) * t128;
t85 = Icges(5,2) * t130 + t195;
t152 = t128 * t85 - t130 * t87;
t84 = Icges(5,5) * t130 - Icges(5,6) * t128;
t225 = t152 * qJD(2) + t84 * qJD(4);
t83 = Icges(5,5) * t128 + Icges(5,6) * t130;
t144 = qJD(4) * t83;
t56 = Icges(5,6) * t129 + t131 * t148;
t201 = t128 * t56;
t88 = Icges(5,1) * t130 - t195;
t58 = Icges(5,5) * t129 + t131 * t88;
t154 = -t130 * t58 + t201;
t191 = Icges(5,3) * t131;
t224 = -t131 * t144 + (-t129 * t84 + t154 + t191) * qJD(2);
t192 = Icges(5,6) * t131;
t55 = Icges(5,4) * t187 - Icges(5,2) * t189 - t192;
t202 = t128 * t55;
t194 = Icges(5,5) * t131;
t99 = Icges(5,4) * t189;
t57 = Icges(5,1) * t187 - t194 - t99;
t155 = -t130 * t57 + t202;
t54 = Icges(5,3) * t129 + t131 * t84;
t190 = qJD(2) * t54;
t223 = qJD(2) * t155 - t129 * t144 + t190;
t53 = Icges(5,5) * t187 - Icges(5,6) * t189 - t191;
t18 = -t129 * t155 - t131 * t53;
t222 = t129 * (-t85 * t131 + t58) - t131 * (-Icges(5,2) * t187 + t57 - t99);
t171 = qJD(2) * qJD(4);
t79 = qJDD(4) * t129 + t131 * t171;
t221 = t79 / 0.2e1;
t80 = -qJDD(4) * t131 + t129 * t171;
t220 = t80 / 0.2e1;
t219 = t129 / 0.2e1;
t218 = -t131 / 0.2e1;
t215 = -t129 * t53 - t57 * t185;
t214 = t129 * t54 + t58 * t185;
t48 = t62 + t226;
t209 = -t85 + t88;
t208 = t87 + t148;
t169 = t129 * t206;
t105 = t129 * t204;
t180 = t131 * rSges(4,3) + t105;
t61 = t169 - t180;
t207 = -t90 - t61;
t205 = rSges(5,1) * t130;
t200 = t129 * t83;
t198 = t131 * t83;
t176 = qJD(2) * t131;
t197 = rSges(5,3) * t176 + qJD(2) * t102;
t28 = -t129 * t152 - t198;
t184 = t28 * qJD(2);
t183 = t84 * qJD(2);
t181 = rSges(4,3) * t176 + qJD(2) * t105;
t172 = qJD(2) * qJD(3);
t179 = qJDD(3) * t129 + t131 * t172;
t110 = qJ(3) * t176;
t178 = t110 + t114;
t177 = qJD(2) * t129;
t173 = m(2) + m(3) + m(4);
t167 = -pkin(2) - t206;
t166 = -t175 / 0.2e1;
t165 = t175 / 0.2e1;
t164 = -t174 / 0.2e1;
t163 = t174 / 0.2e1;
t162 = -t53 + t201;
t44 = t58 * t187;
t161 = t131 * t54 - t44;
t26 = t128 * t58 + t130 * t56;
t145 = qJD(4) * t85;
t32 = -t131 * t145 + (-t129 * t148 + t192) * qJD(2);
t146 = qJD(4) * t87;
t34 = -t131 * t146 + (-t129 * t88 + t194) * qJD(2);
t140 = -qJD(4) * t26 - t128 * t32 + t130 * t34 + t190;
t25 = t128 * t57 + t130 * t55;
t33 = qJD(2) * t56 - t129 * t145;
t35 = qJD(2) * t58 - t129 * t146;
t141 = qJD(2) * t53 - qJD(4) * t25 - t128 * t33 + t130 * t35;
t158 = -(t129 * t223 + t141 * t131) * t131 + t129 * (t129 * t224 + t140 * t131);
t157 = t129 * (t140 * t129 - t131 * t224) - t131 * (t141 * t129 - t131 * t223);
t94 = rSges(3,1) * t131 - rSges(3,2) * t129;
t91 = rSges(3,1) * t129 + rSges(3,2) * t131;
t92 = -rSges(5,2) * t128 + t205;
t153 = t128 * t87 + t130 * t85;
t19 = -t189 * t56 - t161;
t151 = t129 * t19 - t131 * t18;
t20 = -t188 * t55 - t215;
t21 = -t188 * t56 + t214;
t150 = t21 * t129 - t131 * t20;
t147 = -qJDD(3) * t131 + t129 * t172 + qJD(2) * (-pkin(2) * t177 + t178) + qJDD(2) * t226;
t143 = -t56 * t129 + t55 * t131;
t142 = (-t128 * t208 + t130 * t209) * qJD(2);
t74 = t148 * qJD(4);
t75 = t88 * qJD(4);
t139 = qJD(2) * t83 - qJD(4) * t153 - t128 * t74 + t130 * t75;
t138 = -t128 * t222 + t143 * t130;
t82 = qJD(2) * t90;
t76 = t92 * qJD(4);
t70 = qJD(2) * t226 - t115;
t41 = qJD(2) * t48 - t115;
t40 = qJD(2) * t207 + t114;
t37 = -qJD(4) * t71 + (t131 * t92 + t122) * qJD(2);
t36 = -t174 * t203 + (-t128 * t174 - t130 * t177) * rSges(5,1) + t197;
t29 = -t131 * t152 + t200;
t27 = t29 * qJD(2);
t24 = qJD(1) + (t129 * t59 + t131 * t60) * qJD(4);
t17 = qJDD(2) * t62 + qJD(2) * (-qJD(2) * t169 + t181) + t147;
t16 = t207 * qJDD(2) + (-t62 * qJD(2) - t70) * qJD(2) + t179;
t13 = t139 * t129 - t131 * t225;
t12 = t129 * t225 + t139 * t131;
t11 = -qJD(4) * t154 + t128 * t34 + t130 * t32;
t10 = -t155 * qJD(4) + t128 * t35 + t130 * t33;
t9 = t59 * t79 - t60 * t80 + qJDD(1) + (t129 * t37 + t131 * t36) * qJD(4);
t8 = -t76 * t175 - t79 * t89 + t212 * qJDD(2) + (-t110 + t36 + (t217 + t182) * qJD(2)) * qJD(2) + t147;
t7 = -t76 * t174 + t80 * t89 + t170 * qJDD(2) + (-t70 - t37 + (t226 - t160) * qJD(2)) * qJD(2) + t179;
t6 = qJD(4) * t150 + t27;
t5 = qJD(4) * t151 + t184;
t1 = [m(5) * t9 + t173 * qJDD(1) + (-m(5) - t173) * g(3); -m(3) * (-g(1) * t91 + g(2) * t94) + (t27 + ((t19 - t44 + (t54 + t202) * t131 + t215) * t131 + t214 * t129) * qJD(4)) * t163 + (-t152 * qJD(4) + t128 * t75 + t130 * t74) * qJD(2) + (t26 + t29) * t221 + (t25 + t28) * t220 + (-t184 + ((t131 * t162 + t21 - t214) * t131 + (t129 * t162 + t161 + t20) * t129) * qJD(4) + t5) * t166 + (t11 + t12) * t165 + (t10 + t13 + t6) * t164 + (t22 * t115 + t23 * (t114 + t197) + t230 * qJD(4) + ((t22 * (-t127 - t92) - t23 * t136) * t131 + (t22 * (-rSges(5,3) + t136) + t23 * (-t127 - t205)) * t129) * qJD(2) - (qJD(2) * t213 + t159 - t22 - t82) * t23 + (-g(2) + t8) * t228 + (-g(1) + t7) * t227) * m(5) + (t40 * t115 + t41 * (t178 + t181) + (t40 * (t167 + t204) * t131 + (t40 * (-rSges(4,3) - qJ(3)) + t41 * t167) * t129) * qJD(2) - (-qJD(2) * t61 + t114 - t40 - t82) * t41 + (t17 - g(2)) * t48 + (t16 - g(1)) * (t167 * t129 + t117 + t180)) * m(4) + (Icges(4,2) * t135 ^ 2 + (Icges(4,1) * t134 + 0.2e1 * Icges(4,4) * t135) * t134 + m(3) * (t91 ^ 2 + t94 ^ 2) + Icges(3,3) + t153) * qJDD(2); (-m(4) - m(5)) * (g(1) * t129 - g(2) * t131) + 0.2e1 * (t218 * t8 + t219 * t7) * m(5) + 0.2e1 * (t16 * t219 + t17 * t218) * m(4); t6 * t176 / 0.2e1 + (t12 * qJD(2) + qJD(4) * t158 + t29 * qJDD(2) + t20 * t80 + t21 * t79) * t219 + t150 * t221 + ((t20 * t129 + t21 * t131) * qJD(2) + t158) * t165 + t5 * t177 / 0.2e1 + (t13 * qJD(2) + qJD(4) * t157 + t28 * qJDD(2) + t18 * t80 + t19 * t79) * t218 + t151 * t220 + ((t18 * t129 + t19 * t131) * qJD(2) + t157) * t164 + qJDD(2) * (t26 * t129 - t25 * t131) / 0.2e1 + qJD(2) * (-t10 * t131 + t11 * t129 + (t25 * t129 + t131 * t26) * qJD(2)) / 0.2e1 + ((-t175 * t198 + t183) * t129 + (t142 + (t129 * t200 + t138) * qJD(4)) * t131) * t166 + ((-t174 * t200 - t183) * t131 + (t142 + (t131 * t198 + t138) * qJD(4)) * t129) * t163 - qJD(2) * ((t209 * t128 + t208 * t130) * qJD(2) + (t143 * t128 + t130 * t222) * qJD(4)) / 0.2e1 + ((-t22 * t76 + t9 * t60 + t24 * (qJD(2) * t59 + t36) + (-qJD(2) * t23 - t7) * t89) * t131 + (-t23 * t76 + t9 * t59 + t24 * (-qJD(2) * t60 + t37) + (qJD(2) * t22 - t8) * t89) * t129 - t230 * qJD(2) - (t24 * (-t129 * t71 - t131 * t72) + (-t129 * t23 - t131 * t22) * t92) * qJD(4) + g(1) * t72 + g(2) * t71 - g(3) * t92) * m(5);];
tau = t1;
