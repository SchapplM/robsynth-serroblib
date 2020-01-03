% Calculate vector of inverse dynamics joint torques for
% S4PRRR3
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:34
% EndTime: 2019-12-31 16:31:37
% DurationCPUTime: 2.13s
% Computational Cost: add. (5831->281), mult. (4345->378), div. (0->0), fcn. (3366->6), ass. (0->171)
t139 = pkin(7) + qJ(2);
t134 = sin(t139);
t212 = pkin(2) * qJD(2);
t184 = t134 * t212;
t140 = qJD(2) + qJD(3);
t136 = qJ(3) + t139;
t131 = sin(t136);
t132 = cos(t136);
t90 = rSges(4,1) * t131 + rSges(4,2) * t132;
t210 = t140 * t90;
t70 = -t184 - t210;
t141 = sin(qJ(4));
t142 = cos(qJ(4));
t214 = rSges(5,1) * t142;
t119 = -rSges(5,2) * t141 + t214;
t101 = t119 * qJD(4);
t213 = rSges(5,2) * t142;
t118 = rSges(5,1) * t141 + t213;
t138 = qJDD(2) + qJDD(3);
t135 = cos(t139);
t143 = qJD(2) ^ 2;
t156 = (-qJDD(2) * t134 - t135 * t143) * pkin(2);
t189 = qJD(4) * t132;
t198 = t131 * t141;
t110 = rSges(5,2) * t198;
t193 = t132 * rSges(5,3) + t110;
t197 = t131 * t142;
t68 = rSges(5,1) * t197 - t193;
t128 = t132 * pkin(6);
t92 = pkin(3) * t131 - t128;
t215 = -t92 - t68;
t231 = -t132 * pkin(3) - t131 * pkin(6);
t181 = qJD(4) * t213;
t195 = t132 * t141;
t185 = rSges(5,2) * t195;
t187 = qJD(4) * t141;
t182 = t140 * t185 + (rSges(5,1) * t187 + t181) * t131;
t194 = t132 * t142;
t232 = rSges(5,1) * t194 + t131 * rSges(5,3);
t47 = t232 * t140 - t182;
t188 = qJD(4) * t140;
t85 = -qJDD(4) * t132 + t131 * t188;
t12 = -t101 * t189 + t85 * t118 + t215 * t138 + t156 + (t231 * t140 - t47) * t140;
t236 = -g(1) + t12;
t196 = t132 * t140;
t199 = t131 * t140;
t73 = rSges(4,1) * t196 - rSges(4,2) * t199;
t235 = -t138 * t90 - t140 * t73 - g(1) + t156;
t106 = pkin(6) * t196;
t130 = pkin(2) * t135;
t222 = pkin(2) * t134;
t171 = qJDD(2) * t130 - t143 * t222;
t190 = qJD(4) * t131;
t160 = rSges(5,3) * t196 + t140 * t110 - t132 * t181;
t179 = t132 * t187;
t46 = (-t140 * t197 - t179) * rSges(5,1) + t160;
t69 = -t185 + t232;
t54 = t69 - t231;
t84 = qJDD(4) * t131 + t132 * t188;
t13 = -t101 * t190 - t84 * t118 + (-pkin(3) * t199 + t106 + t46) * t140 + t54 * t138 + t171;
t234 = -g(2) + t13;
t91 = t132 * rSges(4,1) - rSges(4,2) * t131;
t233 = t138 * t91 - t140 * t210 - g(2) + t171;
t137 = Icges(5,4) * t142;
t163 = -Icges(5,2) * t141 + t137;
t116 = Icges(5,1) * t141 + t137;
t60 = t140 * t68;
t230 = -rSges(5,1) * t179 + t140 * t92 + t106 + t160 + t60;
t229 = -t118 * t190 + t140 * t54;
t113 = Icges(5,5) * t142 - Icges(5,6) * t141;
t204 = Icges(5,4) * t141;
t114 = Icges(5,2) * t142 + t204;
t161 = t141 * t114 - t142 * t116;
t228 = t113 * qJD(4) + t140 * t161;
t112 = Icges(5,5) * t141 + Icges(5,6) * t142;
t150 = Icges(5,3) * t140 - qJD(4) * t112;
t158 = t163 * t132;
t65 = Icges(5,6) * t131 + t158;
t208 = t141 * t65;
t117 = Icges(5,1) * t142 - t204;
t159 = t117 * t132;
t67 = Icges(5,5) * t131 + t159;
t165 = -t142 * t67 + t208;
t227 = -t113 * t199 + t132 * t150 + t140 * t165;
t157 = t113 * t132;
t64 = Icges(5,4) * t197 - Icges(5,2) * t198 - Icges(5,6) * t132;
t209 = t141 * t64;
t109 = Icges(5,4) * t198;
t66 = Icges(5,1) * t197 - Icges(5,5) * t132 - t109;
t166 = -t142 * t66 + t209;
t226 = t131 * t150 + (t157 + t166) * t140;
t62 = Icges(5,5) * t197 - Icges(5,6) * t198 - Icges(5,3) * t132;
t24 = -t131 * t166 - t132 * t62;
t217 = -Icges(5,2) * t197 - t109 + t66;
t219 = t116 * t131 + t64;
t225 = -t217 * t141 - t219 * t142;
t224 = t84 / 0.2e1;
t223 = t85 / 0.2e1;
t221 = -t131 * t62 - t66 * t194;
t63 = Icges(5,3) * t131 + t157;
t220 = t131 * t63 + t67 * t194;
t218 = -t116 * t132 - t65;
t216 = -t114 * t132 + t67;
t180 = t118 * t189;
t155 = -t180 - t184;
t30 = t140 * t215 + t155;
t207 = t30 * t132;
t201 = t112 * t132;
t49 = -t131 * t161 - t201;
t206 = t49 * t140;
t202 = t112 * t131;
t200 = t113 * t140;
t192 = -t114 + t117;
t191 = t116 + t163;
t186 = m(2) + m(3) + m(4);
t183 = t135 * t212;
t178 = -pkin(3) - t214;
t177 = -t190 / 0.2e1;
t176 = t190 / 0.2e1;
t175 = -t189 / 0.2e1;
t174 = t189 / 0.2e1;
t55 = t67 * t197;
t173 = t132 * t63 - t55;
t172 = -t62 + t208;
t95 = rSges(3,1) * t135 - rSges(3,2) * t134;
t94 = rSges(3,1) * t134 + rSges(3,2) * t135;
t25 = -t198 * t65 - t173;
t170 = t131 * t25 - t132 * t24;
t26 = -t64 * t195 - t221;
t27 = -t195 * t65 + t220;
t169 = t131 * t27 - t132 * t26;
t31 = t183 + t229;
t168 = -t31 * t131 - t207;
t167 = t131 * t68 + t132 * t69;
t35 = t141 * t66 + t142 * t64;
t36 = t141 * t67 + t142 * t65;
t162 = t114 * t142 + t116 * t141;
t154 = -t216 * t141 + t142 * t218;
t53 = t131 * t178 + t128 + t193;
t153 = (-t141 * t191 + t142 * t192) * t140;
t152 = Icges(5,5) * t140 - qJD(4) * t116;
t151 = Icges(5,6) * t140 - qJD(4) * t114;
t50 = -t132 * t161 + t202;
t48 = t50 * t140;
t10 = qJD(4) * t169 + t48;
t100 = t117 * qJD(4);
t43 = t131 * t151 + t140 * t158;
t45 = t131 * t152 + t140 * t159;
t16 = -qJD(4) * t166 + t141 * t45 + t142 * t43;
t42 = t132 * t151 - t163 * t199;
t44 = -t117 * t199 + t132 * t152;
t17 = -qJD(4) * t165 + t141 * t44 + t142 * t42;
t99 = t163 * qJD(4);
t145 = -qJD(4) * t162 + t100 * t142 + t112 * t140 - t141 * t99;
t20 = t228 * t131 + t145 * t132;
t21 = t145 * t131 - t228 * t132;
t9 = qJD(4) * t170 + t206;
t148 = (t48 + ((t25 - t55 + (t63 + t209) * t132 + t221) * t132 + t220 * t131) * qJD(4)) * t174 + (-qJD(4) * t161 + t100 * t141 + t142 * t99) * t140 + (t36 + t50) * t224 + (t35 + t49) * t223 + (-t206 + ((t132 * t172 - t220 + t27) * t132 + (t131 * t172 + t173 + t26) * t131) * qJD(4) + t9) * t177 + (t17 + t20) * t176 + (Icges(4,3) + t162) * t138 + (t16 + t21 + t10) * t175;
t147 = -qJD(4) * t35 + t140 * t62 - t141 * t43 + t142 * t45;
t146 = -qJD(4) * t36 + t140 * t63 - t141 * t42 + t142 * t44;
t144 = (t178 * t207 + (t30 * (-rSges(5,3) - pkin(6)) + t31 * t178) * t131) * t140;
t83 = t118 * t132;
t82 = t118 * t131;
t71 = t140 * t91 + t183;
t34 = qJD(4) * t167 + qJD(1);
t11 = t68 * t84 - t69 * t85 + qJDD(1) + (t131 * t47 + t132 * t46) * qJD(4);
t6 = t146 * t131 - t227 * t132;
t5 = t147 * t131 - t226 * t132;
t4 = t227 * t131 + t146 * t132;
t3 = t226 * t131 + t147 * t132;
t1 = [m(5) * t11 + t186 * qJDD(1) + (-m(5) - t186) * g(3); Icges(3,3) * qJDD(2) + t148 + (t233 * (t130 + t91) + t235 * (-t90 - t222) + (-t73 - t183 + t71) * t70) * m(4) + (g(1) * t94 - g(2) * t95 + (t94 ^ 2 + t95 ^ 2) * qJDD(2)) * m(3) + (t30 * (t182 - t183) + t144 + t234 * (t130 + t54) + t236 * (t53 - t222) + (-t155 + t30 - t184 + t230) * t31) * m(5); t148 + (t144 + t234 * t54 + t236 * t53 + (t180 + t230) * t31 + (t182 + t229) * t30) * m(5) + (-t70 * t73 - t71 * t210 + (t70 * t140 + t233) * t91 + (t71 * t140 - t235) * t90) * m(4); t10 * t196 / 0.2e1 + t131 * (t138 * t50 + t140 * t20 + t26 * t85 + t27 * t84 + (t131 * t4 - t132 * t3) * qJD(4)) / 0.2e1 + t169 * t224 + ((t140 * t27 - t3) * t132 + (t140 * t26 + t4) * t131) * t176 + t9 * t199 / 0.2e1 - t132 * (t138 * t49 + t140 * t21 + t24 * t85 + t25 * t84 + (t131 * t6 - t132 * t5) * qJD(4)) / 0.2e1 + t170 * t223 + ((t140 * t25 - t5) * t132 + (t140 * t24 + t6) * t131) * t175 + t138 * (t131 * t36 - t132 * t35) / 0.2e1 + t140 * ((t140 * t36 - t16) * t132 + (t140 * t35 + t17) * t131) / 0.2e1 + ((-t190 * t201 + t200) * t131 + (t153 + (-t225 * t132 + (t202 + t154) * t131) * qJD(4)) * t132) * t177 + ((-t189 * t202 - t200) * t132 + (t153 + (t154 * t131 + (-t225 + t201) * t132) * qJD(4)) * t131) * t174 - t140 * ((t141 * t192 + t142 * t191) * t140 + ((t131 * t216 - t217 * t132) * t142 + (t218 * t131 + t219 * t132) * t141) * qJD(4)) / 0.2e1 + (t11 * t167 + t34 * ((t46 + t60) * t132 + (-t140 * t69 + t47) * t131) + t168 * t101 + ((-t140 * t31 - t12) * t132 + (t140 * t30 - t13) * t131) * t118 - (t30 * t82 - t31 * t83) * t140 - (t34 * (-t131 * t82 - t132 * t83) + t168 * t119) * qJD(4) + g(1) * t83 + g(2) * t82 - g(3) * t119) * m(5);];
tau = t1;
