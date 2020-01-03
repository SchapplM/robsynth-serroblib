% Calculate vector of inverse dynamics joint torques for
% S4RPRR2
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:11
% DurationCPUTime: 3.01s
% Computational Cost: add. (5930->302), mult. (4561->405), div. (0->0), fcn. (3472->8), ass. (0->184)
t150 = sin(qJ(4));
t152 = cos(qJ(4));
t228 = rSges(5,1) * t152;
t124 = -rSges(5,2) * t150 + t228;
t103 = t124 * qJD(4);
t226 = rSges(5,2) * t152;
t122 = rSges(5,1) * t150 + t226;
t147 = qJDD(1) + qJDD(3);
t148 = qJD(1) + qJD(3);
t149 = qJ(1) + pkin(7);
t142 = sin(t149);
t143 = cos(t149);
t154 = qJD(1) ^ 2;
t151 = sin(qJ(1));
t153 = cos(qJ(1));
t168 = (-qJDD(1) * t151 - t153 * t154) * pkin(1);
t157 = (-qJDD(1) * t142 - t143 * t154) * pkin(2) + t168;
t144 = qJ(3) + t149;
t139 = cos(t144);
t202 = qJD(4) * t139;
t138 = sin(t144);
t212 = t138 * t150;
t114 = rSges(5,2) * t212;
t207 = t139 * rSges(5,3) + t114;
t211 = t138 * t152;
t70 = rSges(5,1) * t211 - t207;
t134 = t139 * pkin(6);
t94 = pkin(3) * t138 - t134;
t229 = -t94 - t70;
t246 = -t139 * pkin(3) - t138 * pkin(6);
t197 = qJD(4) * t226;
t209 = t139 * t150;
t199 = rSges(5,2) * t209;
t200 = qJD(4) * t150;
t198 = t148 * t199 + (rSges(5,1) * t200 + t197) * t138;
t208 = t139 * t152;
t247 = rSges(5,1) * t208 + t138 * rSges(5,3);
t47 = t247 * t148 - t198;
t201 = qJD(4) * t148;
t87 = -qJDD(4) * t139 + t138 * t201;
t12 = -t103 * t202 + t87 * t122 + t229 * t147 + (t246 * t148 - t47) * t148 + t157;
t254 = t12 - g(1);
t210 = t139 * t148;
t110 = pkin(6) * t210;
t137 = pkin(2) * t143;
t146 = t153 * pkin(1);
t140 = qJDD(1) * t146;
t237 = pkin(1) * t151;
t187 = -pkin(2) * t142 - t237;
t166 = qJDD(1) * t137 + t187 * t154 + t140;
t203 = qJD(4) * t138;
t213 = t138 * t148;
t172 = rSges(5,3) * t210 + t148 * t114 - t139 * t197;
t195 = t139 * t200;
t46 = (-t148 * t211 - t195) * rSges(5,1) + t172;
t71 = -t199 + t247;
t54 = t71 - t246;
t86 = qJDD(4) * t138 + t139 * t201;
t13 = -t103 * t203 - t86 * t122 + (-pkin(3) * t213 + t110 + t46) * t148 + t54 * t147 + t166;
t253 = t13 - g(2);
t108 = rSges(4,2) * t213;
t75 = rSges(4,1) * t210 - t108;
t92 = rSges(4,1) * t138 + rSges(4,2) * t139;
t252 = -t147 * t92 - t148 * t75 - g(1) + t157;
t223 = t148 * t92;
t132 = t139 * rSges(4,1);
t93 = -rSges(4,2) * t138 + t132;
t251 = t147 * t93 - t148 * t223 - g(2) + t166;
t62 = t148 * t70;
t250 = -t148 * t94 - t62;
t174 = t187 * qJD(1);
t196 = t122 * t202;
t161 = t174 - t196;
t30 = t229 * t148 + t161;
t249 = t148 * t30;
t97 = t143 * rSges(3,1) - rSges(3,2) * t142;
t90 = t146 + t97;
t245 = t137 + t146;
t145 = Icges(5,4) * t152;
t177 = -Icges(5,2) * t150 + t145;
t120 = Icges(5,1) * t150 + t145;
t117 = Icges(5,5) * t152 - Icges(5,6) * t150;
t116 = Icges(5,5) * t150 + Icges(5,6) * t152;
t162 = Icges(5,3) * t148 - t116 * qJD(4);
t170 = t177 * t139;
t67 = Icges(5,6) * t138 + t170;
t221 = t150 * t67;
t218 = Icges(5,4) * t150;
t121 = Icges(5,1) * t152 - t218;
t171 = t121 * t139;
t69 = Icges(5,5) * t138 + t171;
t179 = -t152 * t69 + t221;
t244 = -t117 * t213 + t162 * t139 + t179 * t148;
t169 = t117 * t139;
t66 = Icges(5,4) * t211 - Icges(5,2) * t212 - Icges(5,6) * t139;
t222 = t150 * t66;
t113 = Icges(5,4) * t212;
t68 = Icges(5,1) * t211 - Icges(5,5) * t139 - t113;
t180 = -t152 * t68 + t222;
t243 = t162 * t138 + (t169 + t180) * t148;
t118 = Icges(5,2) * t152 + t218;
t175 = t150 * t118 - t152 * t120;
t242 = t117 * qJD(4) + t175 * t148;
t64 = Icges(5,5) * t211 - Icges(5,6) * t212 - Icges(5,3) * t139;
t24 = -t180 * t138 - t139 * t64;
t231 = -Icges(5,2) * t211 - t113 + t68;
t233 = t120 * t138 + t66;
t241 = -t231 * t150 - t233 * t152;
t240 = t86 / 0.2e1;
t239 = t87 / 0.2e1;
t238 = m(3) + m(4);
t173 = t245 * qJD(1);
t61 = t148 * t93 + t173;
t236 = t61 * t92;
t235 = -t138 * t64 - t68 * t208;
t65 = Icges(5,3) * t138 + t169;
t234 = t138 * t65 + t69 * t208;
t232 = -t120 * t139 - t67;
t230 = -t118 * t139 + t69;
t225 = t139 * t30;
t215 = t116 * t139;
t49 = -t175 * t138 - t215;
t220 = t49 * t148;
t216 = t116 * t138;
t214 = t117 * t148;
t206 = -t118 + t121;
t205 = t120 + t177;
t194 = -pkin(3) - t228;
t193 = -t203 / 0.2e1;
t192 = t203 / 0.2e1;
t191 = -t202 / 0.2e1;
t190 = t202 / 0.2e1;
t55 = t69 * t211;
t189 = t139 * t65 - t55;
t188 = -t64 + t221;
t125 = rSges(2,1) * t153 - rSges(2,2) * t151;
t123 = rSges(2,1) * t151 + rSges(2,2) * t153;
t96 = rSges(3,1) * t142 + rSges(3,2) * t143;
t25 = -t67 * t212 - t189;
t184 = t138 * t25 - t139 * t24;
t26 = -t66 * t209 - t235;
t27 = -t67 * t209 + t234;
t183 = t138 * t27 - t139 * t26;
t91 = t122 * t203;
t31 = t54 * t148 + t173 - t91;
t182 = -t138 * t31 - t225;
t181 = t138 * t70 + t139 * t71;
t37 = t150 * t68 + t152 * t66;
t38 = t150 * t69 + t152 * t67;
t176 = t118 * t152 + t120 * t150;
t167 = -t230 * t150 + t232 * t152;
t53 = t194 * t138 + t134 + t207;
t60 = t174 - t223;
t165 = (-t205 * t150 + t206 * t152) * t148;
t164 = Icges(5,5) * t148 - qJD(4) * t120;
t163 = Icges(5,6) * t148 - t118 * qJD(4);
t50 = -t175 * t139 + t216;
t48 = t50 * t148;
t10 = t183 * qJD(4) + t48;
t101 = t177 * qJD(4);
t102 = t121 * qJD(4);
t43 = t163 * t138 + t148 * t170;
t45 = t164 * t138 + t148 * t171;
t16 = -t180 * qJD(4) + t150 * t45 + t152 * t43;
t42 = t163 * t139 - t177 * t213;
t44 = -t121 * t213 + t164 * t139;
t17 = -t179 * qJD(4) + t150 * t44 + t152 * t42;
t156 = -t176 * qJD(4) - t101 * t150 + t102 * t152 + t116 * t148;
t20 = t242 * t138 + t156 * t139;
t21 = t156 * t138 - t242 * t139;
t9 = t184 * qJD(4) + t220;
t160 = (t48 + ((t25 - t55 + (t65 + t222) * t139 + t235) * t139 + t234 * t138) * qJD(4)) * t190 + (-t175 * qJD(4) + t101 * t152 + t102 * t150) * t148 + (t38 + t50) * t240 + (t37 + t49) * t239 + (-t220 + ((t188 * t139 - t234 + t27) * t139 + (t188 * t138 + t189 + t26) * t138) * qJD(4) + t9) * t193 + (t17 + t20) * t192 + (Icges(4,3) + t176) * t147 + (t16 + t21 + t10) * t191;
t159 = -t37 * qJD(4) + t148 * t64 - t150 * t43 + t152 * t45;
t158 = -t38 * qJD(4) + t148 * t65 - t150 * t42 + t152 * t44;
t155 = t30 * t198 + t31 * (-rSges(5,1) * t195 + t110 + t172) + (t194 * t225 + (t30 * (-rSges(5,3) - pkin(6)) + t31 * t194) * t138) * t148;
t85 = t122 * t139;
t84 = t122 * t138;
t34 = t181 * qJD(4) + qJD(2);
t11 = t70 * t86 - t71 * t87 + qJDD(2) + (t138 * t47 + t139 * t46) * qJD(4);
t6 = t158 * t138 - t244 * t139;
t5 = t159 * t138 - t243 * t139;
t4 = t244 * t138 + t158 * t139;
t3 = t243 * t138 + t159 * t139;
t1 = [t160 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t60 * t108 + (-t60 * t132 - t236) * t148 + (t61 * t187 - t245 * t60) * qJD(1) + t251 * (t93 + t245) + t252 * (t187 - t92)) * m(4) + ((qJDD(1) * t97 - g(2) + t140) * t90 + (-qJDD(1) * t96 + t168 + (-0.2e1 * t97 + 0.2e1 * t90 - t146) * t154 - g(1)) * (-t96 - t237)) * m(3) + ((t123 ^ 2 + t125 ^ 2) * qJDD(1) + g(1) * t123 - g(2) * t125) * m(2) + ((t31 * t187 - t245 * t30) * qJD(1) + t155 - (-t30 + t161 + t250) * t31 + t253 * (t54 + t245) + t254 * (t187 + t53)) * m(5); m(5) * t11 + t238 * qJDD(2) + (-m(5) - t238) * g(3); t160 + (-t30 * t91 - t31 * (-t196 + t250) + t155 + (t249 + t253) * t54 + t254 * t53) * m(5) + (t236 * t148 - t60 * t75 - t61 * t223 + (t60 * t148 + t251) * t93 - t252 * t92) * m(4); t10 * t210 / 0.2e1 + t138 * (t147 * t50 + t148 * t20 + t26 * t87 + t27 * t86 + (t138 * t4 - t139 * t3) * qJD(4)) / 0.2e1 + t183 * t240 + ((t148 * t27 - t3) * t139 + (t148 * t26 + t4) * t138) * t192 + t9 * t213 / 0.2e1 - t139 * (t147 * t49 + t148 * t21 + t24 * t87 + t25 * t86 + (t138 * t6 - t139 * t5) * qJD(4)) / 0.2e1 + t184 * t239 + ((t148 * t25 - t5) * t139 + (t148 * t24 + t6) * t138) * t191 + t147 * (t138 * t38 - t139 * t37) / 0.2e1 + t148 * ((t148 * t38 - t16) * t139 + (t148 * t37 + t17) * t138) / 0.2e1 + ((-t203 * t215 + t214) * t138 + (t165 + (-t241 * t139 + (t216 + t167) * t138) * qJD(4)) * t139) * t193 + ((-t202 * t216 - t214) * t139 + (t165 + (t167 * t138 + (-t241 + t215) * t139) * qJD(4)) * t138) * t190 - t148 * ((t206 * t150 + t205 * t152) * t148 + ((t230 * t138 - t231 * t139) * t152 + (t232 * t138 + t233 * t139) * t150) * qJD(4)) / 0.2e1 + (t11 * t181 + t34 * ((t46 + t62) * t139 + (-t148 * t71 + t47) * t138) + t182 * t103 + ((-t148 * t31 - t12) * t139 + (-t13 + t249) * t138) * t122 - (t30 * t84 - t31 * t85) * t148 - (t34 * (-t138 * t84 - t139 * t85) + t182 * t124) * qJD(4) + g(1) * t85 + g(2) * t84 - g(3) * t124) * m(5);];
tau = t1;
