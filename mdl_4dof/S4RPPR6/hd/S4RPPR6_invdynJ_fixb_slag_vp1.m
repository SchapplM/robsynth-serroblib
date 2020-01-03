% Calculate vector of inverse dynamics joint torques for
% S4RPPR6
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPPR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR6_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:37
% EndTime: 2019-12-31 16:40:43
% DurationCPUTime: 4.21s
% Computational Cost: add. (3091->395), mult. (8107->512), div. (0->0), fcn. (7853->6), ass. (0->181)
t173 = cos(pkin(6));
t175 = sin(qJ(1));
t234 = t173 * t175;
t176 = cos(qJ(1));
t250 = pkin(5) * t176;
t132 = pkin(3) * t234 + t250;
t174 = sin(qJ(4));
t172 = sin(pkin(6));
t252 = cos(qJ(4));
t210 = t172 * t252;
t116 = t174 * t234 - t175 * t210;
t185 = t172 * t174 + t173 * t252;
t117 = t185 * t175;
t64 = -t117 * rSges(5,1) + t116 * rSges(5,2) - rSges(5,3) * t176;
t271 = -t132 + t64;
t129 = -t173 * t174 + t210;
t118 = t129 * t176;
t119 = t185 * t176;
t104 = Icges(5,4) * t117;
t57 = -Icges(5,2) * t116 + Icges(5,6) * t176 + t104;
t103 = Icges(5,4) * t116;
t61 = -Icges(5,1) * t117 - Icges(5,5) * t176 + t103;
t248 = t118 * t57 - t119 * t61;
t54 = Icges(5,5) * t117 - Icges(5,6) * t116 + Icges(5,3) * t176;
t16 = -t175 * t54 + t248;
t269 = -t116 * t57 - t117 * t61;
t28 = -t129 * t61 - t185 * t57;
t14 = t176 * t54 + t269;
t239 = Icges(5,4) * t119;
t59 = Icges(5,2) * t118 - Icges(5,6) * t175 + t239;
t105 = Icges(5,4) * t118;
t62 = Icges(5,1) * t119 - Icges(5,5) * t175 + t105;
t247 = t118 * t59 + t119 * t62;
t56 = Icges(5,5) * t119 + Icges(5,6) * t118 - Icges(5,3) * t175;
t200 = t175 * t56 - t247;
t267 = t200 - t14;
t92 = Icges(5,5) * t129 - Icges(5,6) * t185;
t238 = Icges(5,4) * t129;
t94 = -Icges(5,2) * t185 + t238;
t127 = Icges(5,4) * t185;
t96 = Icges(5,1) * t129 - t127;
t184 = t116 * t94 - t117 * t96 - t176 * t92;
t266 = t184 * qJD(1);
t265 = -t116 * t59 + t117 * t62;
t233 = t173 * t176;
t154 = pkin(3) * t233;
t264 = pkin(5) * t175 - t154;
t140 = t176 * pkin(1) + t175 * qJ(2);
t235 = t172 * t176;
t112 = rSges(4,1) * t233 + t175 * rSges(4,2) + rSges(4,3) * t235;
t126 = pkin(2) * t233 + qJ(3) * t235;
t201 = t126 + t140;
t262 = t201 + t112;
t261 = qJD(4) * t129;
t113 = rSges(3,1) * t233 - rSges(3,2) * t235 + t175 * rSges(3,3);
t260 = t175 * (-Icges(5,2) * t119 + t105 + t62) - t176 * (-Icges(5,2) * t117 - t103 - t61);
t259 = t173 ^ 2;
t258 = -m(4) - m(5);
t217 = qJD(1) * qJD(4);
t134 = -qJDD(4) * t175 - t176 * t217;
t257 = t134 / 0.2e1;
t135 = qJDD(4) * t176 - t175 * t217;
t256 = t135 / 0.2e1;
t255 = t175 / 0.2e1;
t254 = -t176 / 0.2e1;
t253 = -rSges(5,3) - pkin(5);
t249 = g(1) * t175;
t123 = t185 * qJD(4);
t223 = qJD(1) * t175;
t72 = -t123 * t176 - t129 * t223;
t73 = t176 * t261 - t185 * t223;
t244 = t73 * rSges(5,1) + t72 * rSges(5,2);
t243 = -Icges(5,2) * t129 - t127 + t96;
t242 = -Icges(5,1) * t185 - t238 - t94;
t232 = t119 * rSges(5,1) + t118 * rSges(5,2);
t65 = -rSges(5,3) * t175 + t232;
t240 = t65 - t264;
t237 = qJ(3) * t172;
t236 = t172 * t175;
t100 = t113 + t140;
t194 = pkin(2) * t173 + t237;
t125 = t194 * t175;
t166 = t176 * qJ(2);
t138 = pkin(1) * t175 - t166;
t231 = -t125 - t138;
t215 = rSges(3,1) * t234;
t150 = rSges(3,2) * t236;
t227 = t176 * rSges(3,3) + t150;
t111 = t215 - t227;
t230 = -t138 - t111;
t222 = qJD(1) * t176;
t229 = rSges(3,3) * t222 + qJD(1) * t150;
t221 = qJD(3) * t172;
t148 = t176 * t221;
t163 = qJD(2) * t175;
t228 = t148 + t163;
t218 = qJD(1) * qJD(2);
t226 = qJDD(2) * t175 + t176 * t218;
t225 = qJ(2) * t222 + t163;
t224 = -qJD(1) * t138 + t163;
t220 = qJD(4) * t175;
t219 = qJD(4) * t176;
t216 = qJDD(3) * t172;
t15 = t176 * t56 + t265;
t170 = t176 * rSges(4,2);
t195 = rSges(4,1) * t173 + rSges(4,3) * t172;
t110 = t175 * t195 - t170;
t213 = -t110 + t231;
t212 = t176 * t216 + t226;
t211 = t148 + t225;
t209 = t175 * t221;
t208 = -rSges(3,1) * t173 - pkin(1);
t207 = -t220 / 0.2e1;
t206 = t220 / 0.2e1;
t205 = -t219 / 0.2e1;
t204 = t219 / 0.2e1;
t203 = t231 + t271;
t202 = -qJD(1) * t125 + t148 + t224;
t74 = qJD(1) * t118 - t123 * t175;
t75 = qJD(1) * t119 + t175 * t261;
t199 = -rSges(5,1) * t75 - rSges(5,2) * t74;
t31 = Icges(5,5) * t73 + Icges(5,6) * t72 - Icges(5,3) * t222;
t32 = Icges(5,5) * t75 + Icges(5,6) * t74 - Icges(5,3) * t223;
t33 = Icges(5,4) * t73 + Icges(5,2) * t72 - Icges(5,6) * t222;
t34 = Icges(5,4) * t75 + Icges(5,2) * t74 - Icges(5,6) * t223;
t35 = Icges(5,1) * t73 + Icges(5,4) * t72 - Icges(5,5) * t222;
t36 = Icges(5,1) * t75 + Icges(5,4) * t74 - Icges(5,5) * t223;
t198 = (t118 * t34 + t119 * t36 - t175 * t32 - t222 * t54 + t57 * t72 - t61 * t73) * t176 - t175 * (t118 * t33 + t119 * t35 - t175 * t31 - t222 * t56 + t59 * t72 + t62 * t73);
t197 = -t175 * (-t116 * t33 + t117 * t35 + t176 * t31 - t223 * t56 + t59 * t74 + t62 * t75) + t176 * (-t116 * t34 + t117 * t36 + t176 * t32 - t223 * t54 + t57 * t74 - t61 * t75);
t164 = qJD(2) * t176;
t196 = t164 - t209;
t141 = rSges(2,1) * t176 - rSges(2,2) * t175;
t139 = rSges(2,1) * t175 + rSges(2,2) * t176;
t193 = t14 * t176 - t15 * t175;
t192 = t16 * t176 + t175 * t200;
t191 = t175 * (Icges(5,5) * t118 - Icges(5,6) * t119) - t176 * (-Icges(5,5) * t116 - Icges(5,6) * t117);
t188 = -pkin(1) - t194;
t122 = t140 * qJD(1) - t164;
t187 = -t194 * t222 - t122 - 0.2e1 * t209;
t186 = -qJDD(2) * t176 + qJD(1) * (-pkin(1) * t223 + t225) + qJDD(1) * t140 + t175 * t218;
t183 = -(Icges(5,1) * t118 - t239 - t59) * t175 + (-Icges(5,1) * t116 - t104 - t57) * t176;
t181 = -pkin(3) * t173 + t188;
t179 = -pkin(1) + (-rSges(4,1) - pkin(2)) * t173 + (-rSges(4,3) - qJ(3)) * t172;
t178 = qJDD(1) * t126 + t175 * t216 + t186 + (-t194 * t223 + 0.2e1 * t148) * qJD(1);
t161 = rSges(4,2) * t222;
t98 = rSges(5,1) * t129 - rSges(5,2) * t185;
t97 = -rSges(5,1) * t185 - rSges(5,2) * t129;
t91 = -Icges(5,5) * t185 - Icges(5,6) * t129;
t90 = t98 * t219;
t89 = -rSges(5,1) * t123 - rSges(5,2) * t261;
t88 = -Icges(5,1) * t123 - Icges(5,4) * t261;
t87 = -Icges(5,4) * t123 - Icges(5,2) * t261;
t86 = -Icges(5,5) * t123 - Icges(5,6) * t261;
t85 = qJD(1) * t100 - t164;
t84 = qJD(1) * t230 + t163;
t83 = rSges(5,1) * t118 - rSges(5,2) * t119;
t82 = -rSges(5,1) * t116 - rSges(5,2) * t117;
t52 = t262 * qJD(1) - t196;
t51 = qJD(1) * t213 + t228;
t40 = qJDD(1) * t113 + qJD(1) * (-qJD(1) * t215 + t229) + t186;
t39 = t230 * qJDD(1) + (-t113 * qJD(1) - t122) * qJD(1) + t226;
t38 = -rSges(5,3) * t223 - t199;
t37 = -rSges(5,3) * t222 + t244;
t30 = -qJD(3) * t173 + (t175 * t64 - t176 * t65) * qJD(4);
t29 = t129 * t62 - t185 * t59;
t27 = t118 * t94 + t119 * t96 - t175 * t92;
t25 = t27 * qJD(1);
t24 = -t164 + (qJD(4) * t98 + t221) * t175 + (t201 + t240) * qJD(1);
t23 = qJD(1) * t203 + t228 + t90;
t19 = qJDD(1) * t112 + (-t195 * t223 + t161) * qJD(1) + t178;
t18 = t213 * qJDD(1) + (-t112 * qJD(1) + t187) * qJD(1) + t212;
t13 = -qJDD(3) * t173 - t134 * t64 - t135 * t65 + (-t175 * t38 - t176 * t37) * qJD(4);
t12 = -t116 * t87 + t117 * t88 + t176 * t86 - t223 * t92 + t74 * t94 + t75 * t96;
t11 = t118 * t87 + t119 * t88 - t175 * t86 - t222 * t92 + t72 * t94 + t73 * t96;
t10 = -t123 * t62 + t129 * t35 - t185 * t33 - t261 * t59;
t9 = t123 * t61 + t129 * t36 - t185 * t34 - t261 * t57;
t8 = t89 * t220 - t134 * t98 + t240 * qJDD(1) + (-t132 * qJD(1) + t37) * qJD(1) + t178;
t7 = t89 * t219 + t135 * t98 + t203 * qJDD(1) + (t264 * qJD(1) + t187 - t38) * qJD(1) + t212;
t6 = qJD(4) * t192 + t25;
t5 = qJD(4) * t193 - t266;
t1 = [-m(2) * (-g(1) * t139 + g(2) * t141) + (-t123 * t96 + t129 * t88 - t185 * t87 - t261 * t94) * qJD(1) + (t25 + (t248 * t176 + (t267 + t269) * t175) * qJD(4)) * t205 + (t29 + t27) * t257 + (t28 - t184) * t256 + (t5 + t266 + ((t247 + t267) * t176 + t265 * t175) * qJD(4)) * t206 + (t10 + t11) * t207 + (t12 + t9 + t6) * t204 + (-(t271 * qJD(1) + t202 - t23 + t90) * t24 - (-t237 - pkin(1) + (-pkin(2) - pkin(3)) * t173) * t249 + t23 * (t164 + t199) + t24 * (t211 + t244) + (t7 * t181 - t23 * t221) * t175 + ((t181 * t23 + t24 * t253) * t176 + (t23 * (-qJ(2) - t253) + t24 * t181) * t175) * qJD(1) + (-g(2) + t8) * (t175 * t253 + t154 + t201 + t232) + (-g(1) + t7) * (t166 + t64 - t250)) * m(5) + (-(-qJD(1) * t110 + t202 - t51) * t52 + t51 * t196 + t52 * (t161 + t211) + (t51 * t179 * t176 + (t51 * (-rSges(4,2) - qJ(2)) + t52 * (t188 - t195)) * t175) * qJD(1) + (-g(2) + t19) * t262 + (-g(1) + t18) * (t175 * t179 + t166 + t170)) * m(4) + (-(-qJD(1) * t111 + t224 - t84) * t85 + t84 * t164 + t85 * (t225 + t229) + (t84 * (rSges(3,2) * t172 + t208) * t176 + (t84 * (-rSges(3,3) - qJ(2)) + t85 * t208) * t175) * qJD(1) + (-g(1) + t39) * (t208 * t175 + t166 + t227) + (-g(2) + t40) * t100) * m(3) + (-t185 * t94 + t129 * t96 + m(2) * (t139 ^ 2 + t141 ^ 2) + Icges(2,3) + (Icges(3,2) + Icges(4,3)) * t259 + ((Icges(3,1) + Icges(4,1)) * t172 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t173) * t172) * qJDD(1); (-m(3) + t258) * (-g(2) * t176 + t249) + 0.2e1 * (t254 * t8 + t255 * t7) * m(5) + 0.2e1 * (t18 * t255 + t19 * t254) * m(4) + 0.2e1 * (t254 * t40 + t255 * t39) * m(3); t258 * (-g(3) * t173 + (g(1) * t176 + g(2) * t175) * t172) + m(4) * (qJDD(3) * t259 + t18 * t235 + t19 * t236) + m(5) * (-t13 * t173 + t235 * t7 + t236 * t8); -t6 * t222 / 0.2e1 - t175 * (qJD(1) * t11 + qJD(4) * t198 + qJDD(1) * t27 - t134 * t200 + t135 * t16) / 0.2e1 + t192 * t257 + ((-t16 * t175 + t176 * t200) * qJD(1) + t198) * t207 - t5 * t223 / 0.2e1 + t176 * (qJD(1) * t12 + qJD(4) * t197 - qJDD(1) * t184 + t134 * t15 + t135 * t14) / 0.2e1 + t193 * t256 + ((-t14 * t175 - t15 * t176) * qJD(1) + t197) * t204 + qJDD(1) * (-t175 * t29 + t176 * t28) / 0.2e1 + qJD(1) * (-t10 * t175 + t176 * t9 + (-t28 * t175 - t176 * t29) * qJD(1)) / 0.2e1 + ((t118 * t243 + t119 * t242 - t175 * t91) * qJD(1) + (-t118 * t260 + t119 * t183 + t175 * t191) * qJD(4)) * t206 + ((-t116 * t243 + t117 * t242 + t176 * t91) * qJD(1) + (t260 * t116 + t183 * t117 - t191 * t176) * qJD(4)) * t205 - qJD(1) * ((t242 * t129 - t185 * t243) * qJD(1) + (t129 * t183 + t185 * t260) * qJD(4)) / 0.2e1 + ((t23 * t89 - t13 * t65 + t30 * (qJD(1) * t64 - t37) + (qJD(1) * t24 + t7) * t98) * t176 + (t24 * t89 + t13 * t64 + t30 * (qJD(1) * t65 - t38) + (-qJD(1) * t23 + t8) * t98) * t175 - (-t23 * t82 + t24 * t83) * qJD(1) - (t30 * (-t175 * t82 - t176 * t83) + (t24 * t175 + t23 * t176) * t97) * qJD(4) - g(1) * t83 - g(2) * t82 - g(3) * t97) * m(5);];
tau = t1;
