% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPPR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:12
% DurationCPUTime: 4.58s
% Computational Cost: add. (4765->306), mult. (4383->412), div. (0->0), fcn. (3160->8), ass. (0->173)
t141 = qJ(1) + pkin(7);
t136 = sin(t141);
t138 = cos(t141);
t186 = -rSges(4,2) * t138 + rSges(4,3) * t136;
t146 = cos(qJ(1));
t139 = t146 * pkin(1);
t264 = t138 * pkin(2) + qJ(3) * t136;
t268 = t139 + t264;
t274 = t186 + t268;
t147 = qJD(1) ^ 2;
t230 = qJ(4) * t138;
t142 = sin(pkin(8));
t222 = t136 * t142;
t239 = rSges(5,2) * cos(pkin(8));
t261 = -rSges(5,1) * t222 - rSges(5,3) * t138 - t136 * t239;
t275 = t230 - t261 + t268;
t140 = pkin(8) + qJ(5);
t135 = sin(t140);
t137 = cos(t140);
t179 = rSges(6,1) * t135 + rSges(6,2) * t137;
t223 = t136 * t137;
t102 = Icges(6,4) * t223;
t224 = t135 * t136;
t231 = Icges(6,5) * t138;
t61 = Icges(6,1) * t224 + t102 + t231;
t232 = Icges(6,4) * t137;
t170 = Icges(6,1) * t135 + t232;
t62 = -Icges(6,5) * t136 + t138 * t170;
t89 = -Icges(6,2) * t135 + t232;
t71 = t89 * t138;
t155 = t136 * (t62 + t71) - t138 * (-Icges(6,2) * t224 + t102 + t61);
t233 = Icges(6,4) * t135;
t169 = Icges(6,2) * t137 + t233;
t59 = Icges(6,6) * t138 + t136 * t169;
t60 = -Icges(6,6) * t136 + t138 * t169;
t91 = Icges(6,1) * t137 - t233;
t72 = t91 * t136;
t73 = t91 * t138;
t156 = t136 * (t60 - t73) - t138 * (t59 - t72);
t273 = -t156 * t135 + t155 * t137;
t113 = pkin(4) * t222;
t131 = t138 * rSges(6,3);
t63 = rSges(6,1) * t224 + rSges(6,2) * t223 + t131;
t272 = t268 + t113 + t63;
t244 = t89 + t170;
t245 = -t169 + t91;
t271 = (t135 * t244 - t137 * t245) * qJD(1);
t270 = 2 * qJD(5);
t173 = t135 * t62 + t137 * t60;
t265 = t173 * t138;
t187 = rSges(3,1) * t138 - rSges(3,2) * t136;
t263 = t139 + t187;
t168 = Icges(6,5) * t135 + Icges(6,6) * t137;
t58 = -Icges(6,3) * t136 + t138 * t168;
t227 = qJD(1) * t58;
t30 = t135 * t60 - t137 * t62;
t38 = qJD(1) * t59 - qJD(5) * t71;
t40 = -qJD(5) * t73 + (t136 * t170 + t231) * qJD(1);
t260 = qJD(5) * t30 + t135 * t40 + t137 * t38 + t227;
t82 = t169 * qJD(5);
t83 = t170 * qJD(5);
t87 = Icges(6,5) * t137 - Icges(6,6) * t135;
t259 = qJD(1) * t87 + (t135 * t89 - t137 * t91) * qJD(5) + t135 * t83 + t137 * t82;
t174 = t135 * t59 - t137 * t61;
t57 = Icges(6,3) * t138 + t136 * t168;
t228 = qJD(1) * t57;
t211 = qJD(5) * t136;
t39 = qJD(1) * t60 + t211 * t89;
t41 = qJD(1) * t62 + qJD(5) * t72;
t258 = qJD(5) * t174 - t135 * t41 - t137 * t39 + t228;
t256 = t136 / 0.2e1;
t255 = -t138 / 0.2e1;
t253 = -rSges(6,3) - pkin(2);
t145 = sin(qJ(1));
t252 = pkin(1) * t145;
t251 = pkin(2) * t136;
t250 = -qJD(1) / 0.2e1;
t242 = rSges(6,1) * t137;
t238 = rSges(6,2) * t135;
t236 = rSges(4,3) * t138;
t68 = t136 * t87;
t235 = t138 * t87;
t209 = qJD(1) * qJD(3);
t213 = qJD(1) * t136;
t212 = qJD(1) * t138;
t118 = qJ(3) * t212;
t122 = qJD(3) * t136;
t216 = t118 + t122;
t234 = t136 * t209 + qJD(1) * (-pkin(2) * t213 + t216);
t144 = -pkin(6) - qJ(4);
t123 = qJD(3) * t138;
t185 = qJD(4) * t136 - t123;
t210 = qJD(5) * t138;
t96 = -t238 + t242;
t18 = -t96 * t210 + (-t138 * t144 + t272) * qJD(1) + t185;
t229 = qJD(1) * t18;
t226 = qJD(1) * t168;
t225 = qJD(5) * t96;
t221 = t138 * t142;
t172 = t135 * t91 + t137 * t89;
t35 = t138 * t172 - t68;
t220 = t35 * qJD(1);
t219 = qJ(4) + t144;
t198 = t142 * t212;
t203 = t138 * t239;
t218 = rSges(5,1) * t198 + qJD(1) * t203;
t217 = pkin(4) * t198 + t144 * t213;
t215 = rSges(4,2) * t213 + rSges(4,3) * t212;
t214 = qJD(4) * t138 + t122;
t208 = qJD(1) * qJD(4);
t207 = -rSges(5,3) - pkin(2) - qJ(4);
t19 = t138 * t57 + t223 * t59 + t224 * t61;
t20 = -t138 * t58 - t223 * t60 - t224 * t62;
t206 = t147 * t252;
t205 = t147 * t139;
t204 = t179 * t212 + t211 * t242;
t125 = t138 * qJ(3);
t93 = -t125 + t251;
t85 = qJD(1) * t93;
t202 = -t85 + t214;
t201 = qJD(5) * t238;
t200 = t118 + t214;
t193 = -t211 / 0.2e1;
t191 = -t210 / 0.2e1;
t189 = t125 - t252;
t183 = -t251 - t252;
t95 = rSges(3,1) * t136 + rSges(3,2) * t138;
t180 = rSges(5,1) * t142 + t239;
t177 = -qJ(4) * t136 - t252;
t175 = t135 * t61 + t137 * t59;
t166 = -rSges(5,3) * t136 + t138 * t180 + t177;
t76 = t96 * t138;
t164 = (t136 * t20 + t138 * t19) * qJD(5);
t54 = t136 * t57;
t21 = -t138 * t175 + t54;
t22 = -t136 * t58 + t265;
t163 = (t136 * t22 + t138 * t21) * qJD(5);
t129 = t136 * rSges(6,3);
t64 = t138 * t179 - t129;
t162 = pkin(4) * t221 + t136 * t219 + t177 + t64;
t160 = pkin(4) * t142 + t179;
t154 = -qJD(1) * t173 - qJD(5) * t235 + t228;
t153 = qJD(1) * t175 + qJD(5) * t68 + t227;
t152 = t172 * qJD(1) - qJD(5) * t168;
t117 = t138 * t209;
t151 = -0.2e1 * t136 * t208 + t117 + (-t230 - t139) * t147;
t150 = 0.2e1 * t138 * t208 + t147 * t177 + t234;
t84 = t179 * qJD(5);
t79 = t96 * t211;
t75 = t96 * t136;
t74 = qJD(1) * t264 - t123;
t45 = -t205 + t117 + (-qJD(1) * t186 - t74) * qJD(1);
t44 = qJD(1) * t215 - t206 + t234;
t43 = (-rSges(6,3) * qJD(1) - t201) * t136 + t204;
t42 = -qJD(5) * t76 + (t136 * t179 + t131) * qJD(1);
t34 = t136 * t172 + t235;
t33 = t34 * qJD(1);
t32 = qJD(1) * t275 + t185;
t31 = (t166 - t93) * qJD(1) + t214;
t28 = qJD(2) + (-t136 * t63 - t138 * t64) * qJD(5);
t24 = (qJD(1) * t261 - t74) * qJD(1) + t151;
t23 = qJD(1) * (-rSges(5,3) * t213 + t218) + t150;
t17 = t79 + (t162 - t93) * qJD(1) + t214;
t13 = -t84 * t211 + (-qJD(1) * t113 - t42 - t74 + (qJD(1) * t219 + t225) * t138) * qJD(1) + t151;
t12 = t84 * t210 + (t43 + (qJ(4) * qJD(1) + t225) * t136 + t217) * qJD(1) + t150;
t11 = -t136 * t259 + t138 * t152;
t10 = t136 * t152 + t138 * t259;
t9 = qJD(5) * t173 - t135 * t38 + t137 * t40;
t8 = -qJD(5) * t175 - t135 * t39 + t137 * t41;
t7 = (-t136 * t43 + t138 * t42 + (t136 * t64 - t138 * t63) * qJD(1)) * qJD(5);
t6 = t163 - t220;
t5 = t33 + t164;
t1 = [m(3) * ((-t147 * t95 - t206) * t263 + (-t205 + (-0.2e1 * t187 - t139 + t263) * t147) * (-t95 - t252)) + (-t172 * qJD(5) + t135 * t82 - t137 * t83) * qJD(1) + (t33 + ((-t21 + t54 + t20) * t136 + (t22 - t265 + (-t175 + t58) * t136 + t19) * t138) * qJD(5)) * t193 + (t6 + t220 + (t136 ^ 2 * t58 + (-t54 + t20 + (t175 + t58) * t138) * t138) * qJD(5)) * t191 + (t13 * (t136 * t144 + t125 - t129 + t183) - t17 * t185 + t12 * t272 + t18 * (-t136 * t201 + t200 + t204 + t217) + (-t12 * t144 + t13 * t160 + t17 * t225) * t138 + ((-t145 * t18 - t146 * t17) * pkin(1) + t17 * (t144 + t253) * t138 + (t17 * (-qJ(3) - t160) + t18 * t253) * t136) * qJD(1) - (-t17 + t79 + t202) * t18 - t162 * t229) * m(6) + (t24 * (rSges(5,1) * t221 + t189 + t203) + t31 * t123 + t23 * t275 + t32 * (t200 + t218) + (-t31 * qJD(4) + t24 * t207) * t136 + ((-t145 * t32 - t146 * t31) * pkin(1) + t31 * t207 * t138 + (t31 * (-qJ(3) - t180) + t32 * t207) * t136) * qJD(1) - (qJD(1) * t166 + t202 - t31) * t32) * m(5) + (t45 * (t236 + (rSges(4,2) - pkin(2)) * t136 + t189) + t44 * t274 + (-t122 + t215 + t216 + t85 + (-rSges(4,2) * t136 + t183 - t236 + t252) * qJD(1)) * (qJD(1) * t274 - t123)) * m(4) + (t9 + t10 + t5) * t211 / 0.2e1 + (qJD(1) * t30 + t11 + t8) * t210 / 0.2e1 + (t138 * t35 + (-t174 + t34) * t136) * qJD(5) * t250; m(6) * t7; 0.2e1 * (t12 * t255 + t13 * t256) * m(6) + 0.2e1 * (t23 * t255 + t24 * t256) * m(5) + 0.2e1 * (t255 * t44 + t256 * t45) * m(4); m(5) * (t136 * t23 + t138 * t24) + m(6) * (t12 * t136 + t13 * t138); qJD(1) * (t136 * t9 + t138 * t8 + (t136 * t174 + t30 * t138) * qJD(1)) / 0.2e1 + ((t210 * t68 - t226) * t138 + (-t271 + (-t138 * t235 - t273) * qJD(5)) * t136) * t191 + ((-t211 * t235 - t226) * t136 + (t271 + (t136 * t68 + t273) * qJD(5)) * t138) * t193 + ((-t135 * t245 - t137 * t244) * qJD(1) + (t135 * t155 + t137 * t156) * qJD(5)) * t250 + (qJD(1) * t10 + ((t136 * t153 + t138 * t258) * t138 + t136 * (t136 * t154 - t138 * t260) + (-t21 * t136 + t22 * t138) * qJD(1)) * t270) * t256 + (qJD(1) * t11 + (t136 * (t136 * t260 + t138 * t154) + t138 * (-t136 * t258 + t138 * t153) + (-t19 * t136 + t138 * t20) * qJD(1)) * t270) * t138 / 0.2e1 - (t5 + t164) * t213 / 0.2e1 + (t6 + t163) * t212 / 0.2e1 + ((t18 * t84 - t7 * t64 + t28 * (-qJD(1) * t63 + t42) + (qJD(1) * t17 - t12) * t96) * t138 + (-t17 * t84 - t7 * t63 + t28 * (qJD(1) * t64 - t43) + (t13 + t229) * t96) * t136 - (t17 * t76 + t18 * t75) * qJD(1) - (t28 * (-t136 * t75 - t138 * t76) - (t136 * t17 - t138 * t18) * t179) * qJD(5)) * m(6);];
tauc = t1(:);
