% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:32
% DurationCPUTime: 2.49s
% Computational Cost: add. (4315->280), mult. (4640->372), div. (0->0), fcn. (3476->6), ass. (0->170)
t147 = qJD(1) + qJD(2);
t149 = sin(qJ(4));
t151 = cos(qJ(4));
t148 = qJ(1) + qJ(2);
t143 = sin(t148);
t144 = cos(t148);
t229 = Icges(5,4) * t149;
t184 = Icges(5,2) * t151 + t229;
t70 = Icges(5,6) * t144 + t143 * t184;
t220 = t143 * t151;
t125 = Icges(5,4) * t220;
t221 = t143 * t149;
t72 = Icges(5,1) * t221 + Icges(5,5) * t144 + t125;
t188 = t149 * t72 + t151 * t70;
t71 = -Icges(5,6) * t143 + t144 * t184;
t228 = Icges(5,4) * t151;
t185 = Icges(5,1) * t149 + t228;
t73 = -Icges(5,5) * t143 + t144 * t185;
t40 = t149 * t71 - t151 * t73;
t169 = t184 * t147;
t114 = -Icges(5,2) * t149 + t228;
t256 = -Icges(5,6) * t147 + qJD(4) * t114;
t44 = t143 * t256 + t144 * t169;
t170 = t185 * t147;
t116 = Icges(5,1) * t151 - t229;
t255 = -Icges(5,5) * t147 + qJD(4) * t116;
t46 = t143 * t255 + t144 * t170;
t268 = -qJD(4) * t188 + t147 * t40 - t149 * t44 + t151 * t46;
t216 = t114 + t185;
t217 = -t184 + t116;
t267 = (t149 * t216 - t151 * t217) * t147;
t150 = sin(qJ(1));
t234 = pkin(1) * qJD(1);
t208 = t150 * t234;
t95 = rSges(3,1) * t143 + rSges(3,2) * t144;
t232 = t147 * t95;
t76 = -t208 - t232;
t190 = rSges(5,1) * t149 + rSges(5,2) * t151;
t266 = 0.2e1 * qJD(4);
t237 = rSges(5,1) * t151;
t206 = qJD(4) * t237;
t219 = t144 * t147;
t204 = t143 * t206 + t190 * t219;
t236 = rSges(5,2) * t149;
t205 = qJD(4) * t236;
t110 = qJ(3) * t219;
t128 = qJD(3) * t143;
t218 = t110 + t128;
t75 = -rSges(5,3) * t143 + t144 * t190;
t66 = t147 * t75;
t265 = -t143 * t205 + t204 + t218 - t66;
t222 = t143 * t147;
t215 = rSges(4,2) * t222 + rSges(4,3) * t219;
t94 = rSges(4,2) * t143 + rSges(4,3) * t144;
t264 = -t147 * t94 + t215;
t142 = t144 * pkin(2);
t96 = t143 * qJ(3) + t142;
t177 = -rSges(4,2) * t144 + rSges(4,3) * t143 + t96;
t263 = t147 * t177;
t183 = Icges(5,5) * t149 + Icges(5,6) * t151;
t262 = t183 * t147;
t186 = t149 * t73 + t151 * t71;
t261 = t186 * t144;
t141 = t144 * pkin(6);
t138 = t144 * rSges(5,3);
t74 = rSges(5,1) * t221 + rSges(5,2) * t220 + t138;
t260 = t141 + t74 + t96;
t123 = -t236 + t237;
t212 = qJD(4) * t144;
t91 = t123 * t212;
t259 = t147 * t260 - t91;
t43 = t143 * t169 - t144 * t256;
t45 = t143 * t170 - t144 * t255;
t69 = -Icges(5,3) * t143 + t144 * t183;
t258 = qJD(4) * t40 + t147 * t69 + t149 * t45 + t151 * t43;
t112 = Icges(5,5) * t151 - Icges(5,6) * t149;
t257 = -Icges(5,3) * t147 + qJD(4) * t112;
t100 = t184 * qJD(4);
t101 = t185 * qJD(4);
t254 = (t114 * t149 - t116 * t151) * qJD(4) + t100 * t151 + t101 * t149 + t112 * t147;
t187 = t149 * t70 - t151 * t72;
t68 = Icges(5,3) * t144 + t143 * t183;
t253 = qJD(4) * t187 + t147 * t68 - t149 * t46 - t151 * t44;
t239 = t114 * t144 + t73;
t241 = -t116 * t144 + t71;
t251 = t149 * t241 - t151 * t239;
t240 = -Icges(5,2) * t221 + t125 + t72;
t242 = -t116 * t143 + t70;
t250 = t149 * t242 - t151 * t240;
t249 = -pkin(2) - pkin(6);
t248 = t143 / 0.2e1;
t247 = -t144 / 0.2e1;
t245 = -t147 / 0.2e1;
t244 = pkin(1) * t150;
t243 = pkin(6) * t143;
t152 = cos(qJ(1));
t145 = t152 * pkin(1);
t131 = t144 * qJ(3);
t93 = pkin(2) * t143 - t131;
t178 = -t93 + t94;
t182 = t151 * t114 + t149 * t116;
t80 = t143 * t112;
t53 = t144 * t182 - t80;
t231 = t53 * t147;
t88 = t147 * t93;
t230 = t128 - t88;
t223 = t112 * t144;
t214 = qJD(3) * t147;
t213 = qJD(4) * t143;
t211 = -rSges(5,3) + t249;
t23 = t144 * t68 + t220 * t70 + t221 * t72;
t24 = -t144 * t69 - t220 * t71 - t221 * t73;
t153 = qJD(1) ^ 2;
t210 = t153 * t244;
t209 = t153 * t145;
t207 = t152 * t234;
t201 = -t213 / 0.2e1;
t199 = -t212 / 0.2e1;
t98 = rSges(3,1) * t144 - rSges(3,2) * t143;
t194 = t144 * t214 - t209;
t79 = rSges(3,1) * t219 - rSges(3,2) * t222;
t193 = t128 - t208;
t129 = qJD(3) * t144;
t192 = -t129 + t207;
t90 = t123 * t213;
t179 = t193 + t90;
t29 = (t75 - t93 - t243) * t147 + t179;
t30 = t192 + t259;
t189 = t143 * t29 - t144 * t30;
t180 = t143 * t214 + t147 * (-pkin(2) * t222 + t218) - t210;
t176 = (t205 - t206) * t144;
t173 = (t143 * t24 + t144 * t23) * qJD(4);
t63 = t143 * t68;
t25 = -t144 * t188 + t63;
t26 = -t143 * t69 + t261;
t172 = (t143 * t26 + t144 * t25) * qJD(4);
t36 = (-t143 * t74 - t144 * t75) * qJD(4);
t171 = t129 - t176;
t164 = t143 * t262 - t144 * t257 - t147 * t186;
t163 = t143 * t257 + t144 * t262 + t147 * t188;
t162 = -qJD(4) * t183 + t147 * t182;
t10 = t172 - t231;
t14 = qJD(4) * t186 - t149 * t43 + t151 * t45;
t17 = t143 * t162 + t144 * t254;
t18 = -t143 * t254 + t144 * t162;
t52 = t143 * t182 + t223;
t51 = t52 * t147;
t9 = t51 + t173;
t159 = (t51 + ((-t25 + t63 + t24) * t143 + (t26 - t261 + (-t188 + t69) * t143 + t23) * t144) * qJD(4)) * t201 + (-qJD(4) * t182 + t100 * t149 - t101 * t151) * t147 + (t231 + (t143 ^ 2 * t69 + (-t63 + t24 + (t188 + t69) * t144) * t144) * qJD(4) + t10) * t199 + (t14 + t17 + t9) * t213 / 0.2e1 + (t18 + t268) * t212 / 0.2e1 + (t144 * t53 + (-t187 + t52) * t143) * qJD(4) * t245;
t156 = t143 * t249 + t131 + t75;
t54 = t147 * t178 + t193;
t55 = t192 + t263;
t155 = (-t54 * t142 + (t54 * (-rSges(4,3) - qJ(3)) - t55 * pkin(2)) * t143) * t147;
t154 = (t29 * t211 * t144 + (t29 * (-qJ(3) - t190) + t30 * t211) * t143) * t147;
t146 = t147 ^ 2;
t120 = rSges(4,2) * t219;
t104 = t190 * qJD(4);
t87 = t123 * t144;
t86 = t123 * t143;
t77 = t147 * t98 + t207;
t67 = t147 * t96 - t129;
t62 = -t147 * t79 - t209;
t61 = -t147 * t232 - t210;
t50 = (-rSges(5,3) * t147 - t205) * t143 + t204;
t49 = (t143 * t190 + t138) * t147 + t176;
t32 = (-rSges(4,3) * t222 + t120 - t67) * t147 + t194;
t31 = t147 * t215 + t180;
t20 = -t146 * t141 - t104 * t213 + (-t49 - t67 + t91) * t147 + t194;
t19 = -t146 * t243 + t147 * t50 + (t104 * t144 + t123 * t222) * qJD(4) + t180;
t1 = [t159 + m(3) * (t62 * (-t95 - t244) + t61 * (t145 + t98) + (-t79 - t207 + t77) * t76) + (t20 * (t156 - t244) + t29 * (t171 - t207) + t19 * (t145 + t260) + t154 + (pkin(6) * t222 - t179 - t208 + t265 + t29 + t88) * t30) * m(5) + (t32 * (t178 - t244) + t54 * (t120 - t192) + t31 * (t145 + t177) + t155 + (t110 + t54 + t88 + t264) * t55) * m(4); t159 + (t20 * t156 + t19 * t260 + t154 + (t147 * t243 - t230 + t265 - t90) * t30 + (-t129 + t171 + t259) * t29) * m(5) + (t31 * t177 + t32 * t178 + t155 + (t218 - t230 + t264) * t55 + (t120 + t263) * t54) * m(4) + (t61 * t98 - t62 * t95 - t76 * t79 - t77 * t232 - (-t76 * t98 - t77 * t95) * t147) * m(3); 0.2e1 * (t19 * t247 + t20 * t248) * m(5) + 0.2e1 * (t247 * t31 + t248 * t32) * m(4); t147 * (t268 * t144 + (t147 * t187 + t14) * t143) / 0.2e1 + ((t80 * t212 - t262) * t144 + (-t267 + (t251 * t143 + (-t250 - t223) * t144) * qJD(4)) * t143) * t199 + ((-t213 * t223 - t262) * t143 + (t267 + (t250 * t144 + (-t251 + t80) * t143) * qJD(4)) * t144) * t201 + ((-t149 * t217 - t151 * t216) * t147 + ((t143 * t241 - t144 * t242) * t151 + (t143 * t239 - t144 * t240) * t149) * qJD(4)) * t245 + (t147 * t17 + ((t143 * t163 + t144 * t253 + t147 * t26) * t144 + (t143 * t164 - t144 * t258 - t147 * t25) * t143) * t266) * t248 + (t147 * t18 + ((-t143 * t253 + t144 * t163 + t147 * t24) * t144 + (t143 * t258 + t144 * t164 - t147 * t23) * t143) * t266) * t144 / 0.2e1 - (t9 + t173) * t222 / 0.2e1 + (t10 + t172) * t219 / 0.2e1 + (0.2e1 * t36 * ((-t147 * t74 + t49) * t144 + (-t50 + t66) * t143) - t189 * t104 + ((t147 * t29 - t19) * t144 + (t147 * t30 + t20) * t143) * t123 - (t29 * t87 + t30 * t86) * t147 - (t36 * (-t143 * t86 - t144 * t87) - t189 * t190) * qJD(4)) * m(5);];
tauc = t1(:);
