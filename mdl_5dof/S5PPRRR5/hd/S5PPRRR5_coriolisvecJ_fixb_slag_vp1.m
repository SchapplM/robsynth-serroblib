% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:38
% EndTime: 2019-12-31 17:35:43
% DurationCPUTime: 3.16s
% Computational Cost: add. (7727->259), mult. (8381->375), div. (0->0), fcn. (8858->8), ass. (0->164)
t231 = qJ(3) + qJ(4);
t214 = sin(t231);
t215 = cos(t231);
t245 = sin(pkin(8));
t246 = cos(pkin(8));
t119 = -t214 * t245 - t215 * t246;
t120 = t214 * t246 - t215 * t245;
t149 = sin(qJ(5));
t151 = cos(qJ(5));
t243 = Icges(6,4) * t151;
t186 = -Icges(6,2) * t149 + t243;
t64 = -Icges(6,6) * t120 + t119 * t186;
t244 = Icges(6,4) * t149;
t188 = Icges(6,1) * t151 - t244;
t67 = -Icges(6,5) * t120 + t119 * t188;
t189 = -t149 * t64 + t151 * t67;
t184 = Icges(6,5) * t151 - Icges(6,6) * t149;
t61 = -Icges(6,3) * t120 + t119 * t184;
t27 = t119 * t189 - t120 * t61;
t62 = -Icges(6,3) * t119 - t120 * t184;
t274 = t120 * t62;
t190 = -t149 * t67 - t151 * t64;
t272 = 2 * qJD(5);
t152 = cos(qJ(3));
t271 = pkin(3) * t152;
t220 = qJD(3) + qJD(4);
t87 = rSges(5,1) * t119 + rSges(5,2) * t120;
t270 = t220 * t87;
t232 = t120 * t151;
t233 = t120 * t149;
t71 = -rSges(6,1) * t232 + rSges(6,2) * t233 - rSges(6,3) * t119;
t167 = -pkin(4) * t120 - pkin(7) * t119 + t71;
t183 = Icges(6,5) * t149 + Icges(6,6) * t151;
t77 = t183 * t119;
t76 = t183 * t120;
t174 = t220 * t214;
t175 = t220 * t215;
t105 = t174 * t245 + t175 * t246;
t222 = qJD(5) * t149;
t176 = -t105 * t151 + t120 * t222;
t104 = t174 * t246 - t175 * t245;
t221 = qJD(5) * t151;
t179 = t104 * t149 + t119 * t221;
t225 = qJD(5) * t119;
t211 = t225 / 0.2e1;
t185 = Icges(6,2) * t151 + t244;
t187 = Icges(6,1) * t149 + t243;
t182 = -t149 * t185 + t151 * t187;
t47 = t119 * t182 - t76;
t269 = t220 * t47;
t88 = -rSges(5,1) * t120 + rSges(5,2) * t119;
t201 = t220 * t88;
t74 = -rSges(5,1) * t104 - rSges(5,2) * t105;
t268 = t201 - t74;
t177 = t105 * t149 + t120 * t221;
t42 = rSges(6,1) * t176 + rSges(6,2) * t177 - rSges(6,3) * t104;
t153 = -pkin(4) * t105 - t104 * pkin(7) + t42;
t116 = t120 * pkin(7);
t267 = pkin(4) * t119 - t116;
t102 = t105 * pkin(7);
t195 = rSges(6,1) * t149 + rSges(6,2) * t151;
t223 = qJD(5) * t195;
t180 = -t119 * t223 - t167 * t220;
t219 = t119 * t222;
t181 = -rSges(6,1) * t219 - rSges(6,2) * t179 - t105 * rSges(6,3);
t249 = rSges(6,1) * t151;
t216 = pkin(4) + t249;
t266 = t104 * t216 - t102 - t180 + t181;
t205 = t245 * t152;
t150 = sin(qJ(3));
t207 = t246 * t150;
t168 = t205 - t207;
t122 = t168 * qJD(3);
t199 = pkin(3) * t207;
t171 = (-pkin(3) * t205 + t199) * qJD(3);
t265 = pkin(3) * t122 + t171;
t147 = qJD(2) * t245;
t165 = t245 * t271 - t199;
t230 = qJD(3) * t165 + t147;
t75 = -rSges(5,1) * t105 + rSges(5,2) * t104;
t264 = (t75 - t270) * (t201 + t230);
t224 = qJD(5) * t120;
t28 = -t180 + t230;
t236 = t119 * t149;
t248 = t120 * rSges(6,3);
t204 = -rSges(6,2) * t236 - t248;
t235 = t119 * t151;
t72 = -rSges(6,1) * t235 - t204;
t193 = t119 * t72 + t120 * t71;
t30 = qJD(5) * t193 + qJD(1);
t196 = -rSges(6,2) * t149 + t249;
t70 = t119 * t196 - t248;
t96 = t120 * t223;
t263 = (t153 - t96 - (t70 + t267) * t220) * t28 - t30 * (t70 + t72) * t224;
t68 = -Icges(6,5) * t119 - t120 * t188;
t252 = t120 * t185 + t68;
t65 = -Icges(6,6) * t119 - t120 * t186;
t254 = -t120 * t187 + t65;
t262 = t149 * t252 + t151 * t254;
t256 = t119 * t62 + t232 * t68;
t255 = -t235 * t68 + t274;
t253 = -t119 * t187 - t64;
t251 = t119 * t185 - t67;
t250 = m(4) * qJD(3);
t247 = t149 * t65;
t238 = t104 * t151;
t229 = -t185 + t188;
t228 = -t186 - t187;
t227 = qJD(3) * t171;
t210 = -t224 / 0.2e1;
t206 = t245 * t150;
t203 = qJD(2) * t246;
t200 = t184 * t220;
t125 = -t152 * t246 - t206;
t197 = rSges(4,1) * t125 - rSges(4,2) * t168;
t164 = -pkin(3) * t206 - t246 * t271;
t170 = qJD(3) * t164 - t203;
t29 = t96 + t170 + (-t72 + t267) * t220;
t194 = -t119 * t28 - t120 * t29;
t192 = t149 * t68 + t151 * t65;
t191 = -t151 * t68 + t247;
t25 = t119 * t61 + t232 * t67 - t233 * t64;
t178 = t219 - t238;
t24 = t233 * t65 - t256;
t173 = (-t119 * t24 + t120 * t25) * qJD(5);
t26 = t236 * t65 + t255;
t172 = (-t119 * t26 + t120 * t27) * qJD(5);
t169 = t149 * t251 + t151 * t253;
t123 = t125 * qJD(3);
t163 = t119 * t216 - t116 + t204;
t41 = -rSges(6,1) * t238 - t181;
t162 = t104 * t72 + t105 * t71 + t119 * t41 + t120 * t42;
t158 = (t149 * t228 + t151 * t229) * t220;
t10 = t172 - t269;
t128 = t186 * qJD(5);
t129 = t188 * qJD(5);
t15 = qJD(5) * t191 - t149 * (Icges(6,1) * t176 + Icges(6,4) * t177 - Icges(6,5) * t104) - t151 * (Icges(6,4) * t176 + Icges(6,2) * t177 - Icges(6,6) * t104);
t16 = qJD(5) * t189 - t149 * (Icges(6,1) * t178 + Icges(6,4) * t179 + Icges(6,5) * t105) - t151 * (Icges(6,4) * t178 + Icges(6,2) * t179 + Icges(6,6) * t105);
t127 = t184 * qJD(5);
t154 = -t128 * t149 + t129 * t151 + (-t149 * t187 - t151 * t185) * qJD(5);
t22 = t104 * t182 - t105 * t183 + t119 * t154 - t120 * t127;
t23 = t104 * t183 + t105 * t182 + t119 * t127 + t120 * t154;
t46 = t120 * t182 + t77;
t45 = t46 * t220;
t9 = -t45 + t173;
t155 = -t220 * (-t128 * t151 - t129 * t149) + t9 * t224 / 0.2e1 + (-t45 + t16 + t22) * t210 + (t15 + t23 + t269 + t10) * t211 + (t182 * t220 + (-t192 + t46) * t104 / 0.2e1 - (-t190 + t47) * t105 / 0.2e1 + ((t25 - t26 + t255) * t120 + t256 * t119) * t210 + ((-t24 + (t61 + t247) * t120 - t256) * t120 + (-t25 + (t191 + t61) * t119 + t274) * t119) * t211) * qJD(5);
t130 = t196 * qJD(5);
t112 = qJD(3) ^ 2 * t125 * pkin(3);
t94 = rSges(4,1) * t168 + rSges(4,2) * t125;
t93 = rSges(4,1) * t123 - rSges(4,2) * t122;
t92 = rSges(4,1) * t122 + rSges(4,2) * t123;
t86 = qJD(3) * t197 - t203;
t85 = qJD(3) * t94 + t147;
t83 = t195 * t119;
t82 = t195 * t120;
t57 = t170 + t270;
t55 = -t220 * t74 + t227;
t54 = t220 * t75 + t112;
t36 = Icges(6,5) * t176 + Icges(6,6) * t177 - Icges(6,3) * t104;
t35 = Icges(6,5) * t178 + Icges(6,6) * t179 + Icges(6,3) * t105;
t21 = t227 + (t105 * t195 + t120 * t130) * qJD(5) + (pkin(4) * t104 - t102 - t41) * t220;
t20 = t112 + (t104 * t195 + t119 * t130) * qJD(5) + t153 * t220;
t12 = t162 * qJD(5);
t1 = [m(6) * t12; (t245 * t93 + t246 * t92) * t250 + m(5) * (t245 * t54 - t246 * t55) + m(6) * (t20 * t245 - t21 * t246); t155 - (t197 * t85 - t86 * t94) * t250 + m(4) * (t85 * t93 - t86 * t92 + (-t197 * t92 + t93 * t94) * qJD(3)) + (t20 * (t165 + t167) + t21 * (t163 + t164) + (t265 + t266) * t29 + t263) * m(6) + (t54 * (t165 + t88) + t55 * (t164 + t87) + (t265 + t268) * t57 + t264) * m(5); t155 + (t21 * t163 + t20 * t167 + t266 * t29 + t263) * m(6) + (t268 * t57 + t54 * t88 + t55 * t87 + t264) * m(5); -t220 * (t104 * t192 - t105 * t190 - t15 * t119 + t16 * t120) / 0.2e1 + ((t77 * t224 + t200) * t120 + (-t158 + (-t262 * t119 + (-t76 + t169) * t120) * qJD(5)) * t119) * t210 + ((t76 * t225 - t200) * t119 + (-t158 + (t169 * t120 + (-t262 - t77) * t119) * qJD(5)) * t120) * t211 + t220 * (-(t229 * t149 - t228 * t151) * t220 + ((t119 * t252 - t120 * t251) * t151 + (-t119 * t254 + t120 * t253) * t149) * qJD(5)) / 0.2e1 - (-t220 * t23 + (-t104 * t24 + t105 * t25 - t119 * (-t104 * t62 + t105 * t191 - t119 * t36) + t120 * (t104 * t61 + t105 * t189 - t119 * t35)) * t272) * t119 / 0.2e1 + (-t22 * t220 + (-t104 * t26 + t105 * t27 - t119 * (t104 * t191 + t105 * t62 + t120 * t36) + t120 * (t104 * t189 - t105 * t61 + t120 * t35)) * t272) * t120 / 0.2e1 - (t9 + t173) * t104 / 0.2e1 + (t10 + t172) * t105 / 0.2e1 + (t12 * t193 + t30 * t162 - t194 * t130 - (-t104 * t28 - t105 * t29 - t119 * t20 - t120 * t21) * t195 - (t30 * (t119 * t83 + t120 * t82) - t194 * t196) * qJD(5) - (t28 * t82 - t29 * t83) * t220) * m(6);];
tauc = t1(:);
