% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:43
% EndTime: 2019-12-05 17:04:48
% DurationCPUTime: 2.75s
% Computational Cost: add. (7251->261), mult. (5321->359), div. (0->0), fcn. (4036->8), ass. (0->175)
t139 = sin(qJ(5));
t141 = cos(qJ(5));
t238 = rSges(6,2) * t141;
t113 = rSges(6,1) * t139 + t238;
t138 = qJ(2) + qJ(3);
t134 = qJ(4) + t138;
t128 = cos(t134);
t213 = qJD(5) * t128;
t137 = qJD(2) + qJD(3);
t130 = qJD(4) + t137;
t127 = sin(t134);
t223 = t127 * t139;
t104 = rSges(6,2) * t223;
t222 = t127 * t141;
t67 = rSges(6,1) * t222 - t128 * rSges(6,3) - t104;
t58 = t130 * t67;
t221 = t128 * t130;
t98 = pkin(6) * t221;
t261 = t113 * t213 + t58 - t98;
t140 = sin(qJ(2));
t237 = pkin(2) * qJD(2);
t205 = t140 * t237;
t131 = sin(t138);
t218 = t131 * t137;
t209 = pkin(3) * t218;
t85 = rSges(5,1) * t127 + rSges(5,2) * t128;
t79 = t130 * t85;
t54 = -t205 - t209 - t79;
t260 = t209 + t261;
t132 = cos(t138);
t89 = rSges(4,1) * t131 + rSges(4,2) * t132;
t234 = t137 * t89;
t71 = -t205 - t234;
t259 = 0.2e1 * qJD(5);
t123 = t127 * pkin(6);
t220 = t128 * t139;
t206 = rSges(6,2) * t220;
t219 = t128 * t141;
t105 = rSges(6,1) * t219;
t257 = t127 * rSges(6,3) + t105;
t68 = -t206 + t257;
t258 = t123 + t68;
t28 = -t205 - t260;
t133 = Icges(6,4) * t141;
t180 = -Icges(6,2) * t139 + t133;
t110 = Icges(6,1) * t139 + t133;
t107 = Icges(6,5) * t141 - Icges(6,6) * t139;
t229 = Icges(6,4) * t139;
t108 = Icges(6,2) * t141 + t229;
t179 = t139 * t108 - t141 * t110;
t256 = t107 * qJD(5) + t179 * t130;
t106 = Icges(6,5) * t139 + Icges(6,6) * t141;
t155 = Icges(6,3) * t130 - t106 * qJD(5);
t168 = t180 * t128;
t64 = Icges(6,6) * t127 + t168;
t232 = t139 * t64;
t111 = Icges(6,1) * t141 - t229;
t169 = t111 * t128;
t66 = Icges(6,5) * t127 + t169;
t182 = -t141 * t66 + t232;
t224 = t127 * t130;
t255 = -t107 * t224 + t155 * t128 + t182 * t130;
t167 = t107 * t128;
t63 = Icges(6,4) * t222 - Icges(6,2) * t223 - Icges(6,6) * t128;
t233 = t139 * t63;
t103 = Icges(6,4) * t223;
t65 = Icges(6,1) * t222 - Icges(6,5) * t128 - t103;
t183 = -t141 * t65 + t233;
t254 = t155 * t127 + (t167 + t183) * t130;
t214 = qJD(5) * t127;
t178 = t113 * t214 - t130 * t258;
t61 = Icges(6,5) * t222 - Icges(6,6) * t223 - Icges(6,3) * t128;
t24 = -t183 * t127 - t128 * t61;
t239 = rSges(6,1) * t141;
t142 = cos(qJ(2));
t204 = t142 * t237;
t217 = t132 * t137;
t208 = pkin(3) * t217;
t165 = t204 + t208;
t29 = t165 - t178;
t146 = (-t28 * t105 + (t28 * (-rSges(6,3) - pkin(6)) - t29 * t239) * t127) * t130;
t203 = qJD(5) * t238;
t212 = qJD(5) * t139;
t207 = t130 * t206 + (rSges(6,1) * t212 + t203) * t127;
t253 = t146 + (t207 - t178) * t28;
t241 = -Icges(6,2) * t222 - t103 + t65;
t243 = t110 * t127 + t63;
t252 = -t241 * t139 - t243 * t141;
t249 = t130 / 0.2e1;
t248 = pkin(2) * t140;
t247 = pkin(3) * t131;
t246 = pkin(3) * t137 ^ 2;
t124 = t128 * pkin(6);
t135 = t142 * pkin(2);
t245 = -t127 * t61 - t65 * t219;
t62 = Icges(6,3) * t127 + t167;
t244 = t127 * t62 + t66 * t219;
t242 = -t110 * t128 - t64;
t240 = -t108 * t128 + t66;
t86 = t128 * rSges(5,1) - rSges(5,2) * t127;
t235 = t130 * t86;
t226 = t106 * t128;
t47 = -t179 * t127 - t226;
t231 = t47 * t130;
t227 = t106 * t127;
t225 = t107 * t130;
t216 = -t108 + t111;
t215 = t110 + t180;
t143 = qJD(2) ^ 2;
t211 = t143 * t248;
t210 = t143 * t135;
t202 = t128 * t212;
t198 = -t214 / 0.2e1;
t195 = t213 / 0.2e1;
t51 = t66 * t222;
t193 = t128 * t62 - t51;
t177 = rSges(6,3) * t221 + t130 * t104 - t128 * t203;
t42 = (-t130 * t222 - t202) * rSges(6,1) + t177;
t192 = t42 + t58;
t43 = t130 * t257 - t207;
t191 = -t130 * t68 + t43;
t190 = -t61 + t232;
t90 = t132 * rSges(4,1) - rSges(4,2) * t131;
t70 = rSges(5,1) * t221 - rSges(5,2) * t224;
t83 = rSges(4,1) * t217 - rSges(4,2) * t218;
t126 = pkin(3) * t132;
t187 = t126 + t86;
t184 = -rSges(6,2) * t139 + t239;
t44 = t139 * t65 + t141 * t63;
t45 = t139 * t66 + t141 * t64;
t176 = t124 - t67;
t25 = -t64 * t223 - t193;
t173 = (t127 * t25 - t128 * t24) * qJD(5);
t26 = -t63 * t220 - t245;
t27 = -t64 * t220 + t244;
t172 = (t127 * t27 - t128 * t26) * qJD(5);
t171 = -t131 * t246 - t211;
t170 = -t132 * t246 - t210;
t164 = -t85 - t247;
t163 = t126 + t258;
t162 = -t70 - t208;
t161 = -t240 * t139 + t242 * t141;
t159 = t176 - t247;
t158 = (-t215 * t139 + t216 * t141) * t130;
t157 = Icges(6,5) * t130 - qJD(5) * t110;
t156 = Icges(6,6) * t130 - t108 * qJD(5);
t48 = -t179 * t128 + t227;
t46 = t48 * t130;
t10 = t46 + t172;
t14 = -t183 * qJD(5) + t139 * (t157 * t127 + t130 * t169) + t141 * (t156 * t127 + t130 * t168);
t15 = -t182 * qJD(5) + t139 * (-t111 * t224 + t157 * t128) + t141 * (t156 * t128 - t180 * t224);
t92 = t180 * qJD(5);
t93 = t111 * qJD(5);
t145 = t106 * t130 - t139 * t92 + t141 * t93 + (-t108 * t141 - t110 * t139) * qJD(5);
t18 = t127 * t256 + t145 * t128;
t19 = t145 * t127 - t128 * t256;
t9 = t173 + t231;
t152 = (t46 + ((t25 - t51 + (t62 + t233) * t128 + t245) * t128 + t244 * t127) * qJD(5)) * t195 + (-t179 * qJD(5) + t139 * t93 + t141 * t92) * t130 + (t9 - t231 + ((t190 * t128 - t244 + t27) * t128 + (t190 * t127 + t193 + t26) * t127) * qJD(5)) * t198 + (t15 + t18) * t214 / 0.2e1 - (t10 + t14 + t19) * t213 / 0.2e1 + ((t44 + t47) * t127 + (t45 + t48) * t128) * qJD(5) * t249;
t151 = -rSges(6,1) * t202 + t177 + t98;
t144 = t151 - t209;
t129 = t130 ^ 2;
t96 = t184 * qJD(5);
t81 = t113 * t128;
t80 = t113 * t127;
t72 = t137 * t90 + t204;
t60 = -t137 * t83 - t210;
t59 = -t137 * t234 - t211;
t55 = t165 + t235;
t50 = -t130 * t70 + t170;
t49 = -t130 * t79 + t171;
t32 = qJD(1) + (t127 * t67 + t128 * t68) * qJD(5);
t23 = -t129 * t123 - t130 * t43 + (t113 * t224 - t128 * t96) * qJD(5) + t170;
t22 = t129 * t124 + t130 * t42 + (-t113 * t221 - t127 * t96) * qJD(5) + t171;
t11 = (t191 * t127 + t192 * t128) * qJD(5);
t1 = [m(6) * t11; t152 + m(4) * (t60 * (-t89 - t248) + t59 * (t135 + t90) + (-t83 - t204 + t72) * t71) + m(6) * (t23 * (t159 - t248) + t28 * (-t165 + t207) + t22 * (t135 + t163) + t29 * (t144 - t205) + t146) + m(5) * (t50 * (t164 - t248) + t49 * (t135 + t187) + (t162 - t204 + t55) * t54); t152 + (t23 * t159 + t22 * t163 + (t144 + t260) * t29 + t253) * m(6) + (t50 * t164 + t49 * t187 + (t162 + t208 + t235) * t54) * m(5) + (t59 * t90 - t60 * t89 - t71 * t83 - t72 * t234 - (-t71 * t90 - t72 * t89) * t137) * m(4); t152 + (t22 * t258 + t23 * t176 + (t151 + t261) * t29 + t253) * m(6) + (-(-t54 * t86 - t55 * t85) * t130 + t49 * t86 - t50 * t85 - t54 * t70 - t55 * t79) * m(5); ((t130 * t45 - t14) * t128 + (t130 * t44 + t15) * t127) * t249 + ((-t214 * t226 + t225) * t127 + (t158 + (-t252 * t128 + (t227 + t161) * t127) * qJD(5)) * t128) * t198 + ((-t213 * t227 - t225) * t128 + (t158 + (t161 * t127 + (-t252 + t226) * t128) * qJD(5)) * t127) * t195 - t130 * ((t216 * t139 + t215 * t141) * t130 + ((t240 * t127 - t241 * t128) * t141 + (t242 * t127 + t243 * t128) * t139) * qJD(5)) / 0.2e1 + (t130 * t18 + ((-t127 * t254 + t130 * t27) * t128 + (t127 * t255 + t130 * t26) * t127) * t259) * t127 / 0.2e1 - (t130 * t19 + ((t128 * t254 + t130 * t25) * t128 + (-t128 * t255 + t130 * t24) * t127) * t259) * t128 / 0.2e1 + (t9 + t173) * t224 / 0.2e1 + (t10 + t172) * t221 / 0.2e1 + ((t11 * t68 + t32 * t192 - t28 * t96) * t128 + (t11 * t67 + t32 * t191 - t29 * t96) * t127 + ((-t130 * t29 - t23) * t128 + (t130 * t28 - t22) * t127) * t113 - (t28 * t80 - t29 * t81) * t130 - (t32 * (-t127 * t80 - t128 * t81) + (-t127 * t29 - t128 * t28) * t184) * qJD(5)) * m(6);];
tauc = t1(:);
