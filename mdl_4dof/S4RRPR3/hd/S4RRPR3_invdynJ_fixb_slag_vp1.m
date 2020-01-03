% Calculate vector of inverse dynamics joint torques for
% S4RRPR3
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:32
% DurationCPUTime: 2.65s
% Computational Cost: add. (6123->324), mult. (4671->419), div. (0->0), fcn. (3528->8), ass. (0->196)
t159 = sin(qJ(1));
t161 = cos(qJ(1));
t162 = qJD(1) ^ 2;
t177 = (-qJDD(1) * t159 - t161 * t162) * pkin(1);
t268 = t177 - g(1);
t239 = pkin(1) * qJD(1);
t208 = t159 * t239;
t157 = qJ(1) + qJ(2);
t150 = sin(t157);
t151 = cos(t157);
t103 = rSges(3,1) * t150 + rSges(3,2) * t151;
t156 = qJD(1) + qJD(2);
t230 = t103 * t156;
t76 = -t208 - t230;
t158 = sin(qJ(4));
t160 = cos(qJ(4));
t241 = rSges(5,1) * t160;
t131 = -rSges(5,2) * t158 + t241;
t108 = t131 * qJD(4);
t240 = rSges(5,2) * t160;
t129 = rSges(5,1) * t158 + t240;
t155 = qJDD(1) + qJDD(2);
t248 = pkin(2) * t150;
t149 = pkin(7) + t157;
t143 = sin(t149);
t225 = t143 * t158;
t119 = rSges(5,2) * t225;
t144 = cos(t149);
t218 = t144 * rSges(5,3) + t119;
t224 = t143 * t160;
t72 = rSges(5,1) * t224 - t218;
t140 = t144 * pkin(6);
t99 = pkin(3) * t143 - t140;
t194 = -t72 - t99 - t248;
t214 = qJD(4) * t144;
t154 = t156 ^ 2;
t219 = t151 * t154;
t141 = t144 * pkin(3);
t259 = -t143 * pkin(6) - t141;
t206 = qJD(4) * t240;
t222 = t144 * t158;
t209 = rSges(5,2) * t222;
t212 = qJD(4) * t158;
t205 = t156 * t209 + (rSges(5,1) * t212 + t206) * t143;
t221 = t144 * t160;
t120 = rSges(5,1) * t221;
t260 = t143 * rSges(5,3) + t120;
t47 = t156 * t260 - t205;
t213 = qJD(4) * t156;
t89 = -qJDD(4) * t144 + t143 * t213;
t12 = -pkin(2) * t219 - t108 * t214 + t89 * t129 + t177 + (t156 * t259 - t47) * t156 + t194 * t155;
t267 = -g(1) + t12;
t226 = t143 * t156;
t113 = rSges(4,2) * t226;
t223 = t144 * t156;
t97 = rSges(4,1) * t143 + rSges(4,2) * t144;
t266 = -t155 * t97 - t156 * (rSges(4,1) * t223 - t113) + (-t150 * t155 - t219) * pkin(2) + t268;
t142 = t151 * rSges(3,1);
t220 = t150 * t156;
t92 = -rSges(3,2) * t220 + t156 * t142;
t265 = -t103 * t155 - t156 * t92 + t268;
t115 = pkin(6) * t223;
t145 = pkin(2) * t151;
t153 = t161 * pkin(1);
t249 = pkin(1) * t159;
t192 = qJDD(1) * t153 - t162 * t249;
t174 = t155 * t145 - t154 * t248 + t192;
t215 = qJD(4) * t143;
t73 = -t209 + t260;
t234 = t73 - t259;
t181 = rSges(5,3) * t223 + t156 * t119 - t144 * t206;
t203 = t144 * t212;
t46 = (-t156 * t224 - t203) * rSges(5,1) + t181;
t88 = qJDD(4) * t143 + t144 * t213;
t13 = -t108 * t215 - t88 * t129 + (-pkin(3) * t226 + t115 + t46) * t156 + t234 * t155 + t174;
t264 = -g(2) + t13;
t138 = t144 * rSges(4,1);
t98 = -rSges(4,2) * t143 + t138;
t263 = -t97 * t154 + t155 * t98 - g(2) + t174;
t104 = -rSges(3,2) * t150 + t142;
t262 = t104 * t155 - t156 * t230 - g(2) + t192;
t261 = t145 + t98;
t258 = t145 + t234;
t152 = Icges(5,4) * t160;
t184 = -Icges(5,2) * t158 + t152;
t126 = Icges(5,1) * t158 + t152;
t64 = t156 * t72;
t257 = -rSges(5,1) * t203 + t156 * t99 + t115 + t181 + t64;
t123 = Icges(5,5) * t160 - Icges(5,6) * t158;
t122 = Icges(5,5) * t158 + Icges(5,6) * t160;
t170 = Icges(5,3) * t156 - t122 * qJD(4);
t179 = t184 * t144;
t69 = Icges(5,6) * t143 + t179;
t236 = t158 * t69;
t232 = Icges(5,4) * t158;
t127 = Icges(5,1) * t160 - t232;
t180 = t127 * t144;
t71 = Icges(5,5) * t143 + t180;
t186 = -t160 * t71 + t236;
t256 = -t123 * t226 + t170 * t144 + t186 * t156;
t178 = t123 * t144;
t68 = Icges(5,4) * t224 - Icges(5,2) * t225 - Icges(5,6) * t144;
t237 = t158 * t68;
t118 = Icges(5,4) * t225;
t70 = Icges(5,1) * t224 - Icges(5,5) * t144 - t118;
t187 = -t160 * t70 + t237;
t255 = t170 * t143 + (t178 + t187) * t156;
t254 = -t129 * t215 + t156 * t258;
t124 = Icges(5,2) * t160 + t232;
t182 = t158 * t124 - t160 * t126;
t253 = t123 * qJD(4) + t182 * t156;
t66 = Icges(5,5) * t224 - Icges(5,6) * t225 - Icges(5,3) * t144;
t24 = -t187 * t143 - t144 * t66;
t243 = -Icges(5,2) * t224 - t118 + t70;
t245 = t126 * t143 + t68;
t252 = -t243 * t158 - t245 * t160;
t251 = t88 / 0.2e1;
t250 = t89 / 0.2e1;
t247 = -t143 * t66 - t70 * t221;
t67 = Icges(5,3) * t143 + t178;
t246 = t143 * t67 + t71 * t221;
t244 = -t126 * t144 - t69;
t242 = -t124 * t144 + t71;
t228 = t122 * t144;
t49 = -t182 * t143 - t228;
t235 = t49 * t156;
t229 = t122 * t143;
t227 = t123 * t156;
t217 = -t124 + t127;
t216 = t126 + t184;
t207 = t161 * t239;
t31 = t207 + t254;
t211 = t31 * t248;
t210 = pkin(2) * t220;
t204 = t129 * t214;
t202 = -pkin(3) - t241;
t78 = -t97 - t248;
t200 = -t215 / 0.2e1;
t199 = t215 / 0.2e1;
t198 = -t214 / 0.2e1;
t197 = t214 / 0.2e1;
t57 = t71 * t224;
t196 = t144 * t67 - t57;
t195 = -t66 + t236;
t132 = rSges(2,1) * t161 - rSges(2,2) * t159;
t130 = rSges(2,1) * t159 + rSges(2,2) * t161;
t25 = -t69 * t225 - t196;
t191 = t143 * t25 - t144 * t24;
t26 = -t68 * t222 - t247;
t27 = -t69 * t222 + t246;
t190 = t143 * t27 - t144 * t26;
t176 = -t204 - t208;
t30 = t194 * t156 + t176;
t189 = -t31 * t143 - t30 * t144;
t188 = t143 * t72 + t144 * t73;
t37 = t158 * t70 + t160 * t68;
t38 = t158 * t71 + t160 * t69;
t183 = t124 * t160 + t126 * t158;
t175 = -t242 * t158 + t244 * t160;
t173 = (-t216 * t158 + t217 * t160) * t156;
t172 = Icges(5,5) * t156 - qJD(4) * t126;
t171 = Icges(5,6) * t156 - t124 * qJD(4);
t55 = t202 * t143 + t140 + t218 - t248;
t50 = -t182 * t144 + t229;
t48 = t50 * t156;
t10 = t190 * qJD(4) + t48;
t106 = t184 * qJD(4);
t107 = t127 * qJD(4);
t43 = t171 * t143 + t156 * t179;
t45 = t172 * t143 + t156 * t180;
t16 = -t187 * qJD(4) + t158 * t45 + t160 * t43;
t42 = t171 * t144 - t184 * t226;
t44 = -t127 * t226 + t172 * t144;
t17 = -t186 * qJD(4) + t158 * t44 + t160 * t42;
t165 = -t183 * qJD(4) - t106 * t158 + t107 * t160 + t122 * t156;
t20 = t143 * t253 + t165 * t144;
t21 = t165 * t143 - t144 * t253;
t9 = t191 * qJD(4) + t235;
t168 = (t48 + ((t25 - t57 + (t67 + t237) * t144 + t247) * t144 + t246 * t143) * qJD(4)) * t197 + (-t182 * qJD(4) + t106 * t160 + t107 * t158) * t156 + (t38 + t50) * t251 + (t37 + t49) * t250 + (-t235 + ((t195 * t144 - t246 + t27) * t144 + (t195 * t143 + t196 + t26) * t143) * qJD(4) + t9) * t200 + (t17 + t20) * t199 + (t16 + t21 + t10) * t198 + (Icges(4,3) + Icges(3,3) + t183) * t155;
t167 = -t37 * qJD(4) + t156 * t66 - t158 * t43 + t160 * t45;
t166 = -t38 * qJD(4) + t156 * t67 - t158 * t42 + t160 * t44;
t62 = t78 * t156 - t208;
t63 = t156 * t261 + t207;
t164 = (t62 * (-t138 - t145) + t63 * t78) * t156;
t163 = (t30 * (-t120 - t141 - t145) - t211 + (t30 * (-rSges(5,3) - pkin(6)) + t31 * t202) * t143) * t156;
t90 = t156 * t97;
t87 = t129 * t144;
t86 = t129 * t143;
t77 = t104 * t156 + t207;
t36 = t188 * qJD(4) + qJD(3);
t11 = t72 * t88 - t73 * t89 + qJDD(3) + (t143 * t47 + t144 * t46) * qJD(4);
t6 = t166 * t143 - t144 * t256;
t5 = t167 * t143 - t144 * t255;
t4 = t143 * t256 + t166 * t144;
t3 = t143 * t255 + t167 * t144;
t1 = [Icges(2,3) * qJDD(1) + t168 + (t262 * (t104 + t153) + t265 * (-t103 - t249) + (-t92 - t207 + t77) * t76) * m(3) + (g(1) * t130 - g(2) * t132 + (t130 ^ 2 + t132 ^ 2) * qJDD(1)) * m(2) + (t30 * (t205 - t207) + t163 + t264 * (t153 + t258) + t267 * (t55 - t249) + (-t176 + t30 + t210 - t208 + t257) * t31) * m(5) + (-(-t62 - t90 - t210) * t63 + t62 * (t113 - t207) + t164 + t263 * (t153 + t261) + t266 * (t78 - t249)) * m(4); t168 + (t211 * t156 + t163 + t264 * t258 + t267 * t55 + (t204 + t257) * t31 + (t205 + t254) * t30) * m(5) + (t63 * t90 - (-t63 * t248 - t261 * t62) * t156 + t62 * t113 + t164 + t263 * t261 + t266 * t78) * m(4) + (-t76 * t92 - t77 * t230 + (t76 * t156 + t262) * t104 + (t77 * t156 - t265) * t103) * m(3); (t11 - g(3)) * m(5) + (qJDD(3) - g(3)) * m(4); t10 * t223 / 0.2e1 + t143 * (t155 * t50 + t156 * t20 + t26 * t89 + t27 * t88 + (t143 * t4 - t144 * t3) * qJD(4)) / 0.2e1 + t190 * t251 + ((t156 * t27 - t3) * t144 + (t156 * t26 + t4) * t143) * t199 + t9 * t226 / 0.2e1 - t144 * (t155 * t49 + t156 * t21 + t24 * t89 + t25 * t88 + (t143 * t6 - t144 * t5) * qJD(4)) / 0.2e1 + t191 * t250 + ((t156 * t25 - t5) * t144 + (t156 * t24 + t6) * t143) * t198 + t155 * (t143 * t38 - t144 * t37) / 0.2e1 + t156 * ((t156 * t38 - t16) * t144 + (t156 * t37 + t17) * t143) / 0.2e1 + ((-t215 * t228 + t227) * t143 + (t173 + (-t252 * t144 + (t229 + t175) * t143) * qJD(4)) * t144) * t200 + ((-t214 * t229 - t227) * t144 + (t173 + (t175 * t143 + (-t252 + t228) * t144) * qJD(4)) * t143) * t197 - t156 * ((t217 * t158 + t216 * t160) * t156 + ((t242 * t143 - t243 * t144) * t160 + (t244 * t143 + t245 * t144) * t158) * qJD(4)) / 0.2e1 + (t11 * t188 + t36 * ((t46 + t64) * t144 + (-t156 * t73 + t47) * t143) + t189 * t108 + ((-t156 * t31 - t12) * t144 + (t156 * t30 - t13) * t143) * t129 - (t30 * t86 - t31 * t87) * t156 - (t36 * (-t143 * t86 - t144 * t87) + t189 * t131) * qJD(4) + g(1) * t87 + g(2) * t86 - g(3) * t131) * m(5);];
tau = t1;
