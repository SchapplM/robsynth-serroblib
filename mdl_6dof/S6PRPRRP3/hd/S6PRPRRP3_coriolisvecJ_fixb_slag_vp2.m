% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:01
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:00:41
% EndTime: 2018-11-23 15:00:47
% DurationCPUTime: 5.68s
% Computational Cost: add. (5024->418), mult. (13274->557), div. (0->0), fcn. (10041->10), ass. (0->203)
t297 = Ifges(6,4) + Ifges(7,4);
t150 = sin(pkin(11));
t152 = cos(pkin(11));
t155 = sin(qJ(4));
t158 = cos(qJ(4));
t133 = t150 * t155 - t158 * t152;
t128 = t133 * qJD(4);
t134 = t150 * t158 + t152 * t155;
t129 = t134 * qJD(4);
t156 = sin(qJ(2));
t151 = sin(pkin(6));
t213 = qJD(1) * t151;
t202 = t156 * t213;
t300 = pkin(4) * t129 + pkin(9) * t128 - t202;
t159 = cos(qJ(2));
t215 = t151 * t159;
t165 = t133 * t215;
t242 = pkin(8) + qJ(3);
t138 = t242 * t150;
t139 = t242 * t152;
t269 = -t158 * t138 - t139 * t155;
t279 = qJD(1) * t165 - t133 * qJD(3) + qJD(4) * t269;
t298 = Ifges(6,1) + Ifges(7,1);
t282 = Ifges(7,5) + Ifges(6,5);
t296 = Ifges(6,2) + Ifges(7,2);
t281 = Ifges(7,6) + Ifges(6,6);
t299 = Ifges(5,2) / 0.2e1;
t154 = sin(qJ(5));
t157 = cos(qJ(5));
t295 = -t279 * t154 + t300 * t157;
t209 = qJD(5) * t157;
t145 = -pkin(3) * t152 - pkin(2);
t96 = pkin(4) * t133 - pkin(9) * t134 + t145;
t294 = t300 * t154 + t279 * t157 + t96 * t209;
t127 = t134 * qJD(2);
t107 = qJD(4) * t157 - t127 * t154;
t293 = t297 * t107;
t108 = qJD(4) * t154 + t127 * t157;
t292 = t297 * t108;
t291 = t297 * t157;
t290 = t297 * t154;
t211 = qJD(1) * t159;
t201 = t151 * t211;
t174 = qJD(3) - t201;
t119 = qJD(2) * t145 + t174;
t126 = t133 * qJD(2);
t137 = qJD(2) * qJ(3) + t202;
t153 = cos(pkin(6));
t212 = qJD(1) * t153;
t110 = t152 * t137 + t150 * t212;
t231 = pkin(8) * qJD(2);
t100 = t152 * t231 + t110;
t143 = t152 * t212;
t99 = t143 + (-t137 - t231) * t150;
t53 = t100 * t158 + t155 * t99;
t49 = qJD(4) * pkin(9) + t53;
t56 = pkin(4) * t126 - pkin(9) * t127 + t119;
t17 = t154 * t56 + t157 * t49;
t12 = qJ(6) * t107 + t17;
t16 = -t154 * t49 + t157 * t56;
t173 = t17 * t154 + t16 * t157;
t187 = mrSges(7,1) * t154 + mrSges(7,2) * t157;
t189 = mrSges(6,1) * t154 + mrSges(6,2) * t157;
t247 = t157 / 0.2e1;
t250 = -t154 / 0.2e1;
t274 = qJD(5) + t126;
t252 = -t274 / 0.2e1;
t254 = t108 / 0.2e1;
t257 = -t107 / 0.2e1;
t266 = t157 * t298 - t290;
t267 = -t154 * t296 + t291;
t268 = -t154 * t281 + t157 * t282;
t271 = t108 * t298 + t282 * t274 + t293;
t280 = t107 * t296 + t274 * t281 + t292;
t52 = -t155 * t100 + t158 * t99;
t48 = -qJD(4) * pkin(4) - t52;
t31 = -pkin(5) * t107 + qJD(6) + t48;
t11 = -qJ(6) * t108 + t16;
t7 = pkin(5) * t274 + t11;
t263 = t173 * mrSges(6,3) + (t12 * t154 + t7 * t157) * mrSges(7,3) - t187 * t31 - t189 * t48 + t267 * t257 - t266 * t254 + t268 * t252 - t280 * t250 - t271 * t247;
t288 = t127 * Ifges(5,1) / 0.2e1;
t289 = t119 * mrSges(5,2) - t52 * mrSges(5,3) - Ifges(5,4) * t126 + Ifges(5,5) * qJD(4) - t263 + t288;
t287 = t126 * t299;
t169 = qJ(6) * t128 - qJD(6) * t134;
t102 = -t138 * t155 + t139 * t158;
t98 = t157 * t102;
t286 = pkin(5) * t129 + t169 * t157 + (-t98 + (qJ(6) * t134 - t96) * t154) * qJD(5) + t295;
t199 = t134 * t209;
t285 = -qJ(6) * t199 + (-qJD(5) * t102 + t169) * t154 + t294;
t210 = qJD(5) * t154;
t284 = -t102 * t210 + t294;
t46 = t154 * t96 + t98;
t283 = -qJD(5) * t46 + t295;
t166 = t134 * t215;
t103 = qJD(1) * t166;
t77 = qJD(3) * t134 + qJD(4) * t102;
t278 = t77 - t103;
t277 = t154 * t282 + t157 * t281;
t276 = t157 * t296 + t290;
t275 = t154 * t298 + t291;
t255 = -t108 / 0.2e1;
t118 = qJD(2) * t129;
t117 = qJD(2) * t128;
t69 = qJD(5) * t107 - t117 * t157;
t70 = -qJD(5) * t108 + t117 * t154;
t273 = t118 * t281 + t296 * t70 + t297 * t69;
t272 = t118 * t282 + t297 * t70 + t298 * t69;
t214 = t150 ^ 2 + t152 ^ 2;
t270 = mrSges(4,3) * t214;
t132 = (qJD(3) + t201) * qJD(2);
t26 = qJD(4) * t52 - t133 * t132;
t216 = t151 * t156;
t200 = qJD(2) * t216;
t195 = qJD(1) * t200;
t68 = pkin(4) * t118 + pkin(9) * t117 + t195;
t3 = t154 * t68 + t157 * t26 + t56 * t209 - t210 * t49;
t4 = -qJD(5) * t17 - t154 * t26 + t157 * t68;
t265 = -t154 * t4 + t157 * t3;
t204 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t205 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t206 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t261 = t205 * t107 + t206 * t108 + t204 * t274 - t12 * mrSges(7,2) - t17 * mrSges(6,2) - t53 * mrSges(5,3) - Ifges(5,6) * qJD(4) + t287 - Ifges(5,4) * t127 + t119 * mrSges(5,1) + t16 * mrSges(6,1) + t7 * mrSges(7,1) - t281 * t257 - t282 * t255 - (Ifges(7,3) + Ifges(6,3)) * t252;
t260 = m(4) / 0.2e1;
t259 = t69 / 0.2e1;
t258 = t70 / 0.2e1;
t253 = t118 / 0.2e1;
t124 = -t150 * t216 + t152 * t153;
t125 = t150 * t153 + t152 * t216;
t170 = t158 * t124 - t125 * t155;
t27 = qJD(4) * t53 + t132 * t134;
t244 = t27 * t170;
t241 = -qJ(6) - pkin(9);
t41 = mrSges(7,1) * t118 - mrSges(7,3) * t69;
t42 = mrSges(6,1) * t118 - mrSges(6,3) * t69;
t240 = t41 + t42;
t43 = -mrSges(7,2) * t118 + mrSges(7,3) * t70;
t44 = -mrSges(6,2) * t118 + mrSges(6,3) * t70;
t239 = t43 + t44;
t94 = pkin(4) * t127 + pkin(9) * t126;
t25 = t154 * t94 + t157 * t52;
t72 = -mrSges(7,2) * t274 + mrSges(7,3) * t107;
t73 = -mrSges(6,2) * t274 + mrSges(6,3) * t107;
t238 = t72 + t73;
t74 = mrSges(7,1) * t274 - mrSges(7,3) * t108;
t75 = mrSges(6,1) * t274 - mrSges(6,3) * t108;
t237 = t74 + t75;
t230 = t269 * t27;
t223 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t107 + mrSges(6,2) * t108 + mrSges(5,3) * t127;
t220 = qJ(6) * t157;
t219 = t126 * t154;
t218 = t134 * t154;
t59 = -mrSges(7,1) * t107 + mrSges(7,2) * t108;
t207 = t59 + t223;
t203 = t151 ^ 2 * t211;
t28 = -t70 * mrSges(7,1) + t69 * mrSges(7,2);
t24 = -t154 * t52 + t157 * t94;
t197 = qJD(5) * t241;
t79 = t118 * mrSges(5,1) - t117 * mrSges(5,2);
t45 = -t102 * t154 + t157 * t96;
t196 = t214 * t132;
t1 = pkin(5) * t118 - qJ(6) * t69 - qJD(6) * t108 + t4;
t2 = qJ(6) * t70 + qJD(6) * t107 + t3;
t194 = -t1 * t157 - t2 * t154;
t193 = -t3 * t154 - t4 * t157;
t192 = t12 * t157 - t7 * t154;
t190 = mrSges(6,1) * t157 - mrSges(6,2) * t154;
t188 = mrSges(7,1) * t157 - mrSges(7,2) * t154;
t172 = t16 * t154 - t17 * t157;
t171 = -(-t137 * t150 + t143) * t150 + t110 * t152;
t83 = t124 * t155 + t125 * t158;
t61 = -t154 * t83 - t157 * t215;
t168 = t154 * t215 - t157 * t83;
t167 = t171 * t159;
t164 = t4 * mrSges(6,1) + t1 * mrSges(7,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t160 = qJD(2) ^ 2;
t146 = -pkin(5) * t157 - pkin(4);
t141 = t241 * t157;
t140 = t241 * t154;
t136 = -qJD(2) * pkin(2) + t174;
t123 = -qJD(6) * t154 + t157 * t197;
t122 = qJD(6) * t157 + t154 * t197;
t115 = Ifges(6,3) * t118;
t114 = Ifges(7,3) * t118;
t111 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t126;
t93 = mrSges(5,1) * t126 + mrSges(5,2) * t127;
t78 = pkin(5) * t218 - t269;
t67 = Ifges(6,5) * t69;
t66 = Ifges(7,5) * t69;
t65 = Ifges(6,6) * t70;
t64 = Ifges(7,6) * t70;
t51 = qJD(2) * t166 + qJD(4) * t83;
t50 = -qJD(2) * t165 + qJD(4) * t170;
t34 = (-t154 * t128 + t199) * pkin(5) + t77;
t33 = -pkin(5) * t219 + t53;
t32 = -qJ(6) * t218 + t46;
t30 = pkin(5) * t133 - t134 * t220 + t45;
t29 = -mrSges(6,1) * t70 + mrSges(6,2) * t69;
t18 = qJ(6) * t219 + t25;
t15 = qJD(5) * t168 - t154 * t50 + t157 * t200;
t14 = qJD(5) * t61 + t154 * t200 + t157 * t50;
t13 = pkin(5) * t127 + t126 * t220 + t24;
t10 = -pkin(5) * t70 + t27;
t5 = [-t83 * t118 * mrSges(5,3) + t50 * t111 - t239 * t168 + t240 * t61 + t237 * t15 + t238 * t14 - (-t117 * mrSges(5,3) + t28 + t29) * t170 + t207 * t51 + ((-mrSges(3,1) * t160 + (qJD(2) * (-mrSges(4,1) * t152 + mrSges(4,2) * t150) + t93) * qJD(2)) * t156 + (-t79 + (-mrSges(3,2) + t270) * t160) * t159) * t151 + m(7) * (t1 * t61 - t10 * t170 + t12 * t14 + t15 * t7 - t168 * t2 + t31 * t51) + m(6) * (t14 * t17 + t15 * t16 - t168 * t3 + t4 * t61 + t48 * t51 - t244) + m(5) * (t26 * t83 + t50 * t53 - t51 * t52 - t244) + m(4) * (-t124 * t150 + t125 * t152) * t132 + 0.2e1 * (t151 * t167 * t260 + (m(5) * (t119 * t151 - t203) / 0.2e1 + (t136 * t151 - t203) * t260) * t156) * qJD(2); t283 * t75 + (t16 * t283 + t17 * t284 + t278 * t48 + t3 * t46 + t4 * t45 - t230) * m(6) + t284 * t73 + t285 * t72 + (t1 * t30 + t10 * t78 + t2 * t32 + t286 * t7 + (-t103 + t34) * t31 + t285 * t12) * m(7) + t286 * t74 + (t287 + t261) * t129 + (t10 * t187 + mrSges(5,2) * t195 - Ifges(5,1) * t117 - Ifges(5,4) * t118 + (mrSges(5,3) + t189) * t27 + t194 * mrSges(7,3) + t193 * mrSges(6,3) + (mrSges(6,3) * t172 - mrSges(7,3) * t192 + t188 * t31 + t190 * t48 + t276 * t257 + t275 * t255 + t277 * t252 - t280 * t157 / 0.2e1) * qJD(5) + t266 * t259 + t267 * t258 + t268 * t253 + t272 * t247 + (qJD(5) * t271 + t273) * t250) * t134 + (t102 * t26 - t119 * t202 + t145 * t195 - t278 * t52 + t279 * t53 - t230) * m(5) + t279 * t111 + t223 * t77 - t207 * t103 - t93 * t202 + (-pkin(2) * t195 + qJ(3) * t196 + qJD(3) * t171 - (t136 * t156 + t167) * t213) * m(4) + (mrSges(5,1) * t195 - t26 * mrSges(5,3) + t66 / 0.2e1 + t64 / 0.2e1 + t114 / 0.2e1 + t67 / 0.2e1 + t65 / 0.2e1 + t115 / 0.2e1 + Ifges(5,4) * t117 + t205 * t70 + t206 * t69 + (Ifges(5,2) + t204) * t118 + t164) * t133 + (-t102 * t118 + t117 * t269) * mrSges(5,3) - t269 * t29 - (t288 + t289) * t128 + (qJD(2) * t174 * t214 + t196) * mrSges(4,3) + t30 * t41 + t32 * t43 + t45 * t42 + t46 * t44 + t34 * t59 + t78 * t28 + t145 * t79; t126 * t111 + (m(4) + m(5)) * t195 - t160 * t270 - t207 * t127 + (t238 * t274 + t240) * t157 + (-t237 * t274 + t239) * t154 - m(5) * (-t126 * t53 - t127 * t52) - m(4) * t171 * qJD(2) + t79 + (-t127 * t31 + t192 * t274 - t194) * m(7) + (-t127 * t48 - t172 * t274 - t193) * m(6); t275 * t259 + t276 * t258 + t277 * t253 - t223 * t53 - t10 * t188 + (-mrSges(5,1) - t190) * t27 - ((t299 - Ifges(5,1) / 0.2e1) * t127 - t289) * t126 + m(7) * (t1 * t140 + t10 * t146 + t12 * t122 + t123 * t7 - t141 * t2) - t261 * t127 + ((m(7) * t31 + t59) * t154 * pkin(5) - t263) * qJD(5) + (-t13 + t123) * t74 + t265 * mrSges(6,3) + (m(6) * t265 + (-m(6) * t173 - t154 * t73 - t157 * t75) * qJD(5) - t154 * t42 + t157 * t44) * pkin(9) + (-t1 * t154 + t157 * t2) * mrSges(7,3) - m(7) * (t12 * t18 + t13 * t7 + t31 * t33) + (-t18 + t122) * t72 + t272 * t154 / 0.2e1 + t273 * t247 - t26 * mrSges(5,2) - pkin(4) * t29 + (-pkin(4) * t27 - t16 * t24 - t17 * t25 - t48 * t53) * m(6) - t33 * t59 - t25 * t73 - t24 * t75 - t52 * t111 - Ifges(5,5) * t117 - Ifges(5,6) * t118 + t140 * t41 - t141 * t43 + t146 * t28; (-(t11 - t7) * t12 + (-t108 * t31 + t1) * pkin(5)) * m(7) + (t107 * t7 + t108 * t12) * mrSges(7,3) + (t107 * t16 + t108 * t17) * mrSges(6,3) + t164 + t114 + t115 + t66 + t67 + t65 + t64 - t11 * t72 - t16 * t73 + t12 * t74 + t17 * t75 - t31 * (mrSges(7,1) * t108 + mrSges(7,2) * t107) - t48 * (mrSges(6,1) * t108 + mrSges(6,2) * t107) + (-t108 * t59 + t41) * pkin(5) + (t107 * t298 - t292) * t255 + t280 * t254 + (t107 * t282 - t108 * t281) * t252 + (-t108 * t296 + t271 + t293) * t257; -t107 * t72 + t108 * t74 + 0.2e1 * (t10 / 0.2e1 + t12 * t257 + t7 * t254) * m(7) + t28;];
tauc  = t5(:);
