% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2018-11-23 15:54
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:53:51
% EndTime: 2018-11-23 15:53:56
% DurationCPUTime: 5.86s
% Computational Cost: add. (6168->503), mult. (16649->668), div. (0->0), fcn. (12370->8), ass. (0->219)
t313 = Ifges(5,1) + Ifges(6,1);
t192 = sin(pkin(9));
t265 = pkin(7) + qJ(2);
t182 = t265 * t192;
t194 = cos(pkin(9));
t184 = t265 * t194;
t196 = sin(qJ(3));
t270 = cos(qJ(3));
t135 = -t196 * t182 + t184 * t270;
t176 = t192 * t270 + t196 * t194;
t205 = t176 * qJD(2);
t103 = qJD(3) * t135 + t205;
t193 = cos(pkin(10));
t315 = t193 * t176 * qJD(5) - t103;
t314 = Ifges(4,1) / 0.2e1;
t312 = Ifges(6,4) + Ifges(5,5);
t311 = Ifges(6,5) - Ifges(5,4);
t191 = sin(pkin(10));
t264 = -pkin(8) + qJ(4);
t181 = t264 * t191;
t183 = t264 * t193;
t195 = sin(qJ(6));
t197 = cos(qJ(6));
t132 = t181 * t197 - t183 * t195;
t210 = t191 * t195 + t193 * t197;
t178 = qJD(1) * t184;
t161 = t196 * t178;
t177 = qJD(1) * t182;
t162 = t270 * t177;
t128 = -t161 - t162;
t123 = t191 * t128;
t229 = t270 * t194;
t186 = qJD(1) * t229;
t243 = t196 * t192;
t163 = qJD(1) * t243 - t186;
t164 = t176 * qJD(1);
t126 = pkin(3) * t164 + qJ(4) * t163;
t282 = pkin(4) + pkin(5);
t29 = t123 + (pkin(8) * t163 - t126) * t193 - t282 * t164;
t268 = pkin(8) * t191;
t63 = t191 * t126 + t193 * t128;
t49 = t164 * qJ(5) + t63;
t38 = -t163 * t268 + t49;
t310 = qJD(4) * t210 + qJD(6) * t132 - t195 * t29 - t197 * t38;
t134 = t181 * t195 + t183 * t197;
t175 = t191 * t197 - t193 * t195;
t309 = qJD(4) * t175 - qJD(6) * t134 + t195 * t38 - t197 * t29;
t258 = Ifges(6,5) * t191;
t262 = Ifges(5,4) * t191;
t308 = t193 * t313 + t258 - t262;
t231 = t282 * t191;
t247 = qJ(5) * t193;
t307 = t231 - t247;
t306 = Ifges(5,3) + Ifges(6,2) + Ifges(4,2);
t137 = -t193 * qJD(3) + t164 * t191;
t209 = qJD(3) * t191 + t193 * t164;
t305 = t137 * t195 + t197 * t209;
t83 = t137 * t197 - t195 * t209;
t206 = t229 - t243;
t167 = t206 * qJD(3);
t152 = qJD(1) * t167;
t41 = qJD(6) * t83 + t152 * t210;
t288 = t41 / 0.2e1;
t42 = -qJD(6) * t305 + t152 * t175;
t287 = t42 / 0.2e1;
t281 = -t137 / 0.2e1;
t280 = t137 / 0.2e1;
t304 = t152 / 0.2e1;
t168 = t176 * qJD(3);
t153 = qJD(1) * t168;
t278 = -t153 / 0.2e1;
t303 = t153 / 0.2e1;
t300 = t164 * t314;
t108 = t175 * t163;
t166 = t175 * qJD(6);
t236 = t108 - t166;
t109 = t210 * t163;
t165 = t210 * qJD(6);
t237 = t109 - t165;
t297 = -t270 * t182 - t196 * t184;
t160 = qJD(6) - t163;
t234 = qJD(1) * qJD(2);
t227 = t192 * t234;
t238 = qJD(2) * t186 - qJD(3) * t162;
t87 = -t196 * t227 + (qJD(4) - t161) * qJD(3) + t238;
t79 = t191 * t87;
t81 = pkin(3) * t153 - qJ(4) * t152 - qJD(4) * t164;
t14 = t79 + (-pkin(8) * t152 - t81) * t193 - t282 * t153;
t33 = t191 * t81 + t193 * t87;
t16 = t153 * qJ(5) + t163 * qJD(5) + t33;
t246 = t152 * t191;
t15 = pkin(8) * t246 + t16;
t129 = -t196 * t177 + t178 * t270;
t125 = qJD(3) * qJ(4) + t129;
t230 = -pkin(2) * t194 - pkin(1);
t179 = qJD(1) * t230 + qJD(2);
t99 = pkin(3) * t163 - qJ(4) * t164 + t179;
t53 = -t191 * t125 + t193 * t99;
t220 = qJD(5) - t53;
t19 = -pkin(8) * t209 - t163 * t282 + t220;
t54 = t193 * t125 + t191 * t99;
t44 = t163 * qJ(5) + t54;
t26 = pkin(8) * t137 + t44;
t5 = t19 * t197 - t195 * t26;
t1 = qJD(6) * t5 + t14 * t195 + t15 * t197;
t6 = t19 * t195 + t197 * t26;
t2 = -qJD(6) * t6 + t14 * t197 - t15 * t195;
t295 = -t2 * mrSges(7,1) + t1 * mrSges(7,2) - Ifges(7,5) * t41 - Ifges(7,6) * t42;
t294 = (m(3) * qJ(2) + mrSges(3,3)) * (t192 ^ 2 + t194 ^ 2);
t257 = Ifges(6,5) * t193;
t214 = t191 * Ifges(6,3) + t257;
t261 = Ifges(5,4) * t193;
t215 = -t191 * Ifges(5,2) + t261;
t232 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t245 = t152 * t193;
t92 = qJD(1) * t205 + qJD(3) * t129;
t201 = -qJ(5) * t245 - qJD(5) * t209 + t92;
t30 = pkin(4) * t246 + t201;
t293 = t232 * t153 + t30 * mrSges(6,1) + t92 * mrSges(5,1) + Ifges(6,6) * t303 + t214 * t304 + Ifges(5,6) * t278 - t152 * t215 / 0.2e1 - t16 * mrSges(6,2) - t33 * mrSges(5,3);
t233 = -Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t43 = -pkin(4) * t163 + t220;
t292 = -t137 * t232 - t179 * mrSges(4,1) - t53 * mrSges(5,1) + t43 * mrSges(6,1) + t5 * mrSges(7,1) + t54 * mrSges(5,2) - t6 * mrSges(7,2) + t129 * mrSges(4,3) - t44 * mrSges(6,3) + t164 * Ifges(4,4) + t305 * Ifges(7,5) + Ifges(4,6) * qJD(3) + Ifges(5,6) * t280 + Ifges(6,6) * t281 + t83 * Ifges(7,6) + t160 * Ifges(7,3) - t306 * t163 / 0.2e1 + (t233 - t312 / 0.2e1) * t209;
t208 = qJD(3) * pkin(3) - qJD(4) + t128;
t271 = t193 / 0.2e1;
t272 = t191 / 0.2e1;
t273 = -t191 / 0.2e1;
t204 = qJ(5) * t209 + t208;
t52 = pkin(4) * t137 - t204;
t291 = t179 * mrSges(4,2) + (-t44 * t191 + t43 * t193) * mrSges(6,2) + (-t54 * t191 - t53 * t193) * mrSges(5,3) - Ifges(4,4) * t163 + Ifges(4,5) * qJD(3) + t300 - t128 * mrSges(4,3) + (Ifges(6,5) * t209 + t163 * Ifges(6,6) + t137 * Ifges(6,3)) * t272 + (Ifges(5,4) * t209 - t137 * Ifges(5,2) + t163 * Ifges(5,6)) * t273 + (t137 * t311 + t163 * t312 + t209 * t313) * t271 + t214 * t280 + t215 * t281 + t52 * (mrSges(6,1) * t191 - mrSges(6,3) * t193) + t308 * t209 / 0.2e1 - t208 * (mrSges(5,1) * t191 + mrSges(5,2) * t193);
t290 = Ifges(7,4) * t288 + Ifges(7,2) * t287 + Ifges(7,6) * t278;
t289 = Ifges(7,1) * t288 + Ifges(7,4) * t287 + Ifges(7,5) * t278;
t286 = -t83 / 0.2e1;
t285 = t83 / 0.2e1;
t284 = -t305 / 0.2e1;
t283 = t305 / 0.2e1;
t276 = -t160 / 0.2e1;
t275 = t160 / 0.2e1;
t269 = Ifges(7,4) * t305;
t35 = -mrSges(7,1) * t83 + mrSges(7,2) * t305;
t89 = mrSges(6,1) * t137 - mrSges(6,3) * t209;
t263 = t35 - t89;
t102 = t206 * qJD(2) + qJD(3) * t297;
t93 = pkin(3) * t168 - qJ(4) * t167 - qJD(4) * t176;
t47 = t193 * t102 + t191 * t93;
t260 = Ifges(6,4) * t193;
t259 = Ifges(5,5) * t193;
t256 = Ifges(5,6) * t191;
t255 = Ifges(6,6) * t191;
t254 = t297 * t92;
t250 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t137 - mrSges(5,2) * t209 - t164 * mrSges(4,3);
t104 = -mrSges(6,2) * t137 + mrSges(6,3) * t163;
t105 = -mrSges(5,2) * t163 - mrSges(5,3) * t137;
t242 = t104 + t105;
t106 = mrSges(5,1) * t163 - mrSges(5,3) * t209;
t107 = -mrSges(6,1) * t163 + mrSges(6,2) * t209;
t241 = t106 - t107;
t110 = -mrSges(5,2) * t153 - mrSges(5,3) * t246;
t113 = -mrSges(6,2) * t246 + mrSges(6,3) * t153;
t240 = t110 + t113;
t111 = mrSges(5,1) * t153 - mrSges(5,3) * t245;
t112 = -t153 * mrSges(6,1) + mrSges(6,2) * t245;
t239 = t111 - t112;
t127 = -pkin(3) * t206 - qJ(4) * t176 + t230;
t72 = t191 * t127 + t193 * t135;
t98 = mrSges(5,1) * t246 + mrSges(5,2) * t245;
t57 = -qJ(5) * t206 + t72;
t226 = qJ(5) * t191 + pkin(3);
t32 = t193 * t81 - t79;
t95 = t191 * t102;
t46 = t193 * t93 - t95;
t225 = t153 * mrSges(4,1) + t152 * mrSges(4,2);
t62 = t126 * t193 - t123;
t130 = t191 * t135;
t71 = t127 * t193 - t130;
t25 = t168 * qJ(5) - qJD(5) * t206 + t47;
t223 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1;
t97 = mrSges(6,1) * t246 - mrSges(6,3) * t245;
t13 = -t42 * mrSges(7,1) + t41 * mrSges(7,2);
t213 = pkin(4) * t191 - t247;
t212 = t191 * t53 - t193 * t54;
t36 = t130 + (-pkin(8) * t176 - t127) * t193 + t282 * t206;
t48 = t176 * t268 + t57;
t11 = -t195 * t48 + t197 * t36;
t12 = t195 * t36 + t197 * t48;
t60 = -mrSges(7,2) * t160 + mrSges(7,3) * t83;
t61 = mrSges(7,1) * t160 - mrSges(7,3) * t305;
t211 = -t195 * t61 + t197 * t60;
t117 = t175 * t176;
t21 = -pkin(4) * t153 - t32;
t202 = t92 * mrSges(5,2) + t21 * mrSges(6,2) - t32 * mrSges(5,3) - t30 * mrSges(6,3) - t233 * t153 + t303 * t312 + t304 * t308;
t180 = -pkin(4) * t193 - t226;
t169 = t193 * t282 + t226;
t150 = Ifges(7,3) * t153;
t146 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t163;
t118 = t210 * t176;
t91 = (-qJD(3) * t178 - t227) * t196 + t238;
t82 = Ifges(7,4) * t83;
t78 = t176 * t213 - t297;
t64 = -t163 * t213 + t129;
t59 = -t176 * t307 + t297;
t58 = pkin(4) * t206 - t71;
t56 = -t165 * t176 + t167 * t175;
t55 = qJD(6) * t117 + t167 * t210;
t51 = -pkin(4) * t164 - t62;
t50 = t163 * t307 - t129;
t45 = t167 * t213 - t315;
t37 = -t137 * t282 + t204;
t34 = -pkin(4) * t168 - t46;
t31 = -t167 * t307 + t315;
t28 = mrSges(7,2) * t153 + mrSges(7,3) * t42;
t27 = -mrSges(7,1) * t153 - mrSges(7,3) * t41;
t24 = Ifges(7,1) * t305 + Ifges(7,5) * t160 + t82;
t23 = Ifges(7,2) * t83 + Ifges(7,6) * t160 + t269;
t20 = t152 * t231 + t201;
t18 = t167 * t268 + t25;
t17 = t95 + (-pkin(8) * t167 - t93) * t193 - t282 * t168;
t4 = -qJD(6) * t12 + t17 * t197 - t18 * t195;
t3 = qJD(6) * t11 + t17 * t195 + t18 * t197;
t7 = [(-t135 * t153 - t152 * t297) * mrSges(4,3) - t297 * t98 + m(5) * (-t103 * t208 + t32 * t71 + t33 * t72 + t46 * t53 + t47 * t54 - t254) + m(4) * (t102 * t129 - t103 * t128 + t135 * t91 - t254) + (t163 * t223 - t292) * t168 - t250 * t103 + (t300 + (t255 + t260 - t256 + t259) * t163 / 0.2e1 + t291) * t167 + 0.2e1 * t294 * t234 + t230 * t225 + (t92 * mrSges(4,3) + t152 * Ifges(4,1) - t153 * Ifges(4,4) + ((Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t245 + t202) * t193 + (((Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t191 + t311 * t193) * t152 + t293) * t191) * t176 - (t150 / 0.2e1 - t33 * mrSges(5,2) + t16 * mrSges(6,3) - t91 * mrSges(4,3) + t32 * mrSges(5,1) - t21 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + t306) * t153 + (-Ifges(4,4) + t312 * t193 + (-Ifges(5,6) + Ifges(6,6)) * t191) * t152 + t295) * t206 + (t1 * t117 - t118 * t2 - t5 * t55 + t56 * t6) * mrSges(7,3) + m(7) * (t1 * t12 + t11 * t2 - t20 * t59 + t3 * t6 + t31 * t37 + t4 * t5) + m(6) * (t16 * t57 + t21 * t58 + t25 * t44 + t30 * t78 + t34 * t43 + t45 * t52) + t102 * t146 + t58 * t112 + t57 * t113 - t20 * (-mrSges(7,1) * t117 + mrSges(7,2) * t118) + t25 * t104 + t47 * t105 + t46 * t106 + t34 * t107 + t72 * t110 + t71 * t111 + t78 * t97 + t45 * t89 + t3 * t60 + t4 * t61 + t55 * t24 / 0.2e1 + t37 * (-mrSges(7,1) * t56 + mrSges(7,2) * t55) + t56 * t23 / 0.2e1 + t59 * t13 + t31 * t35 + t11 * t27 + t12 * t28 + (Ifges(7,5) * t55 + Ifges(7,6) * t56) * t275 + (Ifges(7,5) * t118 + Ifges(7,6) * t117) * t278 + (Ifges(7,1) * t55 + Ifges(7,4) * t56) * t283 + (Ifges(7,4) * t55 + Ifges(7,2) * t56) * t285 + (Ifges(7,4) * t118 + Ifges(7,2) * t117) * t287 + (Ifges(7,1) * t118 + Ifges(7,4) * t117) * t288 + t118 * t289 + t117 * t290; -t210 * t27 + t175 * t28 + t236 * t61 + t237 * t60 + t239 * t193 + t240 * t191 + (t250 + t263) * t164 + (-t241 * t191 + t242 * t193 + t146) * t163 - m(4) * (-t128 * t164 - t129 * t163) + t225 - t294 * qJD(1) ^ 2 + (t1 * t175 + t164 * t37 - t2 * t210 + t236 * t5 + t237 * t6) * m(7) + (t16 * t191 - t193 * t21 - t164 * t52 - (-t191 * t43 - t193 * t44) * t163) * m(6) + (-t163 * t212 + t164 * t208 + t191 * t33 + t193 * t32) * m(5); (qJ(4) * t240 + qJD(4) * t242 - t293) * t193 + (-pkin(3) * t92 - t212 * qJD(4) + (-t191 * t32 + t193 * t33) * qJ(4) + t208 * t129 - t53 * t62 - t54 * t63) * m(5) + (-t1 * t210 - t175 * t2 + t236 * t6 - t237 * t5) * mrSges(7,3) - t20 * (mrSges(7,1) * t210 + mrSges(7,2) * t175) + (Ifges(7,5) * t175 - Ifges(7,6) * t210) * t278 + (Ifges(7,4) * t175 - Ifges(7,2) * t210) * t287 + (Ifges(7,1) * t175 - Ifges(7,4) * t210) * t288 - t210 * t290 + t292 * t164 + t250 * t129 + (-t165 / 0.2e1 + t109 / 0.2e1) * t24 + (-Ifges(7,5) * t165 - Ifges(7,6) * t166) * t275 + (-Ifges(7,1) * t165 - Ifges(7,4) * t166) * t283 + (-Ifges(7,4) * t165 - Ifges(7,2) * t166) * t285 + (-t166 / 0.2e1 + t108 / 0.2e1) * t23 + (-mrSges(7,1) * t236 + mrSges(7,2) * t237) * t37 + ((t259 / 0.2e1 - t256 / 0.2e1 + t260 / 0.2e1 + t255 / 0.2e1) * t163 + (t314 - t223) * t164 + t291) * t163 + (Ifges(4,5) + (-Ifges(6,3) * t193 + t258) * t272 + (Ifges(5,2) * t193 + t262) * t273 + (t191 * t313 - t257 + t261) * t271) * t152 + t309 * t61 + (t1 * t134 + t132 * t2 - t169 * t20 + t310 * t6 + t309 * t5 + (qJD(5) * t191 - t50) * t37) * m(7) + t310 * t60 + (-Ifges(7,5) * t109 - Ifges(7,6) * t108) * t276 + (-Ifges(7,1) * t109 - Ifges(7,4) * t108) * t284 + (-Ifges(7,4) * t109 - Ifges(7,2) * t108) * t286 + (-qJ(4) * t239 - qJD(4) * t241 + qJD(5) * t263 + t202) * t191 + t180 * t97 + t169 * t13 - Ifges(4,6) * t153 - t128 * t146 + t132 * t27 + t134 * t28 - t49 * t104 - t63 * t105 - t62 * t106 - t51 * t107 - pkin(3) * t98 - t64 * t89 - t91 * mrSges(4,2) - t92 * mrSges(4,1) - t50 * t35 + (-t43 * t51 - t44 * t49 - t52 * t64 + t180 * t30 + (qJ(4) * t16 + qJD(4) * t44) * t193 + (qJ(4) * t21 + qJD(4) * t43 - qJD(5) * t52) * t191) * m(6) + t175 * t289; t241 * t209 + t242 * t137 + t83 * t60 - t305 * t61 - t13 + t97 + t98 + (-t305 * t5 + t6 * t83 + t20) * m(7) + (t137 * t44 - t209 * t43 + t30) * m(6) + (t137 * t54 + t209 * t53 + t92) * m(5); t195 * t28 + t197 * t27 - t263 * t209 + t211 * qJD(6) + (-t104 - t211) * t163 + t112 + (t1 * t195 - t209 * t37 + t197 * t2 + t160 * (-t195 * t5 + t197 * t6)) * m(7) + (-t163 * t44 + t209 * t52 + t21) * m(6); -t150 - t37 * (mrSges(7,1) * t305 + mrSges(7,2) * t83) + (Ifges(7,1) * t83 - t269) * t284 + t23 * t283 + (Ifges(7,5) * t83 - Ifges(7,6) * t305) * t276 - t5 * t60 + t6 * t61 + (t305 * t6 + t5 * t83) * mrSges(7,3) + (-Ifges(7,2) * t305 + t24 + t82) * t286 - t295;];
tauc  = t7(:);
