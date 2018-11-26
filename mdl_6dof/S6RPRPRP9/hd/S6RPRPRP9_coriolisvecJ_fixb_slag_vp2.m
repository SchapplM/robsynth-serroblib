% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2018-11-23 16:02
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:01:52
% EndTime: 2018-11-23 16:02:00
% DurationCPUTime: 9.10s
% Computational Cost: add. (4921->525), mult. (11270->713), div. (0->0), fcn. (7191->6), ass. (0->228)
t314 = Ifges(6,1) + Ifges(7,1);
t322 = Ifges(6,4) - Ifges(7,5);
t312 = Ifges(7,4) + Ifges(6,5);
t183 = sin(pkin(9));
t187 = cos(qJ(5));
t184 = cos(pkin(9));
t185 = sin(qJ(5));
t249 = t184 * t185;
t156 = t183 * t187 + t249;
t140 = t156 * qJD(1);
t186 = sin(qJ(3));
t119 = t186 * t140;
t301 = t156 * qJD(5);
t244 = t119 + t301;
t155 = t183 * t185 - t187 * t184;
t305 = t155 * t186;
t120 = qJD(1) * t305;
t141 = t155 * qJD(5);
t243 = t120 + t141;
t321 = Ifges(6,6) - Ifges(7,6);
t325 = Ifges(6,3) + Ifges(7,2);
t324 = -Ifges(4,6) / 0.2e1;
t233 = t183 * qJD(3);
t188 = cos(qJ(3));
t240 = qJD(1) * t188;
t153 = t184 * t240 + t233;
t323 = t153 / 0.2e1;
t232 = qJD(1) * qJD(3);
t215 = t188 * t232;
t193 = qJD(3) * t305;
t222 = t183 * t240;
t239 = qJD(3) * t184;
t198 = t222 - t239;
t94 = t185 * t153 + t187 * t198;
t55 = qJD(1) * t193 - qJD(5) * t94;
t192 = t187 * t153 - t185 * t198;
t221 = t186 * t233;
t212 = qJD(1) * t221;
t238 = qJD(3) * t186;
t213 = t238 * t249;
t56 = -qJD(1) * t213 + qJD(5) * t192 - t187 * t212;
t320 = t312 * t215 + t314 * t55 - t322 * t56;
t241 = qJD(1) * t186;
t180 = qJD(5) + t241;
t272 = Ifges(7,5) * t94;
t91 = Ifges(6,4) * t94;
t309 = t312 * t180 + t192 * t314 + t272 - t91;
t189 = -pkin(1) - pkin(7);
t174 = qJD(1) * t189 + qJD(2);
t162 = t186 * t174;
t223 = t183 * t241;
t122 = -pkin(4) * t223 + t162;
t319 = pkin(5) * t244 + qJ(6) * t243 - qJD(6) * t156 - t122;
t318 = -Ifges(6,6) / 0.2e1;
t317 = Ifges(7,6) / 0.2e1;
t295 = t55 / 0.2e1;
t293 = t56 / 0.2e1;
t291 = -t94 / 0.2e1;
t290 = t94 / 0.2e1;
t279 = t180 / 0.2e1;
t288 = t192 / 0.2e1;
t316 = qJD(3) * t324;
t315 = Ifges(5,5) * t323;
t16 = t56 * mrSges(7,1) - t55 * mrSges(7,3);
t17 = t56 * mrSges(6,1) + t55 * mrSges(6,2);
t310 = -t16 - t17;
t130 = t155 * t188;
t308 = -qJD(3) * t130 - t186 * t301 - t140;
t128 = t156 * t188;
t234 = qJD(5) * t187;
t235 = qJD(5) * t185;
t303 = -t183 * t235 + t184 * t234;
t307 = -t155 * qJD(1) + qJD(3) * t128 + t186 * t303;
t306 = (qJ(2) * (m(3) + m(4)));
t163 = pkin(3) * t186 - qJ(4) * t188 + qJ(2);
t145 = t163 * qJD(1);
t148 = qJD(3) * qJ(4) + t162;
t83 = t184 * t145 - t148 * t183;
t57 = pkin(4) * t241 - pkin(8) * t153 + t83;
t84 = t183 * t145 + t184 * t148;
t60 = -pkin(8) * t198 + t84;
t19 = -t185 * t60 + t187 * t57;
t304 = qJD(6) - t19;
t302 = t325 * t215 + t312 * t55 - t321 * t56;
t237 = qJD(3) * t188;
t20 = t185 * t57 + t187 * t60;
t248 = t184 * t186;
t231 = pkin(8) * t248;
t194 = (pkin(4) * t188 + t231) * qJD(1);
t206 = pkin(3) * t188 + qJ(4) * t186;
t131 = qJD(3) * t206 - qJD(4) * t188 + qJD(2);
t111 = t131 * qJD(1);
t245 = t188 * t174;
t134 = (qJD(4) + t245) * qJD(3);
t63 = t184 * t111 - t134 * t183;
t45 = qJD(3) * t194 + t63;
t64 = t183 * t111 + t184 * t134;
t58 = pkin(8) * t212 + t64;
t4 = -qJD(5) * t20 - t185 * t58 + t187 * t45;
t247 = t186 * t189;
t113 = t183 * t163 + t184 * t247;
t251 = t183 * t188;
t102 = -pkin(8) * t251 + t113;
t151 = t184 * t163;
t214 = -t183 * t189 + pkin(4);
t92 = -pkin(8) * t184 * t188 + t186 * t214 + t151;
t265 = t187 * t102 + t185 * t92;
t116 = t184 * t131;
t66 = t116 + (t188 * t214 + t231) * qJD(3);
t236 = qJD(3) * t189;
t220 = t188 * t236;
t99 = t183 * t131 + t184 * t220;
t80 = pkin(8) * t221 + t99;
t9 = -qJD(5) * t265 - t185 * t80 + t187 * t66;
t3 = t185 * t45 + t187 * t58 + t57 * t234 - t235 * t60;
t1 = qJ(6) * t215 + qJD(6) * t180 + t3;
t2 = -pkin(5) * t215 - t4;
t300 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t14 = -pkin(5) * t180 + t304;
t15 = qJ(6) * t180 + t20;
t226 = Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1;
t227 = t317 + t318;
t228 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t299 = t226 * t180 + t228 * t192 + t227 * t94 + t15 * mrSges(7,3) + t19 * mrSges(6,1) + t83 * mrSges(5,1) + t316 - (Ifges(4,4) * t188 - t186 * Ifges(4,2)) * qJD(1) / 0.2e1 + Ifges(6,6) * t291 + Ifges(7,6) * t290 - Ifges(5,6) * t198 / 0.2e1 + Ifges(5,3) * t241 / 0.2e1 + t315 - t14 * mrSges(7,1) - t20 * mrSges(6,2) - t84 * mrSges(5,2) + t312 * t288 + t325 * t279;
t298 = -Ifges(5,1) / 0.2e1;
t297 = Ifges(7,5) * t295 + Ifges(7,3) * t293 + t215 * t317;
t296 = -t55 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t293 + t215 * t318;
t294 = -t56 / 0.2e1;
t289 = -t192 / 0.2e1;
t280 = -t180 / 0.2e1;
t278 = t183 / 0.2e1;
t277 = t184 / 0.2e1;
t275 = mrSges(6,3) * t94;
t274 = mrSges(6,3) * t192;
t273 = Ifges(6,4) * t192;
t270 = pkin(8) + qJ(4);
t39 = -mrSges(7,2) * t56 + mrSges(7,3) * t215;
t42 = -mrSges(6,2) * t215 - mrSges(6,3) * t56;
t269 = t39 + t42;
t40 = mrSges(6,1) * t215 - mrSges(6,3) * t55;
t41 = -mrSges(7,1) * t215 + t55 * mrSges(7,2);
t268 = t41 - t40;
t69 = -mrSges(6,2) * t180 - t275;
t70 = -mrSges(7,2) * t94 + mrSges(7,3) * t180;
t267 = t69 + t70;
t71 = mrSges(6,1) * t180 - t274;
t72 = -mrSges(7,1) * t180 + mrSges(7,2) * t192;
t266 = t72 - t71;
t158 = t206 * qJD(1);
t103 = t184 * t158 - t183 * t245;
t73 = t103 + t194;
t104 = t183 * t158 + t184 * t245;
t82 = pkin(8) * t223 + t104;
t25 = t185 * t73 + t187 * t82;
t264 = Ifges(5,1) * t153;
t263 = Ifges(4,4) * t186;
t262 = Ifges(5,4) * t184;
t260 = Ifges(5,5) * t184;
t259 = Ifges(5,2) * t183;
t258 = Ifges(5,6) * t183;
t257 = qJ(2) * mrSges(4,1);
t256 = qJ(2) * mrSges(4,2);
t255 = Ifges(4,5) * qJD(3);
t253 = qJD(3) * mrSges(4,2);
t117 = -mrSges(5,2) * t241 - mrSges(5,3) * t198;
t250 = t184 * t117;
t242 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t198 - t153 * mrSges(5,2) - mrSges(4,3) * t240;
t225 = Ifges(5,5) * t241;
t224 = Ifges(5,6) * t241;
t182 = -pkin(4) * t184 - pkin(3);
t157 = t174 * t238;
t152 = pkin(4) * t251 - t188 * t189;
t211 = -t153 / 0.2e1 + t233 / 0.2e1;
t210 = mrSges(4,1) * t186 + mrSges(4,2) * t188;
t209 = -mrSges(5,1) * t183 - mrSges(5,2) * t184;
t208 = -Ifges(5,1) * t184 + Ifges(5,4) * t183;
t207 = t259 - t262;
t205 = -t183 * t63 + t184 * t64;
t204 = t183 * t84 + t184 * t83;
t203 = -t183 * t83 + t184 * t84;
t24 = -t185 * t82 + t187 * t73;
t34 = -t102 * t185 + t187 * t92;
t168 = t270 * t183;
t169 = t270 * t184;
t199 = -t187 * t168 - t169 * t185;
t110 = -t168 * t185 + t169 * t187;
t138 = -pkin(4) * t221 + t186 * t236;
t137 = -qJD(3) * pkin(3) + qJD(4) - t245;
t8 = -t102 * t235 + t185 * t66 + t187 * t80 + t92 * t234;
t197 = t209 * qJD(1);
t114 = -pkin(4) * t212 + t157;
t100 = pkin(4) * t198 + t137;
t171 = -mrSges(4,3) * t241 - t253;
t159 = t210 * qJD(1);
t147 = t255 + (Ifges(4,1) * t188 - t263) * qJD(1);
t133 = (mrSges(5,1) * t188 + mrSges(5,3) * t248) * t232;
t132 = (mrSges(5,3) * t183 * t186 - mrSges(5,2) * t188) * t232;
t127 = t156 * t186;
t121 = t197 * t238;
t118 = mrSges(5,1) * t241 - mrSges(5,3) * t153;
t112 = -t183 * t247 + t151;
t108 = (Ifges(5,5) * t188 + t186 * t208) * t232;
t107 = (Ifges(5,6) * t188 + t186 * t207) * t232;
t98 = -t183 * t220 + t116;
t90 = Ifges(7,5) * t192;
t89 = pkin(5) * t155 - qJ(6) * t156 + t182;
t88 = -Ifges(5,4) * t198 + t225 + t264;
t87 = Ifges(5,4) * t153 - Ifges(5,2) * t198 + t224;
t79 = -t187 * t221 + t188 * t303 - t213;
t77 = -t188 * t301 + t193;
t68 = qJD(4) * t156 + qJD(5) * t110;
t67 = -qJD(4) * t155 + qJD(5) * t199;
t59 = pkin(5) * t128 + qJ(6) * t130 + t152;
t38 = mrSges(6,1) * t94 + mrSges(6,2) * t192;
t37 = mrSges(7,1) * t94 - mrSges(7,3) * t192;
t36 = pkin(5) * t192 + qJ(6) * t94;
t33 = -pkin(5) * t186 - t34;
t32 = qJ(6) * t186 + t265;
t29 = -t94 * Ifges(6,2) + t180 * Ifges(6,6) + t273;
t26 = t180 * Ifges(7,6) + t94 * Ifges(7,3) + t90;
t23 = -pkin(5) * t240 - t24;
t22 = qJ(6) * t240 + t25;
t21 = t94 * pkin(5) - qJ(6) * t192 + t100;
t18 = pkin(5) * t79 - qJ(6) * t77 + qJD(6) * t130 + t138;
t7 = -pkin(5) * t237 - t9;
t6 = qJ(6) * t237 + qJD(6) * t186 + t8;
t5 = pkin(5) * t56 - qJ(6) * t55 - qJD(6) * t192 + t114;
t10 = [(Ifges(7,3) * t290 - Ifges(6,2) * t291 - t20 * mrSges(6,3) - t15 * mrSges(7,2) + t26 / 0.2e1 + t21 * mrSges(7,1) - t29 / 0.2e1 + t100 * mrSges(6,1) - t321 * t279 - t322 * t288) * t79 + (-t128 * t322 - t130 * t314) * t295 + (mrSges(6,2) * t100 + mrSges(7,2) * t14 - mrSges(6,3) * t19 - mrSges(7,3) * t21 + Ifges(6,4) * t291 + Ifges(7,5) * t290 + t309 / 0.2e1 + t312 * t279 + t314 * t288) * t77 + ((-m(5) * t189 - t209) * t162 + (Ifges(5,6) * t277 + t324) * qJD(3) + t189 * t171 + t315 + t299) * t237 + t128 * t296 + t128 * t297 + (((2 * mrSges(3,3)) + t210 + (2 * t306)) * qJD(2) + (0.2e1 * t257 + (t260 / 0.2e1 - t258 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t188 - t228 * t130 + t227 * t128) * t237) * qJD(1) + t302 * t186 / 0.2e1 - t320 * t130 / 0.2e1 + (t63 * mrSges(5,1) - t64 * mrSges(5,2) + Ifges(6,6) * t294 + Ifges(7,6) * t293 + t312 * t295 + ((-Ifges(4,5) / 0.2e1 + t207 * t277) * qJD(3) - t147 / 0.2e1 + t137 * t209 + t208 * t323 - t184 * t88 / 0.2e1 + t87 * t278 + t204 * mrSges(5,3) + (m(5) * t137 - t242) * t189 + (-0.2e1 * t256 + (-0.3e1 / 0.2e1 * t260 + 0.3e1 / 0.2e1 * t258 + 0.3e1 / 0.2e1 * Ifges(4,4)) * t186 + (t184 ^ 2 * t298 + 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(4,2) + (0.3e1 / 0.2e1 * t262 - t259) * t183 + t226) * t188) * qJD(1)) * qJD(3) + t300) * t186 + m(5) * (t112 * t63 + t113 * t64 + t83 * t98 + t84 * t99) + m(7) * (t1 * t32 + t14 * t7 + t15 * t6 + t18 * t21 + t2 * t33 + t5 * t59) + t18 * t37 + t32 * t39 + t34 * t40 + t33 * t41 + t59 * t16 + t8 * t69 + t6 * t70 + t9 * t71 + t7 * t72 + m(6) * (t100 * t138 + t114 * t152 + t19 * t9 + t20 * t8 + t265 * t3 + t34 * t4) + t265 * t42 + (-Ifges(6,4) * t130 - Ifges(6,2) * t128) * t294 + (-t1 * t128 - t130 * t2) * mrSges(7,2) + (-Ifges(7,5) * t130 + Ifges(7,3) * t128) * t293 + (-t128 * t3 + t130 * t4) * mrSges(6,3) + t5 * (mrSges(7,1) * t128 + mrSges(7,3) * t130) + t114 * (mrSges(6,1) * t128 - mrSges(6,2) * t130) + t99 * t117 + t98 * t118 + t113 * t132 + t112 * t133 + t138 * t38 + t152 * t17 + qJD(2) * t159 + (-t183 * t107 / 0.2e1 + t108 * t277 - t189 * t121 + (-t183 * t64 - t184 * t63) * mrSges(5,3)) * t188; (t132 * t184 - t133 * t183) * t186 - t269 * t305 + t268 * t127 + (-t121 + t310) * t188 + ((-t118 * t183 + t171 + t250) * t188 + (t37 + t38 - t242) * t186) * qJD(3) + m(5) * (t205 * t186 + (t137 * t186 + (t203 - t162) * t188) * qJD(3)) + t308 * t267 + t307 * t266 + (-t1 * t305 + t127 * t2 + t14 * t307 + t15 * t308 - t188 * t5 + t21 * t238) * m(7) + (t100 * t238 - t114 * t188 - t127 * t4 - t19 * t307 + t20 * t308 - t3 * t305) * m(6) + (-m(5) * t204 - t117 * t183 - t118 * t184 - t159 + (-mrSges(3,3) - t306) * qJD(1)) * qJD(1); ((-t255 / 0.2e1 + t147 / 0.2e1 + (t256 - t263 / 0.2e1) * qJD(1) + (t137 * mrSges(5,2) - t83 * mrSges(5,3) + t264 / 0.2e1 + t88 / 0.2e1 + t225 / 0.2e1) * t184 + (t137 * mrSges(5,1) - t84 * mrSges(5,3) - t87 / 0.2e1 + t239 * t298 - t224 / 0.2e1 + t211 * Ifges(5,4)) * t183) * t186 + ((-t257 + (Ifges(4,4) / 0.2e1 + t258 / 0.2e1) * t188) * qJD(1) + t211 * Ifges(5,5) + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1 - Ifges(5,3) / 0.2e1 + t207 * t278) * t241 + t316 + (-t155 * t321 + t156 * t312) * qJD(3) / 0.2e1 - t299) * t188) * qJD(1) + (-qJD(4) * t118 - t63 * mrSges(5,3) - qJ(4) * t133 + t108 / 0.2e1) * t183 + (Ifges(7,5) * t156 + Ifges(7,3) * t155) * t293 + (Ifges(6,4) * t156 - Ifges(6,2) * t155) * t294 + t155 * t296 + t155 * t297 + (-t100 * t122 + t199 * t4 + t110 * t3 + t114 * t182 + (-t25 + t67) * t20 + (-t24 - t68) * t19) * m(6) - t268 * t199 + (t26 - t29) * (t119 / 0.2e1 + t301 / 0.2e1) + (Ifges(6,4) * t120 - Ifges(7,5) * t141 + Ifges(6,2) * t119 + Ifges(7,3) * t301) * t290 + (-Ifges(6,4) * t141 + Ifges(7,5) * t120 - Ifges(6,2) * t301 - Ifges(7,3) * t119) * t291 + (-pkin(3) * t157 + qJ(4) * t205 + qJD(4) * t203 - t103 * t83 - t104 * t84 - t137 * t162) * m(5) + (t1 * t110 - t199 * t2 + t5 * t89 + t319 * t21 + (-t22 + t67) * t15 + (-t23 + t68) * t14) * m(7) + t319 * t37 + t309 * (-t120 / 0.2e1 - t141 / 0.2e1) + t320 * t156 / 0.2e1 + (t119 * t321 + t120 * t312) * t280 + (-t141 * t312 - t301 * t321) * t279 + (t119 * t322 + t120 * t314) * t289 + (-t155 * t322 + t156 * t314) * t295 + (-t141 * t314 - t301 * t322) * t288 + (qJD(4) * t117 + t64 * mrSges(5,3) + qJ(4) * t132 + t107 / 0.2e1) * t184 - t25 * t69 - t22 * t70 - t24 * t71 - t23 * t72 + t89 * t16 - t104 * t117 - t103 * t118 - pkin(3) * t121 - t122 * t38 + t5 * (mrSges(7,1) * t155 - mrSges(7,3) * t156) + t114 * (mrSges(6,1) * t155 + mrSges(6,2) * t156) + (mrSges(7,1) * t244 + mrSges(7,3) * t243) * t21 + (mrSges(6,1) * t244 - mrSges(6,2) * t243) * t100 + (-t155 * t3 - t156 * t4 + t19 * t243 - t20 * t244) * mrSges(6,3) + (-t1 * t155 - t14 * t243 - t15 * t244 + t156 * t2) * mrSges(7,2) + t182 * t17 + ((-t171 - t253) * t188 + ((-mrSges(5,1) * t184 + mrSges(5,2) * t183 - mrSges(4,1)) * qJD(3) + t242) * t186) * t174 + t266 * t68 + t267 * t67 + t269 * t110; t117 * t222 + t153 * t118 + t267 * t94 - t266 * t192 + (-t250 + (m(5) * t174 + t197) * t186) * qJD(3) - m(5) * (-t83 * t153 - t198 * t84) + (-t14 * t192 + t15 * t94 + t5) * m(7) + (t19 * t192 + t20 * t94 + t114) * m(6) - t310; t29 * t288 + (-t314 * t94 + t26 - t273 + t90) * t289 + (-pkin(5) * t2 + qJ(6) * t1 - t14 * t20 + t15 * t304 - t21 * t36) * m(7) + t300 - t36 * t37 + qJ(6) * t39 - pkin(5) * t41 + qJD(6) * t70 + (Ifges(7,3) * t192 - t272) * t291 + (t14 * t94 + t15 * t192) * mrSges(7,2) - t21 * (mrSges(7,1) * t192 + mrSges(7,3) * t94) - t100 * (mrSges(6,1) * t192 - mrSges(6,2) * t94) + (-t192 * t321 - t312 * t94) * t280 + (-Ifges(6,2) * t192 + t309 - t91) * t290 + (-t266 + t274) * t20 + (-t267 - t275) * t19 + t302; -t180 * t70 + t192 * t37 + 0.2e1 * (t2 / 0.2e1 + t15 * t280 + t21 * t288) * m(7) + t41;];
tauc  = t10(:);
