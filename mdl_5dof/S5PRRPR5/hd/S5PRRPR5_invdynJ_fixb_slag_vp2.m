% Calculate vector of inverse dynamics joint torques for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:20
% EndTime: 2019-12-05 16:25:55
% DurationCPUTime: 12.06s
% Computational Cost: add. (3733->516), mult. (8848->740), div. (0->0), fcn. (6621->14), ass. (0->245)
t167 = sin(qJ(5));
t170 = cos(qJ(5));
t195 = -t170 * mrSges(6,1) + t167 * mrSges(6,2);
t302 = -m(6) * pkin(4) - mrSges(5,1) + t195;
t169 = sin(qJ(2));
t163 = sin(pkin(5));
t242 = qJD(1) * t163;
t216 = t169 * t242;
t168 = sin(qJ(3));
t238 = qJD(3) * t168;
t226 = pkin(3) * t238;
t320 = t226 - t216;
t171 = cos(qJ(3));
t166 = -qJ(4) - pkin(7);
t205 = qJD(3) * t166;
t114 = qJD(4) * t171 + t168 * t205;
t161 = sin(pkin(10));
t164 = cos(pkin(10));
t127 = t161 * t168 - t164 * t171;
t176 = -qJD(4) * t168 + t171 * t205;
t172 = cos(qJ(2));
t215 = t172 * t242;
t303 = t164 * t114 + t127 * t215 + t161 * t176;
t128 = t161 * t171 + t164 * t168;
t121 = t128 * qJD(3);
t122 = t127 * qJD(3);
t319 = pkin(4) * t121 + pkin(8) * t122 + t320;
t233 = t171 * qJD(2);
t234 = t168 * qJD(2);
t120 = -t161 * t233 - t164 * t234;
t269 = mrSges(5,3) * t120;
t90 = qJD(3) * t170 + t120 * t167;
t91 = qJD(3) * t167 - t120 * t170;
t271 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t90 + mrSges(6,2) * t91 - t269;
t135 = qJD(2) * pkin(7) + t216;
t200 = qJ(4) * qJD(2) + t135;
t165 = cos(pkin(5));
t241 = qJD(1) * t165;
t214 = t168 * t241;
t83 = t171 * t200 + t214;
t263 = t161 * t83;
t151 = t171 * t241;
t82 = -t168 * t200 + t151;
t77 = qJD(3) * pkin(3) + t82;
t32 = t164 * t77 - t263;
t318 = m(5) * t32 - t271;
t240 = qJD(2) * t163;
t209 = qJD(1) * t240;
t144 = t172 * t209;
t230 = qJDD(1) * t163;
t110 = t169 * t230 + t144;
t100 = qJDD(2) * pkin(7) + t110;
t316 = qJD(3) * t241 + t100;
t160 = qJ(3) + pkin(10);
t158 = sin(t160);
t159 = cos(t160);
t315 = t302 * t159 + (-m(6) * pkin(8) + mrSges(5,2) - mrSges(6,3)) * t158;
t246 = t163 * t169;
t123 = t165 * t171 - t168 * t246;
t257 = cos(pkin(9));
t202 = t257 * t172;
t162 = sin(pkin(9));
t248 = t162 * t169;
t118 = -t165 * t248 + t202;
t245 = t163 * t171;
t314 = -t118 * t168 + t162 * t245;
t119 = t127 * qJD(2);
t111 = qJD(5) + t119;
t232 = qJD(2) * qJD(3);
t132 = qJDD(2) * t171 - t168 * t232;
t133 = qJDD(2) * t168 + t171 * t232;
t79 = t132 * t161 + t133 * t164;
t36 = qJD(5) * t90 + qJDD(3) * t167 + t170 * t79;
t289 = t36 / 0.2e1;
t37 = -qJD(5) * t91 + qJDD(3) * t170 - t167 * t79;
t288 = t37 / 0.2e1;
t78 = t132 * t164 - t133 * t161;
t73 = qJDD(5) - t78;
t287 = t73 / 0.2e1;
t313 = m(5) + m(6);
t312 = -t121 / 0.2e1;
t311 = -t122 / 0.2e1;
t310 = t132 / 0.2e1;
t156 = pkin(3) * t171 + pkin(2);
t66 = pkin(4) * t127 - pkin(8) * t128 - t156;
t140 = t166 * t171;
t210 = t166 * t168;
t88 = -t164 * t140 + t161 * t210;
t30 = -t167 * t88 + t170 * t66;
t309 = qJD(5) * t30 + t319 * t167 + t170 * t303;
t75 = t164 * t83;
t33 = t161 * t77 + t75;
t29 = qJD(3) * pkin(8) + t33;
t107 = -qJD(2) * t156 + qJD(4) - t215;
t46 = pkin(4) * t119 + pkin(8) * t120 + t107;
t14 = -t167 * t29 + t170 * t46;
t308 = t14 * mrSges(6,1);
t15 = t167 * t46 + t170 * t29;
t307 = t15 * mrSges(6,2);
t11 = -mrSges(6,1) * t37 + mrSges(6,2) * t36;
t62 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t79;
t306 = t11 - t62;
t31 = t167 * t66 + t170 * t88;
t305 = -qJD(5) * t31 - t167 * t303 + t319 * t170;
t260 = qJDD(3) / 0.2e1;
t194 = mrSges(6,1) * t167 + mrSges(6,2) * t170;
t301 = -mrSges(5,3) - t194;
t197 = -mrSges(4,1) * t171 + mrSges(4,2) * t168;
t63 = mrSges(5,1) * t119 - mrSges(5,2) * t120;
t300 = t197 * qJD(2) + t63;
t203 = t257 * t169;
t247 = t162 * t172;
t116 = t165 * t203 + t247;
t204 = t163 * t257;
t299 = -t116 * t168 - t171 * t204;
t235 = qJD(5) * t170;
t183 = -t122 * t167 + t128 * t235;
t229 = qJDD(1) * t165;
t43 = -t135 * t238 + t168 * t229 + t171 * t316;
t148 = t171 * t229;
t94 = t135 * t171 + t214;
t44 = -t94 * qJD(3) - t100 * t168 + t148;
t298 = -t168 * t44 + t171 * t43;
t143 = t169 * t209;
t109 = t172 * t230 - t143;
t99 = -qJDD(2) * pkin(2) - t109;
t64 = -pkin(3) * t132 + qJDD(4) + t99;
t20 = -pkin(4) * t78 - pkin(8) * t79 + t64;
t231 = qJD(2) * qJD(4);
t237 = qJD(3) * t171;
t23 = -t135 * t237 + qJDD(3) * pkin(3) - qJ(4) * t133 + t148 + (-t231 - t316) * t168;
t24 = qJ(4) * t132 + t171 * t231 + t43;
t8 = t161 * t23 + t164 * t24;
t6 = qJDD(3) * pkin(8) + t8;
t1 = qJD(5) * t14 + t167 * t20 + t170 * t6;
t2 = -qJD(5) * t15 - t167 * t6 + t170 * t20;
t297 = t1 * t170 - t167 * t2;
t296 = 0.2e1 * t260;
t225 = m(4) * pkin(7) + mrSges(4,3);
t295 = mrSges(3,2) - t225 + t301;
t179 = m(4) * pkin(2) - t197;
t294 = mrSges(3,1) + t179 - t315;
t293 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t173 = qJD(2) ^ 2;
t292 = Ifges(6,1) * t289 + Ifges(6,4) * t288 + Ifges(6,5) * t287;
t291 = m(5) * pkin(3);
t286 = -t90 / 0.2e1;
t285 = -t91 / 0.2e1;
t284 = t91 / 0.2e1;
t283 = -t111 / 0.2e1;
t281 = -t120 / 0.2e1;
t279 = t170 / 0.2e1;
t278 = pkin(3) * t161;
t277 = pkin(3) * t164;
t276 = g(3) * t163;
t273 = t91 * Ifges(6,4);
t270 = mrSges(5,3) * t119;
t268 = Ifges(4,4) * t168;
t267 = Ifges(4,4) * t171;
t266 = Ifges(5,4) * t120;
t265 = Ifges(6,4) * t167;
t264 = Ifges(6,4) * t170;
t254 = t119 * t167;
t253 = t119 * t170;
t251 = t128 * t167;
t250 = t128 * t170;
t249 = t162 * t163;
t244 = t163 * t172;
t239 = qJD(2) * t169;
t236 = qJD(5) * t167;
t228 = Ifges(6,5) * t36 + Ifges(6,6) * t37 + Ifges(6,3) * t73;
t227 = pkin(3) * t234;
t223 = mrSges(4,3) * t234;
t222 = mrSges(4,3) * t233;
t220 = t167 * t244;
t218 = t170 * t244;
t89 = Ifges(6,4) * t90;
t27 = Ifges(6,1) * t91 + Ifges(6,5) * t111 + t89;
t217 = t27 * t279;
t213 = t163 * t239;
t212 = t172 * t240;
t38 = -t78 * mrSges(5,1) + t79 * mrSges(5,2);
t207 = -t236 / 0.2e1;
t199 = t314 * pkin(3);
t193 = Ifges(6,1) * t170 - t265;
t192 = t171 * Ifges(4,2) + t268;
t191 = -Ifges(6,2) * t167 + t264;
t190 = Ifges(4,5) * t171 - Ifges(4,6) * t168;
t189 = Ifges(6,5) * t170 - Ifges(6,6) * t167;
t7 = -t161 * t24 + t164 * t23;
t49 = -mrSges(6,2) * t111 + mrSges(6,3) * t90;
t50 = mrSges(6,1) * t111 - mrSges(6,3) * t91;
t187 = -t167 * t50 + t170 * t49;
t186 = t123 * pkin(3);
t124 = t165 * t168 + t169 * t245;
t58 = t123 * t161 + t124 * t164;
t47 = -t167 * t58 - t218;
t185 = -t170 * t58 + t220;
t28 = -qJD(3) * pkin(4) - t32;
t184 = t28 * t194;
t182 = t122 * t170 + t128 * t236;
t136 = -qJD(2) * pkin(2) - t215;
t181 = t136 * (mrSges(4,1) * t168 + mrSges(4,2) * t171);
t180 = t168 * (Ifges(4,1) * t171 - t268);
t115 = -t165 * t202 + t248;
t117 = t165 * t247 + t203;
t175 = -g(1) * t117 - g(2) * t115 + g(3) * t244;
t157 = Ifges(4,4) * t233;
t155 = -pkin(4) - t277;
t139 = -qJD(3) * mrSges(4,2) + t222;
t138 = qJD(3) * mrSges(4,1) - t223;
t126 = Ifges(4,1) * t234 + Ifges(4,5) * qJD(3) + t157;
t125 = Ifges(4,6) * qJD(3) + qJD(2) * t192;
t113 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t133;
t112 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t132;
t108 = Ifges(5,4) * t119;
t103 = t158 * t165 + t159 * t246;
t97 = -qJD(3) * mrSges(5,2) - t270;
t93 = -t135 * t168 + t151;
t87 = -t140 * t161 - t164 * t210;
t84 = -mrSges(4,1) * t132 + mrSges(4,2) * t133;
t81 = qJD(3) * t123 + t171 * t212;
t80 = -qJD(3) * t124 - t168 * t212;
t71 = t118 * t159 + t158 * t249;
t69 = t116 * t159 - t158 * t204;
t61 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t78;
t60 = -t120 * Ifges(5,1) + Ifges(5,5) * qJD(3) - t108;
t59 = -t119 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t266;
t57 = -t164 * t123 + t124 * t161;
t53 = -pkin(4) * t120 + pkin(8) * t119 + t227;
t42 = t164 * t82 - t263;
t41 = t161 * t80 + t164 * t81;
t40 = t161 * t82 + t75;
t39 = t161 * t81 - t164 * t80;
t26 = t90 * Ifges(6,2) + t111 * Ifges(6,6) + t273;
t25 = t91 * Ifges(6,5) + t90 * Ifges(6,6) + t111 * Ifges(6,3);
t19 = -mrSges(6,2) * t73 + mrSges(6,3) * t37;
t18 = mrSges(6,1) * t73 - mrSges(6,3) * t36;
t17 = t167 * t53 + t170 * t42;
t16 = -t167 * t42 + t170 * t53;
t13 = qJD(5) * t185 - t167 * t41 + t170 * t213;
t12 = qJD(5) * t47 + t167 * t213 + t170 * t41;
t5 = -qJDD(3) * pkin(4) - t7;
t3 = t36 * Ifges(6,4) + t37 * Ifges(6,2) + t73 * Ifges(6,6);
t4 = [m(2) * qJDD(1) + t124 * t112 + t123 * t113 + t12 * t49 + t13 * t50 + t80 * t138 + t81 * t139 + t47 * t18 - t185 * t19 + t41 * t97 + t58 * t61 + t306 * t57 + t271 * t39 + (-m(2) - m(3) - m(4) - t313) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t173 - t38 - t84) * t172 + (-mrSges(3,1) * t173 - mrSges(3,2) * qJDD(2) + qJD(2) * t300) * t169) * t163 + m(3) * (qJDD(1) * t165 ^ 2 + (t109 * t172 + t110 * t169) * t163) + m(4) * (t123 * t44 + t124 * t43 + t80 * t93 + t81 * t94 + (t136 * t239 - t172 * t99) * t163) + m(5) * (-t32 * t39 + t33 * t41 - t57 * t7 + t58 * t8 + (t107 * t239 - t172 * t64) * t163) + m(6) * (-t1 * t185 + t12 * t15 + t13 * t14 + t2 * t47 + t28 * t39 + t5 * t57); t133 * t267 / 0.2e1 + (-t1 * t251 + t14 * t182 - t15 * t183 - t2 * t250) * mrSges(6,3) + (-(t136 * t169 + (-t168 * t93 + t171 * t94) * t172) * t242 - pkin(2) * t99) * m(4) + t99 * t197 + t28 * (mrSges(6,1) * t183 - mrSges(6,2) * t182) + t111 * (-Ifges(6,5) * t182 - Ifges(6,6) * t183 + Ifges(6,3) * t121) / 0.2e1 + t90 * (-Ifges(6,4) * t182 - Ifges(6,2) * t183 + Ifges(6,6) * t121) / 0.2e1 + ((t315 * t172 + (t166 * t313 + t301) * t169) * t163 - t313 * t156 * t244) * g(3) + (-m(6) * t28 + t318) * (-t114 * t161 + t128 * t215 + t164 * t176) + (t1 * t31 + t14 * t305 + t15 * t309 + t2 * t30 + t5 * t87) * m(6) + (t320 * t107 - t156 * t64 + t303 * t33 - t7 * t87 + t8 * t88) * m(5) - t300 * t216 - t183 * t26 / 0.2e1 + (t171 * t112 - t168 * t113 + m(4) * ((-t94 * t168 - t93 * t171) * qJD(3) + t298) - t138 * t237 - t139 * t238) * pkin(7) + (-t237 * t93 - t238 * t94 + t298) * mrSges(4,3) + (t64 * mrSges(5,2) - t7 * mrSges(5,3) + Ifges(5,1) * t79 + Ifges(5,4) * t78 + Ifges(5,5) * t296 + t189 * t287 + t191 * t288 + t193 * t289 + t5 * t194 + t27 * t207) * t128 + (-Ifges(5,2) * t78 + t64 * mrSges(5,1) - Ifges(5,4) * t79 - t8 * mrSges(5,3) + Ifges(6,3) * t287 + Ifges(6,6) * t288 + Ifges(6,5) * t289 + t228 / 0.2e1 - t296 * Ifges(5,6) + t293) * t127 + (t180 + t171 * (-Ifges(4,2) * t168 + t267)) * t232 / 0.2e1 + (-t121 * t33 + t122 * t32) * mrSges(5,3) - t119 * (-Ifges(5,4) * t122 - Ifges(5,2) * t121) / 0.2e1 + t107 * (mrSges(5,1) * t121 - mrSges(5,2) * t122) + (-Ifges(5,1) * t122 - Ifges(5,4) * t121) * t281 - t156 * t38 + t121 * t25 / 0.2e1 - pkin(2) * t84 + t88 * t61 + t250 * t292 + (-Ifges(6,1) * t182 - Ifges(6,4) * t183 + Ifges(6,5) * t121) * t284 + (t169 * t276 - t110 + t144) * mrSges(3,2) + t30 * t18 + t31 * t19 + t63 * t226 - t122 * t217 + (-t172 * t276 + t109 + t143) * mrSges(3,1) + Ifges(3,3) * qJDD(2) - (t169 * t225 + t172 * t179) * t276 - t171 * t139 * t215 + Ifges(4,6) * t171 * t260 + t303 * t97 + t305 * t50 + t306 * t87 - t121 * t307 + t309 * t49 + (Ifges(4,1) * t133 + Ifges(4,4) * t310 + t296 * Ifges(4,5) + t138 * t215) * t168 + (Ifges(5,5) * t311 + Ifges(5,6) * t312 + t181 + t190 * qJD(3) / 0.2e1) * qJD(3) + (-t313 * (-t117 * t156 - t118 * t166) + t295 * t118 + t294 * t117) * g(1) + (-t313 * (-t115 * t156 - t116 * t166) + t295 * t116 + t294 * t115) * g(2) + t171 * (Ifges(4,4) * t133 + Ifges(4,2) * t132 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t121 * t308 + t192 * t310 + t60 * t311 + t59 * t312 + t126 * t237 / 0.2e1 - t125 * t238 / 0.2e1 - t3 * t251 / 0.2e1; (-t14 * t16 - t15 * t17 + t155 * t5 - t28 * t40) * m(6) + t5 * t195 + t119 * t184 - qJD(2) * t181 + (-m(6) * (pkin(8) * t71 + t199) - m(5) * t199 + t71 * mrSges(5,2) - t314 * mrSges(4,1) - (-t118 * t171 - t168 * t249) * mrSges(4,2) + t302 * (-t118 * t158 + t159 * t249)) * g(1) + (t217 + t184) * qJD(5) + t318 * t40 - m(5) * (t107 * t227 + t33 * t42) + (Ifges(5,2) * t120 - t108 + t60) * t119 / 0.2e1 + (-m(6) * (pkin(3) * t299 + t69 * pkin(8)) + t69 * mrSges(5,2) - (-t116 * t171 + t168 * t204) * mrSges(4,2) + t302 * (-t116 * t158 - t159 * t204) + (-t291 - mrSges(4,1)) * t299) * g(2) + (-m(5) * t186 + t103 * mrSges(5,2) - t123 * mrSges(4,1) + t124 * mrSges(4,2) - m(6) * (pkin(8) * t103 + t186) + t302 * (-t158 * t246 + t159 * t165)) * g(3) + (-g(1) * t71 - g(2) * t69 - g(3) * t103 + (-t236 - t254) * t15 + (-t235 - t253) * t14 + t297) * mrSges(6,3) + (t170 * t19 - t167 * t18 + m(6) * ((-t14 * t170 - t15 * t167) * qJD(5) + t297) - t50 * t235 - t49 * t236) * (pkin(8) + t278) - (-Ifges(4,2) * t234 + t126 + t157) * t233 / 0.2e1 + (t111 * t189 + t191 * t90 + t193 * t91) * qJD(5) / 0.2e1 + (-Ifges(5,1) * t119 + t25 + t266) * t120 / 0.2e1 - t107 * (-mrSges(5,1) * t120 - mrSges(5,2) * t119) - qJD(3) * (-Ifges(5,5) * t119 + Ifges(5,6) * t120) / 0.2e1 + (-Ifges(6,5) * t120 - t119 * t193) * t285 + (-Ifges(6,6) * t120 - t119 * t191) * t286 + (-Ifges(6,3) * t120 - t119 * t189) * t283 + (t207 - t254 / 0.2e1) * t26 + (t138 + t223) * t94 + (-t139 + t222) * t93 + t155 * t11 + Ifges(4,6) * t132 + Ifges(4,5) * t133 - t42 * t97 + Ifges(5,6) * t78 + Ifges(5,5) * t79 - t17 * t49 - t16 * t50 - t43 * mrSges(4,2) + t44 * mrSges(4,1) + (Ifges(6,5) * t167 + Ifges(6,6) * t170) * t287 + (Ifges(6,2) * t170 + t265) * t288 + (Ifges(6,1) * t167 + t264) * t289 + (t161 * t8 + t164 * t7) * t291 + t167 * t292 + t59 * t281 - t32 * t270 + t62 * t277 + t61 * t278 + t3 * t279 - t8 * mrSges(5,2) + t7 * mrSges(5,1) + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) - t173 * t180 / 0.2e1 - t33 * t269 - t120 * t307 + t120 * t308 - t63 * t227 - t190 * t232 / 0.2e1 + t125 * t234 / 0.2e1 + t27 * t253 / 0.2e1; t167 * t19 + t170 * t18 + t271 * t120 + t187 * qJD(5) - (-t187 - t97) * t119 + t38 + (t1 * t167 + t120 * t28 + t2 * t170 + t175 + t111 * (-t14 * t167 + t15 * t170)) * m(6) + (t119 * t33 - t120 * t32 + t175 + t64) * m(5); -t28 * (mrSges(6,1) * t91 + mrSges(6,2) * t90) + (Ifges(6,1) * t90 - t273) * t285 + t26 * t284 + (Ifges(6,5) * t90 - Ifges(6,6) * t91) * t283 - t14 * t49 + t15 * t50 - g(1) * ((t117 * t170 - t167 * t71) * mrSges(6,1) + (-t117 * t167 - t170 * t71) * mrSges(6,2)) - g(2) * ((t115 * t170 - t167 * t69) * mrSges(6,1) + (-t115 * t167 - t170 * t69) * mrSges(6,2)) - g(3) * ((-t103 * t167 - t218) * mrSges(6,1) + (-t103 * t170 + t220) * mrSges(6,2)) + (t14 * t90 + t15 * t91) * mrSges(6,3) + t228 + (-Ifges(6,2) * t91 + t27 + t89) * t286 + t293;];
tau = t4;
