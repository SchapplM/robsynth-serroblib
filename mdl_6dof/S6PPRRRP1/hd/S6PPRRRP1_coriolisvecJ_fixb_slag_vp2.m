% Calculate vector of centrifugal and coriolis load on the joints for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2018-11-23 14:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PPRRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:50:44
% EndTime: 2018-11-23 14:50:50
% DurationCPUTime: 5.69s
% Computational Cost: add. (4253->434), mult. (11566->598), div. (0->0), fcn. (9217->12), ass. (0->210)
t288 = Ifges(6,4) + Ifges(7,4);
t289 = Ifges(6,1) + Ifges(7,1);
t279 = Ifges(6,5) + Ifges(7,5);
t287 = Ifges(6,2) + Ifges(7,2);
t278 = Ifges(6,6) + Ifges(7,6);
t286 = qJD(4) / 0.2e1;
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t198 = qJD(4) * t145;
t143 = sin(qJ(4));
t203 = qJD(3) * t143;
t115 = -t142 * t203 + t198;
t285 = t288 * t115;
t199 = qJD(4) * t142;
t116 = t145 * t203 + t199;
t284 = t288 * t116;
t283 = t288 * t145;
t282 = t288 * t142;
t146 = cos(qJ(4));
t201 = qJD(3) * t146;
t135 = Ifges(5,4) * t201;
t182 = Ifges(5,5) * t286;
t141 = cos(pkin(6));
t128 = qJD(1) * t141 + qJD(2);
t137 = sin(pkin(7));
t140 = cos(pkin(7));
t139 = cos(pkin(12));
t138 = sin(pkin(6));
t204 = qJD(1) * t138;
t185 = t139 * t204;
t94 = t128 * t140 - t137 * t185;
t220 = t143 * t94;
t136 = sin(pkin(12));
t144 = sin(qJ(3));
t147 = cos(qJ(3));
t210 = t139 * t140;
t154 = (t136 * t147 + t144 * t210) * t138;
t212 = t137 * t144;
t59 = qJD(1) * t154 + t128 * t212;
t57 = qJD(3) * pkin(9) + t59;
t36 = t146 * t57 + t220;
t33 = qJD(4) * pkin(10) + t36;
t122 = -pkin(4) * t146 - pkin(10) * t143 - pkin(3);
t213 = t136 * t144;
t58 = t147 * (t128 * t137 + t140 * t185) - t204 * t213;
t50 = qJD(3) * t122 - t58;
t11 = -t142 * t33 + t145 * t50;
t12 = t142 * t50 + t145 * t33;
t162 = t11 * t145 + t12 * t142;
t173 = mrSges(7,1) * t142 + mrSges(7,2) * t145;
t175 = mrSges(6,1) * t142 + mrSges(6,2) * t145;
t35 = -t143 * t57 + t146 * t94;
t32 = -qJD(4) * pkin(4) - t35;
t24 = -pkin(5) * t115 + qJD(6) + t32;
t242 = t145 / 0.2e1;
t245 = -t142 / 0.2e1;
t248 = t116 / 0.2e1;
t131 = qJD(5) - t201;
t263 = t289 * t116 + t279 * t131 + t285;
t264 = t287 * t115 + t278 * t131 + t284;
t265 = t131 / 0.2e1;
t266 = t115 / 0.2e1;
t268 = t289 * t145 - t282;
t270 = -t287 * t142 + t283;
t8 = -qJ(6) * t116 + t11;
t7 = pkin(5) * t131 + t8;
t9 = qJ(6) * t115 + t12;
t259 = t162 * mrSges(6,3) + (t9 * t142 + t7 * t145) * mrSges(7,3) - t24 * t173 - t32 * t175 - t270 * t266 - t268 * t248 - (-t278 * t142 + t279 * t145) * t265 - t264 * t245 - t263 * t242;
t280 = t203 / 0.2e1;
t56 = -qJD(3) * pkin(3) - t58;
t281 = t56 * mrSges(5,2) + Ifges(5,1) * t280 + t135 / 0.2e1 + t182 - t35 * mrSges(5,3) - t259;
t207 = t145 * t146;
t160 = pkin(5) * t143 - qJ(6) * t207;
t233 = -qJ(6) - pkin(10);
t179 = qJD(5) * t233;
t178 = pkin(4) * t143 - pkin(10) * t146;
t118 = t178 * qJD(3);
t22 = t145 * t118 - t142 * t35;
t277 = -qJD(3) * t160 - qJD(6) * t142 + t145 * t179 - t22;
t194 = qJD(6) * t145;
t23 = t142 * t118 + t145 * t35;
t276 = t194 - t23 + (qJ(6) * t201 + t179) * t142;
t193 = qJD(3) * qJD(4);
t180 = t143 * t193;
t192 = qJD(4) * qJD(5);
t196 = qJD(5) * t142;
t197 = qJD(4) * t146;
t87 = t145 * t192 + (-t143 * t196 + t145 * t197) * qJD(3);
t195 = qJD(5) * t145;
t156 = t142 * t197 + t143 * t195;
t88 = -qJD(3) * t156 - t142 * t192;
t275 = t278 * t180 + t287 * t88 + t288 * t87;
t274 = t279 * t180 + t288 * t88 + t289 * t87;
t273 = qJD(4) * t36;
t155 = (t147 * t210 - t213) * t138;
t211 = t137 * t147;
t272 = t141 * t211 + t155;
t271 = t287 * t145 + t282;
t269 = t289 * t142 + t283;
t267 = -m(5) * t35 + m(6) * t32;
t181 = -Ifges(5,6) * qJD(4) / 0.2e1;
t132 = pkin(9) * t207;
t96 = t142 * t122 + t132;
t262 = t279 * t142 + t278 * t145;
t54 = (qJD(1) * t155 + t128 * t211) * qJD(3);
t14 = qJD(4) * t35 + t146 * t54;
t119 = t178 * qJD(4);
t44 = (t119 + t59) * qJD(3);
t3 = t145 * t14 + t142 * t44 + t50 * t195 - t196 * t33;
t4 = -qJD(5) * t12 - t14 * t142 + t145 * t44;
t261 = -t142 * t4 + t145 * t3;
t188 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t189 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t190 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t258 = -t189 * t115 - t190 * t116 - t188 * t131 - t11 * mrSges(6,1) - t56 * mrSges(5,1) - t7 * mrSges(7,1) - t181 + (Ifges(5,4) * t143 + t146 * Ifges(5,2)) * qJD(3) / 0.2e1 + t12 * mrSges(6,2) + t36 * mrSges(5,3) + t9 * mrSges(7,2) - t278 * t266 - (Ifges(7,3) + Ifges(6,3)) * t265 - t279 * t248;
t257 = 0.2e1 * m(5);
t256 = t87 / 0.2e1;
t255 = t88 / 0.2e1;
t252 = m(7) * t24;
t251 = -t115 / 0.2e1;
t249 = -t116 / 0.2e1;
t247 = -t131 / 0.2e1;
t241 = pkin(5) * t142;
t15 = t143 * t54 + t273;
t100 = -t137 * t138 * t139 + t140 * t141;
t74 = t141 * t212 + t154;
t161 = t146 * t100 - t143 * t74;
t238 = t15 * t161;
t55 = qJD(3) * t59;
t236 = t55 * t272;
t48 = -t88 * mrSges(7,1) + t87 * mrSges(7,2);
t49 = -mrSges(6,1) * t88 + mrSges(6,2) * t87;
t232 = t48 + t49;
t68 = mrSges(7,1) * t180 - mrSges(7,3) * t87;
t69 = mrSges(6,1) * t180 - mrSges(6,3) * t87;
t231 = t68 + t69;
t70 = -mrSges(7,2) * t180 + mrSges(7,3) * t88;
t71 = -mrSges(6,2) * t180 + mrSges(6,3) * t88;
t230 = t70 + t71;
t89 = -mrSges(7,2) * t131 + t115 * mrSges(7,3);
t90 = -mrSges(6,2) * t131 + mrSges(6,3) * t115;
t229 = t90 + t89;
t91 = mrSges(7,1) * t131 - mrSges(7,3) * t116;
t92 = mrSges(6,1) * t131 - mrSges(6,3) * t116;
t228 = t91 + t92;
t101 = -t146 * t140 + t143 * t212;
t223 = t101 * t15;
t218 = t147 * t55;
t217 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t115 - mrSges(6,2) * t116 - mrSges(5,3) * t203;
t214 = qJ(6) * t143;
t209 = t142 * t146;
t208 = t143 * t145;
t206 = t142 * t119 + t122 * t195;
t205 = t143 * pkin(9) * t199 + t145 * t119;
t202 = qJD(3) * t144;
t200 = qJD(3) * t147;
t75 = -mrSges(7,1) * t115 + mrSges(7,2) * t116;
t191 = t75 - t217;
t187 = mrSges(5,3) * t201;
t184 = t137 * t202;
t183 = t137 * t200;
t176 = mrSges(6,1) * t145 - mrSges(6,2) * t142;
t174 = mrSges(7,1) * t145 - mrSges(7,2) * t142;
t47 = t100 * t143 + t146 * t74;
t21 = -t142 * t272 + t145 * t47;
t20 = -t142 * t47 - t145 * t272;
t102 = t140 * t143 + t146 * t212;
t79 = -t102 * t142 - t145 * t211;
t159 = -t102 * t145 + t142 * t211;
t1 = pkin(5) * t180 - qJ(6) * t87 - qJD(6) * t116 + t4;
t2 = qJ(6) * t88 + qJD(6) * t115 + t3;
t153 = -t4 * mrSges(6,1) - t1 * mrSges(7,1) + t3 * mrSges(6,2) + t2 * mrSges(7,2);
t148 = qJD(3) ^ 2;
t134 = -pkin(5) * t145 - pkin(4);
t130 = Ifges(6,3) * t180;
t129 = Ifges(7,3) * t180;
t126 = t233 * t145;
t125 = t233 * t142;
t124 = -qJD(4) * mrSges(5,2) + t187;
t120 = (pkin(9) + t241) * t143;
t117 = (-mrSges(5,1) * t146 + mrSges(5,2) * t143) * qJD(3);
t114 = t145 * t122;
t109 = (mrSges(5,1) * t143 + mrSges(5,2) * t146) * t193;
t95 = -pkin(9) * t209 + t114;
t93 = pkin(5) * t156 + pkin(9) * t197;
t86 = Ifges(6,5) * t87;
t85 = Ifges(7,5) * t87;
t84 = Ifges(6,6) * t88;
t83 = Ifges(7,6) * t88;
t81 = -t142 * t214 + t96;
t78 = qJD(4) * t102 + t143 * t183;
t77 = -qJD(4) * t101 + t146 * t183;
t72 = -qJ(6) * t208 + t114 + (-pkin(9) * t142 - pkin(5)) * t146;
t67 = t74 * qJD(3);
t66 = t272 * qJD(3);
t52 = -qJD(5) * t96 + t205;
t51 = (-t143 * t198 - t146 * t196) * pkin(9) + t206;
t34 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t208 + (-qJD(6) * t143 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t146) * t142 + t206;
t31 = qJD(5) * t159 - t142 * t77 + t145 * t184;
t30 = qJD(5) * t79 + t142 * t184 + t145 * t77;
t28 = -t143 * t194 + t160 * qJD(4) + (-t132 + (-t122 + t214) * t142) * qJD(5) + t205;
t27 = t220 + (qJD(3) * t241 + t57) * t146;
t26 = t142 * t59 + t207 * t58;
t25 = t145 * t59 - t209 * t58;
t19 = qJD(4) * t47 + t143 * t66;
t18 = qJD(4) * t161 + t146 * t66;
t10 = -pkin(5) * t88 + t15;
t6 = -qJD(5) * t21 - t142 * t18 + t145 * t67;
t5 = qJD(5) * t20 + t142 * t67 + t145 * t18;
t13 = [-t272 * t109 + t67 * t117 + t18 * t124 + t228 * t6 + t229 * t5 - t232 * t161 + t230 * t21 + t231 * t20 + t191 * t19 + (-t67 * mrSges(4,1) - t66 * mrSges(4,2) + (-t143 * t47 - t146 * t161) * qJD(4) * mrSges(5,3)) * qJD(3) + m(7) * (t1 * t20 - t10 * t161 + t19 * t24 + t2 * t21 + t5 * t9 + t6 * t7) + m(4) * (t54 * t74 - t58 * t67 + t59 * t66 - t236) + m(6) * (t11 * t6 + t12 * t5 + t19 * t32 + t20 * t4 + t21 * t3 - t238) + m(5) * (t14 * t47 + t18 * t36 - t19 * t35 + t56 * t67 - t236 - t238); -t102 * mrSges(5,3) * t180 + t77 * t124 - t230 * t159 + t231 * t79 + t228 * t31 + t229 * t30 + t191 * t78 + (qJD(4) * t187 + t232) * t101 + m(7) * (t1 * t79 + t10 * t101 - t159 * t2 + t24 * t78 + t30 * t9 + t31 * t7) + m(5) * (t102 * t14 - t35 * t78 + t36 * t77 + t223) + m(6) * (t11 * t31 + t12 * t30 - t159 * t3 + t32 * t78 + t4 * t79 + t223) + ((-t148 * mrSges(4,2) - t109) * t147 + (-t148 * mrSges(4,1) + qJD(3) * t117) * t144 + m(4) * (t144 * t54 + t200 * t59 - t202 * t58 - t218) + m(5) * (t202 * t56 - t218)) * t137; -pkin(3) * t109 - t59 * t117 + t120 * t48 + t28 * t91 + t34 * t89 + t51 * t90 + t52 * t92 + t72 * t68 + t95 * t69 + t81 * t70 + t96 * t71 + t93 * t75 - t229 * t26 - t228 * t25 + (qJD(3) * t58 - t54) * mrSges(4,2) - m(6) * (t11 * t25 + t12 * t26) + m(6) * (t11 * t52 + t12 * t51 + t3 * t96 + t4 * t95) + (-pkin(3) * t55 / 0.2e1 - t56 * t59 / 0.2e1) * t257 + (t14 * mrSges(5,3) - t58 * t124 - t55 * mrSges(5,1) - t129 / 0.2e1 - t130 / 0.2e1 - t85 / 0.2e1 - t86 / 0.2e1 - t83 / 0.2e1 - t84 / 0.2e1 - t189 * t88 - t190 * t87 + (pkin(9) * t14 / 0.2e1 - t36 * t58 / 0.2e1) * t257 + ((-t217 + t267) * pkin(9) + t182 + 0.3e1 / 0.2e1 * t135 + t281) * qJD(4) + t153) * t146 + (t10 * t173 + t55 * mrSges(5,2) + (mrSges(5,3) + t175) * t15 + (-t1 * t145 - t2 * t142) * mrSges(7,3) + (-t3 * t142 - t4 * t145) * mrSges(6,3) + (((-0.3e1 / 0.2e1 * Ifges(5,4) + t190 * t145 - t189 * t142) * t143 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2) - t188) * t146) * qJD(3) + t181 - t258) * qJD(4) + (t24 * t174 + t32 * t176 + (t7 * t142 - t9 * t145) * mrSges(7,3) + (t11 * t142 - t12 * t145) * mrSges(6,3) + t271 * t251 + t269 * t249 + t262 * t247 - t264 * t145 / 0.2e1) * qJD(5) + t268 * t256 + t270 * t255 + (qJD(5) * t263 + t275) * t245 + t274 * t242 + (m(6) * t15 - qJD(4) * t124 + t49 + (t15 - t273) * m(5)) * pkin(9) + (-t191 - t252 - t267) * t58) * t143 + (t1 * t72 + t10 * t120 + t2 * t81 + t24 * t93 + (-t26 + t34) * t9 + (-t25 + t28) * t7) * m(7); m(6) * (-pkin(4) * t15 + pkin(10) * t261) + t261 * mrSges(6,3) + ((Ifges(5,4) * t280 + t262 * t286 + t181 + t258) * t143 + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t203 - t135 / 0.2e1 + t182 - t281) * t146) * qJD(3) - t10 * t174 - m(6) * (t11 * t22 + t12 * t23 + t32 * t36) + (-t259 + (t75 + t252) * t241 + (-m(6) * t162 - t142 * t90 - t145 * t92) * pkin(10)) * qJD(5) + (-t1 * t142 + t145 * t2) * mrSges(7,3) + (-t142 * t69 + t145 * t71) * pkin(10) + (-mrSges(5,1) - t176) * t15 - t14 * mrSges(5,2) - pkin(4) * t49 - t27 * t75 - t23 * t90 - t22 * t92 - t35 * t124 + t125 * t68 - t126 * t70 + t134 * t48 + t217 * t36 + t269 * t256 + t271 * t255 + t274 * t142 / 0.2e1 + t275 * t242 + t276 * t89 + t277 * t91 + (t1 * t125 + t10 * t134 - t126 * t2 - t24 * t27 + t276 * t9 + t277 * t7) * m(7); -t153 + (-t116 * t75 + t68) * pkin(5) + (t115 * t7 + t116 * t9) * mrSges(7,3) + (t11 * t115 + t116 * t12) * mrSges(6,3) + t129 + t130 + t85 + t86 + t83 + t84 + (-(-t7 + t8) * t9 + (-t116 * t24 + t1) * pkin(5)) * m(7) - t8 * t89 - t11 * t90 + t9 * t91 + t12 * t92 - t24 * (t116 * mrSges(7,1) + t115 * mrSges(7,2)) - t32 * (mrSges(6,1) * t116 + mrSges(6,2) * t115) + (t289 * t115 - t284) * t249 + t264 * t248 + (t279 * t115 - t278 * t116) * t247 + (-t287 * t116 + t263 + t285) * t251; -t115 * t89 + t116 * t91 + 0.2e1 * (t10 / 0.2e1 + t9 * t251 + t7 * t248) * m(7) + t48;];
tauc  = t13(:);
