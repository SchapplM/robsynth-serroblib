% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:23:21
% EndTime: 2019-12-31 19:23:47
% DurationCPUTime: 12.76s
% Computational Cost: add. (2611->476), mult. (8281->628), div. (0->0), fcn. (5750->6), ass. (0->215)
t171 = sin(pkin(8));
t175 = sin(qJ(2));
t219 = qJD(1) * t175;
t209 = t171 * t219;
t174 = cos(pkin(5));
t172 = sin(pkin(5));
t217 = qJD(2) * t172;
t176 = cos(qJ(2));
t218 = qJD(1) * t176;
t149 = t174 * t218 + t217;
t173 = cos(pkin(8));
t266 = t149 * t173;
t97 = t209 - t266;
t259 = -t97 / 0.2e1;
t258 = t97 / 0.2e1;
t221 = t175 * t173;
t222 = t174 * t176;
t136 = t171 * t222 + t221;
t98 = qJD(1) * t136 + t171 * t217;
t257 = -t98 / 0.2e1;
t256 = t98 / 0.2e1;
t212 = t174 * qJD(2);
t148 = t172 * t218 - t212;
t244 = -t148 / 0.2e1;
t283 = t148 / 0.2e1;
t279 = Ifges(6,4) + Ifges(5,5);
t278 = Ifges(4,5) + Ifges(6,5);
t271 = Ifges(5,1) + Ifges(6,1) + Ifges(4,3);
t201 = qJ(3) * t174 + pkin(7);
t189 = qJD(1) * t201;
t144 = t175 * t189;
t131 = qJD(2) * pkin(2) - t144;
t265 = t172 * t175;
t150 = -pkin(2) * t176 - qJ(3) * t265 - pkin(1);
t132 = t150 * qJD(1);
t67 = -t131 * t172 + t132 * t174 + qJD(3);
t180 = -qJ(4) * t98 + t67;
t16 = pkin(3) * t97 + t180;
t118 = pkin(7) * t218 + qJ(3) * t149;
t30 = -t171 * t118 + t173 * (t131 * t174 + t132 * t172);
t177 = qJD(4) - t30;
t22 = pkin(3) * t148 + t177;
t280 = -Ifges(4,4) + Ifges(6,6);
t281 = Ifges(4,1) + Ifges(6,3);
t240 = pkin(3) + qJ(5);
t6 = pkin(4) * t98 + t148 * t240 + t177;
t7 = t240 * t97 + t180;
t285 = -t22 * mrSges(5,1) - t6 * mrSges(6,1) - t67 * mrSges(4,2) + t7 * mrSges(6,2) + t30 * mrSges(4,3) + t16 * mrSges(5,3) + Ifges(5,4) * t244 + Ifges(5,6) * t258 + t280 * t259 + t278 * t283 + (Ifges(5,2) + t281) * t257;
t227 = t171 * t174;
t228 = t171 * t172;
t31 = t118 * t173 + t131 * t227 + t132 * t228;
t23 = qJ(4) * t148 - t31;
t11 = -pkin(4) * t97 + qJD(5) - t23;
t276 = -Ifges(5,6) + Ifges(6,6);
t277 = Ifges(6,2) + Ifges(5,3);
t284 = t67 * mrSges(4,1) + t23 * mrSges(5,1) - t11 * mrSges(6,1) - t16 * mrSges(5,2) - t31 * mrSges(4,3) + t7 * mrSges(6,3) + Ifges(4,4) * t257 + Ifges(4,6) * t283 + t279 * t244 + t276 * t256 + (Ifges(4,2) + t277) * t258;
t282 = -t30 * mrSges(4,1) + t31 * mrSges(4,2) - t22 * mrSges(5,2) - t11 * mrSges(6,2) + t23 * mrSges(5,3) + t6 * mrSges(6,3) + Ifges(5,4) * t256 + Ifges(4,6) * t258 + t257 * t278 + t259 * t279 + t271 * t283;
t184 = t171 * t176 + t174 * t221;
t123 = t184 * qJD(1);
t200 = t174 * t209;
t124 = t173 * t218 - t200;
t224 = t172 * t176;
t185 = pkin(2) * t175 - qJ(3) * t224;
t143 = t185 * qJD(1);
t145 = t176 * t189;
t76 = t143 * t174 + t145 * t172;
t183 = -qJ(4) * t124 + t76;
t273 = -t123 * t240 + (-qJD(4) * t171 - qJD(5) * t173) * t172 - t183;
t128 = t171 * t144;
t223 = t173 * t174;
t190 = t145 * t223 - t128;
t204 = t240 * t175;
t214 = qJD(3) * t172;
t229 = t143 * t173;
t272 = -qJD(5) * t174 + t171 * t214 - pkin(4) * t124 - (-qJD(1) * t204 - t229) * t172 - t190;
t270 = -Ifges(5,4) + t278;
t269 = -Ifges(4,6) + t279;
t268 = -Ifges(5,4) / 0.2e1;
t125 = t184 * qJD(2);
t113 = qJD(1) * t125;
t267 = -t113 / 0.2e1;
t215 = qJD(2) * t176;
t206 = t173 * t215;
t114 = qJD(1) * t206 - qJD(2) * t200;
t251 = t114 / 0.2e1;
t151 = t201 * t175;
t264 = t173 * (t150 * t172 - t151 * t174);
t263 = t172 ^ 2 + t174 ^ 2;
t216 = qJD(2) * t175;
t213 = qJD(3) * t175;
t110 = qJD(2) * t185 - t172 * t213;
t188 = qJD(2) * t201;
t111 = qJD(3) * t222 - t175 * t188;
t112 = -t174 * t213 - t176 * t188;
t27 = t110 * t228 + t111 * t173 + t112 * t227;
t17 = -t172 * (qJ(4) * t216 - qJD(4) * t176) - t27;
t211 = qJD(1) * qJD(2);
t203 = t175 * t211;
t199 = t172 * t203;
t260 = Ifges(5,2) * t251 + Ifges(5,6) * t267 + t199 * t268;
t95 = t110 * qJD(1);
t96 = t112 * qJD(1);
t48 = -t172 * t96 + t174 * t95;
t255 = m(4) * t48;
t254 = pkin(1) * mrSges(3,1);
t253 = pkin(1) * mrSges(3,2);
t241 = 0.3e1 / 0.2e1 * t176;
t71 = mrSges(5,1) * t97 + mrSges(5,3) * t148;
t72 = -mrSges(6,1) * t97 - mrSges(6,2) * t148;
t239 = -t71 + t72;
t69 = -mrSges(4,1) * t148 - mrSges(4,3) * t98;
t73 = mrSges(5,1) * t98 - mrSges(5,2) * t148;
t238 = t73 - t69;
t237 = Ifges(3,4) * t175;
t236 = Ifges(3,2) * t176;
t235 = t173 * t95;
t234 = Ifges(3,5) * qJD(2);
t233 = Ifges(3,6) * qJD(2);
t232 = qJD(2) * mrSges(3,1);
t231 = qJD(2) * mrSges(3,2);
t230 = t110 * t173;
t226 = t171 * t175;
t225 = t172 * t173;
t152 = t201 * t176;
t139 = t171 * t152;
t220 = pkin(3) * t224 + t139;
t138 = pkin(2) * t227 + qJ(3) * t225;
t88 = qJD(1) * t111 + qJD(2) * t214;
t15 = t173 * t88 + t227 * t96 + t228 * t95;
t46 = t143 * t228 - t144 * t173 - t145 * t227;
t62 = t150 * t228 - t151 * t227 + t152 * t173;
t210 = -pkin(2) * t173 - pkin(3);
t202 = -qJ(4) * t171 - pkin(2);
t82 = t114 * mrSges(5,1) + mrSges(5,2) * t199;
t81 = -t113 * mrSges(6,1) + mrSges(6,2) * t199;
t64 = -t114 * mrSges(6,2) + mrSges(6,3) * t113;
t63 = t110 * t174 - t112 * t172;
t85 = t150 * t174 + t151 * t172;
t197 = Ifges(6,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1;
t196 = -Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1 + Ifges(4,6) / 0.2e1;
t195 = Ifges(6,5) / 0.2e1 + Ifges(4,5) / 0.2e1 + t268;
t194 = -Ifges(6,6) / 0.2e1 + Ifges(4,4) / 0.2e1 + Ifges(5,6) / 0.2e1;
t77 = t171 * t88;
t193 = -t223 * t96 + t77;
t192 = qJD(2) * t204;
t99 = t171 * t111;
t191 = -t112 * t223 + t99;
t119 = -qJ(4) * t174 - t138;
t182 = -qJ(4) * t136 + t85;
t44 = qJ(4) * t224 - t62;
t79 = mrSges(6,1) * t114 - mrSges(6,3) * t199;
t41 = -qJ(4) * t172 * t219 - t46;
t179 = -qJ(4) * t114 - qJD(4) * t98 + t48;
t126 = -t212 * t226 + t206;
t178 = -qJ(4) * t126 - qJD(4) * t136 + t63;
t9 = -qJ(4) * t199 + qJD(4) * t148 - t15;
t169 = Ifges(3,4) * t218;
t165 = qJ(3) * t228;
t162 = mrSges(3,3) * t218 - t231;
t161 = -mrSges(3,3) * t219 + t232;
t147 = qJD(4) * t174 + t173 * t214;
t142 = Ifges(3,1) * t219 + t169 + t234;
t141 = t233 + (t236 + t237) * qJD(1);
t137 = pkin(2) * t223 - t165;
t135 = -t173 * t222 + t226;
t121 = (-pkin(3) * t173 + t202) * t172;
t120 = t174 * t210 + t165;
t109 = t114 * mrSges(5,3);
t107 = t114 * mrSges(4,2);
t94 = (-t173 * t240 + t202) * t172;
t93 = pkin(4) * t225 - t119;
t89 = pkin(4) * t228 + t165 + (-qJ(5) + t210) * t174;
t84 = mrSges(4,1) * t199 - mrSges(4,3) * t114;
t83 = -mrSges(4,2) * t199 - mrSges(4,3) * t113;
t80 = mrSges(5,1) * t113 - mrSges(5,3) * t199;
t75 = -t209 * t263 + t266;
t74 = t173 * t219 * t263 + t149 * t171;
t70 = mrSges(6,1) * t98 + mrSges(6,3) * t148;
t68 = mrSges(4,2) * t148 - mrSges(4,3) * t97;
t66 = -t113 * mrSges(5,2) - t109;
t65 = t113 * mrSges(4,1) + t107;
t61 = -t139 + t264;
t60 = t114 * Ifges(4,1) - Ifges(4,4) * t113 + Ifges(4,5) * t199;
t59 = t114 * Ifges(4,4) - t113 * Ifges(4,2) + Ifges(4,6) * t199;
t58 = t114 * Ifges(4,5) - t113 * Ifges(4,6) + Ifges(4,3) * t199;
t57 = Ifges(5,1) * t199 - t114 * Ifges(5,4) + t113 * Ifges(5,5);
t56 = Ifges(6,1) * t199 + t113 * Ifges(6,4) + t114 * Ifges(6,5);
t54 = Ifges(6,4) * t199 + t113 * Ifges(6,2) + t114 * Ifges(6,6);
t53 = Ifges(5,5) * t199 - t114 * Ifges(5,6) + t113 * Ifges(5,3);
t52 = Ifges(6,5) * t199 + Ifges(6,6) * t113 + t114 * Ifges(6,3);
t51 = -mrSges(5,2) * t97 - mrSges(5,3) * t98;
t50 = mrSges(4,1) * t97 + mrSges(4,2) * t98;
t49 = -mrSges(6,2) * t98 + mrSges(6,3) * t97;
t47 = t220 - t264;
t45 = t128 + (t143 * t172 - t145 * t174) * t173;
t43 = pkin(3) * t135 + t182;
t42 = (-pkin(3) * t219 - t229) * t172 + t190;
t29 = pkin(3) * t123 + t183;
t28 = -pkin(4) * t135 - t44;
t26 = -t99 + (t110 * t172 + t112 * t174) * t173;
t25 = t135 * t240 + t182;
t24 = t151 * t223 + t136 * pkin(4) + (qJ(5) * t176 - t150 * t173) * t172 + t220;
t21 = -pkin(4) * t123 - t41;
t20 = (-pkin(3) * t216 - t230) * t172 + t191;
t14 = -t77 + (t172 * t95 + t174 * t96) * t173;
t13 = pkin(3) * t125 + t178;
t12 = (-pkin(3) * t203 - t235) * t172 + t193;
t10 = -pkin(4) * t125 - t17;
t8 = t126 * pkin(4) + (qJD(5) * t176 - t192 - t230) * t172 + t191;
t5 = pkin(3) * t113 + t179;
t4 = qJD(5) * t135 + t125 * t240 + t178;
t3 = -pkin(4) * t113 - t9;
t2 = pkin(4) * t114 + qJD(5) * t148 + (-qJD(1) * t192 - t235) * t172 + t193;
t1 = qJD(5) * t97 + t113 * t240 + t179;
t18 = [-(t58 + t57 + t56) * t224 / 0.2e1 + (t60 + t52) * t136 / 0.2e1 + (t54 + t53) * t135 / 0.2e1 + qJD(2) ^ 2 * (Ifges(3,5) * t176 - Ifges(3,6) * t175) / 0.2e1 - t135 * t59 / 0.2e1 + t48 * (mrSges(4,1) * t135 + mrSges(4,2) * t136) + t1 * (-mrSges(6,2) * t136 + mrSges(6,3) * t135) + t5 * (-mrSges(5,2) * t135 - mrSges(5,3) * t136) + t24 * t79 + t44 * t80 + t28 * t81 + t47 * t82 + t62 * t83 + t61 * t84 + t85 * t65 + t25 * t64 + t43 * t66 + t27 * t68 + t26 * t69 + t8 * t70 + t17 * t71 + t10 * t72 + t20 * t73 + t4 * t49 + t13 * t51 + t63 * t50 + (t135 * t277 + t136 * t276 - t224 * t279) * t113 / 0.2e1 + (t135 * t280 + t136 * t281 - t224 * t278) * t251 + (-t161 * t176 - t162 * t175) * qJD(2) * pkin(7) + ((Ifges(5,4) * t257 + Ifges(4,6) * t259 + t244 * t271 + t256 * t278 + t258 * t279 - t282) * t172 - t141 / 0.2e1) * t216 + ((Ifges(3,4) * t241 - 0.2e1 * t253) * t176 + (Ifges(3,1) * t241 - 0.3e1 / 0.2e1 * t236 - 0.2e1 * t254 - 0.3e1 / 0.2e1 * t237) * t175 + (t135 * t269 + t136 * t270 - t224 * t271) * t265 / 0.2e1) * t211 + m(4) * (t14 * t61 + t15 * t62 + t26 * t30 + t27 * t31 + t48 * t85 + t63 * t67) + m(5) * (t12 * t47 + t13 * t16 + t17 * t23 + t20 * t22 + t43 * t5 + t44 * t9) + m(6) * (t1 * t25 + t10 * t11 + t2 * t24 + t28 * t3 + t4 * t7 + t6 * t8) + (Ifges(4,4) * t136 - Ifges(4,2) * t135 - Ifges(4,6) * t224) * t267 + t136 * t260 + (-Ifges(4,2) * t259 + Ifges(5,6) * t257 + t269 * t244 + t280 * t256 + t277 * t258 + t284) * t125 + (Ifges(4,4) * t259 - Ifges(5,2) * t257 + t270 * t244 + t281 * t256 + t276 * t258 - t285) * t126 + t142 * t215 / 0.2e1 - t114 * (-Ifges(5,4) * t224 - Ifges(5,2) * t136 + Ifges(5,6) * t135) / 0.2e1 + t3 * (-mrSges(6,1) * t135 - mrSges(6,2) * t224) + t12 * (mrSges(5,1) * t136 - mrSges(5,2) * t224) + t14 * (-mrSges(4,1) * t224 - mrSges(4,3) * t136) + t2 * (mrSges(6,1) * t136 + mrSges(6,3) * t224) + t9 * (mrSges(5,1) * t135 + mrSges(5,3) * t224) + t15 * (mrSges(4,2) * t224 - mrSges(4,3) * t135); ((t234 / 0.2e1 + qJD(1) * t253 - t169 / 0.2e1 - t142 / 0.2e1 + (t161 - t232) * pkin(7)) * t176 + (t141 / 0.2e1 - t233 / 0.2e1 + (t237 / 0.2e1 + t254 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t176) * qJD(1) + (t162 + t231) * pkin(7) + (-t195 * t98 + t196 * t97 + t197 * t148 + ((t171 * t195 + t173 * t196) * t172 + t197 * t174) * qJD(2) + t282) * t172) * t175) * qJD(1) + t137 * t84 + t138 * t83 + t119 * t80 + t120 * t82 + t121 * t66 + m(4) * (t137 * t14 + t138 * t15) + m(5) * (t119 * t9 + t12 * t120 + t121 * t5 - t147 * t23) + t93 * t81 + t94 * t64 - t76 * t50 + t89 * t79 - t46 * t68 - t45 * t69 - t41 * t71 - t21 * t72 - t42 * t73 - t29 * t51 + t272 * t70 + (t1 * t94 + t2 * t89 + t3 * t93 + t273 * t7 + t272 * t6 + (t147 - t21) * t11) * m(6) + t273 * t49 + (-t15 * mrSges(4,2) - t9 * mrSges(5,3) + t3 * mrSges(6,2) - t2 * mrSges(6,3) + t12 * mrSges(5,2) + t14 * mrSges(4,1) + t56 / 0.2e1 + t57 / 0.2e1 + t58 / 0.2e1 + t195 * t114 - t196 * t113) * t174 - m(4) * (t30 * t45 + t31 * t46 + t67 * t76) - m(5) * (t16 * t29 + t22 * t42 + t23 * t41) + ((-t65 - t255) * pkin(2) + (t15 * mrSges(4,3) - t9 * mrSges(5,1) + t3 * mrSges(6,1) - t53 / 0.2e1 - t54 / 0.2e1 + t59 / 0.2e1 - t48 * mrSges(4,1) + t5 * mrSges(5,2) - t1 * mrSges(6,3) + (m(4) * t31 + t68) * qJD(3) + t194 * t114 + (-Ifges(4,2) / 0.2e1 - Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1) * t113) * t173 + (-t14 * mrSges(4,3) + t2 * mrSges(6,1) + t12 * mrSges(5,1) + t52 / 0.2e1 + t260 + t60 / 0.2e1 + t48 * mrSges(4,2) - t5 * mrSges(5,3) - t1 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t114 - t194 * t113 + (-m(5) * t16 - t51) * qJD(4) + (-m(4) * t30 + m(5) * t22 + t238) * qJD(3)) * t171) * t172 + (-Ifges(4,2) * t258 + Ifges(5,6) * t256 + t280 * t257 + t277 * t259 + t269 * t283 - t284) * t123 + (Ifges(4,4) * t258 - Ifges(5,2) * t256 + t281 * t257 + t276 * t259 + t270 * t283 + t285) * t124 + t239 * t147; t107 - t109 + (mrSges(4,1) - mrSges(5,2)) * t113 + (-t68 - t239) * t75 + (-t70 - t238) * t74 + t255 - m(4) * (-t30 * t74 + t31 * t75) + t64 + (-t11 * t75 - t6 * t74 + t1) * m(6) + (-t22 * t74 + t23 * t75 + t5) * m(5); (t49 + t51) * t98 + t239 * t148 + t79 + t82 + (t11 * t148 + t7 * t98 + t2) * m(6) + (-t148 * t23 + t16 * t98 + t12) * m(5); -t148 * t70 - t97 * t49 + 0.2e1 * (t3 / 0.2e1 + t6 * t244 + t7 * t259) * m(6) + t81;];
tauc = t18(:);
