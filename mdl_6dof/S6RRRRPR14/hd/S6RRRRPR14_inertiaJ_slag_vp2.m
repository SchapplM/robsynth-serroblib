% Calculate joint inertia matrix for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:23:42
% EndTime: 2018-11-23 18:23:44
% DurationCPUTime: 2.92s
% Computational Cost: add. (8343->577), mult. (21577->833), div. (0->0), fcn. (24687->14), ass. (0->211)
t278 = 2 * pkin(11);
t219 = cos(pkin(13));
t257 = t219 / 0.2e1;
t225 = sin(qJ(2));
t218 = sin(pkin(6));
t229 = cos(qJ(2));
t243 = t218 * t229;
t221 = cos(pkin(6));
t256 = pkin(1) * t221;
t164 = pkin(9) * t243 + t225 * t256;
t217 = sin(pkin(7));
t220 = cos(pkin(7));
t240 = t220 * t229;
t238 = t218 * t240;
t108 = (t217 * t221 + t238) * pkin(10) + t164;
t203 = t229 * t256;
t244 = t218 * t225;
t119 = pkin(2) * t221 + t203 + (-pkin(10) * t220 - pkin(9)) * t244;
t134 = (-pkin(10) * t217 * t225 - pkin(2) * t229 - pkin(1)) * t218;
t224 = sin(qJ(3));
t228 = cos(qJ(3));
t53 = -t224 * t108 + (t119 * t220 + t134 * t217) * t228;
t245 = t217 * t228;
t113 = -t221 * t245 + t224 * t244 - t228 * t238;
t216 = sin(pkin(13));
t246 = t217 * t224;
t114 = t221 * t246 + (t224 * t240 + t225 * t228) * t218;
t156 = -t217 * t243 + t220 * t221;
t223 = sin(qJ(4));
t227 = cos(qJ(4));
t87 = t114 * t227 + t156 * t223;
t56 = t113 * t219 - t216 * t87;
t57 = t113 * t216 + t219 * t87;
t86 = t114 * t223 - t227 * t156;
t22 = Ifges(6,4) * t57 + Ifges(6,2) * t56 + Ifges(6,6) * t86;
t277 = t22 / 0.2e1;
t23 = Ifges(6,1) * t57 + Ifges(6,4) * t56 + Ifges(6,5) * t86;
t276 = t23 / 0.2e1;
t42 = Ifges(5,1) * t87 - Ifges(5,4) * t86 + Ifges(5,5) * t113;
t275 = t42 / 0.2e1;
t158 = t220 * t223 + t227 * t246;
t120 = -t158 * t216 - t219 * t245;
t121 = t158 * t219 - t216 * t245;
t157 = -t227 * t220 + t223 * t246;
t67 = Ifges(6,4) * t121 + Ifges(6,2) * t120 + Ifges(6,6) * t157;
t274 = t67 / 0.2e1;
t68 = Ifges(6,1) * t121 + Ifges(6,4) * t120 + Ifges(6,5) * t157;
t273 = t68 / 0.2e1;
t222 = sin(qJ(6));
t226 = cos(qJ(6));
t171 = t216 * t226 + t219 * t222;
t150 = t171 * t223;
t170 = -t216 * t222 + t219 * t226;
t151 = t170 * t223;
t97 = Ifges(7,4) * t151 - Ifges(7,2) * t150 - Ifges(7,6) * t227;
t272 = t97 / 0.2e1;
t98 = Ifges(7,1) * t151 - Ifges(7,4) * t150 - Ifges(7,5) * t227;
t271 = t98 / 0.2e1;
t101 = Ifges(5,1) * t158 - Ifges(5,4) * t157 - Ifges(5,5) * t245;
t270 = t101 / 0.2e1;
t117 = Ifges(7,4) * t171 + Ifges(7,2) * t170;
t269 = t117 / 0.2e1;
t118 = Ifges(7,1) * t171 + Ifges(7,4) * t170;
t268 = t118 / 0.2e1;
t250 = Ifges(6,4) * t219;
t145 = -Ifges(6,6) * t227 + (-Ifges(6,2) * t216 + t250) * t223;
t267 = t145 / 0.2e1;
t251 = Ifges(6,4) * t216;
t146 = -Ifges(6,5) * t227 + (Ifges(6,1) * t219 - t251) * t223;
t266 = t146 / 0.2e1;
t265 = -t150 / 0.2e1;
t264 = t151 / 0.2e1;
t263 = t170 / 0.2e1;
t262 = t171 / 0.2e1;
t183 = Ifges(6,2) * t219 + t251;
t261 = t183 / 0.2e1;
t184 = Ifges(6,1) * t216 + t250;
t260 = t184 / 0.2e1;
t188 = Ifges(5,1) * t223 + Ifges(5,4) * t227;
t259 = t188 / 0.2e1;
t258 = -t216 / 0.2e1;
t255 = pkin(2) * t228;
t254 = pkin(11) * t223;
t253 = pkin(11) * t227;
t252 = pkin(12) + qJ(5);
t78 = -t119 * t217 + t220 * t134;
t45 = pkin(3) * t113 - pkin(11) * t114 + t78;
t241 = t220 * t224;
t54 = t228 * t108 + t119 * t241 + t134 * t246;
t49 = pkin(11) * t156 + t54;
t20 = t223 * t45 + t227 * t49;
t15 = qJ(5) * t113 + t20;
t48 = -pkin(3) * t156 - t53;
t26 = pkin(4) * t86 - qJ(5) * t87 + t48;
t6 = t219 * t15 + t216 * t26;
t196 = pkin(10) * t246;
t147 = t196 + (-pkin(3) - t255) * t220;
t85 = pkin(4) * t157 - qJ(5) * t158 + t147;
t163 = pkin(2) * t241 + pkin(10) * t245;
t148 = pkin(11) * t220 + t163;
t149 = (-pkin(3) * t228 - pkin(11) * t224 - pkin(2)) * t217;
t95 = t227 * t148 + t223 * t149;
t88 = -qJ(5) * t245 + t95;
t51 = t216 * t85 + t219 * t88;
t162 = -pkin(9) * t244 + t203;
t249 = t162 * mrSges(3,1);
t248 = t164 * mrSges(3,2);
t247 = t216 * t223;
t242 = t219 * t223;
t116 = Ifges(7,5) * t171 + Ifges(7,6) * t170;
t160 = mrSges(6,1) * t247 + mrSges(6,2) * t242;
t178 = -pkin(4) * t227 - qJ(5) * t223 - pkin(3);
t136 = t216 * t178 + t219 * t253;
t186 = Ifges(5,5) * t223 + Ifges(5,6) * t227;
t239 = t216 ^ 2 + t219 ^ 2;
t30 = -t222 * t57 + t226 * t56;
t31 = t222 * t56 + t226 * t57;
t7 = Ifges(7,5) * t31 + Ifges(7,6) * t30 + Ifges(7,3) * t86;
t74 = t120 * t226 - t121 * t222;
t75 = t120 * t222 + t121 * t226;
t34 = Ifges(7,5) * t75 + Ifges(7,6) * t74 + Ifges(7,3) * t157;
t40 = Ifges(5,5) * t87 - Ifges(5,6) * t86 + Ifges(5,3) * t113;
t63 = Ifges(4,5) * t114 - Ifges(4,6) * t113 + Ifges(4,3) * t156;
t140 = Ifges(4,5) * t246 + Ifges(4,6) * t245 + Ifges(4,3) * t220;
t237 = Ifges(3,5) * t244 + Ifges(3,6) * t243 + Ifges(3,3) * t221;
t236 = t116 / 0.2e1 + Ifges(6,5) * t216 / 0.2e1 + Ifges(6,6) * t257;
t32 = -t56 * mrSges(6,1) + t57 * mrSges(6,2);
t10 = -t30 * mrSges(7,1) + t31 * mrSges(7,2);
t43 = -t74 * mrSges(7,1) + t75 * mrSges(7,2);
t5 = -t15 * t216 + t219 * t26;
t50 = -t216 * t88 + t219 * t85;
t19 = -t223 * t49 + t227 * t45;
t77 = -t120 * mrSges(6,1) + t121 * mrSges(6,2);
t180 = -t219 * mrSges(6,1) + t216 * mrSges(6,2);
t102 = t150 * mrSges(7,1) + t151 * mrSges(7,2);
t115 = -t170 * mrSges(7,1) + t171 * mrSges(7,2);
t94 = -t223 * t148 + t149 * t227;
t21 = Ifges(6,5) * t57 + Ifges(6,6) * t56 + Ifges(6,3) * t86;
t41 = Ifges(5,4) * t87 - Ifges(5,2) * t86 + Ifges(5,6) * t113;
t235 = t7 / 0.2e1 + t21 / 0.2e1 - t41 / 0.2e1;
t100 = Ifges(5,4) * t158 - Ifges(5,2) * t157 - Ifges(5,6) * t245;
t66 = Ifges(6,5) * t121 + Ifges(6,6) * t120 + Ifges(6,3) * t157;
t234 = t34 / 0.2e1 + t66 / 0.2e1 - t100 / 0.2e1;
t96 = Ifges(7,5) * t151 - Ifges(7,6) * t150 - Ifges(7,3) * t227;
t144 = -Ifges(6,3) * t227 + (Ifges(6,5) * t219 - Ifges(6,6) * t216) * t223;
t187 = Ifges(5,4) * t223 + Ifges(5,2) * t227;
t233 = t96 / 0.2e1 + t144 / 0.2e1 - t187 / 0.2e1;
t89 = pkin(4) * t245 - t94;
t99 = Ifges(5,5) * t158 - Ifges(5,6) * t157 - Ifges(5,3) * t245;
t16 = -pkin(4) * t113 - t19;
t231 = pkin(11) ^ 2;
t215 = t227 ^ 2;
t214 = t223 ^ 2;
t211 = t214 * t231;
t205 = -pkin(5) * t219 - pkin(4);
t185 = -mrSges(5,1) * t227 + mrSges(5,2) * t223;
t181 = t252 * t219;
t179 = t252 * t216;
t176 = (pkin(5) * t216 + pkin(11)) * t223;
t175 = -mrSges(6,1) * t227 - mrSges(6,3) * t242;
t174 = mrSges(6,2) * t227 - mrSges(6,3) * t247;
t173 = -mrSges(4,2) * t220 + mrSges(4,3) * t245;
t172 = mrSges(4,1) * t220 - mrSges(4,3) * t246;
t169 = t219 * t178;
t161 = t220 * t255 - t196;
t159 = (-mrSges(4,1) * t228 + mrSges(4,2) * t224) * t217;
t142 = Ifges(4,5) * t220 + (Ifges(4,1) * t224 + Ifges(4,4) * t228) * t217;
t141 = Ifges(4,6) * t220 + (Ifges(4,4) * t224 + Ifges(4,2) * t228) * t217;
t135 = -t216 * t253 + t169;
t130 = -mrSges(7,1) * t227 - mrSges(7,3) * t151;
t129 = mrSges(7,2) * t227 - mrSges(7,3) * t150;
t128 = -mrSges(5,1) * t245 - mrSges(5,3) * t158;
t127 = mrSges(5,2) * t245 - mrSges(5,3) * t157;
t125 = -t179 * t222 + t181 * t226;
t124 = -t179 * t226 - t181 * t222;
t122 = -pkin(12) * t247 + t136;
t107 = -pkin(12) * t242 + t169 + (-pkin(11) * t216 - pkin(5)) * t227;
t103 = mrSges(5,1) * t157 + mrSges(5,2) * t158;
t93 = mrSges(6,1) * t157 - mrSges(6,3) * t121;
t92 = -mrSges(6,2) * t157 + mrSges(6,3) * t120;
t91 = mrSges(4,1) * t156 - mrSges(4,3) * t114;
t90 = -mrSges(4,2) * t156 - mrSges(4,3) * t113;
t76 = mrSges(4,1) * t113 + mrSges(4,2) * t114;
t70 = t107 * t222 + t122 * t226;
t69 = t107 * t226 - t122 * t222;
t65 = Ifges(4,1) * t114 - Ifges(4,4) * t113 + Ifges(4,5) * t156;
t64 = Ifges(4,4) * t114 - Ifges(4,2) * t113 + Ifges(4,6) * t156;
t62 = -pkin(5) * t120 + t89;
t61 = mrSges(7,1) * t157 - mrSges(7,3) * t75;
t60 = -mrSges(7,2) * t157 + mrSges(7,3) * t74;
t59 = mrSges(5,1) * t113 - mrSges(5,3) * t87;
t58 = -mrSges(5,2) * t113 - mrSges(5,3) * t86;
t52 = mrSges(5,1) * t86 + mrSges(5,2) * t87;
t39 = pkin(12) * t120 + t51;
t38 = mrSges(6,1) * t86 - mrSges(6,3) * t57;
t37 = -mrSges(6,2) * t86 + mrSges(6,3) * t56;
t36 = Ifges(7,1) * t75 + Ifges(7,4) * t74 + Ifges(7,5) * t157;
t35 = Ifges(7,4) * t75 + Ifges(7,2) * t74 + Ifges(7,6) * t157;
t33 = pkin(5) * t157 - pkin(12) * t121 + t50;
t18 = mrSges(7,1) * t86 - mrSges(7,3) * t31;
t17 = -mrSges(7,2) * t86 + mrSges(7,3) * t30;
t13 = t222 * t33 + t226 * t39;
t12 = -t222 * t39 + t226 * t33;
t11 = -pkin(5) * t56 + t16;
t9 = Ifges(7,1) * t31 + Ifges(7,4) * t30 + Ifges(7,5) * t86;
t8 = Ifges(7,4) * t31 + Ifges(7,2) * t30 + Ifges(7,6) * t86;
t4 = pkin(12) * t56 + t6;
t3 = pkin(5) * t86 - pkin(12) * t57 + t5;
t2 = t222 * t3 + t226 * t4;
t1 = -t222 * t4 + t226 * t3;
t14 = [m(4) * (t53 ^ 2 + t54 ^ 2 + t78 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t48 ^ 2) + m(6) * (t16 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + (t7 + t21 - t41) * t86 + m(3) * (t162 ^ 2 + t164 ^ 2) + (t40 - t64) * t113 + (t237 - 0.2e1 * t248 + 0.2e1 * t249) * t221 + Ifges(2,3) + 0.2e1 * t11 * t10 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + t30 * t8 + t31 * t9 + 0.2e1 * t16 * t32 + 0.2e1 * t6 * t37 + 0.2e1 * t5 * t38 + 0.2e1 * t48 * t52 + ((Ifges(3,5) * t225 + Ifges(3,6) * t229) * t221 + 0.2e1 * (-t162 * t225 + t164 * t229) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t229 + mrSges(3,2) * t225) + t229 * (Ifges(3,4) * t225 + Ifges(3,2) * t229) + t225 * (Ifges(3,1) * t225 + Ifges(3,4) * t229) + m(3) * pkin(1) ^ 2) * t218) * t218 + t56 * t22 + t57 * t23 + 0.2e1 * t20 * t58 + 0.2e1 * t19 * t59 + 0.2e1 * t78 * t76 + t87 * t42 + 0.2e1 * t54 * t90 + 0.2e1 * t53 * t91 + t114 * t65 + t156 * t63; (-pkin(2) * t76 + t224 * t65 / 0.2e1 + (t64 / 0.2e1 - t40 / 0.2e1) * t228) * t217 + m(5) * (t147 * t48 + t19 * t94 + t20 * t95) + m(6) * (t16 * t89 + t5 * t50 + t51 * t6) + m(7) * (t1 * t12 + t11 * t62 + t13 * t2) + m(4) * (-pkin(2) * t217 * t78 + t161 * t53 + t163 * t54) + (t99 / 0.2e1 - t141 / 0.2e1) * t113 + t158 * t275 + t121 * t276 + t120 * t277 + t87 * t270 + t57 * t273 + t56 * t274 + t237 - t248 + t249 + t234 * t86 + t235 * t157 + t13 * t17 + t12 * t18 + t30 * t35 / 0.2e1 + t31 * t36 / 0.2e1 + t11 * t43 + t50 * t38 + t51 * t37 + t2 * t60 + t1 * t61 + t62 * t10 + t74 * t8 / 0.2e1 + t75 * t9 / 0.2e1 + t16 * t77 + t89 * t32 + t6 * t92 + t5 * t93 + t94 * t59 + t95 * t58 + t48 * t103 + t20 * t127 + t19 * t128 + t114 * t142 / 0.2e1 + t147 * t52 + t156 * t140 / 0.2e1 + t78 * t159 + t161 * t91 + t163 * t90 + t53 * t172 + t54 * t173 + t220 * t63 / 0.2e1; Ifges(3,3) + 0.2e1 * t13 * t60 + 0.2e1 * t12 * t61 + 0.2e1 * t62 * t43 + t74 * t35 + t75 * t36 + 0.2e1 * t89 * t77 + 0.2e1 * t51 * t92 + 0.2e1 * t50 * t93 + t120 * t67 + t121 * t68 + 0.2e1 * t95 * t127 + 0.2e1 * t94 * t128 + 0.2e1 * t147 * t103 + t158 * t101 + 0.2e1 * t161 * t172 + 0.2e1 * t163 * t173 + t220 * t140 + (t34 + t66 - t100) * t157 + (-0.2e1 * pkin(2) * t159 + t224 * t142 + (t141 - t99) * t228) * t217 + m(5) * (t147 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2 + t89 ^ 2) + m(7) * (t12 ^ 2 + t13 ^ 2 + t62 ^ 2) + m(4) * (pkin(2) ^ 2 * t217 ^ 2 + t161 ^ 2 + t163 ^ 2); m(7) * (t1 * t69 + t11 * t176 + t2 * t70) + m(6) * (t135 * t5 + t136 * t6 + t16 * t254) + t63 + (-t19 * mrSges(5,3) + t22 * t258 + t23 * t257 + t275 + (-t59 + t32) * pkin(11)) * t223 + t87 * t259 + t9 * t264 + t8 * t265 + t57 * t266 + t56 * t267 + t31 * t271 + t30 * t272 + m(5) * (-pkin(3) * t48 + (-t19 * t223 + t20 * t227) * pkin(11)) + t233 * t86 + (t20 * mrSges(5,3) + pkin(11) * t58 - t235) * t227 - pkin(3) * t52 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + t69 * t18 + t70 * t17 + t11 * t102 + t2 * t129 + t1 * t130 + t135 * t38 + t136 * t37 + t16 * t160 + t6 * t174 + t5 * t175 + t176 * t10 + t48 * t185 + t113 * t186 / 0.2e1; -t186 * t245 / 0.2e1 + m(7) * (t12 * t69 + t13 * t70 + t176 * t62) + m(6) * (t135 * t50 + t136 * t51 + t254 * t89) + t158 * t259 + t36 * t264 + t35 * t265 + t121 * t266 + t120 * t267 + t75 * t271 + t74 * t272 + (-t94 * mrSges(5,3) + t67 * t258 + t68 * t257 + t270 + (-t128 + t77) * pkin(11)) * t223 + m(5) * (-pkin(3) * t147 + (-t223 * t94 + t227 * t95) * pkin(11)) + t233 * t157 + (t95 * mrSges(5,3) + pkin(11) * t127 - t234) * t227 + t140 + t69 * t61 + t70 * t60 + t62 * t102 - pkin(3) * t103 + t13 * t129 + t12 * t130 + t135 * t93 + t136 * t92 + t89 * t160 + t161 * mrSges(4,1) - t163 * mrSges(4,2) + t51 * t174 + t50 * t175 + t176 * t43 + t147 * t185; -0.2e1 * pkin(3) * t185 + 0.2e1 * t176 * t102 + 0.2e1 * t70 * t129 + 0.2e1 * t69 * t130 + 0.2e1 * t135 * t175 + 0.2e1 * t136 * t174 - t150 * t97 + t151 * t98 + Ifges(4,3) + (t214 + t215) * mrSges(5,3) * t278 + (-t144 + t187 - t96) * t227 + m(5) * (pkin(3) ^ 2 + t215 * t231 + t211) + m(7) * (t176 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t135 ^ 2 + t136 ^ 2 + t211) + (-t145 * t216 + t146 * t219 + t160 * t278 + t188) * t223; (t6 * mrSges(6,3) + qJ(5) * t37 + t277) * t219 + (-t5 * mrSges(6,3) - qJ(5) * t38 + t276) * t216 + m(7) * (t1 * t124 + t11 * t205 + t125 * t2) + t40 + (-t1 * t171 + t170 * t2) * mrSges(7,3) + m(6) * (-pkin(4) * t16 + (-t216 * t5 + t219 * t6) * qJ(5)) + t236 * t86 + t19 * mrSges(5,1) - t20 * mrSges(5,2) - pkin(4) * t32 + t11 * t115 + t30 * t269 + t31 * t268 + t124 * t18 + t125 * t17 + t8 * t263 + t9 * t262 + t16 * t180 + t56 * t261 + t57 * t260 + t205 * t10; t236 * t157 + t99 + (-t12 * t171 + t13 * t170) * mrSges(7,3) + m(6) * (-pkin(4) * t89 + (-t216 * t50 + t219 * t51) * qJ(5)) + m(7) * (t12 * t124 + t125 * t13 + t205 * t62) + (-t50 * mrSges(6,3) - qJ(5) * t93 + t273) * t216 + (t51 * mrSges(6,3) + qJ(5) * t92 + t274) * t219 - pkin(4) * t77 + t94 * mrSges(5,1) - t95 * mrSges(5,2) + t62 * t115 + t74 * t269 + t75 * t268 + t124 * t61 + t125 * t60 + t35 * t263 + t36 * t262 + t89 * t180 + t120 * t261 + t121 * t260 + t205 * t43; t125 * t129 + t124 * t130 + t117 * t265 + t118 * t264 - pkin(4) * t160 + t97 * t263 + t98 * t262 + t176 * t115 + t205 * t102 + (-mrSges(5,1) + t180) * t254 + (t170 * t70 - t171 * t69) * mrSges(7,3) + (t136 * mrSges(6,3) + qJ(5) * t174 + t223 * t260 + t267) * t219 + (-t135 * mrSges(6,3) - qJ(5) * t175 - t223 * t183 / 0.2e1 + t266) * t216 + m(7) * (t124 * t69 + t125 * t70 + t176 * t205) + m(6) * (-pkin(4) * t254 + (-t135 * t216 + t136 * t219) * qJ(5)) + (-pkin(11) * mrSges(5,2) - t236) * t227 + t186; -0.2e1 * pkin(4) * t180 + 0.2e1 * t205 * t115 + t170 * t117 + t171 * t118 + t219 * t183 + t216 * t184 + Ifges(5,3) + m(7) * (t124 ^ 2 + t125 ^ 2 + t205 ^ 2) + m(6) * (qJ(5) ^ 2 * t239 + pkin(4) ^ 2) + 0.2e1 * (-t124 * t171 + t125 * t170) * mrSges(7,3) + 0.2e1 * t239 * qJ(5) * mrSges(6,3); m(6) * t16 + m(7) * t11 + t10 + t32; m(6) * t89 + m(7) * t62 + t43 + t77; m(6) * t254 + m(7) * t176 + t102 + t160; -m(6) * pkin(4) + m(7) * t205 + t115 + t180; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t12 - mrSges(7,2) * t13 + t34; mrSges(7,1) * t69 - mrSges(7,2) * t70 + t96; mrSges(7,1) * t124 - mrSges(7,2) * t125 + t116; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
