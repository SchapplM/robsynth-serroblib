% Calculate joint inertia matrix for
% S6RRRRPR12
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
% Datum: 2018-11-23 18:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:21:34
% EndTime: 2018-11-23 18:21:37
% DurationCPUTime: 2.80s
% Computational Cost: add. (8443->539), mult. (21649->787), div. (0->0), fcn. (24920->14), ass. (0->202)
t207 = cos(pkin(6));
t210 = sin(qJ(3));
t214 = cos(qJ(3));
t204 = sin(pkin(6));
t206 = cos(pkin(7));
t215 = cos(qJ(2));
t234 = t206 * t215;
t227 = t204 * t234;
t211 = sin(qJ(2));
t237 = t204 * t211;
t203 = sin(pkin(7));
t238 = t203 * t214;
t120 = -t207 * t238 + t210 * t237 - t214 * t227;
t202 = sin(pkin(13));
t205 = cos(pkin(13));
t239 = t203 * t210;
t121 = t207 * t239 + (t210 * t234 + t211 * t214) * t204;
t236 = t204 * t215;
t149 = -t203 * t236 + t206 * t207;
t209 = sin(qJ(4));
t213 = cos(qJ(4));
t82 = -t121 * t209 + t149 * t213;
t83 = t121 * t213 + t149 * t209;
t54 = t202 * t83 - t205 * t82;
t55 = t202 * t82 + t205 * t83;
t20 = Ifges(6,5) * t55 - Ifges(6,6) * t54 + Ifges(6,3) * t120;
t38 = Ifges(5,5) * t83 + Ifges(5,6) * t82 + Ifges(5,3) * t120;
t271 = t20 + t38;
t251 = pkin(1) * t207;
t156 = pkin(9) * t236 + t211 * t251;
t115 = (t203 * t207 + t227) * pkin(10) + t156;
t186 = t215 * t251;
t126 = pkin(2) * t207 + t186 + (-pkin(10) * t206 - pkin(9)) * t237;
t137 = (-pkin(10) * t203 * t211 - pkin(2) * t215 - pkin(1)) * t204;
t58 = -t210 * t115 + t214 * (t126 * t206 + t137 * t203);
t248 = -qJ(5) - pkin(11);
t169 = t248 * t213;
t223 = t248 * t209;
t128 = -t169 * t202 - t205 * t223;
t270 = t128 ^ 2;
t269 = 0.2e1 * t128;
t22 = Ifges(6,1) * t55 - Ifges(6,4) * t54 + Ifges(6,5) * t120;
t268 = t22 / 0.2e1;
t39 = Ifges(5,4) * t83 + Ifges(5,2) * t82 + Ifges(5,6) * t120;
t267 = t39 / 0.2e1;
t40 = Ifges(5,1) * t83 + Ifges(5,4) * t82 + Ifges(5,5) * t120;
t266 = t40 / 0.2e1;
t150 = t206 * t213 - t209 * t239;
t151 = t206 * t209 + t213 * t239;
t104 = -t205 * t150 + t151 * t202;
t105 = t150 * t202 + t151 * t205;
t64 = Ifges(6,1) * t105 - Ifges(6,4) * t104 - Ifges(6,5) * t238;
t265 = t64 / 0.2e1;
t161 = t202 * t209 - t205 * t213;
t162 = t202 * t213 + t205 * t209;
t208 = sin(qJ(6));
t212 = cos(qJ(6));
t246 = Ifges(7,4) * t212;
t91 = Ifges(7,6) * t161 + (-Ifges(7,2) * t208 + t246) * t162;
t264 = t91 / 0.2e1;
t247 = Ifges(7,4) * t208;
t92 = Ifges(7,5) * t161 + (Ifges(7,1) * t212 - t247) * t162;
t263 = t92 / 0.2e1;
t98 = Ifges(5,4) * t151 + Ifges(5,2) * t150 - Ifges(5,6) * t238;
t262 = t98 / 0.2e1;
t99 = Ifges(5,1) * t151 + Ifges(5,4) * t150 - Ifges(5,5) * t238;
t261 = t99 / 0.2e1;
t125 = Ifges(6,1) * t162 - Ifges(6,4) * t161;
t260 = t125 / 0.2e1;
t170 = Ifges(7,5) * t208 + Ifges(7,6) * t212;
t259 = t170 / 0.2e1;
t172 = Ifges(7,2) * t212 + t247;
t258 = t172 / 0.2e1;
t173 = Ifges(5,4) * t209 + Ifges(5,2) * t213;
t257 = t173 / 0.2e1;
t174 = Ifges(7,1) * t208 + t246;
t256 = t174 / 0.2e1;
t175 = Ifges(5,1) * t209 + Ifges(5,4) * t213;
t255 = t175 / 0.2e1;
t254 = -t208 / 0.2e1;
t253 = t208 / 0.2e1;
t252 = t212 / 0.2e1;
t250 = pkin(2) * t214;
t249 = -Ifges(5,3) - Ifges(6,3);
t78 = -t126 * t203 + t206 * t137;
t45 = pkin(3) * t120 - pkin(11) * t121 + t78;
t235 = t206 * t210;
t59 = t214 * t115 + t126 * t235 + t137 * t239;
t49 = pkin(11) * t149 + t59;
t23 = -t209 * t49 + t213 * t45;
t12 = pkin(4) * t120 - qJ(5) * t83 + t23;
t24 = t209 * t45 + t213 * t49;
t15 = qJ(5) * t82 + t24;
t6 = t202 * t12 + t205 * t15;
t155 = pkin(2) * t235 + pkin(10) * t238;
t142 = pkin(11) * t206 + t155;
t143 = (-pkin(3) * t214 - pkin(11) * t210 - pkin(2)) * t203;
t95 = -t142 * t209 + t213 * t143;
t72 = -pkin(4) * t238 - qJ(5) * t151 + t95;
t96 = t213 * t142 + t209 * t143;
t77 = qJ(5) * t150 + t96;
t42 = t202 * t72 + t205 * t77;
t154 = -pkin(9) * t237 + t186;
t245 = t154 * mrSges(3,1);
t244 = t156 * mrSges(3,2);
t243 = t162 * t208;
t242 = t162 * t212;
t189 = pkin(4) * t202 + pkin(12);
t241 = t189 * t208;
t240 = t189 * t212;
t233 = Ifges(6,5) * t105 - Ifges(6,6) * t104;
t232 = Ifges(5,5) * t151 + Ifges(5,6) * t150;
t123 = Ifges(6,5) * t162 - Ifges(6,6) * t161;
t171 = Ifges(5,5) * t209 + Ifges(5,6) * t213;
t231 = t208 ^ 2 + t212 ^ 2;
t230 = t209 ^ 2 + t213 ^ 2;
t29 = t120 * t212 - t208 * t55;
t30 = t120 * t208 + t212 * t55;
t7 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t54;
t21 = Ifges(6,4) * t55 - Ifges(6,2) * t54 + Ifges(6,6) * t120;
t229 = t7 / 0.2e1 - t21 / 0.2e1;
t86 = -t105 * t208 - t212 * t238;
t87 = t105 * t212 - t208 * t238;
t33 = Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t104;
t63 = Ifges(6,4) * t105 - Ifges(6,2) * t104 - Ifges(6,6) * t238;
t228 = t33 / 0.2e1 - t63 / 0.2e1;
t124 = Ifges(6,4) * t162 - Ifges(6,2) * t161;
t90 = Ifges(7,5) * t242 - Ifges(7,6) * t243 + Ifges(7,3) * t161;
t226 = t90 / 0.2e1 - t124 / 0.2e1;
t67 = Ifges(4,5) * t121 - Ifges(4,6) * t120 + Ifges(4,3) * t149;
t138 = Ifges(4,5) * t239 + Ifges(4,6) * t238 + Ifges(4,3) * t206;
t225 = Ifges(3,5) * t237 + Ifges(3,6) * t236 + Ifges(3,3) * t207;
t191 = -pkin(4) * t213 - pkin(3);
t224 = t123 / 0.2e1 + t171 / 0.2e1;
t25 = t54 * mrSges(6,1) + t55 * mrSges(6,2);
t70 = t104 * mrSges(6,1) + t105 * mrSges(6,2);
t122 = t161 * mrSges(6,1) + t162 * mrSges(6,2);
t48 = -pkin(3) * t149 - t58;
t26 = -pkin(4) * t82 + t48;
t10 = pkin(5) * t54 - pkin(12) * t55 + t26;
t4 = pkin(12) * t120 + t6;
t1 = t10 * t212 - t208 * t4;
t2 = t10 * t208 + t212 * t4;
t222 = -t1 * t208 + t2 * t212;
t221 = mrSges(7,1) * t208 + mrSges(7,2) * t212;
t5 = t12 * t205 - t15 * t202;
t37 = -pkin(12) * t238 + t42;
t181 = pkin(10) * t239;
t141 = t181 + (-pkin(3) - t250) * t206;
t106 = -pkin(4) * t150 + t141;
t47 = pkin(5) * t104 - pkin(12) * t105 + t106;
t16 = -t208 * t37 + t212 * t47;
t17 = t208 * t47 + t212 * t37;
t220 = -t16 * t208 + t17 * t212;
t41 = -t202 * t77 + t205 * t72;
t112 = pkin(5) * t161 - pkin(12) * t162 + t191;
t130 = -t205 * t169 + t202 * t223;
t73 = t112 * t212 - t130 * t208;
t74 = t112 * t208 + t130 * t212;
t219 = -t208 * t73 + t212 * t74;
t190 = -pkin(4) * t205 - pkin(5);
t168 = -mrSges(5,1) * t213 + mrSges(5,2) * t209;
t167 = -mrSges(7,1) * t212 + mrSges(7,2) * t208;
t165 = -mrSges(4,2) * t206 + mrSges(4,3) * t238;
t164 = mrSges(4,1) * t206 - mrSges(4,3) * t239;
t153 = t206 * t250 - t181;
t152 = (-mrSges(4,1) * t214 + mrSges(4,2) * t210) * t203;
t140 = Ifges(4,5) * t206 + (Ifges(4,1) * t210 + Ifges(4,4) * t214) * t203;
t139 = Ifges(4,6) * t206 + (Ifges(4,4) * t210 + Ifges(4,2) * t214) * t203;
t133 = -mrSges(5,1) * t238 - mrSges(5,3) * t151;
t132 = mrSges(5,2) * t238 + mrSges(5,3) * t150;
t114 = mrSges(7,1) * t161 - mrSges(7,3) * t242;
t113 = -mrSges(7,2) * t161 - mrSges(7,3) * t243;
t108 = -mrSges(5,1) * t150 + mrSges(5,2) * t151;
t107 = t221 * t162;
t97 = -Ifges(5,3) * t238 + t232;
t94 = -mrSges(6,1) * t238 - mrSges(6,3) * t105;
t93 = mrSges(6,2) * t238 - mrSges(6,3) * t104;
t89 = mrSges(4,1) * t149 - mrSges(4,3) * t121;
t88 = -mrSges(4,2) * t149 - mrSges(4,3) * t120;
t76 = mrSges(4,1) * t120 + mrSges(4,2) * t121;
t69 = Ifges(4,1) * t121 - Ifges(4,4) * t120 + Ifges(4,5) * t149;
t68 = Ifges(4,4) * t121 - Ifges(4,2) * t120 + Ifges(4,6) * t149;
t66 = mrSges(5,1) * t120 - mrSges(5,3) * t83;
t65 = -mrSges(5,2) * t120 + mrSges(5,3) * t82;
t62 = -Ifges(6,3) * t238 + t233;
t61 = mrSges(7,1) * t104 - mrSges(7,3) * t87;
t60 = -mrSges(7,2) * t104 + mrSges(7,3) * t86;
t57 = -mrSges(7,1) * t86 + mrSges(7,2) * t87;
t56 = -mrSges(5,1) * t82 + mrSges(5,2) * t83;
t36 = pkin(5) * t238 - t41;
t35 = Ifges(7,1) * t87 + Ifges(7,4) * t86 + Ifges(7,5) * t104;
t34 = Ifges(7,4) * t87 + Ifges(7,2) * t86 + Ifges(7,6) * t104;
t32 = mrSges(6,1) * t120 - mrSges(6,3) * t55;
t31 = -mrSges(6,2) * t120 - mrSges(6,3) * t54;
t19 = mrSges(7,1) * t54 - mrSges(7,3) * t30;
t18 = -mrSges(7,2) * t54 + mrSges(7,3) * t29;
t13 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t9 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t54;
t8 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t54;
t3 = -pkin(5) * t120 - t5;
t11 = [m(4) * (t58 ^ 2 + t59 ^ 2 + t78 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2 + t48 ^ 2) + m(6) * (t26 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(3) * (t154 ^ 2 + t156 ^ 2) + ((Ifges(3,5) * t211 + Ifges(3,6) * t215) * t207 + 0.2e1 * (-t154 * t211 + t156 * t215) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t215 + mrSges(3,2) * t211) + t211 * (Ifges(3,1) * t211 + Ifges(3,4) * t215) + t215 * (Ifges(3,4) * t211 + Ifges(3,2) * t215)) * t204) * t204 + (t7 - t21) * t54 + (t225 - 0.2e1 * t244 + 0.2e1 * t245) * t207 + Ifges(2,3) + (-t68 + t271) * t120 + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t26 * t25 + t29 * t8 + t30 * t9 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + t55 * t22 + 0.2e1 * t48 * t56 + 0.2e1 * t24 * t65 + 0.2e1 * t23 * t66 + 0.2e1 * t78 * t76 + t82 * t39 + t83 * t40 + 0.2e1 * t59 * t88 + 0.2e1 * t58 * t89 + t121 * t69 + t149 * t67; m(7) * (t1 * t16 + t17 * t2 + t3 * t36) + m(5) * (t141 * t48 + t23 * t95 + t24 * t96) + m(4) * (-pkin(2) * t203 * t78 + t153 * t58 + t155 * t59) + m(6) * (t106 * t26 + t41 * t5 + t42 * t6) + (t62 / 0.2e1 + t97 / 0.2e1 - t139 / 0.2e1) * t120 - t244 + t245 + (-pkin(2) * t76 + t210 * t69 / 0.2e1 + (t68 / 0.2e1 - t38 / 0.2e1 - t20 / 0.2e1) * t214) * t203 + t225 + t83 * t261 + t82 * t262 + t55 * t265 + t151 * t266 + t150 * t267 + t105 * t268 + t229 * t104 + t17 * t18 + t16 * t19 + t228 * t54 + t29 * t34 / 0.2e1 + t30 * t35 / 0.2e1 + t36 * t13 + t41 * t32 + t42 * t31 + t3 * t57 + t2 * t60 + t1 * t61 + t26 * t70 + t86 * t8 / 0.2e1 + t87 * t9 / 0.2e1 + t6 * t93 + t5 * t94 + t95 * t66 + t96 * t65 + t106 * t25 + t48 * t108 + t24 * t132 + t23 * t133 + t121 * t140 / 0.2e1 + t141 * t56 + t149 * t138 / 0.2e1 + t78 * t152 + t153 * t89 + t155 * t88 + t58 * t164 + t59 * t165 + t206 * t67 / 0.2e1; Ifges(3,3) + 0.2e1 * t36 * t57 + 0.2e1 * t17 * t60 + 0.2e1 * t16 * t61 + t86 * t34 + t87 * t35 + 0.2e1 * t42 * t93 + 0.2e1 * t41 * t94 + t105 * t64 + 0.2e1 * t106 * t70 + 0.2e1 * t96 * t132 + 0.2e1 * t95 * t133 + 0.2e1 * t141 * t108 + t150 * t98 + t151 * t99 + 0.2e1 * t153 * t164 + 0.2e1 * t155 * t165 + t206 * t138 + (t33 - t63) * t104 + (-0.2e1 * pkin(2) * t152 + t210 * t140 + (t139 - t62 - t97) * t214) * t203 + m(6) * (t106 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(7) * (t16 ^ 2 + t17 ^ 2 + t36 ^ 2) + m(5) * (t141 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(4) * (pkin(2) ^ 2 * t203 ^ 2 + t153 ^ 2 + t155 ^ 2); t224 * t120 + (t13 - t32) * t128 + m(5) * (-pkin(3) * t48 + (-t209 * t23 + t213 * t24) * pkin(11)) + m(7) * (t1 * t73 + t128 * t3 + t2 * t74) + m(6) * (-t128 * t5 + t130 * t6 + t191 * t26) + t67 + (-t5 * mrSges(6,3) + t252 * t9 + t254 * t8 + t268) * t162 + (-t23 * mrSges(5,3) - pkin(11) * t66 + t266) * t209 + (t24 * mrSges(5,3) + pkin(11) * t65 + t267) * t213 + t82 * t257 + t55 * t260 + t30 * t263 + t29 * t264 + t83 * t255 + (-t6 * mrSges(6,3) + t229) * t161 + t226 * t54 - pkin(3) * t56 + t58 * mrSges(4,1) - t59 * mrSges(4,2) + t73 * t19 + t74 * t18 + t3 * t107 + t2 * t113 + t1 * t114 + t26 * t122 + t130 * t31 + t48 * t168 + t191 * t25; m(5) * (-pkin(3) * t141 + (-t209 * t95 + t213 * t96) * pkin(11)) + (t57 - t94) * t128 + m(7) * (t128 * t36 + t16 * t73 + t17 * t74) + m(6) * (t106 * t191 - t128 * t41 + t130 * t42) + (-t41 * mrSges(6,3) + t252 * t35 + t254 * t34 + t265) * t162 + (-t95 * mrSges(5,3) - pkin(11) * t133 + t261) * t209 + (t96 * mrSges(5,3) + pkin(11) * t132 + t262) * t213 + t138 + t150 * t257 + t105 * t260 + t87 * t263 + t86 * t264 + t151 * t255 - t224 * t238 + t226 * t104 + (-t42 * mrSges(6,3) + t228) * t161 + t73 * t61 + t74 * t60 + t36 * t107 - pkin(3) * t108 + t17 * t113 + t16 * t114 + t106 * t122 + t130 * t93 + t153 * mrSges(4,1) - t155 * mrSges(4,2) + t141 * t168 + t191 * t70; -0.2e1 * pkin(3) * t168 + t107 * t269 + 0.2e1 * t74 * t113 + 0.2e1 * t73 * t114 + 0.2e1 * t191 * t122 + t213 * t173 + t209 * t175 + Ifges(4,3) + 0.2e1 * t230 * pkin(11) * mrSges(5,3) + (-0.2e1 * mrSges(6,3) * t130 - t124 + t90) * t161 + m(7) * (t73 ^ 2 + t74 ^ 2 + t270) + m(6) * (t130 ^ 2 + t191 ^ 2 + t270) + m(5) * (pkin(11) ^ 2 * t230 + pkin(3) ^ 2) + (mrSges(6,3) * t269 - t208 * t91 + t212 * t92 + t125) * t162; m(7) * (t189 * t222 + t190 * t3) + t222 * mrSges(7,3) + (t202 * t31 + t205 * t32 + m(6) * (t202 * t6 + t205 * t5)) * pkin(4) - t19 * t241 + t18 * t240 + t5 * mrSges(6,1) - t6 * mrSges(6,2) + t23 * mrSges(5,1) - t24 * mrSges(5,2) + t3 * t167 + t54 * t259 + t29 * t258 + t30 * t256 + t190 * t13 + t9 * t253 + t8 * t252 + t271; -t61 * t241 + t60 * t240 + m(7) * (t189 * t220 + t190 * t36) + t41 * mrSges(6,1) - t42 * mrSges(6,2) + t95 * mrSges(5,1) - t96 * mrSges(5,2) + t36 * t167 + t104 * t259 + t86 * t258 + t87 * t256 + t190 * t57 + t35 * t253 + t34 * t252 + t249 * t238 + t220 * mrSges(7,3) + (t202 * t93 + t205 * t94 + m(6) * (t202 * t42 + t205 * t41)) * pkin(4) + t232 + t233; -t114 * t241 + t113 * t240 + t190 * t107 + m(7) * (t128 * t190 + t189 * t219) + t128 * t167 + t161 * t259 + t91 * t252 + t92 * t253 - t128 * mrSges(6,1) - t130 * mrSges(6,2) + (t172 * t254 + t174 * t252) * t162 + (-mrSges(5,1) * t209 - mrSges(5,2) * t213) * pkin(11) + t219 * mrSges(7,3) + (m(6) * (-t128 * t205 + t130 * t202) + (-t161 * t202 - t162 * t205) * mrSges(6,3)) * pkin(4) + t171 + t123; 0.2e1 * t190 * t167 + t212 * t172 + t208 * t174 + m(7) * (t189 ^ 2 * t231 + t190 ^ 2) + m(6) * (t202 ^ 2 + t205 ^ 2) * pkin(4) ^ 2 - t249 + 0.2e1 * (mrSges(6,1) * t205 - mrSges(6,2) * t202) * pkin(4) + 0.2e1 * t231 * t189 * mrSges(7,3); t208 * t18 + t212 * t19 + m(7) * (t1 * t212 + t2 * t208) + m(6) * t26 + t25; t208 * t60 + t212 * t61 + m(7) * (t16 * t212 + t17 * t208) + m(6) * t106 + t70; t208 * t113 + t212 * t114 + m(7) * (t208 * t74 + t212 * t73) + m(6) * t191 + t122; 0; m(7) * t231 + m(6); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t16 - mrSges(7,2) * t17 + t33; mrSges(7,1) * t73 - mrSges(7,2) * t74 + t90; -t189 * t221 + t170; -t167; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
