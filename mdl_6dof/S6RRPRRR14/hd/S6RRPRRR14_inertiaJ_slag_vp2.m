% Calculate joint inertia matrix for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:57:32
% EndTime: 2019-03-09 14:57:40
% DurationCPUTime: 3.30s
% Computational Cost: add. (12853->531), mult. (35560->792), div. (0->0), fcn. (41453->16), ass. (0->210)
t274 = 2 * pkin(12);
t199 = sin(pkin(8));
t203 = cos(pkin(8));
t208 = sin(qJ(4));
t212 = cos(qJ(4));
t200 = sin(pkin(7));
t204 = cos(pkin(7));
t205 = cos(pkin(6));
t201 = sin(pkin(6));
t213 = cos(qJ(2));
t237 = t201 * t213;
t150 = -t200 * t237 + t204 * t205;
t198 = sin(pkin(14));
t202 = cos(pkin(14));
t209 = sin(qJ(2));
t234 = t204 * t213;
t239 = t200 * t205;
t123 = t202 * t239 + (-t198 * t209 + t202 * t234) * t201;
t246 = t123 * t203;
t221 = t150 * t199 + t246;
t259 = pkin(1) * t205;
t158 = pkin(10) * t237 + t209 * t259;
t216 = t201 * t234 + t239;
t119 = qJ(3) * t216 + t158;
t183 = t213 * t259;
t238 = t201 * t209;
t128 = pkin(2) * t205 + t183 + (-qJ(3) * t204 - pkin(10)) * t238;
t137 = (-qJ(3) * t200 * t209 - pkin(2) * t213 - pkin(1)) * t201;
t243 = t198 * t204;
t244 = t198 * t200;
t76 = t202 * t119 + t128 * t243 + t137 * t244;
t53 = pkin(11) * t221 + t76;
t124 = t198 * t216 + t202 * t238;
t257 = pkin(11) * t203;
t236 = t202 * t204;
t240 = t200 * t202;
t75 = -t119 * t198 + t128 * t236 + t137 * t240;
t56 = pkin(3) * t150 - t124 * t257 + t75;
t258 = pkin(11) * t199;
t97 = -t128 * t200 + t204 * t137;
t67 = -pkin(3) * t123 - t124 * t258 + t97;
t19 = -t208 * t53 + t212 * (t199 * t67 + t203 * t56);
t153 = pkin(2) * t243 + qJ(3) * t240;
t228 = t203 * t240;
t117 = (t199 * t204 + t228) * pkin(11) + t153;
t178 = pkin(2) * t236;
t127 = pkin(3) * t204 + t178 + (-qJ(3) - t257) * t244;
t136 = (-pkin(3) * t202 - t198 * t258 - pkin(2)) * t200;
t73 = -t208 * t117 + (t127 * t203 + t136 * t199) * t212;
t206 = sin(qJ(6));
t210 = cos(qJ(6));
t166 = -mrSges(7,1) * t210 + mrSges(7,2) * t206;
t273 = -m(7) * pkin(5) - mrSges(6,1) + t166;
t207 = sin(qJ(5));
t211 = cos(qJ(5));
t242 = t199 * t208;
t154 = -t211 * t203 + t207 * t242;
t272 = t154 ^ 2;
t103 = -t123 * t199 + t150 * t203;
t81 = t124 * t212 + t208 * t221;
t54 = -t103 * t211 + t207 * t81;
t55 = t103 * t207 + t211 * t81;
t241 = t199 * t212;
t80 = t124 * t208 - t150 * t241 - t212 * t246;
t25 = Ifges(6,1) * t55 - Ifges(6,4) * t54 + Ifges(6,5) * t80;
t271 = t25 / 0.2e1;
t235 = t203 * t208;
t126 = t204 * t242 + (t198 * t212 + t202 * t235) * t200;
t149 = -t199 * t240 + t203 * t204;
t104 = t126 * t207 - t211 * t149;
t105 = t126 * t211 + t149 * t207;
t125 = -t204 * t241 + t208 * t244 - t212 * t228;
t63 = Ifges(6,1) * t105 - Ifges(6,4) * t104 + Ifges(6,5) * t125;
t270 = t63 / 0.2e1;
t249 = Ifges(7,4) * t210;
t146 = -Ifges(7,6) * t211 + (-Ifges(7,2) * t206 + t249) * t207;
t269 = t146 / 0.2e1;
t250 = Ifges(7,4) * t206;
t147 = -Ifges(7,5) * t211 + (Ifges(7,1) * t210 - t250) * t207;
t268 = t147 / 0.2e1;
t168 = Ifges(7,5) * t206 + Ifges(7,6) * t210;
t267 = t168 / 0.2e1;
t169 = Ifges(6,5) * t207 + Ifges(6,6) * t211;
t266 = t169 / 0.2e1;
t170 = Ifges(7,2) * t210 + t250;
t265 = t170 / 0.2e1;
t172 = Ifges(7,1) * t206 + t249;
t264 = t172 / 0.2e1;
t173 = Ifges(6,1) * t207 + Ifges(6,4) * t211;
t263 = t173 / 0.2e1;
t262 = -t206 / 0.2e1;
t261 = t206 / 0.2e1;
t260 = t210 / 0.2e1;
t256 = pkin(12) * t207;
t255 = pkin(12) * t211;
t254 = pkin(13) * t206;
t253 = pkin(13) * t210;
t32 = -t206 * t55 + t210 * t80;
t33 = t206 * t80 + t210 * t55;
t13 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t35 = mrSges(6,1) * t80 - mrSges(6,3) * t55;
t252 = t13 - t35;
t20 = t212 * t53 + t56 * t235 + t67 * t242;
t16 = pkin(12) * t103 + t20;
t27 = -t199 * t56 + t203 * t67;
t18 = pkin(4) * t80 - pkin(12) * t81 + t27;
t6 = t211 * t16 + t207 * t18;
t84 = -t105 * t206 + t125 * t210;
t85 = t105 * t210 + t125 * t206;
t46 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t87 = mrSges(6,1) * t125 - mrSges(6,3) * t105;
t251 = t46 - t87;
t96 = -t127 * t199 + t203 * t136;
t68 = pkin(4) * t125 - pkin(12) * t126 + t96;
t74 = t212 * t117 + t127 * t235 + t136 * t242;
t71 = pkin(12) * t149 + t74;
t37 = t207 * t68 + t211 * t71;
t157 = -pkin(10) * t238 + t183;
t248 = t157 * mrSges(3,1);
t247 = t158 * mrSges(3,2);
t245 = t154 * t207;
t233 = t206 * t207;
t232 = t207 * t210;
t231 = t206 ^ 2 + t210 ^ 2;
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t54;
t23 = Ifges(6,5) * t55 - Ifges(6,6) * t54 + Ifges(6,3) * t80;
t38 = Ifges(5,5) * t81 - Ifges(5,6) * t80 + Ifges(5,3) * t103;
t24 = Ifges(6,4) * t55 - Ifges(6,2) * t54 + Ifges(6,6) * t80;
t230 = t8 / 0.2e1 - t24 / 0.2e1;
t41 = Ifges(7,5) * t85 + Ifges(7,6) * t84 + Ifges(7,3) * t104;
t62 = Ifges(6,4) * t105 - Ifges(6,2) * t104 + Ifges(6,6) * t125;
t229 = t41 / 0.2e1 - t62 / 0.2e1;
t61 = Ifges(6,5) * t105 - Ifges(6,6) * t104 + Ifges(6,3) * t125;
t91 = Ifges(5,5) * t126 - Ifges(5,6) * t125 + Ifges(5,3) * t149;
t227 = Ifges(3,5) * t238 + Ifges(3,6) * t237 + Ifges(3,3) * t205;
t145 = Ifges(7,5) * t232 - Ifges(7,6) * t233 - Ifges(7,3) * t211;
t171 = Ifges(6,4) * t207 + Ifges(6,2) * t211;
t226 = t145 / 0.2e1 - t171 / 0.2e1;
t94 = -t123 * mrSges(4,1) + t124 * mrSges(4,2);
t4 = pkin(13) * t80 + t6;
t15 = -pkin(4) * t103 - t19;
t7 = pkin(5) * t54 - pkin(13) * t55 + t15;
t1 = -t206 * t4 + t210 * t7;
t2 = t206 * t7 + t210 * t4;
t225 = -t1 * t206 + t2 * t210;
t224 = mrSges(7,1) * t206 + mrSges(7,2) * t210;
t29 = pkin(13) * t125 + t37;
t70 = -pkin(4) * t149 - t73;
t44 = pkin(5) * t104 - pkin(13) * t105 + t70;
t11 = -t206 * t29 + t210 * t44;
t12 = t206 * t44 + t210 * t29;
t223 = -t11 * t206 + t12 * t210;
t5 = -t16 * t207 + t18 * t211;
t36 = -t207 * t71 + t211 * t68;
t165 = -pkin(5) * t211 - pkin(13) * t207 - pkin(4);
t142 = t165 * t210 - t206 * t255;
t143 = t165 * t206 + t210 * t255;
t218 = -t142 * t206 + t143 * t210;
t156 = t203 * t207 + t211 * t242;
t217 = t156 * t211 + t245;
t215 = pkin(12) ^ 2;
t197 = t211 ^ 2;
t195 = t207 ^ 2;
t193 = t199 ^ 2;
t192 = t195 * t215;
t185 = t193 * t212 ^ 2;
t176 = mrSges(4,2) * t244;
t167 = -mrSges(6,1) * t211 + mrSges(6,2) * t207;
t164 = -mrSges(7,1) * t211 - mrSges(7,3) * t232;
t163 = mrSges(7,2) * t211 - mrSges(7,3) * t233;
t161 = -mrSges(4,2) * t204 + mrSges(4,3) * t240;
t160 = mrSges(4,1) * t204 - mrSges(4,3) * t244;
t159 = t224 * t207;
t152 = -mrSges(4,1) * t240 + t176;
t151 = -qJ(3) * t244 + t178;
t141 = Ifges(4,5) * t204 + (Ifges(4,1) * t198 + Ifges(4,4) * t202) * t200;
t140 = Ifges(4,6) * t204 + (Ifges(4,4) * t198 + Ifges(4,2) * t202) * t200;
t139 = Ifges(4,3) * t204 + (Ifges(4,5) * t198 + Ifges(4,6) * t202) * t200;
t133 = t156 * t210 - t206 * t241;
t132 = -t156 * t206 - t210 * t241;
t109 = mrSges(5,1) * t149 - mrSges(5,3) * t126;
t108 = mrSges(4,1) * t150 - mrSges(4,3) * t124;
t107 = -mrSges(4,2) * t150 + mrSges(4,3) * t123;
t106 = -mrSges(5,2) * t149 - mrSges(5,3) * t125;
t95 = mrSges(5,1) * t125 + mrSges(5,2) * t126;
t93 = Ifges(5,1) * t126 - Ifges(5,4) * t125 + Ifges(5,5) * t149;
t92 = Ifges(5,4) * t126 - Ifges(5,2) * t125 + Ifges(5,6) * t149;
t90 = Ifges(4,1) * t124 + Ifges(4,4) * t123 + Ifges(4,5) * t150;
t89 = Ifges(4,4) * t124 + Ifges(4,2) * t123 + Ifges(4,6) * t150;
t88 = Ifges(4,5) * t124 + Ifges(4,6) * t123 + Ifges(4,3) * t150;
t86 = -mrSges(6,2) * t125 - mrSges(6,3) * t104;
t72 = mrSges(6,1) * t104 + mrSges(6,2) * t105;
t60 = mrSges(7,1) * t104 - mrSges(7,3) * t85;
t59 = -mrSges(7,2) * t104 + mrSges(7,3) * t84;
t58 = mrSges(5,1) * t103 - mrSges(5,3) * t81;
t57 = -mrSges(5,2) * t103 - mrSges(5,3) * t80;
t45 = mrSges(5,1) * t80 + mrSges(5,2) * t81;
t43 = Ifges(7,1) * t85 + Ifges(7,4) * t84 + Ifges(7,5) * t104;
t42 = Ifges(7,4) * t85 + Ifges(7,2) * t84 + Ifges(7,6) * t104;
t40 = Ifges(5,1) * t81 - Ifges(5,4) * t80 + Ifges(5,5) * t103;
t39 = Ifges(5,4) * t81 - Ifges(5,2) * t80 + Ifges(5,6) * t103;
t34 = -mrSges(6,2) * t80 - mrSges(6,3) * t54;
t28 = -pkin(5) * t125 - t36;
t26 = mrSges(6,1) * t54 + mrSges(6,2) * t55;
t22 = mrSges(7,1) * t54 - mrSges(7,3) * t33;
t21 = -mrSges(7,2) * t54 + mrSges(7,3) * t32;
t10 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t54;
t9 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t54;
t3 = -pkin(5) * t80 - t5;
t14 = [m(4) * (t75 ^ 2 + t76 ^ 2 + t97 ^ 2) + m(5) * (t19 ^ 2 + t20 ^ 2 + t27 ^ 2) + m(6) * (t15 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (t227 - 0.2e1 * t247 + 0.2e1 * t248) * t205 + 0.2e1 * t97 * t94 + t103 * t38 + t81 * t40 + m(3) * (t157 ^ 2 + t158 ^ 2) + t124 * t90 + t123 * t89 + 0.2e1 * t76 * t107 + 0.2e1 * t75 * t108 + (t23 - t39) * t80 + (t8 - t24) * t54 + ((Ifges(3,5) * t209 + Ifges(3,6) * t213) * t205 + 0.2e1 * (-t157 * t209 + t158 * t213) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t213 + mrSges(3,2) * t209) + t209 * (Ifges(3,1) * t209 + Ifges(3,4) * t213) + t213 * (Ifges(3,4) * t209 + Ifges(3,2) * t213)) * t201) * t201 + Ifges(2,3) + t150 * t88 + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t15 * t26 + t32 * t9 + t33 * t10 + 0.2e1 * t6 * t34 + 0.2e1 * t5 * t35 + 0.2e1 * t27 * t45 + t55 * t25 + 0.2e1 * t20 * t57 + 0.2e1 * t19 * t58; t227 - t247 + t248 + t123 * t140 / 0.2e1 + t124 * t141 / 0.2e1 + (-t92 / 0.2e1 + t61 / 0.2e1) * t80 + (-pkin(2) * t94 + t198 * t90 / 0.2e1 + t202 * t89 / 0.2e1) * t200 + m(5) * (t19 * t73 + t20 * t74 + t27 * t96) + m(6) * (t15 * t70 + t36 * t5 + t37 * t6) + m(7) * (t1 * t11 + t12 * t2 + t28 * t3) + m(4) * (-pkin(2) * t200 * t97 + t151 * t75 + t153 * t76) + t96 * t45 + t103 * t91 / 0.2e1 + t85 * t10 / 0.2e1 + t6 * t86 + t5 * t87 + t81 * t93 / 0.2e1 + t27 * t95 + t84 * t9 / 0.2e1 + t55 * t270 + t105 * t271 + (t23 / 0.2e1 - t39 / 0.2e1) * t125 + t126 * t40 / 0.2e1 + t20 * t106 + t19 * t109 + t70 * t26 + t15 * t72 + t73 * t58 + t74 * t57 + t149 * t38 / 0.2e1 + t150 * t139 / 0.2e1 + t151 * t108 + t97 * t152 + t153 * t107 + t229 * t54 + t230 * t104 + t75 * t160 + t76 * t161 + t12 * t21 + t11 * t22 + t204 * t88 / 0.2e1 + t28 * t13 + t36 * t35 + t37 * t34 + t32 * t42 / 0.2e1 + t33 * t43 / 0.2e1 + t3 * t46 + t2 * t59 + t1 * t60; (t61 - t92) * t125 + (t41 - t62) * t104 + 0.2e1 * t96 * t95 + t85 * t43 + 0.2e1 * t37 * t86 + 0.2e1 * t36 * t87 + t84 * t42 + m(5) * (t73 ^ 2 + t74 ^ 2 + t96 ^ 2) + m(6) * (t36 ^ 2 + t37 ^ 2 + t70 ^ 2) + m(7) * (t11 ^ 2 + t12 ^ 2 + t28 ^ 2) + m(4) * (pkin(2) ^ 2 * t200 ^ 2 + t151 ^ 2 + t153 ^ 2) + t126 * t93 + t105 * t63 + 0.2e1 * t74 * t106 + 0.2e1 * t73 * t109 + 0.2e1 * t70 * t72 + (-0.2e1 * pkin(2) * t152 + t202 * t140 + t198 * t141) * t200 + Ifges(3,3) + t149 * t91 + 0.2e1 * t151 * t160 + 0.2e1 * t153 * t161 + t204 * t139 + 0.2e1 * t28 * t46 + 0.2e1 * t12 * t59 + 0.2e1 * t11 * t60; t132 * t22 + t133 * t21 + t156 * t34 + t203 * t45 + t252 * t154 + (t208 * t57 + (-t26 + t58) * t212) * t199 + m(7) * (t1 * t132 + t133 * t2 + t154 * t3) + m(6) * (-t15 * t241 - t154 * t5 + t156 * t6) + m(5) * (t203 * t27 + (t19 * t212 + t20 * t208) * t199) + m(4) * t97 + t94; t132 * t60 + t133 * t59 + t156 * t86 + t203 * t95 + t176 + (-m(4) * pkin(2) - mrSges(4,1) * t202) * t200 + t251 * t154 + (t208 * t106 + (t109 - t72) * t212) * t199 + m(7) * (t11 * t132 + t12 * t133 + t154 * t28) + m(6) * (-t154 * t36 + t156 * t37 - t241 * t70) + m(5) * (t203 * t96 + (t208 * t74 + t212 * t73) * t199); m(4) + m(6) * (t156 ^ 2 + t185 + t272) + m(7) * (t132 ^ 2 + t133 ^ 2 + t272) + m(5) * (t193 * t208 ^ 2 + t203 ^ 2 + t185); t142 * t22 + t143 * t21 + t32 * t269 + t38 + m(6) * (-pkin(4) * t15 + (-t207 * t5 + t211 * t6) * pkin(12)) + (-t5 * mrSges(6,3) + pkin(12) * t252 + t10 * t260 + t262 * t9 + t271) * t207 + m(7) * (t1 * t142 + t143 * t2 + t256 * t3) + t33 * t268 + t226 * t54 + (t6 * mrSges(6,3) + pkin(12) * t34 - t230) * t211 + t3 * t159 + t2 * t163 + t1 * t164 + t15 * t167 + t80 * t266 + t19 * mrSges(5,1) - t20 * mrSges(5,2) - pkin(4) * t26 + t55 * t263; t142 * t60 + t143 * t59 + t84 * t269 + m(6) * (-pkin(4) * t70 + (-t207 * t36 + t211 * t37) * pkin(12)) + t91 + (-t36 * mrSges(6,3) + pkin(12) * t251 + t260 * t43 + t262 * t42 + t270) * t207 + (t37 * mrSges(6,3) + pkin(12) * t86 - t229) * t211 + m(7) * (t11 * t142 + t12 * t143 + t256 * t28) - pkin(4) * t72 + t73 * mrSges(5,1) - t74 * mrSges(5,2) + t85 * t268 + t226 * t104 + t28 * t159 + t12 * t163 + t11 * t164 + t70 * t167 + t125 * t266 + t105 * t263; t132 * t164 + t133 * t163 + t154 * t159 + t217 * mrSges(6,3) + (-t208 * mrSges(5,2) + (mrSges(5,1) - t167) * t212) * t199 + m(6) * (pkin(4) * t241 + pkin(12) * t217) + m(7) * (pkin(12) * t245 + t132 * t142 + t133 * t143); -0.2e1 * pkin(4) * t167 + 0.2e1 * t142 * t164 + 0.2e1 * t143 * t163 + Ifges(5,3) + (-t145 + t171) * t211 + (t195 + t197) * mrSges(6,3) * t274 + m(7) * (t142 ^ 2 + t143 ^ 2 + t192) + m(6) * (pkin(4) ^ 2 + t197 * t215 + t192) + (-t146 * t206 + t147 * t210 + t159 * t274 + t173) * t207; m(7) * (-pkin(5) * t3 + pkin(13) * t225) + t21 * t253 - t22 * t254 + t9 * t260 + t3 * t166 + t54 * t267 + t33 * t264 + t32 * t265 - pkin(5) * t13 + t10 * t261 - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t225 * mrSges(7,3) + t23; t42 * t260 + t43 * t261 + t59 * t253 - pkin(5) * t46 + m(7) * (-pkin(5) * t28 + pkin(13) * t223) + t84 * t265 + t104 * t267 + t85 * t264 + t28 * t166 - t60 * t254 - t37 * mrSges(6,2) + t36 * mrSges(6,1) + t223 * mrSges(7,3) + t61; -t156 * mrSges(6,2) + (m(7) * pkin(13) + mrSges(7,3)) * (-t132 * t206 + t133 * t210) + t273 * t154; t147 * t261 + t146 * t260 - pkin(5) * t159 + (m(7) * t218 + t210 * t163 - t206 * t164) * pkin(13) + (t273 * pkin(12) + t170 * t262 + t172 * t260) * t207 + (-pkin(12) * mrSges(6,2) - t168 / 0.2e1) * t211 + t218 * mrSges(7,3) + t169; Ifges(6,3) - 0.2e1 * pkin(5) * t166 + m(7) * (pkin(13) ^ 2 * t231 + pkin(5) ^ 2) + t206 * t172 + t210 * t170 + 0.2e1 * t231 * pkin(13) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t11 - mrSges(7,2) * t12 + t41; mrSges(7,1) * t132 - mrSges(7,2) * t133; mrSges(7,1) * t142 - mrSges(7,2) * t143 + t145; -pkin(13) * t224 + t168; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
