% Calculate joint inertia matrix for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 18:42
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:41:30
% EndTime: 2018-11-23 18:41:33
% DurationCPUTime: 2.44s
% Computational Cost: add. (5765->459), mult. (12704->660), div. (0->0), fcn. (14481->12), ass. (0->176)
t203 = sin(qJ(5));
t208 = cos(qJ(5));
t204 = sin(qJ(4));
t205 = sin(qJ(3));
t209 = cos(qJ(4));
t210 = cos(qJ(3));
t163 = t204 * t205 - t209 * t210;
t165 = t204 * t210 + t205 * t209;
t189 = -pkin(3) * t210 - pkin(2);
t114 = pkin(4) * t163 - pkin(11) * t165 + t189;
t268 = -pkin(10) - pkin(9);
t177 = t268 * t210;
t233 = t268 * t205;
t135 = -t209 * t177 + t204 * t233;
t68 = t208 * t114 - t135 * t203;
t69 = t203 * t114 + t208 * t135;
t226 = -t203 * t68 + t208 * t69;
t200 = sin(pkin(6));
t211 = cos(qJ(2));
t240 = t200 * t211;
t201 = cos(pkin(6));
t206 = sin(qJ(2));
t241 = t200 * t206;
t149 = t201 * t205 + t210 * t241;
t151 = t201 * t206 * pkin(1) + pkin(8) * t240;
t141 = pkin(9) * t201 + t151;
t142 = (-pkin(2) * t211 - pkin(9) * t206 - pkin(1)) * t200;
t87 = -t141 * t205 + t210 * t142;
t62 = -pkin(3) * t240 - pkin(10) * t149 + t87;
t148 = t201 * t210 - t205 * t241;
t88 = t210 * t141 + t205 * t142;
t70 = pkin(10) * t148 + t88;
t37 = t204 * t62 + t209 * t70;
t29 = -pkin(11) * t240 + t37;
t101 = -t209 * t148 + t149 * t204;
t102 = t148 * t204 + t149 * t209;
t180 = pkin(8) * t241;
t257 = pkin(1) * t211;
t140 = t180 + (-pkin(2) - t257) * t201;
t107 = -pkin(3) * t148 + t140;
t40 = pkin(4) * t101 - pkin(11) * t102 + t107;
t10 = -t203 * t29 + t208 * t40;
t11 = t203 * t40 + t208 * t29;
t228 = -t10 * t203 + t11 * t208;
t202 = sin(qJ(6));
t207 = cos(qJ(6));
t162 = -t202 * t203 + t207 * t208;
t164 = t202 * t208 + t203 * t207;
t119 = Ifges(7,5) * t164 + Ifges(7,6) * t162;
t170 = Ifges(6,5) * t203 + Ifges(6,6) * t208;
t275 = t119 / 0.2e1 + t170 / 0.2e1;
t247 = t162 * mrSges(7,3);
t274 = t202 * pkin(5) * t247 + t170;
t273 = -Ifges(4,5) * t149 - Ifges(4,6) * t148;
t132 = -t177 * t204 - t209 * t233;
t272 = t132 ^ 2;
t117 = -mrSges(7,1) * t162 + mrSges(7,2) * t164;
t271 = 0.2e1 * t117;
t270 = 0.2e1 * t132;
t79 = -t102 * t203 - t208 * t240;
t269 = t79 / 0.2e1;
t120 = Ifges(7,4) * t164 + Ifges(7,2) * t162;
t266 = t120 / 0.2e1;
t122 = Ifges(7,1) * t164 + Ifges(7,4) * t162;
t265 = t122 / 0.2e1;
t264 = t162 / 0.2e1;
t263 = t164 / 0.2e1;
t252 = Ifges(6,4) * t208;
t173 = Ifges(6,1) * t203 + t252;
t261 = t173 / 0.2e1;
t260 = t203 / 0.2e1;
t259 = t208 / 0.2e1;
t256 = pkin(3) * t209;
t255 = Ifges(4,3) + Ifges(5,3);
t254 = -Ifges(5,5) * t102 + Ifges(5,6) * t101;
t253 = Ifges(6,4) * t203;
t150 = t201 * t257 - t180;
t249 = t150 * mrSges(3,1);
t248 = t151 * mrSges(3,2);
t246 = t164 * mrSges(7,3);
t243 = t165 * t203;
t242 = t165 * t208;
t239 = Ifges(5,5) * t165 - Ifges(5,6) * t163;
t238 = Ifges(4,5) * t205 + Ifges(4,6) * t210;
t237 = t203 ^ 2 + t208 ^ 2;
t236 = t205 ^ 2 + t210 ^ 2;
t235 = 0.2e1 * mrSges(7,3);
t80 = t102 * t208 - t203 * t240;
t43 = -t202 * t80 + t207 * t79;
t44 = t202 * t79 + t207 * t80;
t12 = Ifges(7,5) * t44 + Ifges(7,6) * t43 + Ifges(7,3) * t101;
t30 = Ifges(6,5) * t80 + Ifges(6,6) * t79 + Ifges(6,3) * t101;
t234 = t207 * t246;
t103 = t164 * t165;
t104 = t162 * t165;
t51 = Ifges(7,5) * t104 - Ifges(7,6) * t103 + Ifges(7,3) * t163;
t232 = Ifges(3,5) * t241 + Ifges(3,6) * t240 + Ifges(3,3) * t201;
t188 = -pkin(5) * t208 - pkin(4);
t36 = -t204 * t70 + t209 * t62;
t186 = pkin(3) * t204 + pkin(11);
t231 = t237 * t186;
t28 = pkin(4) * t240 - t36;
t171 = Ifges(6,2) * t208 + t253;
t230 = t162 * t120 + t164 * t122 + t208 * t171 + t203 * t173 + Ifges(5,3);
t229 = mrSges(6,1) * t203 + mrSges(6,2) * t208;
t49 = -mrSges(6,2) * t101 + mrSges(6,3) * t79;
t50 = mrSges(6,1) * t101 - mrSges(6,3) * t80;
t227 = -t203 * t50 + t208 * t49;
t225 = 0.2e1 * mrSges(6,3) * t237;
t115 = -mrSges(6,2) * t163 - mrSges(6,3) * t243;
t116 = mrSges(6,1) * t163 - mrSges(6,3) * t242;
t224 = t208 * t115 - t203 * t116;
t82 = Ifges(6,5) * t242 - Ifges(6,6) * t243 + Ifges(6,3) * t163;
t158 = (-pkin(12) - t186) * t203;
t195 = t208 * pkin(12);
t159 = t186 * t208 + t195;
t111 = t158 * t207 - t159 * t202;
t112 = t158 * t202 + t159 * t207;
t223 = t111 * mrSges(7,1) - t112 * mrSges(7,2) + t119;
t175 = (-pkin(12) - pkin(11)) * t203;
t176 = pkin(11) * t208 + t195;
t131 = t175 * t207 - t176 * t202;
t134 = t175 * t202 + t176 * t207;
t222 = t131 * mrSges(7,1) - t134 * mrSges(7,2) + t119;
t5 = pkin(5) * t101 - pkin(12) * t80 + t10;
t6 = pkin(12) * t79 + t11;
t3 = -t202 * t6 + t207 * t5;
t4 = t202 * t5 + t207 * t6;
t221 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t12;
t48 = pkin(5) * t163 - pkin(12) * t242 + t68;
t56 = -pkin(12) * t243 + t69;
t20 = -t202 * t56 + t207 * t48;
t21 = t202 * t48 + t207 * t56;
t220 = t20 * mrSges(7,1) - t21 * mrSges(7,2) + t51;
t219 = (mrSges(5,1) * t209 - mrSges(5,2) * t204) * pkin(3);
t218 = (mrSges(7,1) * t207 - mrSges(7,2) * t202) * pkin(5);
t13 = Ifges(7,4) * t44 + Ifges(7,2) * t43 + Ifges(7,6) * t101;
t14 = Ifges(7,1) * t44 + Ifges(7,4) * t43 + Ifges(7,5) * t101;
t168 = -mrSges(6,1) * t208 + mrSges(6,2) * t203;
t17 = -pkin(5) * t79 + t28;
t31 = Ifges(6,4) * t80 + Ifges(6,2) * t79 + Ifges(6,6) * t101;
t32 = Ifges(6,1) * t80 + Ifges(6,4) * t79 + Ifges(6,5) * t101;
t217 = t36 * mrSges(5,1) - t37 * mrSges(5,2) + t17 * t117 + t13 * t264 + t14 * t263 + t28 * t168 + t171 * t269 - t246 * t3 + t4 * t247 + t31 * t259 + t32 * t260 + t80 * t261 + t44 * t265 + t43 * t266 - t254 + t275 * t101 + t228 * mrSges(6,3);
t52 = Ifges(7,4) * t104 - Ifges(7,2) * t103 + Ifges(7,6) * t163;
t53 = Ifges(7,1) * t104 - Ifges(7,4) * t103 + Ifges(7,5) * t163;
t83 = Ifges(6,6) * t163 + (-Ifges(6,2) * t203 + t252) * t165;
t84 = Ifges(6,5) * t163 + (Ifges(6,1) * t208 - t253) * t165;
t94 = pkin(5) * t243 + t132;
t216 = -t135 * mrSges(5,2) - t103 * t266 + t104 * t265 + t94 * t117 - t20 * t246 + t21 * t247 + t83 * t259 + t84 * t260 + t53 * t263 + t52 * t264 - t171 * t243 / 0.2e1 + t242 * t261 + t239 + t275 * t163 + (t168 - mrSges(5,1)) * t132 + t226 * mrSges(6,3);
t187 = -pkin(4) - t256;
t174 = Ifges(4,1) * t205 + Ifges(4,4) * t210;
t172 = Ifges(4,4) * t205 + Ifges(4,2) * t210;
t169 = -mrSges(4,1) * t210 + mrSges(4,2) * t205;
t167 = t188 - t256;
t130 = -mrSges(4,1) * t240 - mrSges(4,3) * t149;
t129 = mrSges(4,2) * t240 + mrSges(4,3) * t148;
t123 = Ifges(5,1) * t165 - Ifges(5,4) * t163;
t121 = Ifges(5,4) * t165 - Ifges(5,2) * t163;
t118 = mrSges(5,1) * t163 + mrSges(5,2) * t165;
t113 = t229 * t165;
t109 = -mrSges(4,1) * t148 + mrSges(4,2) * t149;
t91 = Ifges(4,1) * t149 + Ifges(4,4) * t148 - Ifges(4,5) * t240;
t90 = Ifges(4,4) * t149 + Ifges(4,2) * t148 - Ifges(4,6) * t240;
t86 = -mrSges(5,1) * t240 - mrSges(5,3) * t102;
t85 = mrSges(5,2) * t240 - mrSges(5,3) * t101;
t74 = mrSges(7,1) * t163 - mrSges(7,3) * t104;
t73 = -mrSges(7,2) * t163 - mrSges(7,3) * t103;
t58 = mrSges(7,1) * t103 + mrSges(7,2) * t104;
t57 = mrSges(5,1) * t101 + mrSges(5,2) * t102;
t55 = Ifges(5,1) * t102 - Ifges(5,4) * t101 - Ifges(5,5) * t240;
t54 = Ifges(5,4) * t102 - Ifges(5,2) * t101 - Ifges(5,6) * t240;
t47 = -mrSges(6,1) * t79 + mrSges(6,2) * t80;
t24 = mrSges(7,1) * t101 - mrSges(7,3) * t44;
t23 = -mrSges(7,2) * t101 + mrSges(7,3) * t43;
t16 = -mrSges(7,1) * t43 + mrSges(7,2) * t44;
t1 = [m(7) * (t17 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t28 ^ 2) + m(5) * (t107 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(4) * (t140 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(3) * (pkin(1) ^ 2 * t200 ^ 2 + t150 ^ 2 + t151 ^ 2) + (t12 + t30 - t54) * t101 + Ifges(2,3) + ((-0.2e1 * t150 * mrSges(3,3) + Ifges(3,5) * t201 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t206) * t200) * t206 + (0.2e1 * t151 * mrSges(3,3) + Ifges(3,6) * t201 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t206 + (Ifges(3,2) + t255) * t211) * t200 + t254 + t273) * t211) * t200 + (t232 - 0.2e1 * t248 + 0.2e1 * t249) * t201 + 0.2e1 * t17 * t16 + 0.2e1 * t4 * t23 + 0.2e1 * t3 * t24 + t43 * t13 + t44 * t14 + 0.2e1 * t28 * t47 + 0.2e1 * t11 * t49 + 0.2e1 * t10 * t50 + t79 * t31 + t80 * t32 + 0.2e1 * t37 * t85 + 0.2e1 * t36 * t86 + t102 * t55 + 0.2e1 * t107 * t57 + 0.2e1 * t88 * t129 + 0.2e1 * t87 * t130 + 0.2e1 * t140 * t109 + t148 * t90 + t149 * t91; t249 - t248 + t232 + (t88 * mrSges(4,3) + pkin(9) * t129 + t90 / 0.2e1) * t210 + (-t87 * mrSges(4,3) - pkin(9) * t130 + t91 / 0.2e1) * t205 + (t47 - t86) * t132 + (t51 / 0.2e1 + t82 / 0.2e1 - t121 / 0.2e1) * t101 + m(7) * (t17 * t94 + t20 * t3 + t21 * t4) + m(6) * (t10 * t68 + t11 * t69 + t132 * t28) + m(5) * (t107 * t189 - t132 * t36 + t135 * t37) + (-t37 * mrSges(5,3) + t12 / 0.2e1 + t30 / 0.2e1 - t54 / 0.2e1) * t163 + m(4) * (-pkin(2) * t140 + (-t205 * t87 + t210 * t88) * pkin(9)) + t83 * t269 - (t238 + t239) * t240 / 0.2e1 + (-t36 * mrSges(5,3) - t203 * t31 / 0.2e1 + t32 * t259 + t55 / 0.2e1) * t165 + t21 * t23 + t20 * t24 + t43 * t52 / 0.2e1 + t44 * t53 / 0.2e1 + t17 * t58 + t68 * t50 + t69 * t49 + t4 * t73 + t3 * t74 + t80 * t84 / 0.2e1 + t94 * t16 - t103 * t13 / 0.2e1 + t104 * t14 / 0.2e1 - pkin(2) * t109 + t28 * t113 + t11 * t115 + t10 * t116 + t107 * t118 + t102 * t123 / 0.2e1 + t135 * t85 + t140 * t169 + t148 * t172 / 0.2e1 + t149 * t174 / 0.2e1 + t189 * t57; -0.2e1 * pkin(2) * t169 - t103 * t52 + t104 * t53 + t113 * t270 + 0.2e1 * t69 * t115 + 0.2e1 * t68 * t116 + 0.2e1 * t189 * t118 + t210 * t172 + t205 * t174 + 0.2e1 * t20 * t74 + 0.2e1 * t21 * t73 + 0.2e1 * t94 * t58 + Ifges(3,3) + 0.2e1 * t236 * pkin(9) * mrSges(4,3) + (mrSges(5,3) * t270 - t203 * t83 + t208 * t84 + t123) * t165 + (-0.2e1 * mrSges(5,3) * t135 - t121 + t51 + t82) * t163 + m(7) * (t20 ^ 2 + t21 ^ 2 + t94 ^ 2) + m(6) * (t68 ^ 2 + t69 ^ 2 + t272) + m(5) * (t135 ^ 2 + t189 ^ 2 + t272) + m(4) * (pkin(9) ^ 2 * t236 + pkin(2) ^ 2); m(7) * (t111 * t3 + t112 * t4 + t167 * t17) - t255 * t240 + (t204 * t85 + t209 * t86 + m(5) * (t204 * t37 + t209 * t36)) * pkin(3) + t217 + m(6) * (t186 * t228 + t187 * t28) + t227 * t186 + t87 * mrSges(4,1) - t88 * mrSges(4,2) + t111 * t24 + t112 * t23 + t167 * t16 + t187 * t47 - t273; t224 * t186 + m(7) * (t111 * t20 + t112 * t21 + t167 * t94) + (-mrSges(4,1) * t205 - mrSges(4,2) * t210) * pkin(9) + (m(5) * (-t132 * t209 + t135 * t204) + (-t163 * t204 - t165 * t209) * mrSges(5,3)) * pkin(3) + t216 + m(6) * (t132 * t187 + t186 * t226) + t111 * t74 + t112 * t73 + t167 * t58 + t187 * t113 + t238; t167 * t271 + 0.2e1 * t187 * t168 + Ifges(4,3) + 0.2e1 * t219 + (-t111 * t164 + t112 * t162) * t235 + t186 * t225 + m(7) * (t111 ^ 2 + t112 ^ 2 + t167 ^ 2) + m(6) * (t186 ^ 2 * t237 + t187 ^ 2) + m(5) * (t204 ^ 2 + t209 ^ 2) * pkin(3) ^ 2 + t230; t227 * pkin(11) + m(7) * (t131 * t3 + t134 * t4 + t17 * t188) + t217 + m(6) * (-pkin(4) * t28 + pkin(11) * t228) - Ifges(5,3) * t240 - pkin(4) * t47 + t131 * t24 + t134 * t23 + t188 * t16; t216 + m(6) * (-pkin(4) * t132 + pkin(11) * t226) - pkin(4) * t113 + t131 * t74 + t134 * t73 + t188 * t58 + m(7) * (t131 * t20 + t134 * t21 + t188 * t94) + t224 * pkin(11); (t187 - pkin(4)) * t168 + (t167 + t188) * t117 + t219 + m(7) * (t111 * t131 + t112 * t134 + t167 * t188) + m(6) * (-pkin(4) * t187 + pkin(11) * t231) + ((-t111 - t131) * t164 + (t112 + t134) * t162) * mrSges(7,3) + (pkin(11) * t237 + t231) * mrSges(6,3) + t230; -0.2e1 * pkin(4) * t168 + t188 * t271 + (-t131 * t164 + t134 * t162) * t235 + pkin(11) * t225 + m(7) * (t131 ^ 2 + t134 ^ 2 + t188 ^ 2) + m(6) * (pkin(11) ^ 2 * t237 + pkin(4) ^ 2) + t230; t10 * mrSges(6,1) - t11 * mrSges(6,2) + (m(7) * (t202 * t4 + t207 * t3) + t202 * t23 + t207 * t24) * pkin(5) + t221 + t30; t68 * mrSges(6,1) - t69 * mrSges(6,2) + (t202 * t73 + m(7) * (t20 * t207 + t202 * t21) + t207 * t74) * pkin(5) + t220 + t82; -t229 * t186 + (-t234 + m(7) * (t111 * t207 + t112 * t202)) * pkin(5) + t223 + t274; -t229 * pkin(11) + (-t234 + m(7) * (t131 * t207 + t134 * t202)) * pkin(5) + t222 + t274; Ifges(6,3) + Ifges(7,3) + m(7) * (t202 ^ 2 + t207 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t218; t221; t220; t223; t222; Ifges(7,3) + t218; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
