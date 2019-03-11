% Calculate joint inertia matrix for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:33
% EndTime: 2019-03-10 01:34:38
% DurationCPUTime: 2.59s
% Computational Cost: add. (3933->439), mult. (8702->602), div. (0->0), fcn. (9599->10), ass. (0->160)
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t230 = Ifges(7,4) * t190;
t232 = Ifges(6,4) * t190;
t254 = t230 + t232 + (Ifges(6,1) + Ifges(7,1)) * t186;
t260 = t254 / 0.2e1;
t184 = sin(pkin(6));
t193 = cos(qJ(2));
t220 = t184 * t193;
t185 = cos(pkin(6));
t188 = sin(qJ(3));
t192 = cos(qJ(3));
t189 = sin(qJ(2));
t221 = t184 * t189;
t127 = t185 * t192 - t188 * t221;
t128 = t185 * t188 + t192 * t221;
t187 = sin(qJ(4));
t191 = cos(qJ(4));
t85 = t127 * t187 + t128 * t191;
t63 = -t186 * t85 - t190 * t220;
t64 = -t186 * t220 + t190 * t85;
t84 = -t127 * t191 + t128 * t187;
t19 = Ifges(7,4) * t64 + Ifges(7,2) * t63 + Ifges(7,6) * t84;
t20 = Ifges(6,4) * t64 + Ifges(6,2) * t63 + Ifges(6,6) * t84;
t259 = t19 + t20;
t21 = Ifges(7,1) * t64 + Ifges(7,4) * t63 + Ifges(7,5) * t84;
t22 = Ifges(6,1) * t64 + Ifges(6,4) * t63 + Ifges(6,5) * t84;
t258 = t21 + t22;
t141 = t187 * t188 - t191 * t192;
t142 = t187 * t192 + t188 * t191;
t69 = Ifges(7,6) * t141 + (-Ifges(7,2) * t186 + t230) * t142;
t70 = Ifges(6,6) * t141 + (-Ifges(6,2) * t186 + t232) * t142;
t257 = t69 + t70;
t231 = Ifges(7,4) * t186;
t71 = Ifges(7,5) * t141 + (Ifges(7,1) * t190 - t231) * t142;
t233 = Ifges(6,4) * t186;
t72 = Ifges(6,5) * t141 + (Ifges(6,1) * t190 - t233) * t142;
t256 = t71 + t72;
t255 = t231 + t233 + (Ifges(6,2) + Ifges(7,2)) * t190;
t150 = Ifges(7,5) * t186 + Ifges(7,6) * t190;
t151 = Ifges(6,5) * t186 + Ifges(6,6) * t190;
t253 = t151 / 0.2e1 + t150 / 0.2e1;
t17 = Ifges(7,5) * t64 + Ifges(7,6) * t63 + Ifges(7,3) * t84;
t18 = Ifges(6,5) * t64 + Ifges(6,6) * t63 + Ifges(6,3) * t84;
t252 = t17 + t18;
t251 = -Ifges(4,5) * t128 - Ifges(4,6) * t127;
t245 = -pkin(10) - pkin(9);
t158 = t245 * t192;
t212 = t245 * t188;
t109 = -t158 * t187 - t191 * t212;
t250 = t109 ^ 2;
t249 = 0.2e1 * t109;
t146 = -mrSges(7,1) * t190 + mrSges(7,2) * t186;
t248 = 0.2e1 * t146;
t242 = t186 / 0.2e1;
t241 = t190 / 0.2e1;
t239 = pkin(1) * t193;
t238 = pkin(3) * t191;
t237 = pkin(11) * t190;
t130 = pkin(1) * t185 * t189 + pkin(8) * t220;
t118 = pkin(9) * t185 + t130;
t119 = (-pkin(2) * t193 - pkin(9) * t189 - pkin(1)) * t184;
t75 = -t118 * t188 + t119 * t192;
t43 = -pkin(3) * t220 - pkin(10) * t128 + t75;
t76 = t118 * t192 + t119 * t188;
t48 = pkin(10) * t127 + t76;
t25 = t187 * t43 + t191 * t48;
t16 = -pkin(11) * t220 + t25;
t161 = pkin(8) * t221;
t117 = t161 + (-pkin(2) - t239) * t185;
t88 = -pkin(3) * t127 + t117;
t28 = pkin(4) * t84 - pkin(11) * t85 + t88;
t6 = t16 * t190 + t186 * t28;
t236 = t190 * t6;
t235 = Ifges(4,3) + Ifges(5,3);
t234 = -Ifges(5,5) * t85 + Ifges(5,6) * t84;
t129 = t185 * t239 - t161;
t229 = t129 * mrSges(3,1);
t228 = t130 * mrSges(3,2);
t227 = t186 * mrSges(7,3);
t226 = t190 * mrSges(7,3);
t111 = -t158 * t191 + t187 * t212;
t170 = -pkin(3) * t192 - pkin(2);
t93 = pkin(4) * t141 - pkin(11) * t142 + t170;
t47 = t111 * t190 + t186 * t93;
t225 = t190 * t47;
t224 = t142 * t186;
t223 = t142 * t190;
t167 = pkin(3) * t187 + pkin(11);
t222 = t167 * t190;
t172 = t190 * qJ(6);
t91 = mrSges(7,1) * t224 + mrSges(7,2) * t223;
t219 = Ifges(7,5) * t223 + Ifges(7,3) * t141;
t218 = Ifges(6,5) * t223 + Ifges(6,3) * t141;
t217 = Ifges(5,5) * t142 - Ifges(5,6) * t141;
t216 = Ifges(4,5) * t188 + Ifges(4,6) * t192;
t215 = t186 ^ 2 + t190 ^ 2;
t214 = t188 ^ 2 + t192 ^ 2;
t213 = 0.2e1 * mrSges(7,3);
t211 = Ifges(3,5) * t221 + Ifges(3,6) * t220 + Ifges(3,3) * t185;
t169 = -pkin(5) * t190 - pkin(4);
t29 = -mrSges(7,1) * t63 + mrSges(7,2) * t64;
t5 = -t16 * t186 + t190 * t28;
t24 = -t187 * t48 + t191 * t43;
t46 = -t111 * t186 + t190 * t93;
t208 = t215 * t167;
t207 = t151 + t150;
t1 = pkin(5) * t84 - qJ(6) * t64 + t5;
t206 = -mrSges(6,3) * t5 - mrSges(7,3) * t1;
t15 = pkin(4) * t220 - t24;
t31 = pkin(5) * t141 - t142 * t172 + t46;
t205 = -t46 * mrSges(6,3) - t31 * mrSges(7,3);
t204 = -t186 * t5 + t236;
t203 = t186 * t254 + t190 * t255 + Ifges(5,3);
t202 = mrSges(6,1) * t186 + mrSges(6,2) * t190;
t201 = -t186 * t46 + t225;
t200 = 0.2e1 * t215 * mrSges(6,3);
t199 = (mrSges(5,1) * t191 - mrSges(5,2) * t187) * pkin(3);
t147 = -mrSges(6,1) * t190 + mrSges(6,2) * t186;
t3 = qJ(6) * t63 + t6;
t8 = -pkin(5) * t63 + t15;
t198 = t24 * mrSges(5,1) - t25 * mrSges(5,2) + mrSges(6,3) * t236 + t8 * t146 + t15 * t147 + t3 * t226 - t234 + t253 * t84 + t255 * t63 / 0.2e1 + t64 * t260 + t258 * t242 + t259 * t241;
t39 = -qJ(6) * t224 + t47;
t79 = pkin(5) * t224 + t109;
t197 = -t111 * mrSges(5,2) + mrSges(6,3) * t225 + t79 * t146 + t39 * t226 + t217 + t256 * t242 + t257 * t241 - t255 * t224 / 0.2e1 + t223 * t260 + t253 * t141 + (-mrSges(5,1) + t147) * t109;
t168 = -pkin(4) - t238;
t157 = Ifges(4,1) * t188 + Ifges(4,4) * t192;
t154 = Ifges(4,4) * t188 + Ifges(4,2) * t192;
t149 = t172 + t237;
t148 = -mrSges(4,1) * t192 + mrSges(4,2) * t188;
t145 = (-qJ(6) - pkin(11)) * t186;
t144 = t169 - t238;
t132 = t172 + t222;
t131 = (-qJ(6) - t167) * t186;
t108 = -mrSges(4,1) * t220 - mrSges(4,3) * t128;
t107 = mrSges(4,2) * t220 + mrSges(4,3) * t127;
t100 = Ifges(5,1) * t142 - Ifges(5,4) * t141;
t99 = Ifges(5,4) * t142 - Ifges(5,2) * t141;
t98 = mrSges(5,1) * t141 + mrSges(5,2) * t142;
t97 = mrSges(6,1) * t141 - mrSges(6,3) * t223;
t96 = mrSges(7,1) * t141 - mrSges(7,3) * t223;
t95 = -mrSges(6,2) * t141 - mrSges(6,3) * t224;
t94 = -mrSges(7,2) * t141 - mrSges(7,3) * t224;
t92 = t202 * t142;
t89 = -mrSges(4,1) * t127 + mrSges(4,2) * t128;
t78 = Ifges(4,1) * t128 + Ifges(4,4) * t127 - Ifges(4,5) * t220;
t77 = Ifges(4,4) * t128 + Ifges(4,2) * t127 - Ifges(4,6) * t220;
t74 = -mrSges(5,1) * t220 - mrSges(5,3) * t85;
t73 = mrSges(5,2) * t220 - mrSges(5,3) * t84;
t68 = -Ifges(6,6) * t224 + t218;
t67 = -Ifges(7,6) * t224 + t219;
t40 = mrSges(5,1) * t84 + mrSges(5,2) * t85;
t37 = Ifges(5,1) * t85 - Ifges(5,4) * t84 - Ifges(5,5) * t220;
t36 = Ifges(5,4) * t85 - Ifges(5,2) * t84 - Ifges(5,6) * t220;
t35 = mrSges(6,1) * t84 - mrSges(6,3) * t64;
t34 = mrSges(7,1) * t84 - mrSges(7,3) * t64;
t33 = -mrSges(6,2) * t84 + mrSges(6,3) * t63;
t32 = -mrSges(7,2) * t84 + mrSges(7,3) * t63;
t30 = -mrSges(6,1) * t63 + mrSges(6,2) * t64;
t2 = [(-t36 + t252) * t84 + ((-0.2e1 * t129 * mrSges(3,3) + Ifges(3,5) * t185 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t189) * t184) * t189 + (0.2e1 * t130 * mrSges(3,3) + Ifges(3,6) * t185 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t189 + (Ifges(3,2) + t235) * t193) * t184 + t234 + t251) * t193) * t184 + (t211 - 0.2e1 * t228 + 0.2e1 * t229) * t185 + t258 * t64 + t259 * t63 + 0.2e1 * t8 * t29 + 0.2e1 * t15 * t30 + 0.2e1 * t3 * t32 + 0.2e1 * t6 * t33 + 0.2e1 * t1 * t34 + 0.2e1 * t5 * t35 + m(6) * (t15 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t3 ^ 2 + t8 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t88 ^ 2) + m(4) * (t117 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(3) * (pkin(1) ^ 2 * t184 ^ 2 + t129 ^ 2 + t130 ^ 2) + 0.2e1 * t25 * t73 + 0.2e1 * t24 * t74 + t85 * t37 + 0.2e1 * t88 * t40 + Ifges(2,3) + 0.2e1 * t76 * t107 + 0.2e1 * t75 * t108 + 0.2e1 * t117 * t89 + t127 * t77 + t128 * t78; -(t216 + t217) * t220 / 0.2e1 - t228 + t229 + t211 + (t30 - t74) * t109 + (-t25 * mrSges(5,3) + t17 / 0.2e1 + t18 / 0.2e1 - t36 / 0.2e1) * t141 + t39 * t32 + t46 * t35 + t47 * t33 + t31 * t34 + (-t24 * mrSges(5,3) + t37 / 0.2e1 + (t21 / 0.2e1 + t22 / 0.2e1) * t190 + (-t19 / 0.2e1 - t20 / 0.2e1) * t186) * t142 + m(4) * (-pkin(2) * t117 + (-t188 * t75 + t192 * t76) * pkin(9)) + t79 * t29 - pkin(2) * t89 + t8 * t91 + t15 * t92 + t3 * t94 + t6 * t95 + t1 * t96 + t5 * t97 + t88 * t98 + t85 * t100 / 0.2e1 + (t76 * mrSges(4,3) + pkin(9) * t107 + t77 / 0.2e1) * t192 + (-t75 * mrSges(4,3) - pkin(9) * t108 + t78 / 0.2e1) * t188 + t111 * t73 + t117 * t148 + t127 * t154 / 0.2e1 + t128 * t157 / 0.2e1 + t170 * t40 + (t69 / 0.2e1 + t70 / 0.2e1) * t63 + (t71 / 0.2e1 + t72 / 0.2e1) * t64 + (t67 / 0.2e1 + t68 / 0.2e1 - t99 / 0.2e1) * t84 + m(7) * (t1 * t31 + t3 * t39 + t79 * t8) + m(6) * (t109 * t15 + t46 * t5 + t47 * t6) + m(5) * (-t109 * t24 + t111 * t25 + t170 * t88); -0.2e1 * pkin(2) * t148 + t92 * t249 + t192 * t154 + t188 * t157 + 0.2e1 * t170 * t98 + 0.2e1 * t31 * t96 + 0.2e1 * t39 * t94 + 0.2e1 * t46 * t97 + 0.2e1 * t47 * t95 + 0.2e1 * t79 * t91 + Ifges(3,3) + 0.2e1 * t214 * pkin(9) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t111 + t67 + t68 - t99) * t141 + m(7) * (t31 ^ 2 + t39 ^ 2 + t79 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2 + t250) + m(5) * (t111 ^ 2 + t170 ^ 2 + t250) + m(4) * (pkin(9) ^ 2 * t214 + pkin(2) ^ 2) + (mrSges(5,3) * t249 - t186 * t257 + t190 * t256 + t100) * t142; t198 + (t187 * t73 + t191 * t74 + m(5) * (t187 * t25 + t191 * t24)) * pkin(3) - t235 * t220 + (-t167 * t35 + t206) * t186 + m(7) * (t1 * t131 + t132 * t3 + t144 * t8) + m(6) * (t15 * t168 + t167 * t204) + t75 * mrSges(4,1) - t76 * mrSges(4,2) + t33 * t222 + t131 * t34 + t132 * t32 + t144 * t29 + t168 * t30 - t251; t197 + (-mrSges(4,1) * t188 - mrSges(4,2) * t192) * pkin(9) + (m(5) * (-t109 * t191 + t111 * t187) + (-t141 * t187 - t142 * t191) * mrSges(5,3)) * pkin(3) + m(6) * (t109 * t168 + t167 * t201) + (-t167 * t97 + t205) * t186 + t95 * t222 + t131 * t96 + t132 * t94 + t144 * t91 + t168 * t92 + m(7) * (t131 * t31 + t132 * t39 + t144 * t79) + t216; t144 * t248 + 0.2e1 * t168 * t147 + Ifges(4,3) + 0.2e1 * t199 + (-t131 * t186 + t132 * t190) * t213 + t167 * t200 + m(6) * (t167 ^ 2 * t215 + t168 ^ 2) + m(7) * (t131 ^ 2 + t132 ^ 2 + t144 ^ 2) + m(5) * (t187 ^ 2 + t191 ^ 2) * pkin(3) ^ 2 + t203; t198 + (-pkin(11) * t35 + t206) * t186 - pkin(4) * t30 + m(7) * (t1 * t145 + t149 * t3 + t169 * t8) + m(6) * (-pkin(4) * t15 + pkin(11) * t204) + t33 * t237 - Ifges(5,3) * t220 + t145 * t34 + t149 * t32 + t169 * t29; t197 + (-pkin(11) * t97 + t205) * t186 + m(7) * (t145 * t31 + t149 * t39 + t169 * t79) + m(6) * (-pkin(4) * t109 + pkin(11) * t201) - pkin(4) * t92 + t95 * t237 + t145 * t96 + t149 * t94 + t169 * t91; (t168 - pkin(4)) * t147 + (t144 + t169) * t146 + t199 + m(6) * (-pkin(4) * t168 + pkin(11) * t208) + m(7) * (t131 * t145 + t132 * t149 + t144 * t169) + ((t132 + t149) * t190 + (-t131 - t145) * t186) * mrSges(7,3) + (pkin(11) * t215 + t208) * mrSges(6,3) + t203; -0.2e1 * pkin(4) * t147 + t169 * t248 + (-t145 * t186 + t149 * t190) * t213 + pkin(11) * t200 + m(7) * (t145 ^ 2 + t149 ^ 2 + t169 ^ 2) + m(6) * (pkin(11) ^ 2 * t215 + pkin(4) ^ 2) + t203; mrSges(6,1) * t5 + mrSges(7,1) * t1 - mrSges(6,2) * t6 - mrSges(7,2) * t3 + (m(7) * t1 + t34) * pkin(5) + t252; mrSges(6,1) * t46 + mrSges(7,1) * t31 - mrSges(6,2) * t47 - mrSges(7,2) * t39 + (-Ifges(6,6) - Ifges(7,6)) * t224 + (m(7) * t31 + t96) * pkin(5) + t218 + t219; mrSges(7,1) * t131 - mrSges(7,2) * t132 - t202 * t167 + (m(7) * t131 - t227) * pkin(5) + t207; mrSges(7,1) * t145 - mrSges(7,2) * t149 - t202 * pkin(11) + (m(7) * t145 - t227) * pkin(5) + t207; Ifges(6,3) + Ifges(7,3) + (m(7) * pkin(5) + 0.2e1 * mrSges(7,1)) * pkin(5); m(7) * t8 + t29; m(7) * t79 + t91; m(7) * t144 + t146; m(7) * t169 + t146; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t2(1) t2(2) t2(4) t2(7) t2(11) t2(16); t2(2) t2(3) t2(5) t2(8) t2(12) t2(17); t2(4) t2(5) t2(6) t2(9) t2(13) t2(18); t2(7) t2(8) t2(9) t2(10) t2(14) t2(19); t2(11) t2(12) t2(13) t2(14) t2(15) t2(20); t2(16) t2(17) t2(18) t2(19) t2(20) t2(21);];
Mq  = res;
