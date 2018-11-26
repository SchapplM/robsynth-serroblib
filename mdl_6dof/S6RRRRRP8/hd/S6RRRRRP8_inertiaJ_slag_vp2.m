% Calculate joint inertia matrix for
% S6RRRRRP8
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:32
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:31:51
% EndTime: 2018-11-23 18:31:53
% DurationCPUTime: 2.40s
% Computational Cost: add. (3969->416), mult. (8790->575), div. (0->0), fcn. (9666->10), ass. (0->160)
t186 = sin(qJ(5));
t190 = cos(qJ(5));
t236 = Ifges(7,5) * t190;
t238 = Ifges(6,4) * t190;
t269 = -t236 + t238 + (Ifges(6,1) + Ifges(7,1)) * t186;
t272 = t269 / 0.2e1;
t184 = sin(pkin(6));
t193 = cos(qJ(2));
t225 = t184 * t193;
t185 = cos(pkin(6));
t188 = sin(qJ(3));
t192 = cos(qJ(3));
t189 = sin(qJ(2));
t226 = t184 * t189;
t127 = t185 * t192 - t188 * t226;
t128 = t185 * t188 + t192 * t226;
t187 = sin(qJ(4));
t191 = cos(qJ(4));
t86 = t127 * t187 + t128 * t191;
t65 = t186 * t86 + t190 * t225;
t66 = -t186 * t225 + t190 * t86;
t85 = -t191 * t127 + t128 * t187;
t22 = Ifges(7,1) * t66 + Ifges(7,4) * t85 + Ifges(7,5) * t65;
t23 = Ifges(6,1) * t66 - Ifges(6,4) * t65 + Ifges(6,5) * t85;
t271 = t22 + t23;
t139 = t187 * t188 - t191 * t192;
t140 = t187 * t192 + t188 * t191;
t237 = Ifges(7,5) * t186;
t73 = Ifges(7,4) * t139 + (Ifges(7,1) * t190 + t237) * t140;
t239 = Ifges(6,4) * t186;
t74 = Ifges(6,5) * t139 + (Ifges(6,1) * t190 - t239) * t140;
t270 = t73 + t74;
t268 = t186 ^ 2 + t190 ^ 2;
t257 = -pkin(10) - pkin(9);
t156 = t257 * t192;
t217 = t257 * t188;
t111 = -t191 * t156 + t187 * t217;
t169 = -pkin(3) * t192 - pkin(2);
t93 = pkin(4) * t139 - pkin(11) * t140 + t169;
t48 = -t111 * t186 + t190 * t93;
t49 = t190 * t111 + t186 * t93;
t207 = -t186 * t48 + t190 * t49;
t40 = qJ(6) * t139 + t49;
t41 = -pkin(5) * t139 - t48;
t208 = t186 * t41 + t190 * t40;
t130 = t185 * t189 * pkin(1) + pkin(8) * t225;
t118 = pkin(9) * t185 + t130;
t119 = (-pkin(2) * t193 - pkin(9) * t189 - pkin(1)) * t184;
t77 = -t118 * t188 + t192 * t119;
t44 = -pkin(3) * t225 - pkin(10) * t128 + t77;
t78 = t192 * t118 + t188 * t119;
t50 = pkin(10) * t127 + t78;
t26 = t187 * t44 + t191 * t50;
t17 = -pkin(11) * t225 + t26;
t161 = pkin(8) * t226;
t247 = pkin(1) * t193;
t117 = t161 + (-pkin(2) - t247) * t185;
t88 = -pkin(3) * t127 + t117;
t28 = pkin(4) * t85 - pkin(11) * t86 + t88;
t6 = -t17 * t186 + t190 * t28;
t7 = t190 * t17 + t186 * t28;
t212 = -t186 * t6 + t190 * t7;
t148 = -Ifges(7,3) * t190 + t237;
t151 = Ifges(6,2) * t190 + t239;
t267 = -t151 / 0.2e1 + t148 / 0.2e1;
t149 = Ifges(6,5) * t186 + Ifges(6,6) * t190;
t150 = Ifges(7,4) * t186 - Ifges(7,6) * t190;
t266 = t149 / 0.2e1 + t150 / 0.2e1;
t19 = Ifges(6,5) * t66 - Ifges(6,6) * t65 + Ifges(6,3) * t85;
t20 = Ifges(7,4) * t66 + Ifges(7,2) * t85 + Ifges(7,6) * t65;
t265 = t19 + t20;
t227 = t140 * t190;
t228 = t140 * t186;
t70 = Ifges(6,5) * t227 - Ifges(6,6) * t228 + Ifges(6,3) * t139;
t71 = Ifges(7,4) * t227 + Ifges(7,2) * t139 + Ifges(7,6) * t228;
t264 = t70 + t71;
t263 = -Ifges(4,5) * t128 - Ifges(4,6) * t127;
t262 = (mrSges(6,3) + mrSges(7,2)) * t268;
t109 = -t156 * t187 - t191 * t217;
t261 = t109 ^ 2;
t260 = 0.2e1 * t109;
t145 = -mrSges(7,1) * t190 - mrSges(7,3) * t186;
t259 = 0.2e1 * t145;
t252 = t186 / 0.2e1;
t251 = -t190 / 0.2e1;
t250 = t190 / 0.2e1;
t248 = m(7) * t186;
t246 = pkin(3) * t191;
t3 = qJ(6) * t85 + t7;
t243 = t190 * t3;
t241 = Ifges(4,3) + Ifges(5,3);
t240 = -Ifges(5,5) * t86 + Ifges(5,6) * t85;
t129 = t185 * t247 - t161;
t235 = t129 * mrSges(3,1);
t234 = t130 * mrSges(3,2);
t172 = t186 * mrSges(7,2);
t229 = qJ(6) * t190;
t224 = Ifges(5,5) * t140 - Ifges(5,6) * t139;
t167 = pkin(3) * t187 + pkin(11);
t223 = t268 * pkin(11) * t167;
t222 = t268 * t167 ^ 2;
t221 = Ifges(4,5) * t188 + Ifges(4,6) * t192;
t220 = t268 * pkin(11) ^ 2;
t218 = t188 ^ 2 + t192 ^ 2;
t216 = Ifges(3,5) * t226 + Ifges(3,6) * t225 + Ifges(3,3) * t185;
t34 = -t85 * mrSges(7,1) + t66 * mrSges(7,2);
t25 = -t187 * t50 + t191 * t44;
t96 = -t139 * mrSges(7,1) + mrSges(7,2) * t227;
t16 = pkin(4) * t225 - t25;
t4 = -pkin(5) * t85 - t6;
t213 = t186 * t4 + t243;
t211 = t186 * mrSges(6,1) + t190 * mrSges(6,2);
t210 = t186 * mrSges(7,1) - t190 * mrSges(7,3);
t209 = pkin(5) * t186 - t229;
t142 = -pkin(5) * t190 - qJ(6) * t186 - pkin(4);
t206 = (mrSges(5,1) * t191 - mrSges(5,2) * t187) * pkin(3);
t205 = Ifges(5,3) + (-t148 + t151) * t190 + t269 * t186;
t31 = -mrSges(6,2) * t85 - mrSges(6,3) * t65;
t32 = -mrSges(7,2) * t65 + mrSges(7,3) * t85;
t33 = mrSges(6,1) * t85 - mrSges(6,3) * t66;
t204 = (t31 + t32) * t190 + (-t33 + t34) * t186;
t94 = -mrSges(6,2) * t139 - mrSges(6,3) * t228;
t95 = mrSges(6,1) * t139 - mrSges(6,3) * t227;
t97 = -mrSges(7,2) * t228 + mrSges(7,3) * t139;
t203 = (t94 + t97) * t190 + (-t95 + t96) * t186;
t201 = mrSges(7,2) * t229 - pkin(5) * t172 + t149 + t150;
t200 = 0.2e1 * t262;
t199 = -m(7) * t209 - t210 - t211;
t146 = -mrSges(6,1) * t190 + mrSges(6,2) * t186;
t18 = Ifges(7,5) * t66 + Ifges(7,6) * t85 + Ifges(7,3) * t65;
t21 = Ifges(6,4) * t66 - Ifges(6,2) * t65 + Ifges(6,6) * t85;
t9 = pkin(5) * t65 - qJ(6) * t66 + t16;
t198 = t25 * mrSges(5,1) - t26 * mrSges(5,2) + mrSges(7,2) * t243 + t212 * mrSges(6,3) + t9 * t145 + t16 * t146 + t4 * t172 + t18 * t251 + t21 * t250 + t271 * t252 + t266 * t85 + t267 * t65 + t66 * t272 - t240;
t55 = t209 * t140 + t109;
t69 = Ifges(7,6) * t139 + (Ifges(7,3) * t186 + t236) * t140;
t72 = Ifges(6,6) * t139 + (-Ifges(6,2) * t186 + t238) * t140;
t197 = -t111 * mrSges(5,2) + t55 * t145 + t72 * t250 + t69 * t251 + t224 + t270 * t252 + t267 * t228 + t227 * t272 + t266 * t139 + (-mrSges(5,1) + t146) * t109 + t207 * mrSges(6,3) + t208 * mrSges(7,2);
t168 = -pkin(4) - t246;
t155 = Ifges(4,1) * t188 + Ifges(4,4) * t192;
t152 = Ifges(4,4) * t188 + Ifges(4,2) * t192;
t147 = -mrSges(4,1) * t192 + mrSges(4,2) * t188;
t131 = t142 - t246;
t108 = -mrSges(4,1) * t225 - mrSges(4,3) * t128;
t107 = mrSges(4,2) * t225 + mrSges(4,3) * t127;
t100 = Ifges(5,1) * t140 - Ifges(5,4) * t139;
t99 = Ifges(5,4) * t140 - Ifges(5,2) * t139;
t98 = mrSges(5,1) * t139 + mrSges(5,2) * t140;
t92 = t211 * t140;
t91 = t210 * t140;
t89 = -mrSges(4,1) * t127 + mrSges(4,2) * t128;
t80 = Ifges(4,1) * t128 + Ifges(4,4) * t127 - Ifges(4,5) * t225;
t79 = Ifges(4,4) * t128 + Ifges(4,2) * t127 - Ifges(4,6) * t225;
t76 = -mrSges(5,1) * t225 - mrSges(5,3) * t86;
t75 = mrSges(5,2) * t225 - mrSges(5,3) * t85;
t39 = mrSges(5,1) * t85 + mrSges(5,2) * t86;
t36 = Ifges(5,1) * t86 - Ifges(5,4) * t85 - Ifges(5,5) * t225;
t35 = Ifges(5,4) * t86 - Ifges(5,2) * t85 - Ifges(5,6) * t225;
t30 = mrSges(6,1) * t65 + mrSges(6,2) * t66;
t29 = mrSges(7,1) * t65 - mrSges(7,3) * t66;
t1 = [(t216 - 0.2e1 * t234 + 0.2e1 * t235) * t185 + t271 * t66 + 0.2e1 * t117 * t89 + 0.2e1 * t78 * t107 + 0.2e1 * t77 * t108 + t86 * t36 + 0.2e1 * t88 * t39 + 0.2e1 * t26 * t75 + 0.2e1 * t25 * t76 + 0.2e1 * t7 * t31 + 0.2e1 * t3 * t32 + 0.2e1 * t6 * t33 + 0.2e1 * t4 * t34 + 0.2e1 * t9 * t29 + 0.2e1 * t16 * t30 + (-t35 + t265) * t85 + ((-0.2e1 * t129 * mrSges(3,3) + Ifges(3,5) * t185 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t189) * t184) * t189 + (0.2e1 * t130 * mrSges(3,3) + Ifges(3,6) * t185 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t189 + (Ifges(3,2) + t241) * t193) * t184 + t240 + t263) * t193) * t184 + m(7) * (t3 ^ 2 + t4 ^ 2 + t9 ^ 2) + m(6) * (t16 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t25 ^ 2 + t26 ^ 2 + t88 ^ 2) + m(4) * (t117 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(3) * (pkin(1) ^ 2 * t184 ^ 2 + t129 ^ 2 + t130 ^ 2) + (t18 - t21) * t65 + t127 * t79 + t128 * t80 + Ifges(2,3); -t234 + t235 + t216 + t111 * t75 + t7 * t94 + t6 * t95 + t4 * t96 + t3 * t97 + t88 * t98 + t86 * t100 / 0.2e1 - pkin(2) * t89 + t9 * t91 + t16 * t92 + t48 * t33 + t49 * t31 + t55 * t29 + t40 * t32 + t41 * t34 + (t19 / 0.2e1 + t20 / 0.2e1 - t35 / 0.2e1 - t26 * mrSges(5,3)) * t139 + (t78 * mrSges(4,3) + pkin(9) * t107 + t79 / 0.2e1) * t192 + (-t77 * mrSges(4,3) - pkin(9) * t108 + t80 / 0.2e1) * t188 + (t30 - t76) * t109 - (t221 + t224) * t225 / 0.2e1 + m(4) * (-pkin(2) * t117 + (-t77 * t188 + t78 * t192) * pkin(9)) + t117 * t147 + t127 * t152 / 0.2e1 + t128 * t155 / 0.2e1 + t169 * t39 + (t36 / 0.2e1 - t25 * mrSges(5,3) + (t22 / 0.2e1 + t23 / 0.2e1) * t190 + (t18 / 0.2e1 - t21 / 0.2e1) * t186) * t140 + (t69 / 0.2e1 - t72 / 0.2e1) * t65 + (t73 / 0.2e1 + t74 / 0.2e1) * t66 + (-t99 / 0.2e1 + t70 / 0.2e1 + t71 / 0.2e1) * t85 + m(7) * (t3 * t40 + t4 * t41 + t55 * t9) + m(6) * (t109 * t16 + t48 * t6 + t49 * t7) + m(5) * (-t109 * t25 + t111 * t26 + t169 * t88); -0.2e1 * pkin(2) * t147 + t92 * t260 + t192 * t152 + t188 * t155 + 0.2e1 * t169 * t98 + 0.2e1 * t40 * t97 + 0.2e1 * t41 * t96 + 0.2e1 * t48 * t95 + 0.2e1 * t49 * t94 + 0.2e1 * t55 * t91 + Ifges(3,3) + 0.2e1 * t218 * pkin(9) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t111 + t264 - t99) * t139 + m(7) * (t40 ^ 2 + t41 ^ 2 + t55 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2 + t261) + m(5) * (t111 ^ 2 + t169 ^ 2 + t261) + m(4) * (t218 * pkin(9) ^ 2 + pkin(2) ^ 2) + (mrSges(5,3) * t260 + t100 + t270 * t190 + (t69 - t72) * t186) * t140; t198 + t77 * mrSges(4,1) - t78 * mrSges(4,2) + (t187 * t75 + t191 * t76 + m(5) * (t187 * t26 + t191 * t25)) * pkin(3) + t204 * t167 + m(6) * (t16 * t168 + t212 * t167) + m(7) * (t131 * t9 + t213 * t167) - t241 * t225 + t131 * t29 + t168 * t30 - t263; t197 + (m(5) * (-t109 * t191 + t111 * t187) + (-t139 * t187 - t140 * t191) * mrSges(5,3)) * pkin(3) + m(7) * (t131 * t55 + t208 * t167) + m(6) * (t109 * t168 + t207 * t167) + t203 * t167 + (-mrSges(4,1) * t188 - mrSges(4,2) * t192) * pkin(9) + t131 * t91 + t168 * t92 + t221; t131 * t259 + 0.2e1 * t168 * t146 + Ifges(4,3) + 0.2e1 * t206 + m(7) * (t131 ^ 2 + t222) + m(6) * (t168 ^ 2 + t222) + m(5) * (t187 ^ 2 + t191 ^ 2) * pkin(3) ^ 2 + t200 * t167 + t205; t198 - pkin(4) * t30 + m(7) * (t213 * pkin(11) + t142 * t9) + m(6) * (-pkin(4) * t16 + t212 * pkin(11)) + t204 * pkin(11) + t142 * t29 - Ifges(5,3) * t225; t197 - pkin(4) * t92 + t203 * pkin(11) + m(7) * (t208 * pkin(11) + t142 * t55) + m(6) * (-pkin(4) * t109 + t207 * pkin(11)) + t142 * t91; (t168 - pkin(4)) * t146 + (t131 + t142) * t145 + t206 + m(7) * (t131 * t142 + t223) + m(6) * (-pkin(4) * t168 + t223) + t205 + (pkin(11) + t167) * t262; -0.2e1 * pkin(4) * t146 + t142 * t259 + m(7) * (t142 ^ 2 + t220) + m(6) * (pkin(4) ^ 2 + t220) + t200 * pkin(11) + t205; -pkin(5) * t34 + m(7) * (-pkin(5) * t4 + qJ(6) * t3) + qJ(6) * t32 + t3 * mrSges(7,3) - t7 * mrSges(6,2) + t6 * mrSges(6,1) - t4 * mrSges(7,1) + t265; -t49 * mrSges(6,2) - pkin(5) * t96 + qJ(6) * t97 + m(7) * (-pkin(5) * t41 + qJ(6) * t40) + t40 * mrSges(7,3) + t48 * mrSges(6,1) - t41 * mrSges(7,1) + t264; t199 * t167 + t201; t199 * pkin(11) + t201; Ifges(7,2) + Ifges(6,3) + 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2); m(7) * t4 + t34; m(7) * t41 + t96; t167 * t248 + t172; pkin(11) * t248 + t172; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
