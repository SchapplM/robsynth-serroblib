% Calculate joint inertia matrix for
% S6RRRRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-10 00:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR15_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:34:18
% EndTime: 2019-03-10 00:34:26
% DurationCPUTime: 3.01s
% Computational Cost: add. (5286->521), mult. (13572->725), div. (0->0), fcn. (15082->12), ass. (0->186)
t184 = sin(qJ(4));
t188 = cos(qJ(4));
t247 = t184 ^ 2 + t188 ^ 2;
t183 = sin(qJ(6));
t235 = -t183 / 0.2e1;
t179 = sin(pkin(7));
t181 = cos(pkin(7));
t182 = cos(pkin(6));
t180 = sin(pkin(6));
t190 = cos(qJ(2));
t222 = t180 * t190;
t123 = -t179 * t222 + t181 * t182;
t185 = sin(qJ(3));
t186 = sin(qJ(2));
t189 = cos(qJ(3));
t220 = t181 * t190;
t225 = t179 * t185;
t91 = t182 * t225 + (t185 * t220 + t186 * t189) * t180;
t61 = -t123 * t188 + t184 * t91;
t62 = t123 * t184 + t188 * t91;
t211 = t180 * t220;
t223 = t180 * t186;
t224 = t179 * t189;
t90 = -t182 * t224 + t185 * t223 - t189 * t211;
t20 = Ifges(5,5) * t62 - Ifges(5,6) * t61 + Ifges(5,3) * t90;
t23 = Ifges(6,1) * t90 - Ifges(6,4) * t62 + Ifges(6,5) * t61;
t246 = t20 + t23;
t245 = -m(6) * pkin(4) + mrSges(6,2);
t106 = (-pkin(10) * t179 * t186 - pkin(2) * t190 - pkin(1)) * t180;
t233 = pkin(1) * t182;
t130 = pkin(9) * t222 + t186 * t233;
t82 = (t179 * t182 + t211) * pkin(10) + t130;
t164 = t190 * t233;
t94 = pkin(2) * t182 + t164 + (-pkin(10) * t181 - pkin(9)) * t223;
t33 = -t185 * t82 + (t106 * t179 + t181 * t94) * t189;
t187 = cos(qJ(6));
t37 = -t183 * t90 + t187 * t61;
t38 = t183 * t61 + t187 * t90;
t12 = Ifges(7,1) * t38 + Ifges(7,4) * t37 + Ifges(7,5) * t62;
t244 = t12 / 0.2e1;
t125 = t181 * t184 + t188 * t225;
t124 = -t188 * t181 + t184 * t225;
t95 = t124 * t187 + t183 * t224;
t96 = t124 * t183 - t187 * t224;
t48 = Ifges(7,1) * t96 + Ifges(7,4) * t95 + Ifges(7,5) * t125;
t243 = t48 / 0.2e1;
t242 = pkin(4) + pkin(12);
t241 = pkin(5) + pkin(11);
t229 = Ifges(7,4) * t183;
t115 = Ifges(7,6) * t184 + (-Ifges(7,2) * t187 - t229) * t188;
t240 = t115 / 0.2e1;
t228 = Ifges(7,4) * t187;
t116 = Ifges(7,5) * t184 + (-Ifges(7,1) * t183 - t228) * t188;
t239 = t116 / 0.2e1;
t170 = Ifges(7,5) * t187;
t238 = Ifges(7,6) * t235 + t170 / 0.2e1;
t146 = -Ifges(7,2) * t183 + t228;
t237 = t146 / 0.2e1;
t149 = Ifges(7,1) * t187 - t229;
t236 = t149 / 0.2e1;
t234 = -t187 / 0.2e1;
t232 = pkin(2) * t189;
t231 = -Ifges(6,1) - Ifges(5,3);
t54 = t181 * t106 - t179 * t94;
t26 = pkin(3) * t90 - pkin(11) * t91 + t54;
t221 = t181 * t185;
t34 = t106 * t225 + t189 * t82 + t94 * t221;
t30 = pkin(11) * t123 + t34;
t9 = t184 * t26 + t188 * t30;
t230 = mrSges(7,3) * t188;
t128 = -pkin(9) * t223 + t164;
t227 = t128 * mrSges(3,1);
t226 = t130 * mrSges(3,2);
t136 = -mrSges(7,2) * t184 - t187 * t230;
t219 = t183 * t136;
t218 = t187 * t242;
t129 = pkin(2) * t221 + pkin(10) * t224;
t112 = pkin(11) * t181 + t129;
t113 = (-pkin(3) * t189 - pkin(11) * t185 - pkin(2)) * t179;
t70 = t188 * t112 + t184 * t113;
t217 = -Ifges(6,4) * t125 + Ifges(6,5) * t124;
t216 = Ifges(5,5) * t125 - Ifges(5,6) * t124;
t145 = Ifges(5,5) * t184 + Ifges(5,6) * t188;
t215 = t247 * pkin(11) ^ 2;
t214 = -t183 ^ 2 - t187 ^ 2;
t10 = Ifges(7,5) * t38 + Ifges(7,6) * t37 + Ifges(7,3) * t62;
t43 = Ifges(4,5) * t91 - Ifges(4,6) * t90 + Ifges(4,3) * t123;
t46 = Ifges(7,5) * t96 + Ifges(7,6) * t95 + Ifges(7,3) * t125;
t19 = Ifges(6,5) * t90 - Ifges(6,6) * t62 + Ifges(6,3) * t61;
t22 = Ifges(5,4) * t62 - Ifges(5,2) * t61 + Ifges(5,6) * t90;
t213 = t19 / 0.2e1 - t22 / 0.2e1;
t71 = -Ifges(6,5) * t224 - Ifges(6,6) * t125 + Ifges(6,3) * t124;
t75 = Ifges(5,4) * t125 - Ifges(5,2) * t124 - Ifges(5,6) * t224;
t212 = t71 / 0.2e1 - t75 / 0.2e1;
t107 = Ifges(4,5) * t225 + Ifges(4,6) * t224 + Ifges(4,3) * t181;
t210 = Ifges(3,5) * t223 + Ifges(3,6) * t222 + Ifges(3,3) * t182;
t142 = -Ifges(6,6) * t184 - Ifges(6,3) * t188;
t148 = Ifges(5,4) * t184 + Ifges(5,2) * t188;
t209 = t142 / 0.2e1 - t148 / 0.2e1;
t208 = t145 / 0.2e1 - Ifges(6,4) * t184 / 0.2e1 - Ifges(6,5) * t188 / 0.2e1;
t40 = t62 * mrSges(6,1) + t90 * mrSges(6,2);
t207 = m(7) * t214;
t206 = -qJ(5) * t184 - pkin(3);
t8 = -t184 * t30 + t188 * t26;
t205 = t214 * mrSges(7,3);
t69 = -t184 * t112 + t113 * t188;
t6 = -qJ(5) * t90 - t9;
t21 = Ifges(6,4) * t90 - Ifges(6,2) * t62 + Ifges(6,6) * t61;
t24 = Ifges(5,1) * t62 - Ifges(5,4) * t61 + Ifges(5,5) * t90;
t203 = t10 / 0.2e1 + t24 / 0.2e1 - t21 / 0.2e1;
t72 = -Ifges(6,4) * t224 - Ifges(6,2) * t125 + Ifges(6,6) * t124;
t76 = Ifges(5,1) * t125 - Ifges(5,4) * t124 - Ifges(5,5) * t224;
t202 = t46 / 0.2e1 - t72 / 0.2e1 + t76 / 0.2e1;
t3 = pkin(5) * t62 - t242 * t90 - t8;
t29 = -pkin(3) * t123 - t33;
t194 = -qJ(5) * t62 + t29;
t5 = t242 * t61 + t194;
t1 = -t183 * t5 + t187 * t3;
t2 = t183 * t3 + t187 * t5;
t201 = t1 * t187 + t183 * t2;
t64 = pkin(4) * t224 - t69;
t114 = Ifges(7,3) * t184 + (-Ifges(7,5) * t183 - Ifges(7,6) * t187) * t188;
t143 = -Ifges(6,2) * t184 - Ifges(6,6) * t188;
t150 = Ifges(5,1) * t184 + Ifges(5,4) * t188;
t200 = t114 / 0.2e1 - t143 / 0.2e1 + t150 / 0.2e1;
t199 = mrSges(7,1) * t187 - mrSges(7,2) * t183;
t49 = pkin(5) * t125 + pkin(12) * t224 + t64;
t158 = pkin(10) * t225;
t111 = t158 + (-pkin(3) - t232) * t181;
t195 = -qJ(5) * t125 + t111;
t50 = t242 * t124 + t195;
t15 = -t183 * t50 + t187 * t49;
t16 = t183 * t49 + t187 * t50;
t198 = t15 * t187 + t16 * t183;
t132 = -t242 * t188 + t206;
t152 = t241 * t184;
t88 = -t132 * t183 + t152 * t187;
t89 = t132 * t187 + t152 * t183;
t197 = t183 * t89 + t187 * t88;
t63 = qJ(5) * t224 - t70;
t100 = t125 * mrSges(6,1) - mrSges(6,2) * t224;
t192 = qJ(5) ^ 2;
t153 = t241 * t188;
t141 = mrSges(7,1) * t183 + mrSges(7,2) * t187;
t140 = -mrSges(5,1) * t188 + mrSges(5,2) * t184;
t139 = mrSges(6,2) * t188 - mrSges(6,3) * t184;
t138 = -pkin(4) * t188 + t206;
t135 = mrSges(7,1) * t184 + t183 * t230;
t134 = -mrSges(4,2) * t181 + mrSges(4,3) * t224;
t133 = mrSges(4,1) * t181 - mrSges(4,3) * t225;
t131 = t199 * t188;
t127 = t181 * t232 - t158;
t126 = (-mrSges(4,1) * t189 + mrSges(4,2) * t185) * t179;
t109 = Ifges(4,5) * t181 + (Ifges(4,1) * t185 + Ifges(4,4) * t189) * t179;
t108 = Ifges(4,6) * t181 + (Ifges(4,4) * t185 + Ifges(4,2) * t189) * t179;
t102 = -mrSges(5,1) * t224 - mrSges(5,3) * t125;
t101 = mrSges(5,2) * t224 - mrSges(5,3) * t124;
t99 = mrSges(6,1) * t124 + mrSges(6,3) * t224;
t78 = -mrSges(6,2) * t124 - mrSges(6,3) * t125;
t77 = mrSges(5,1) * t124 + mrSges(5,2) * t125;
t74 = -Ifges(5,3) * t224 + t216;
t73 = -Ifges(6,1) * t224 + t217;
t68 = mrSges(7,1) * t125 - mrSges(7,3) * t96;
t67 = -mrSges(7,2) * t125 + mrSges(7,3) * t95;
t66 = mrSges(4,1) * t123 - mrSges(4,3) * t91;
t65 = -mrSges(4,2) * t123 - mrSges(4,3) * t90;
t60 = pkin(4) * t124 + t195;
t53 = -mrSges(7,1) * t95 + mrSges(7,2) * t96;
t52 = mrSges(4,1) * t90 + mrSges(4,2) * t91;
t51 = -pkin(5) * t124 - t63;
t47 = Ifges(7,4) * t96 + Ifges(7,2) * t95 + Ifges(7,6) * t125;
t45 = Ifges(4,1) * t91 - Ifges(4,4) * t90 + Ifges(4,5) * t123;
t44 = Ifges(4,4) * t91 - Ifges(4,2) * t90 + Ifges(4,6) * t123;
t42 = mrSges(5,1) * t90 - mrSges(5,3) * t62;
t41 = -mrSges(5,2) * t90 - mrSges(5,3) * t61;
t39 = mrSges(6,1) * t61 - mrSges(6,3) * t90;
t32 = -mrSges(6,2) * t61 - mrSges(6,3) * t62;
t31 = mrSges(5,1) * t61 + mrSges(5,2) * t62;
t18 = mrSges(7,1) * t62 - mrSges(7,3) * t38;
t17 = -mrSges(7,2) * t62 + mrSges(7,3) * t37;
t14 = -mrSges(7,1) * t37 + mrSges(7,2) * t38;
t13 = pkin(4) * t61 + t194;
t11 = Ifges(7,4) * t38 + Ifges(7,2) * t37 + Ifges(7,6) * t62;
t7 = -pkin(4) * t90 - t8;
t4 = -pkin(5) * t61 - t6;
t25 = [m(5) * (t29 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(6) * (t13 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2 + t54 ^ 2) + (t10 + t24 - t21) * t62 + (t19 - t22) * t61 + ((Ifges(3,5) * t186 + Ifges(3,6) * t190) * t182 + 0.2e1 * (-t128 * t186 + t130 * t190) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t190 + mrSges(3,2) * t186) + t186 * (Ifges(3,1) * t186 + Ifges(3,4) * t190) + t190 * (Ifges(3,4) * t186 + Ifges(3,2) * t190)) * t180) * t180 + (t210 - 0.2e1 * t226 + 0.2e1 * t227) * t182 + (-t44 + t246) * t90 + Ifges(2,3) + 0.2e1 * t4 * t14 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t29 * t31 + 0.2e1 * t13 * t32 + t37 * t11 + t38 * t12 + 0.2e1 * t6 * t39 + 0.2e1 * t7 * t40 + 0.2e1 * t9 * t41 + 0.2e1 * t8 * t42 + 0.2e1 * t54 * t52 + 0.2e1 * t34 * t65 + 0.2e1 * t33 * t66 + t91 * t45 + t123 * t43 + m(3) * (t128 ^ 2 + t130 ^ 2); t202 * t62 + t203 * t125 + t210 + (-pkin(2) * t52 + t185 * t45 / 0.2e1 + (t44 / 0.2e1 - t20 / 0.2e1 - t23 / 0.2e1) * t189) * t179 + m(7) * (t1 * t15 + t16 * t2 + t4 * t51) + m(6) * (t13 * t60 + t6 * t63 + t64 * t7) + m(5) * (t111 * t29 + t69 * t8 + t70 * t9) + m(4) * (-pkin(2) * t179 * t54 + t127 * t33 + t129 * t34) + (t73 / 0.2e1 + t74 / 0.2e1 - t108 / 0.2e1) * t90 + t212 * t61 + t213 * t124 + t227 - t226 + t38 * t243 + t96 * t244 + t16 * t17 + t15 * t18 + t37 * t47 / 0.2e1 + t51 * t14 + t4 * t53 + t60 * t32 + t63 * t39 + t64 * t40 + t2 * t67 + t1 * t68 + t69 * t42 + t70 * t41 + t29 * t77 + t13 * t78 + t95 * t11 / 0.2e1 + t6 * t99 + t7 * t100 + t9 * t101 + t8 * t102 + t91 * t109 / 0.2e1 + t111 * t31 + t123 * t107 / 0.2e1 + t54 * t126 + t127 * t66 + t129 * t65 + t33 * t133 + t34 * t134 + t181 * t43 / 0.2e1; 0.2e1 * t64 * t100 + 0.2e1 * t70 * t101 + 0.2e1 * t69 * t102 + t181 * t107 + 0.2e1 * t111 * t77 + 0.2e1 * t127 * t133 + 0.2e1 * t129 * t134 + 0.2e1 * t15 * t68 + 0.2e1 * t16 * t67 + t95 * t47 + t96 * t48 + 0.2e1 * t51 * t53 + 0.2e1 * t60 * t78 + 0.2e1 * t63 * t99 + Ifges(3,3) + (t71 - t75) * t124 + (t46 + t76 - t72) * t125 + (-0.2e1 * pkin(2) * t126 + t185 * t109 + (t108 - t73 - t74) * t189) * t179 + m(7) * (t15 ^ 2 + t16 ^ 2 + t51 ^ 2) + m(6) * (t60 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t111 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(4) * (pkin(2) ^ 2 * t179 ^ 2 + t127 ^ 2 + t129 ^ 2); t200 * t62 + (t7 * mrSges(6,1) - t8 * mrSges(5,3) + (t40 - t42) * pkin(11) + t203) * t184 + t43 + (-t6 * mrSges(6,1) + t9 * mrSges(5,3) + t12 * t235 + t11 * t234 + (-t39 + t41) * pkin(11) - t213) * t188 + t208 * t90 + t209 * t61 + m(6) * (t13 * t138 + (t184 * t7 - t188 * t6) * pkin(11)) + m(5) * (-pkin(3) * t29 + (-t184 * t8 + t188 * t9) * pkin(11)) + m(7) * (t1 * t88 + t153 * t4 + t2 * t89) - pkin(3) * t31 + t33 * mrSges(4,1) - t34 * mrSges(4,2) + t88 * t18 + t89 * t17 + t37 * t240 + t38 * t239 + t4 * t131 + t1 * t135 + t2 * t136 + t138 * t32 + t13 * t139 + t29 * t140 + t153 * t14; m(7) * (t15 * t88 + t153 * t51 + t16 * t89) + t200 * t125 + (t64 * mrSges(6,1) - t69 * mrSges(5,3) + (t100 - t102) * pkin(11) + t202) * t184 - t208 * t224 + t107 + (-t63 * mrSges(6,1) + t70 * mrSges(5,3) + t48 * t235 + t47 * t234 + (-t99 + t101) * pkin(11) - t212) * t188 + t209 * t124 + m(6) * (t138 * t60 + (t184 * t64 - t188 * t63) * pkin(11)) + m(5) * (-pkin(3) * t111 + (-t184 * t69 + t188 * t70) * pkin(11)) - pkin(3) * t77 + t88 * t68 + t89 * t67 + t95 * t240 + t96 * t239 + t127 * mrSges(4,1) - t129 * mrSges(4,2) + t51 * t131 + t15 * t135 + t16 * t136 + t138 * t78 + t60 * t139 + t111 * t140 + t153 * t53; -0.2e1 * pkin(3) * t140 + 0.2e1 * t153 * t131 + 0.2e1 * t88 * t135 + 0.2e1 * t89 * t136 + 0.2e1 * t138 * t139 + Ifges(4,3) + (t114 - t143 + t150) * t184 + m(7) * (t153 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(6) * (t138 ^ 2 + t215) + m(5) * (pkin(3) ^ 2 + t215) + (-t115 * t187 - t116 * t183 - t142 + t148) * t188 + 0.2e1 * (mrSges(6,1) + mrSges(5,3)) * pkin(11) * t247; -t6 * mrSges(6,3) + t7 * mrSges(6,2) + t8 * mrSges(5,1) - t9 * mrSges(5,2) - pkin(4) * t40 + t4 * t141 + t62 * t238 + t37 * t237 + t38 * t236 + (t14 - t39) * qJ(5) + (-t1 * mrSges(7,3) - t18 * t242 + t244) * t187 + (-t2 * mrSges(7,3) - t242 * t17 - t11 / 0.2e1) * t183 + m(6) * (-pkin(4) * t7 - qJ(5) * t6) + m(7) * (qJ(5) * t4 - t201 * t242) + t246; -t63 * mrSges(6,3) + t64 * mrSges(6,2) + t69 * mrSges(5,1) - t70 * mrSges(5,2) - pkin(4) * t100 + t51 * t141 + t125 * t238 + t95 * t237 + t96 * t236 + t231 * t224 + (t53 - t99) * qJ(5) + (-t15 * mrSges(7,3) - t242 * t68 + t243) * t187 + (-t16 * mrSges(7,3) - t242 * t67 - t47 / 0.2e1) * t183 + m(6) * (-pkin(4) * t64 - qJ(5) * t63) + m(7) * (qJ(5) * t51 - t198 * t242) + t216 + t217; -t242 * t219 - t135 * t218 + m(7) * (qJ(5) * t153 - t197 * t242) + t153 * t141 + t187 * t239 + t115 * t235 + qJ(5) * t131 - t197 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) - Ifges(6,4) + t238) * t184 + (qJ(5) * mrSges(6,1) + t146 * t234 + t149 * t235 - Ifges(6,5)) * t188 + ((m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t188 + (-mrSges(5,1) + t245) * t184) * pkin(11) + t145; -0.2e1 * pkin(4) * mrSges(6,2) - t183 * t146 + t187 * t149 + m(7) * (-t214 * t242 ^ 2 + t192) + m(6) * (pkin(4) ^ 2 + t192) - t231 + 0.2e1 * (t141 + mrSges(6,3)) * qJ(5) - 0.2e1 * t242 * t205; m(6) * t7 + m(7) * t201 + t183 * t17 + t187 * t18 + t40; m(6) * t64 + m(7) * t198 + t183 * t67 + t187 * t68 + t100; m(7) * t197 + t219 + t187 * t135 + (m(6) * pkin(11) + mrSges(6,1)) * t184; t207 * t242 + t205 + t245; m(6) - t207; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t10; mrSges(7,1) * t15 - mrSges(7,2) * t16 + t46; mrSges(7,1) * t88 - mrSges(7,2) * t89 + t114; -mrSges(7,1) * t218 + t170 + (mrSges(7,2) * t242 - Ifges(7,6)) * t183; t199; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
