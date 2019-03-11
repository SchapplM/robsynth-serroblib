% Calculate joint inertia matrix for
% S6RRRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:41
% EndTime: 2019-03-09 18:56:49
% DurationCPUTime: 2.75s
% Computational Cost: add. (7895->496), mult. (21016->726), div. (0->0), fcn. (24018->14), ass. (0->191)
t189 = sin(pkin(13));
t176 = pkin(3) * t189 + pkin(11);
t256 = 0.2e1 * t176;
t195 = sin(qJ(6));
t199 = cos(qJ(6));
t190 = sin(pkin(7));
t192 = cos(pkin(13));
t197 = sin(qJ(3));
t201 = cos(qJ(3));
t134 = (t189 * t201 + t192 * t197) * t190;
t193 = cos(pkin(7));
t196 = sin(qJ(5));
t200 = cos(qJ(5));
t113 = t134 * t196 - t200 * t193;
t114 = t134 * t200 + t193 * t196;
t224 = t190 * t201;
t225 = t190 * t197;
t133 = t189 * t225 - t192 * t224;
t86 = -t114 * t195 + t133 * t199;
t61 = -mrSges(7,2) * t113 + mrSges(7,3) * t86;
t87 = t114 * t199 + t133 * t195;
t62 = mrSges(7,1) * t113 - mrSges(7,3) * t87;
t255 = -t195 * t62 + t199 * t61;
t194 = cos(pkin(6));
t191 = sin(pkin(6));
t202 = cos(qJ(2));
t222 = t191 * t202;
t140 = -t190 * t222 + t193 * t194;
t198 = sin(qJ(2));
t219 = t193 * t202;
t107 = t194 * t224 + (-t197 * t198 + t201 * t219) * t191;
t108 = t194 * t225 + (t197 * t219 + t198 * t201) * t191;
t74 = t107 * t189 + t108 * t192;
t57 = t140 * t196 + t200 * t74;
t73 = -t192 * t107 + t108 * t189;
t29 = -t195 * t57 + t199 * t73;
t56 = -t140 * t200 + t196 * t74;
t19 = -mrSges(7,2) * t56 + mrSges(7,3) * t29;
t30 = t195 * t73 + t199 * t57;
t20 = mrSges(7,1) * t56 - mrSges(7,3) * t30;
t254 = t199 * t19 - t195 * t20;
t23 = Ifges(6,1) * t57 - Ifges(6,4) * t56 + Ifges(6,5) * t73;
t253 = t23 / 0.2e1;
t65 = Ifges(6,1) * t114 - Ifges(6,4) * t113 + Ifges(6,5) * t133;
t252 = t65 / 0.2e1;
t235 = Ifges(7,4) * t199;
t136 = -Ifges(7,6) * t200 + (-Ifges(7,2) * t195 + t235) * t196;
t251 = t136 / 0.2e1;
t236 = Ifges(7,4) * t195;
t137 = -Ifges(7,5) * t200 + (Ifges(7,1) * t199 - t236) * t196;
t250 = t137 / 0.2e1;
t155 = Ifges(7,5) * t195 + Ifges(7,6) * t199;
t249 = t155 / 0.2e1;
t156 = Ifges(6,5) * t196 + Ifges(6,6) * t200;
t248 = t156 / 0.2e1;
t157 = Ifges(7,2) * t199 + t236;
t247 = t157 / 0.2e1;
t159 = Ifges(7,1) * t195 + t235;
t246 = t159 / 0.2e1;
t160 = Ifges(6,1) * t196 + Ifges(6,4) * t200;
t245 = t160 / 0.2e1;
t244 = -t195 / 0.2e1;
t243 = t195 / 0.2e1;
t242 = t199 / 0.2e1;
t241 = pkin(1) * t194;
t240 = pkin(5) * t200;
t145 = pkin(9) * t222 + t198 * t241;
t101 = (t190 * t194 + t191 * t219) * pkin(10) + t145;
t172 = t202 * t241;
t223 = t191 * t198;
t112 = pkin(2) * t194 + t172 + (-pkin(10) * t193 - pkin(9)) * t223;
t124 = (-pkin(10) * t190 * t198 - pkin(2) * t202 - pkin(1)) * t191;
t220 = t193 * t201;
t51 = -t101 * t197 + t112 * t220 + t124 * t224;
t37 = pkin(3) * t140 - qJ(4) * t108 + t51;
t221 = t193 * t197;
t52 = t201 * t101 + t112 * t221 + t124 * t225;
t46 = qJ(4) * t107 + t52;
t16 = t189 * t37 + t192 * t46;
t14 = pkin(11) * t140 + t16;
t80 = -t112 * t190 + t193 * t124;
t58 = -pkin(3) * t107 + t80;
t25 = pkin(4) * t73 - pkin(11) * t74 + t58;
t6 = t200 * t14 + t196 * t25;
t11 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t32 = mrSges(6,1) * t73 - mrSges(6,3) * t57;
t239 = -t32 + t11;
t33 = Ifges(5,5) * t74 - Ifges(5,6) * t73 + Ifges(5,3) * t140;
t66 = Ifges(4,5) * t108 + Ifges(4,6) * t107 + Ifges(4,3) * t140;
t238 = t33 + t66;
t50 = -mrSges(7,1) * t86 + mrSges(7,2) * t87;
t90 = mrSges(6,1) * t133 - mrSges(6,3) * t114;
t237 = t50 - t90;
t170 = pkin(2) * t220;
t117 = pkin(3) * t193 + t170 + (-pkin(10) - qJ(4)) * t225;
t144 = pkin(2) * t221 + pkin(10) * t224;
t123 = qJ(4) * t224 + t144;
t84 = t189 * t117 + t192 * t123;
t79 = pkin(11) * t193 + t84;
t150 = (-pkin(3) * t201 - pkin(2)) * t190;
t88 = pkin(4) * t133 - pkin(11) * t134 + t150;
t49 = t196 * t88 + t200 * t79;
t143 = -pkin(9) * t223 + t172;
t234 = t143 * mrSges(3,1);
t233 = t145 * mrSges(3,2);
t129 = Ifges(4,5) * t225 + Ifges(4,6) * t224 + Ifges(4,3) * t193;
t93 = Ifges(5,5) * t134 - Ifges(5,6) * t133 + Ifges(5,3) * t193;
t228 = t93 + t129;
t227 = t176 * t196;
t226 = t176 * t200;
t218 = t195 * t196;
t217 = t196 * t199;
t216 = t195 ^ 2 + t199 ^ 2;
t186 = t196 ^ 2;
t188 = t200 ^ 2;
t215 = t186 + t188;
t8 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t56;
t21 = Ifges(6,5) * t57 - Ifges(6,6) * t56 + Ifges(6,3) * t73;
t22 = Ifges(6,4) * t57 - Ifges(6,2) * t56 + Ifges(6,6) * t73;
t214 = t8 / 0.2e1 - t22 / 0.2e1;
t38 = Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t113;
t64 = Ifges(6,4) * t114 - Ifges(6,2) * t113 + Ifges(6,6) * t133;
t213 = t38 / 0.2e1 - t64 / 0.2e1;
t63 = Ifges(6,5) * t114 - Ifges(6,6) * t113 + Ifges(6,3) * t133;
t212 = Ifges(3,5) * t223 + Ifges(3,6) * t222 + Ifges(3,3) * t194;
t177 = -pkin(3) * t192 - pkin(4);
t135 = Ifges(7,5) * t217 - Ifges(7,6) * t218 - Ifges(7,3) * t200;
t158 = Ifges(6,4) * t196 + Ifges(6,2) * t200;
t211 = t135 / 0.2e1 - t158 / 0.2e1;
t43 = t73 * mrSges(5,1) + t74 * mrSges(5,2);
t15 = -t189 * t46 + t192 * t37;
t96 = t133 * mrSges(5,1) + t134 * mrSges(5,2);
t210 = t196 * t216;
t83 = t117 * t192 - t189 * t123;
t4 = pkin(12) * t73 + t6;
t13 = -pkin(4) * t140 - t15;
t7 = pkin(5) * t56 - pkin(12) * t57 + t13;
t1 = -t195 * t4 + t199 * t7;
t2 = t195 * t7 + t199 * t4;
t209 = -t1 * t195 + t199 * t2;
t154 = -t200 * mrSges(6,1) + t196 * mrSges(6,2);
t208 = mrSges(7,1) * t195 + mrSges(7,2) * t199;
t5 = -t14 * t196 + t200 * t25;
t42 = pkin(12) * t133 + t49;
t78 = -pkin(4) * t193 - t83;
t47 = pkin(5) * t113 - pkin(12) * t114 + t78;
t17 = -t195 * t42 + t199 * t47;
t18 = t195 * t47 + t199 * t42;
t207 = -t17 * t195 + t18 * t199;
t48 = -t196 * t79 + t200 * t88;
t147 = -pkin(12) * t196 + t177 - t240;
t110 = t147 * t199 - t195 * t226;
t111 = t147 * t195 + t199 * t226;
t206 = -t110 * t195 + t111 * t199;
t151 = mrSges(7,2) * t200 - mrSges(7,3) * t218;
t152 = -mrSges(7,1) * t200 - mrSges(7,3) * t217;
t205 = t199 * t151 - t195 * t152;
t174 = t176 ^ 2;
t161 = t186 * t174;
t153 = -mrSges(7,1) * t199 + mrSges(7,2) * t195;
t149 = -mrSges(4,2) * t193 + mrSges(4,3) * t224;
t148 = mrSges(4,1) * t193 - mrSges(4,3) * t225;
t146 = t208 * t196;
t142 = -pkin(10) * t225 + t170;
t141 = (-mrSges(4,1) * t201 + mrSges(4,2) * t197) * t190;
t131 = Ifges(4,5) * t193 + (Ifges(4,1) * t197 + Ifges(4,4) * t201) * t190;
t130 = Ifges(4,6) * t193 + (Ifges(4,4) * t197 + Ifges(4,2) * t201) * t190;
t119 = mrSges(5,1) * t193 - mrSges(5,3) * t134;
t118 = -mrSges(5,2) * t193 - mrSges(5,3) * t133;
t95 = Ifges(5,1) * t134 - Ifges(5,4) * t133 + Ifges(5,5) * t193;
t94 = Ifges(5,4) * t134 - Ifges(5,2) * t133 + Ifges(5,6) * t193;
t92 = mrSges(4,1) * t140 - mrSges(4,3) * t108;
t91 = -mrSges(4,2) * t140 + mrSges(4,3) * t107;
t89 = -mrSges(6,2) * t133 - mrSges(6,3) * t113;
t77 = mrSges(6,1) * t113 + mrSges(6,2) * t114;
t75 = -mrSges(4,1) * t107 + mrSges(4,2) * t108;
t68 = Ifges(4,1) * t108 + Ifges(4,4) * t107 + Ifges(4,5) * t140;
t67 = Ifges(4,4) * t108 + Ifges(4,2) * t107 + Ifges(4,6) * t140;
t60 = mrSges(5,1) * t140 - mrSges(5,3) * t74;
t59 = -mrSges(5,2) * t140 - mrSges(5,3) * t73;
t41 = -pkin(5) * t133 - t48;
t40 = Ifges(7,1) * t87 + Ifges(7,4) * t86 + Ifges(7,5) * t113;
t39 = Ifges(7,4) * t87 + Ifges(7,2) * t86 + Ifges(7,6) * t113;
t35 = Ifges(5,1) * t74 - Ifges(5,4) * t73 + Ifges(5,5) * t140;
t34 = Ifges(5,4) * t74 - Ifges(5,2) * t73 + Ifges(5,6) * t140;
t31 = -mrSges(6,2) * t73 - mrSges(6,3) * t56;
t26 = mrSges(6,1) * t56 + mrSges(6,2) * t57;
t10 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t56;
t9 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t56;
t3 = -pkin(5) * t73 - t5;
t12 = [m(3) * (t143 ^ 2 + t145 ^ 2) + m(6) * (t13 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t58 ^ 2) + m(4) * (t51 ^ 2 + t52 ^ 2 + t80 ^ 2) + (t21 - t34) * t73 + (t8 - t22) * t56 + ((Ifges(3,5) * t198 + Ifges(3,6) * t202) * t194 + 0.2e1 * (-t143 * t198 + t145 * t202) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t202 + mrSges(3,2) * t198) + t198 * (Ifges(3,1) * t198 + Ifges(3,4) * t202) + t202 * (Ifges(3,4) * t198 + Ifges(3,2) * t202) + m(3) * pkin(1) ^ 2) * t191) * t191 + Ifges(2,3) + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t19 + 0.2e1 * t1 * t20 + 0.2e1 * t13 * t26 + t29 * t9 + t30 * t10 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + t57 * t23 + 0.2e1 * t58 * t43 + 0.2e1 * t16 * t59 + 0.2e1 * t15 * t60 + t74 * t35 + 0.2e1 * t80 * t75 + 0.2e1 * t52 * t91 + 0.2e1 * t51 * t92 + t107 * t67 + t108 * t68 + (t212 - 0.2e1 * t233 + 0.2e1 * t234) * t194 + t238 * t140; t234 - t233 + t57 * t252 + t114 * t253 + (t63 / 0.2e1 - t94 / 0.2e1) * t73 + t212 + (-pkin(2) * t75 + t197 * t68 / 0.2e1 + t201 * t67 / 0.2e1) * t190 + m(7) * (t1 * t17 + t18 * t2 + t3 * t41) + m(6) * (t13 * t78 + t48 * t5 + t49 * t6) + m(5) * (t15 * t83 + t150 * t58 + t16 * t84) + m(4) * (-pkin(2) * t190 * t80 + t142 * t51 + t144 * t52) + (t93 / 0.2e1 + t129 / 0.2e1) * t140 + (t21 / 0.2e1 - t34 / 0.2e1) * t133 + (t33 / 0.2e1 + t66 / 0.2e1) * t193 + t18 * t19 + t17 * t20 + t29 * t39 / 0.2e1 + t30 * t40 / 0.2e1 + t41 * t11 + t48 * t32 + t49 * t31 + t3 * t50 + t2 * t61 + t1 * t62 + t13 * t77 + t78 * t26 + t83 * t60 + t84 * t59 + t86 * t9 / 0.2e1 + t87 * t10 / 0.2e1 + t6 * t89 + t5 * t90 + t74 * t95 / 0.2e1 + t58 * t96 + t16 * t118 + t15 * t119 + t107 * t130 / 0.2e1 + t108 * t131 / 0.2e1 + t213 * t56 + t134 * t35 / 0.2e1 + t214 * t113 + t80 * t141 + t142 * t92 + t144 * t91 + t51 * t148 + t52 * t149 + t150 * t43; t114 * t65 + 0.2e1 * t84 * t118 + 0.2e1 * t83 * t119 + t134 * t95 + 0.2e1 * t142 * t148 + 0.2e1 * t144 * t149 + 0.2e1 * t150 * t96 + 0.2e1 * t17 * t62 + 0.2e1 * t18 * t61 + t86 * t39 + t87 * t40 + 0.2e1 * t41 * t50 + 0.2e1 * t48 * t90 + 0.2e1 * t49 * t89 + 0.2e1 * t78 * t77 + Ifges(3,3) + t228 * t193 + (t63 - t94) * t133 + (t38 - t64) * t113 + (-0.2e1 * pkin(2) * t141 + t201 * t130 + t197 * t131) * t190 + m(7) * (t17 ^ 2 + t18 ^ 2 + t41 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2 + t78 ^ 2) + m(5) * (t150 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (pkin(2) ^ 2 * t190 ^ 2 + t142 ^ 2 + t144 ^ 2); (t189 * t59 + t192 * t60 + m(5) * (t15 * t192 + t16 * t189)) * pkin(3) + m(6) * (t13 * t177 + (-t5 * t196 + t6 * t200) * t176) + t238 + t15 * mrSges(5,1) - t16 * mrSges(5,2) + t51 * mrSges(4,1) - t52 * mrSges(4,2) + t211 * t56 + t110 * t20 + t111 * t19 + t29 * t251 + t30 * t250 + (t6 * mrSges(6,3) + t176 * t31 - t214) * t200 + t3 * t146 + t2 * t151 + t1 * t152 + t13 * t154 + t73 * t248 + t57 * t245 + t177 * t26 + m(7) * (t1 * t110 + t111 * t2 + t3 * t227) + (-t5 * mrSges(6,3) + t10 * t242 + t239 * t176 + t9 * t244 + t253) * t196; (-t48 * mrSges(6,3) + t237 * t176 + t40 * t242 + t39 * t244 + t252) * t196 + (t189 * t118 + t192 * t119 + m(5) * (t189 * t84 + t192 * t83)) * pkin(3) + t228 + m(6) * (t177 * t78 + (-t48 * t196 + t49 * t200) * t176) + t83 * mrSges(5,1) - t84 * mrSges(5,2) + t211 * t113 + t110 * t62 + t111 * t61 + (t49 * mrSges(6,3) + t176 * t89 - t213) * t200 + t86 * t251 + t87 * t250 + t142 * mrSges(4,1) - t144 * mrSges(4,2) + t41 * t146 + t18 * t151 + t17 * t152 + t78 * t154 + t133 * t248 + t114 * t245 + t177 * t77 + m(7) * (t110 * t17 + t111 * t18 + t41 * t227); 0.2e1 * t110 * t152 + 0.2e1 * t111 * t151 + 0.2e1 * t177 * t154 + Ifges(4,3) + Ifges(5,3) + (-t135 + t158) * t200 + m(7) * (t110 ^ 2 + t111 ^ 2 + t161) + m(6) * (t174 * t188 + t177 ^ 2 + t161) + m(5) * (t189 ^ 2 + t192 ^ 2) * pkin(3) ^ 2 + (-t136 * t195 + t137 * t199 + t146 * t256 + t160) * t196 + 0.2e1 * (mrSges(5,1) * t192 - mrSges(5,2) * t189) * pkin(3) + t215 * mrSges(6,3) * t256; -t239 * t200 + (t31 + t254) * t196 + m(7) * (t209 * t196 - t200 * t3) + m(6) * (t196 * t6 + t200 * t5) + m(5) * t58 + t43; -t237 * t200 + (t89 + t255) * t196 + m(7) * (t207 * t196 - t200 * t41) + m(6) * (t196 * t49 + t200 * t48) + m(5) * t150 + t96; -t200 * t146 + (m(7) * (t206 - t226) + t205) * t196; m(5) + m(6) * t215 + m(7) * (t216 * t186 + t188); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t209 * mrSges(7,3) + t10 * t243 + t3 * t153 + t9 * t242 + t30 * t246 + t29 * t247 + t56 * t249 + t21 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t209 + t254) * pkin(12); t48 * mrSges(6,1) - t49 * mrSges(6,2) + t207 * mrSges(7,3) + t113 * t249 + t41 * t153 + t39 * t242 + t40 * t243 + t87 * t246 + t86 * t247 + t63 + (-m(7) * t41 - t50) * pkin(5) + (m(7) * t207 + t255) * pkin(12); -pkin(5) * t146 + t137 * t243 + t136 * t242 + (m(7) * t206 + t205) * pkin(12) + (t159 * t242 + t157 * t244 + (-m(7) * pkin(5) - mrSges(6,1) + t153) * t176) * t196 + (-t176 * mrSges(6,2) - t155 / 0.2e1) * t200 + t206 * mrSges(7,3) + t156; m(7) * (pkin(12) * t210 + t240) - t200 * t153 + mrSges(7,3) * t210 - t154; Ifges(6,3) - 0.2e1 * pkin(5) * t153 + t195 * t159 + t199 * t157 + m(7) * (t216 * pkin(12) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t216 * pkin(12) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t17 - mrSges(7,2) * t18 + t38; mrSges(7,1) * t110 - mrSges(7,2) * t111 + t135; -t146; -t208 * pkin(12) + t155; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
