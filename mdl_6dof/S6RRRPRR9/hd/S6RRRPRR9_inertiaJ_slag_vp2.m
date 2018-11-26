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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:58
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:57:03
% EndTime: 2018-11-23 17:57:05
% DurationCPUTime: 2.53s
% Computational Cost: add. (7895->496), mult. (21016->728), div. (0->0), fcn. (24018->14), ass. (0->189)
t197 = sin(qJ(6));
t201 = cos(qJ(6));
t191 = sin(pkin(13));
t192 = sin(pkin(7));
t194 = cos(pkin(13));
t199 = sin(qJ(3));
t203 = cos(qJ(3));
t134 = (t191 * t203 + t194 * t199) * t192;
t195 = cos(pkin(7));
t198 = sin(qJ(5));
t202 = cos(qJ(5));
t113 = t134 * t198 - t202 * t195;
t114 = t134 * t202 + t195 * t198;
t225 = t192 * t203;
t226 = t192 * t199;
t133 = t191 * t226 - t194 * t225;
t86 = -t114 * t197 + t133 * t201;
t61 = -mrSges(7,2) * t113 + mrSges(7,3) * t86;
t87 = t114 * t201 + t133 * t197;
t62 = mrSges(7,1) * t113 - mrSges(7,3) * t87;
t256 = -t197 * t62 + t201 * t61;
t196 = cos(pkin(6));
t193 = sin(pkin(6));
t204 = cos(qJ(2));
t223 = t193 * t204;
t140 = -t192 * t223 + t195 * t196;
t200 = sin(qJ(2));
t220 = t195 * t204;
t107 = t196 * t225 + (-t199 * t200 + t203 * t220) * t193;
t108 = t196 * t226 + (t199 * t220 + t200 * t203) * t193;
t74 = t107 * t191 + t108 * t194;
t57 = t140 * t198 + t202 * t74;
t73 = -t194 * t107 + t108 * t191;
t29 = -t197 * t57 + t201 * t73;
t56 = -t202 * t140 + t198 * t74;
t19 = -mrSges(7,2) * t56 + mrSges(7,3) * t29;
t30 = t197 * t73 + t201 * t57;
t20 = mrSges(7,1) * t56 - mrSges(7,3) * t30;
t255 = t201 * t19 - t197 * t20;
t23 = Ifges(6,1) * t57 - Ifges(6,4) * t56 + Ifges(6,5) * t73;
t254 = t23 / 0.2e1;
t65 = Ifges(6,1) * t114 - Ifges(6,4) * t113 + Ifges(6,5) * t133;
t253 = t65 / 0.2e1;
t236 = Ifges(7,4) * t201;
t136 = -Ifges(7,6) * t202 + (-Ifges(7,2) * t197 + t236) * t198;
t252 = t136 / 0.2e1;
t237 = Ifges(7,4) * t197;
t137 = -Ifges(7,5) * t202 + (Ifges(7,1) * t201 - t237) * t198;
t251 = t137 / 0.2e1;
t155 = Ifges(7,5) * t197 + Ifges(7,6) * t201;
t250 = t155 / 0.2e1;
t156 = Ifges(6,5) * t198 + Ifges(6,6) * t202;
t249 = t156 / 0.2e1;
t157 = Ifges(7,2) * t201 + t237;
t248 = t157 / 0.2e1;
t159 = Ifges(7,1) * t197 + t236;
t247 = t159 / 0.2e1;
t160 = Ifges(6,1) * t198 + Ifges(6,4) * t202;
t246 = t160 / 0.2e1;
t245 = -t197 / 0.2e1;
t244 = t197 / 0.2e1;
t243 = t201 / 0.2e1;
t242 = pkin(1) * t196;
t241 = pkin(5) * t202;
t145 = pkin(9) * t223 + t200 * t242;
t101 = (t192 * t196 + t193 * t220) * pkin(10) + t145;
t172 = t204 * t242;
t224 = t193 * t200;
t112 = pkin(2) * t196 + t172 + (-pkin(10) * t195 - pkin(9)) * t224;
t124 = (-pkin(10) * t192 * t200 - pkin(2) * t204 - pkin(1)) * t193;
t221 = t195 * t203;
t51 = -t101 * t199 + t112 * t221 + t124 * t225;
t37 = pkin(3) * t140 - qJ(4) * t108 + t51;
t222 = t195 * t199;
t52 = t203 * t101 + t112 * t222 + t124 * t226;
t46 = qJ(4) * t107 + t52;
t16 = t191 * t37 + t194 * t46;
t14 = pkin(11) * t140 + t16;
t80 = -t112 * t192 + t195 * t124;
t58 = -pkin(3) * t107 + t80;
t25 = pkin(4) * t73 - pkin(11) * t74 + t58;
t6 = t202 * t14 + t198 * t25;
t11 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t32 = mrSges(6,1) * t73 - mrSges(6,3) * t57;
t240 = -t32 + t11;
t33 = Ifges(5,5) * t74 - Ifges(5,6) * t73 + Ifges(5,3) * t140;
t66 = Ifges(4,5) * t108 + Ifges(4,6) * t107 + Ifges(4,3) * t140;
t239 = t33 + t66;
t50 = -mrSges(7,1) * t86 + mrSges(7,2) * t87;
t90 = mrSges(6,1) * t133 - mrSges(6,3) * t114;
t238 = t50 - t90;
t170 = pkin(2) * t221;
t117 = pkin(3) * t195 + t170 + (-pkin(10) - qJ(4)) * t226;
t144 = pkin(2) * t222 + pkin(10) * t225;
t123 = qJ(4) * t225 + t144;
t84 = t191 * t117 + t194 * t123;
t79 = pkin(11) * t195 + t84;
t150 = (-pkin(3) * t203 - pkin(2)) * t192;
t88 = pkin(4) * t133 - pkin(11) * t134 + t150;
t49 = t198 * t88 + t202 * t79;
t143 = -pkin(9) * t224 + t172;
t235 = t143 * mrSges(3,1);
t234 = t145 * mrSges(3,2);
t129 = Ifges(4,5) * t226 + Ifges(4,6) * t225 + Ifges(4,3) * t195;
t93 = Ifges(5,5) * t134 - Ifges(5,6) * t133 + Ifges(5,3) * t195;
t229 = t93 + t129;
t178 = pkin(3) * t191 + pkin(11);
t228 = t178 * t198;
t227 = t178 * t202;
t219 = t197 * t198;
t218 = t198 * t201;
t217 = t197 ^ 2 + t201 ^ 2;
t188 = t198 ^ 2;
t190 = t202 ^ 2;
t216 = t188 + t190;
t8 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t56;
t21 = Ifges(6,5) * t57 - Ifges(6,6) * t56 + Ifges(6,3) * t73;
t22 = Ifges(6,4) * t57 - Ifges(6,2) * t56 + Ifges(6,6) * t73;
t215 = t8 / 0.2e1 - t22 / 0.2e1;
t38 = Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t113;
t64 = Ifges(6,4) * t114 - Ifges(6,2) * t113 + Ifges(6,6) * t133;
t214 = t38 / 0.2e1 - t64 / 0.2e1;
t63 = Ifges(6,5) * t114 - Ifges(6,6) * t113 + Ifges(6,3) * t133;
t213 = Ifges(3,5) * t224 + Ifges(3,6) * t223 + Ifges(3,3) * t196;
t179 = -pkin(3) * t194 - pkin(4);
t135 = Ifges(7,5) * t218 - Ifges(7,6) * t219 - Ifges(7,3) * t202;
t158 = Ifges(6,4) * t198 + Ifges(6,2) * t202;
t212 = t135 / 0.2e1 - t158 / 0.2e1;
t43 = t73 * mrSges(5,1) + t74 * mrSges(5,2);
t15 = -t191 * t46 + t194 * t37;
t96 = t133 * mrSges(5,1) + t134 * mrSges(5,2);
t211 = t217 * t198;
t83 = t117 * t194 - t191 * t123;
t4 = pkin(12) * t73 + t6;
t13 = -pkin(4) * t140 - t15;
t7 = pkin(5) * t56 - pkin(12) * t57 + t13;
t1 = -t197 * t4 + t201 * t7;
t2 = t197 * t7 + t201 * t4;
t210 = -t1 * t197 + t2 * t201;
t154 = -t202 * mrSges(6,1) + t198 * mrSges(6,2);
t5 = -t14 * t198 + t202 * t25;
t42 = pkin(12) * t133 + t49;
t78 = -pkin(4) * t195 - t83;
t47 = pkin(5) * t113 - pkin(12) * t114 + t78;
t17 = -t197 * t42 + t201 * t47;
t18 = t197 * t47 + t201 * t42;
t209 = -t17 * t197 + t18 * t201;
t48 = -t198 * t79 + t202 * t88;
t147 = -pkin(12) * t198 + t179 - t241;
t110 = t147 * t201 - t197 * t227;
t111 = t147 * t197 + t201 * t227;
t208 = -t110 * t197 + t111 * t201;
t151 = mrSges(7,2) * t202 - mrSges(7,3) * t219;
t152 = -mrSges(7,1) * t202 - mrSges(7,3) * t218;
t207 = t201 * t151 - t197 * t152;
t176 = t178 ^ 2;
t161 = t188 * t176;
t153 = -mrSges(7,1) * t201 + mrSges(7,2) * t197;
t149 = -mrSges(4,2) * t195 + mrSges(4,3) * t225;
t148 = mrSges(4,1) * t195 - mrSges(4,3) * t226;
t146 = -mrSges(7,1) * t219 - mrSges(7,2) * t218;
t142 = -pkin(10) * t226 + t170;
t141 = (-mrSges(4,1) * t203 + mrSges(4,2) * t199) * t192;
t131 = Ifges(4,5) * t195 + (Ifges(4,1) * t199 + Ifges(4,4) * t203) * t192;
t130 = Ifges(4,6) * t195 + (Ifges(4,4) * t199 + Ifges(4,2) * t203) * t192;
t119 = mrSges(5,1) * t195 - mrSges(5,3) * t134;
t118 = -mrSges(5,2) * t195 - mrSges(5,3) * t133;
t95 = Ifges(5,1) * t134 - Ifges(5,4) * t133 + Ifges(5,5) * t195;
t94 = Ifges(5,4) * t134 - Ifges(5,2) * t133 + Ifges(5,6) * t195;
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
t12 = [m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t13 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(4) * (t51 ^ 2 + t52 ^ 2 + t80 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t58 ^ 2) + (t21 - t34) * t73 + (t8 - t22) * t56 + m(3) * (t143 ^ 2 + t145 ^ 2) + ((Ifges(3,5) * t200 + Ifges(3,6) * t204) * t196 + 0.2e1 * (-t143 * t200 + t145 * t204) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t204 + mrSges(3,2) * t200) + t204 * (Ifges(3,4) * t200 + Ifges(3,2) * t204) + t200 * (Ifges(3,1) * t200 + Ifges(3,4) * t204) + m(3) * pkin(1) ^ 2) * t193) * t193 + Ifges(2,3) + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t19 + 0.2e1 * t1 * t20 + 0.2e1 * t13 * t26 + t29 * t9 + t30 * t10 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + t57 * t23 + 0.2e1 * t58 * t43 + 0.2e1 * t16 * t59 + 0.2e1 * t15 * t60 + t74 * t35 + 0.2e1 * t80 * t75 + 0.2e1 * t52 * t91 + 0.2e1 * t51 * t92 + t107 * t67 + t108 * t68 + (t213 - 0.2e1 * t234 + 0.2e1 * t235) * t196 + t239 * t140; -t234 + t235 + (t93 / 0.2e1 + t129 / 0.2e1) * t140 + (t21 / 0.2e1 - t34 / 0.2e1) * t133 + (t63 / 0.2e1 - t94 / 0.2e1) * t73 + t57 * t253 + t114 * t254 + t213 + (-pkin(2) * t75 + t199 * t68 / 0.2e1 + t203 * t67 / 0.2e1) * t192 + (t33 / 0.2e1 + t66 / 0.2e1) * t195 + m(7) * (t1 * t17 + t18 * t2 + t3 * t41) + m(6) * (t13 * t78 + t48 * t5 + t49 * t6) + m(5) * (t15 * t83 + t150 * t58 + t16 * t84) + m(4) * (-pkin(2) * t192 * t80 + t142 * t51 + t144 * t52) + t18 * t19 + t17 * t20 + t29 * t39 / 0.2e1 + t30 * t40 / 0.2e1 + t41 * t11 + t48 * t32 + t49 * t31 + t3 * t50 + t2 * t61 + t1 * t62 + t13 * t77 + t78 * t26 + t83 * t60 + t84 * t59 + t86 * t9 / 0.2e1 + t87 * t10 / 0.2e1 + t6 * t89 + t5 * t90 + t74 * t95 / 0.2e1 + t58 * t96 + t16 * t118 + t15 * t119 + t107 * t130 / 0.2e1 + t108 * t131 / 0.2e1 + t134 * t35 / 0.2e1 + t80 * t141 + t142 * t92 + t144 * t91 + t51 * t148 + t52 * t149 + t150 * t43 + t214 * t56 + t215 * t113; t114 * t65 + 0.2e1 * t84 * t118 + 0.2e1 * t83 * t119 + t134 * t95 + 0.2e1 * t142 * t148 + 0.2e1 * t144 * t149 + 0.2e1 * t150 * t96 + 0.2e1 * t17 * t62 + 0.2e1 * t18 * t61 + t86 * t39 + t87 * t40 + 0.2e1 * t41 * t50 + 0.2e1 * t48 * t90 + 0.2e1 * t49 * t89 + 0.2e1 * t78 * t77 + Ifges(3,3) + t229 * t195 + (t63 - t94) * t133 + (t38 - t64) * t113 + (-0.2e1 * pkin(2) * t141 + t203 * t130 + t199 * t131) * t192 + m(7) * (t17 ^ 2 + t18 ^ 2 + t41 ^ 2) + m(6) * (t48 ^ 2 + t49 ^ 2 + t78 ^ 2) + m(5) * (t150 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(4) * (pkin(2) ^ 2 * t192 ^ 2 + t142 ^ 2 + t144 ^ 2); (t191 * t59 + m(5) * (t15 * t194 + t16 * t191) + t194 * t60) * pkin(3) + m(6) * (t13 * t179 + (-t5 * t198 + t6 * t202) * t178) + t239 + t15 * mrSges(5,1) - t16 * mrSges(5,2) + t51 * mrSges(4,1) - t52 * mrSges(4,2) + t110 * t20 + t111 * t19 + t29 * t252 + t30 * t251 - t3 * t146 + t2 * t151 + t1 * t152 + t13 * t154 + t73 * t249 + t57 * t246 + t179 * t26 + t212 * t56 + (-t5 * mrSges(6,3) + t10 * t243 + t240 * t178 + t9 * t245 + t254) * t198 + (t6 * mrSges(6,3) + t178 * t31 - t215) * t202 + m(7) * (t1 * t110 + t111 * t2 + t3 * t228); (t191 * t118 + t194 * t119 + m(5) * (t191 * t84 + t194 * t83)) * pkin(3) + (-t48 * mrSges(6,3) + t238 * t178 + t40 * t243 + t39 * t245 + t253) * t198 + m(6) * (t179 * t78 + (-t48 * t198 + t49 * t202) * t178) + t229 + t83 * mrSges(5,1) - t84 * mrSges(5,2) + t110 * t62 + t111 * t61 + t86 * t252 + t87 * t251 + t142 * mrSges(4,1) - t144 * mrSges(4,2) - t41 * t146 + t18 * t151 + t17 * t152 + t78 * t154 + t133 * t249 + t114 * t246 + t179 * t77 + t212 * t113 + (t49 * mrSges(6,3) + t178 * t89 - t214) * t202 + m(7) * (t110 * t17 + t111 * t18 + t41 * t228); 0.2e1 * t110 * t152 + 0.2e1 * t111 * t151 + 0.2e1 * t179 * t154 + Ifges(4,3) + Ifges(5,3) + (-t135 + t158) * t202 + m(7) * (t110 ^ 2 + t111 ^ 2 + t161) + m(6) * (t176 * t190 + t179 ^ 2 + t161) + m(5) * (t191 ^ 2 + t194 ^ 2) * pkin(3) ^ 2 + (-t136 * t197 + t137 * t201 - 0.2e1 * t146 * t178 + t160) * t198 + 0.2e1 * (mrSges(5,1) * t194 - mrSges(5,2) * t191) * pkin(3) + 0.2e1 * t216 * t178 * mrSges(6,3); -t240 * t202 + (t31 + t255) * t198 + m(7) * (t210 * t198 - t202 * t3) + m(6) * (t198 * t6 + t202 * t5) + m(5) * t58 + t43; -t238 * t202 + (t89 + t256) * t198 + m(7) * (t209 * t198 - t202 * t41) + m(6) * (t198 * t49 + t202 * t48) + m(5) * t150 + t96; t202 * t146 + (m(7) * (t208 - t227) + t207) * t198; m(5) + m(7) * (t217 * t188 + t190) + m(6) * t216; t5 * mrSges(6,1) - t6 * mrSges(6,2) + t210 * mrSges(7,3) + t10 * t244 + t3 * t153 + t9 * t243 + t30 * t247 + t29 * t248 + t56 * t250 + t21 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t210 + t255) * pkin(12); t48 * mrSges(6,1) - t49 * mrSges(6,2) + t209 * mrSges(7,3) + t113 * t250 + t41 * t153 + t39 * t243 + t40 * t244 + t87 * t247 + t86 * t248 + t63 + (-m(7) * t41 - t50) * pkin(5) + (m(7) * t209 + t256) * pkin(12); pkin(5) * t146 + t136 * t243 + t137 * t244 + (m(7) * t208 + t207) * pkin(12) + (t157 * t245 + t159 * t243 + (-m(7) * pkin(5) - mrSges(6,1) + t153) * t178) * t198 + (-t178 * mrSges(6,2) - t155 / 0.2e1) * t202 + t208 * mrSges(7,3) + t156; m(7) * (pkin(12) * t211 + t241) - t202 * t153 + mrSges(7,3) * t211 - t154; Ifges(6,3) + m(7) * (t217 * pkin(12) ^ 2 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t153 + t201 * t157 + t197 * t159 + 0.2e1 * t217 * pkin(12) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t17 - mrSges(7,2) * t18 + t38; mrSges(7,1) * t110 - mrSges(7,2) * t111 + t135; t146; (-mrSges(7,1) * t197 - mrSges(7,2) * t201) * pkin(12) + t155; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
