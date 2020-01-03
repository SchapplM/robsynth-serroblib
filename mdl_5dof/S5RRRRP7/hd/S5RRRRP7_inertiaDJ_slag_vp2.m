% Calculate time derivative of joint inertia matrix for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:56:09
% EndTime: 2019-12-31 21:56:17
% DurationCPUTime: 2.75s
% Computational Cost: add. (2723->279), mult. (6294->397), div. (0->0), fcn. (5321->6), ass. (0->140)
t244 = Ifges(5,1) + Ifges(6,1);
t206 = Ifges(6,4) + Ifges(5,5);
t132 = sin(qJ(3));
t133 = sin(qJ(2));
t135 = cos(qJ(3));
t136 = cos(qJ(2));
t98 = t132 * t133 - t135 * t136;
t243 = t206 * t98;
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t242 = t131 ^ 2 + t134 ^ 2;
t197 = Ifges(6,5) * t131;
t199 = Ifges(5,4) * t131;
t241 = t244 * t134 + t197 - t199;
t196 = Ifges(6,5) * t134;
t198 = Ifges(5,4) * t134;
t240 = t244 * t131 - t196 + t198;
t239 = Ifges(5,6) * t134 + t131 * t206;
t238 = Ifges(6,2) + Ifges(5,3);
t179 = qJD(4) * t131;
t227 = qJD(2) + qJD(3);
t70 = t227 * t98;
t189 = t134 * t70;
t99 = t132 * t136 + t135 * t133;
t142 = t99 * t179 + t189;
t178 = qJD(4) * t134;
t168 = t99 * t178;
t192 = t131 * t70;
t143 = t168 - t192;
t71 = t227 * t99;
t237 = t206 * t71 + (-Ifges(5,4) + Ifges(6,5)) * t143 - t244 * t142;
t204 = t241 * t99 + t243;
t110 = -t134 * mrSges(5,1) + t131 * mrSges(5,2);
t236 = -mrSges(4,1) + t110;
t235 = t241 * qJD(4);
t234 = Ifges(6,6) * t179 + t206 * t178;
t193 = pkin(2) * qJD(3);
t172 = t135 * t193;
t233 = t242 * t172;
t174 = pkin(2) * qJD(2) * t133;
t29 = pkin(3) * t71 + pkin(8) * t70 + t174;
t217 = -pkin(7) - pkin(6);
t165 = qJD(2) * t217;
t107 = t133 * t165;
t157 = t136 * t165;
t115 = t217 * t133;
t116 = t217 * t136;
t228 = t135 * t115 + t116 * t132;
t36 = t228 * qJD(3) + t135 * t107 + t132 * t157;
t122 = -pkin(2) * t136 - pkin(1);
t61 = t98 * pkin(3) - t99 * pkin(8) + t122;
t84 = t115 * t132 - t116 * t135;
t6 = t131 * t29 + t134 * t36 + t61 * t178 - t179 * t84;
t201 = t131 * t61 + t134 * t84;
t7 = -qJD(4) * t201 - t131 * t36 + t134 * t29;
t231 = -t7 * t131 + t134 * t6;
t2 = qJ(5) * t71 + qJD(5) * t98 + t6;
t25 = qJ(5) * t98 + t201;
t31 = -t131 * t84 + t134 * t61;
t26 = -pkin(4) * t98 - t31;
t4 = -pkin(4) * t71 - t7;
t230 = t131 * t4 + t134 * t2 + t26 * t178 - t25 * t179;
t19 = mrSges(5,1) * t143 - mrSges(5,2) * t142;
t37 = qJD(3) * t84 + t107 * t132 - t135 * t157;
t229 = m(5) * t37 + t19;
t226 = (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t242) * t172;
t225 = 0.2e1 * m(5);
t224 = 2 * m(6);
t223 = -2 * Ifges(4,4);
t222 = 0.2e1 * t37;
t221 = -0.2e1 * t228;
t109 = -t134 * mrSges(6,1) - t131 * mrSges(6,3);
t220 = 0.2e1 * t109;
t219 = t71 / 0.2e1;
t112 = Ifges(5,2) * t134 + t199;
t215 = -t112 / 0.2e1;
t213 = pkin(2) * t135;
t209 = t37 * t228;
t207 = t98 * Ifges(5,6);
t150 = Ifges(6,3) * t131 + t196;
t42 = Ifges(6,6) * t98 + t150 * t99;
t151 = -Ifges(5,2) * t131 + t198;
t43 = t151 * t99 + t207;
t205 = t42 - t43;
t191 = t131 * t99;
t62 = -mrSges(5,2) * t98 - mrSges(5,3) * t191;
t65 = -mrSges(6,2) * t191 + mrSges(6,3) * t98;
t203 = t62 + t65;
t188 = t134 * t99;
t63 = mrSges(5,1) * t98 - mrSges(5,3) * t188;
t64 = -mrSges(6,1) * t98 + mrSges(6,2) * t188;
t202 = -t63 + t64;
t120 = pkin(2) * t132 + pkin(8);
t200 = t233 * t120;
t195 = Ifges(5,6) * t131;
t190 = t132 * t228;
t187 = qJD(4) * t99;
t185 = t131 * t135;
t183 = t134 * t135;
t182 = t233 * pkin(8);
t177 = qJD(5) * t134;
t175 = m(6) * t177;
t173 = t132 * t193;
t164 = -t179 / 0.2e1;
t162 = 0.2e1 * t174;
t159 = mrSges(6,2) * t177 + t234;
t155 = t131 * mrSges(5,1) + t134 * mrSges(5,2);
t154 = t131 * mrSges(6,1) - t134 * mrSges(6,3);
t149 = -pkin(4) * t134 - qJ(5) * t131;
t148 = pkin(4) * t131 - qJ(5) * t134;
t145 = t143 * Ifges(6,6) - t206 * t189 + t238 * t71;
t144 = t236 * t173;
t108 = -pkin(3) + t149;
t87 = pkin(4) * t179 - qJ(5) * t178 - qJD(5) * t131;
t22 = -t71 * mrSges(6,1) - t142 * mrSges(6,2);
t141 = mrSges(6,2) * t149 - t195;
t103 = t150 * qJD(4);
t104 = t151 * qJD(4);
t111 = -Ifges(6,3) * t134 + t197;
t140 = (t111 - t112) * t179 + t240 * t178 + (-t103 + t104) * t134 + t235 * t131;
t139 = m(6) * t149 + t109 + t110;
t21 = mrSges(5,1) * t71 + mrSges(5,3) * t142;
t23 = -mrSges(5,2) * t71 - mrSges(5,3) * t143;
t24 = -mrSges(6,2) * t143 + mrSges(6,3) * t71;
t138 = (t23 + t24) * t134 + (-t21 + t22) * t131 + (-t131 * t203 + t134 * t202) * qJD(4) + m(6) * t230 + m(5) * (-t178 * t31 - t179 * t201 + t231);
t101 = t154 * qJD(4);
t102 = t155 * qJD(4);
t14 = -Ifges(6,5) * t142 + Ifges(6,6) * t71 + Ifges(6,3) * t143;
t15 = -Ifges(5,4) * t142 - Ifges(5,2) * t143 + Ifges(5,6) * t71;
t35 = t148 * t99 - t228;
t9 = (-pkin(4) * t70 + qJ(5) * t187) * t131 + (qJ(5) * t70 + (pkin(4) * qJD(4) - qJD(5)) * t99) * t134 + t37;
t137 = t42 * t179 / 0.2e1 + t43 * t164 + t168 * t215 - t36 * mrSges(4,2) - Ifges(4,5) * t70 - Ifges(4,6) * t71 + t35 * t101 - t228 * t102 + t9 * t109 + t236 * t37 + t239 * t219 + (-Ifges(5,6) * t179 + t234) * t98 / 0.2e1 + t237 * t131 / 0.2e1 - (t111 / 0.2e1 + t215) * t192 + (-t104 / 0.2e1 + t103 / 0.2e1) * t191 + t235 * t188 / 0.2e1 + (-Ifges(6,6) * t219 + t15 / 0.2e1 - t14 / 0.2e1) * t134 + ((-t131 * t201 - t134 * t31) * qJD(4) + t231) * mrSges(5,3) + (t99 * t111 + t204) * t178 / 0.2e1 + t230 * mrSges(6,2) + t240 * (t99 * t164 - t189 / 0.2e1);
t124 = mrSges(6,2) * t178;
t121 = -pkin(3) - t213;
t92 = t108 - t213;
t85 = t87 + t173;
t56 = t155 * t99;
t55 = t154 * t99;
t18 = mrSges(6,1) * t143 + mrSges(6,3) * t142;
t1 = [0.2e1 * t35 * t18 + t19 * t221 + 0.2e1 * t2 * t65 + 0.2e1 * t31 * t21 + 0.2e1 * t26 * t22 + 0.2e1 * t201 * t23 + 0.2e1 * t25 * t24 + t56 * t222 + 0.2e1 * t4 * t64 + 0.2e1 * t9 * t55 + 0.2e1 * t6 * t62 + 0.2e1 * t7 * t63 + 0.2e1 * (mrSges(4,1) * t122 - mrSges(4,3) * t84) * t71 + 0.2e1 * m(4) * (t122 * t174 + t36 * t84 - t209) + (t201 * t6 + t31 * t7 - t209) * t225 + (t2 * t25 + t26 * t4 + t35 * t9) * t224 - (0.2e1 * t122 * mrSges(4,2) + mrSges(4,3) * t221 + t131 * t205 + t134 * t204) * t70 + (mrSges(4,1) * t162 - 0.2e1 * t36 * mrSges(4,3) - (t223 - t195) * t70 + ((2 * Ifges(4,2)) + t238) * t71 + t145) * t98 + (mrSges(4,2) * t162 + mrSges(4,3) * t222 - 0.2e1 * Ifges(4,1) * t70 + t237 * t134 + (t14 - t15) * t131 + (t223 + t206 * t134 + (-Ifges(5,6) + Ifges(6,6)) * t131) * t71 + ((t205 - t207) * t134 + (-t204 - t243) * t131) * qJD(4)) * t99 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t133 + mrSges(3,2) * t136) + (-Ifges(3,2) + Ifges(3,1)) * t133 * t136 + (-t133 ^ 2 + t136 ^ 2) * Ifges(3,4)) * qJD(2); t137 + t138 * t120 + (m(4) * (t132 * t36 - t135 * t37) + (-t132 * t71 + t135 * t70) * mrSges(4,3) + ((t99 * mrSges(4,3) + t56) * t132 + (-t98 * mrSges(4,3) + t131 * t202 + t134 * t203) * t135 + m(6) * (t183 * t25 + t185 * t26) + m(5) * (t183 * t201 - t185 * t31 - t190) + m(4) * (t135 * t84 - t190)) * qJD(3)) * pkin(2) + (Ifges(3,5) * t136 - Ifges(3,6) * t133 + (-mrSges(3,1) * t136 + mrSges(3,2) * t133) * pkin(6)) * qJD(2) + m(6) * (t35 * t85 + t9 * t92) + t85 * t55 + t92 * t18 + t229 * t121; 0.2e1 * t92 * t101 + 0.2e1 * t121 * t102 + t85 * t220 + 0.2e1 * t144 + (t85 * t92 + t200) * t224 + (t121 * t173 + t200) * t225 + 0.2e1 * t226 + t140; t137 + t138 * pkin(8) + t87 * t55 + t108 * t18 + m(6) * (t108 * t9 + t35 * t87) - t229 * pkin(3); (t87 + t85) * t109 + (t121 - pkin(3)) * t102 + (t92 + t108) * t101 + t144 + m(6) * (t108 * t85 + t87 * t92 + t182) + m(5) * (-pkin(3) * t173 + t182) + t226 + t140; -0.2e1 * pkin(3) * t102 + t87 * t220 + 0.2e1 * (m(6) * t87 + t101) * t108 + t140; Ifges(5,6) * t192 - pkin(4) * t22 + t2 * mrSges(6,3) + m(6) * (-pkin(4) * t4 + qJ(5) * t2 + qJD(5) * t25) + qJD(5) * t65 + qJ(5) * t24 - t4 * mrSges(6,1) + t7 * mrSges(5,1) - t6 * mrSges(5,2) - t239 * t187 + t145; t120 * t175 + (-m(6) * t148 - t154 - t155) * t172 + (t120 * t139 + t141) * qJD(4) + t159; t141 * qJD(4) + (qJD(4) * t139 + t175) * pkin(8) + t159; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); m(6) * t4 + t22; t124 + (t120 * t178 + t131 * t172) * m(6); m(6) * pkin(8) * t178 + t124; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
