% Calculate time derivative of joint inertia matrix for
% S5RRRRP6
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:52:51
% EndTime: 2019-12-31 21:53:00
% DurationCPUTime: 3.56s
% Computational Cost: add. (2726->310), mult. (6310->440), div. (0->0), fcn. (5316->6), ass. (0->140)
t250 = Ifges(5,4) + Ifges(6,4);
t249 = Ifges(5,1) + Ifges(6,1);
t248 = Ifges(5,2) + Ifges(6,2);
t206 = Ifges(5,5) + Ifges(6,5);
t135 = sin(qJ(3));
t136 = sin(qJ(2));
t138 = cos(qJ(3));
t139 = cos(qJ(2));
t99 = t135 * t136 - t138 * t139;
t247 = t206 * t99;
t236 = Ifges(5,6) + Ifges(6,6);
t246 = t236 * t99;
t137 = cos(qJ(4));
t245 = t250 * t137;
t134 = sin(qJ(4));
t244 = t250 * t134;
t242 = -t248 * t134 + t245;
t241 = t249 * t137 - t244;
t240 = t206 * t134 + t236 * t137;
t238 = t248 * t137 + t244;
t237 = t249 * t134 + t245;
t226 = (t134 ^ 2 + t137 ^ 2) * t138;
t235 = Ifges(5,3) + Ifges(6,3);
t100 = t135 * t139 + t138 * t136;
t181 = qJD(4) * t134;
t225 = qJD(2) + qJD(3);
t72 = t225 * t99;
t193 = t137 * t72;
t142 = t100 * t181 + t193;
t180 = qJD(4) * t137;
t168 = t100 * t180;
t196 = t134 * t72;
t143 = t168 - t196;
t73 = t225 * t100;
t234 = -t142 * t250 - t248 * t143 + t236 * t73;
t233 = -t249 * t142 - t143 * t250 + t206 * t73;
t232 = t100 * t242 + t246;
t202 = t100 * t241 + t247;
t231 = t242 * qJD(4);
t230 = t241 * qJD(4);
t182 = t206 * t180;
t214 = -pkin(7) - pkin(6);
t119 = t214 * t136;
t120 = t214 * t139;
t227 = t138 * t119 + t120 * t135;
t141 = t230 * t134 + t231 * t137 + t237 * t180 - t238 * t181;
t224 = 2 * m(5);
t223 = 2 * m(6);
t170 = qJD(2) * t214;
t109 = t136 * t170;
t157 = t139 * t170;
t86 = t119 * t135 - t120 * t138;
t34 = t86 * qJD(3) + t109 * t135 - t138 * t157;
t222 = 0.2e1 * t34;
t221 = -0.2e1 * t227;
t102 = mrSges(6,1) * t181 + mrSges(6,2) * t180;
t220 = 0.2e1 * t102;
t112 = -mrSges(6,1) * t137 + mrSges(6,2) * t134;
t219 = 0.2e1 * t112;
t218 = -0.2e1 * t134;
t217 = 0.2e1 * t137;
t213 = mrSges(6,3) * pkin(4);
t210 = pkin(2) * t138;
t176 = pkin(2) * qJD(2) * t136;
t27 = pkin(3) * t73 + pkin(8) * t72 + t176;
t33 = qJD(3) * t227 + t138 * t109 + t135 * t157;
t125 = -pkin(2) * t139 - pkin(1);
t61 = t99 * pkin(3) - t100 * pkin(8) + t125;
t179 = t134 * t27 + t137 * t33 + t61 * t180;
t5 = -t86 * t181 + t179;
t209 = t137 * t5;
t208 = t34 * t227;
t162 = -t134 * t33 + t137 * t27;
t78 = t137 * t86;
t30 = t134 * t61 + t78;
t6 = -t30 * qJD(4) + t162;
t207 = t6 * t134;
t204 = -qJ(5) - pkin(8);
t197 = pkin(2) * qJD(3);
t195 = t135 * mrSges(4,1);
t192 = t138 * mrSges(4,2);
t191 = t100 * t134;
t190 = t100 * t137;
t122 = pkin(2) * t135 + pkin(8);
t188 = t122 * t137;
t113 = -mrSges(5,1) * t137 + mrSges(5,2) * t134;
t185 = t135 * t113;
t131 = t137 * qJ(5);
t183 = -qJ(5) - t122;
t177 = 0.2e1 * qJD(4);
t175 = t138 * t197;
t174 = pkin(4) * t181;
t173 = m(6) * pkin(4) + mrSges(6,1);
t124 = -pkin(4) * t137 - pkin(3);
t165 = t236 * t134;
t164 = -t181 / 0.2e1;
t29 = -t134 * t86 + t137 * t61;
t161 = qJD(4) * t204;
t160 = 0.2e1 * t176;
t159 = -t206 * t193 + t235 * t73;
t158 = qJD(4) * t183;
t154 = mrSges(5,3) * t226;
t153 = -(2 * Ifges(4,4)) - t165;
t152 = mrSges(5,1) * t134 + mrSges(5,2) * t137;
t19 = pkin(4) * t99 - t100 * t131 + t29;
t24 = -qJ(5) * t191 + t30;
t147 = -t134 * t24 - t137 * t19;
t146 = -t134 * t30 - t137 * t29;
t144 = qJ(5) * t72 - qJD(5) * t100;
t17 = t143 * mrSges(6,1) - t142 * mrSges(6,2);
t103 = t152 * qJD(4);
t16 = t143 * pkin(4) + t34;
t3 = -qJ(5) * t168 + (-qJD(4) * t86 + t144) * t134 + t179;
t48 = pkin(4) * t191 - t227;
t140 = -t33 * mrSges(4,2) + mrSges(5,3) * t209 - Ifges(4,5) * t72 + t48 * t102 - t227 * t103 + t16 * t112 + (-mrSges(4,1) + t113) * t34 + (-t236 * t181 + t182) * t99 / 0.2e1 + t233 * t134 / 0.2e1 - t231 * t191 / 0.2e1 + t230 * t190 / 0.2e1 + t232 * t164 + t202 * t180 / 0.2e1 + t237 * (t100 * t164 - t193 / 0.2e1) + t238 * (t196 / 0.2e1 - t168 / 0.2e1) + (-Ifges(4,6) + t240 / 0.2e1) * t73 + (t3 * mrSges(6,3) + t234 / 0.2e1) * t137;
t130 = t137 * qJD(5);
t123 = -pkin(3) - t210;
t114 = pkin(8) * t137 + t131;
t111 = t204 * t134;
t110 = t124 - t210;
t108 = t135 * t197 + t174;
t96 = t131 + t188;
t95 = t183 * t134;
t89 = -qJD(5) * t134 + t137 * t161;
t88 = t134 * t161 + t130;
t71 = (-qJD(5) - t175) * t134 + t137 * t158;
t70 = t134 * t158 + t137 * t175 + t130;
t65 = mrSges(5,1) * t99 - mrSges(5,3) * t190;
t64 = mrSges(6,1) * t99 - mrSges(6,3) * t190;
t63 = -mrSges(5,2) * t99 - mrSges(5,3) * t191;
t62 = -mrSges(6,2) * t99 - mrSges(6,3) * t191;
t56 = t152 * t100;
t55 = (mrSges(6,1) * t134 + mrSges(6,2) * t137) * t100;
t23 = -mrSges(5,2) * t73 - t143 * mrSges(5,3);
t22 = -mrSges(6,2) * t73 - t143 * mrSges(6,3);
t21 = mrSges(5,1) * t73 + t142 * mrSges(5,3);
t20 = mrSges(6,1) * t73 + t142 * mrSges(6,3);
t18 = t143 * mrSges(5,1) - t142 * mrSges(5,2);
t1 = pkin(4) * t73 + t144 * t137 + (-t78 + (qJ(5) * t100 - t61) * t134) * qJD(4) + t162;
t2 = [0.2e1 * t1 * t64 + 0.2e1 * t16 * t55 + 0.2e1 * t48 * t17 + t18 * t221 + 0.2e1 * t19 * t20 + 0.2e1 * t29 * t21 + 0.2e1 * t24 * t22 + 0.2e1 * t30 * t23 + 0.2e1 * t3 * t62 + t56 * t222 + 0.2e1 * t5 * t63 + 0.2e1 * t6 * t65 + 0.2e1 * (t125 * mrSges(4,1) - t86 * mrSges(4,3)) * t73 + (t29 * t6 + t30 * t5 - t208) * t224 + (t1 * t19 + t16 * t48 + t24 * t3) * t223 + 0.2e1 * m(4) * (t125 * t176 + t33 * t86 - t208) - (0.2e1 * mrSges(4,2) * t125 + mrSges(4,3) * t221 - t134 * t232 + t202 * t137) * t72 + (mrSges(4,1) * t160 - 0.2e1 * mrSges(4,3) * t33 + ((2 * Ifges(4,2)) + t235) * t73 - t153 * t72 + t159) * t99 + (mrSges(4,2) * t160 + mrSges(4,3) * t222 - 0.2e1 * Ifges(4,1) * t72 + t233 * t137 - t234 * t134 + (t206 * t137 + t153) * t73 + ((-t232 - t246) * t137 + (-t202 - t247) * t134) * qJD(4)) * t100 + 0.2e1 * (-pkin(1) * (mrSges(3,1) * t136 + mrSges(3,2) * t139) + (Ifges(3,1) - Ifges(3,2)) * t136 * t139 + (-t136 ^ 2 + t139 ^ 2) * Ifges(3,4)) * qJD(2); (t147 * mrSges(6,3) + t146 * mrSges(5,3) + (m(5) * t146 - t134 * t63 - t137 * t65) * t122) * qJD(4) + (m(4) * (t135 * t33 - t138 * t34) + (-t135 * t73 + t138 * t72) * mrSges(4,3) + ((-t99 * mrSges(4,3) - t134 * t65 + t137 * t63 + m(4) * t86 + m(5) * (-t134 * t29 + t137 * t30)) * t138 + (t100 * mrSges(4,3) + t56 - (m(4) + m(5)) * t227) * t135) * qJD(3)) * pkin(2) + (Ifges(3,5) * t139 - Ifges(3,6) * t136 + (-mrSges(3,1) * t139 + mrSges(3,2) * t136) * pkin(6)) * qJD(2) + (-t6 * mrSges(5,3) - t1 * mrSges(6,3) - t122 * t21) * t134 + m(5) * (-t122 * t207 + t123 * t34 + t5 * t188) + t140 + m(6) * (t1 * t95 + t108 * t48 + t110 * t16 + t19 * t71 + t24 * t70 + t3 * t96) + t23 * t188 + t70 * t62 + t71 * t64 + t95 * t20 + t96 * t22 + t108 * t55 + t110 * t17 + t123 * t18; 0.2e1 * t123 * t103 + t108 * t219 + t110 * t220 + (t108 * t110 + t70 * t96 + t71 * t95) * t223 + (t71 * t218 + t70 * t217 + (-t134 * t96 - t137 * t95) * t177) * mrSges(6,3) + (-0.2e1 * t192 - 0.2e1 * t195 + (t122 * t226 + t123 * t135) * t224 + 0.2e1 * t185 + 0.2e1 * t154) * t197 + t141; (-t134 * t21 + t137 * t23 + m(5) * (-t29 * t180 - t30 * t181 - t207 + t209) - t63 * t181 - t65 * t180) * pkin(8) + t140 + (t146 * qJD(4) - t207) * mrSges(5,3) + (t147 * qJD(4) - t1 * t134) * mrSges(6,3) + m(6) * (t1 * t111 + t114 * t3 + t124 * t16 + t48 * t174 + t19 * t89 + t24 * t88) + t55 * t174 + t88 * t62 + t89 * t64 + t111 * t20 + t114 * t22 + t124 * t17 + (-m(5) * t34 - t18) * pkin(3); m(6) * (t108 * t124 + t110 * t174 + t111 * t71 + t114 * t70 + t88 * t96 + t89 * t95) + (t108 + t174) * t112 + (-pkin(3) + t123) * t103 + (t110 + t124) * t102 + (m(5) * (-pkin(3) * t135 + pkin(8) * t226) - t195 + t185 - t192 + t154) * t197 + ((t70 + t88) * t137 + (-t71 - t89) * t134 + ((-t111 - t95) * t137 + (-t114 - t96) * t134) * qJD(4)) * mrSges(6,3) + t141; -0.2e1 * pkin(3) * t103 + t174 * t219 + t124 * t220 + (t111 * t89 + t114 * t88 + t124 * t174) * t223 + (t89 * t218 + t88 * t217 + (-t111 * t137 - t114 * t134) * t177) * mrSges(6,3) + t141; mrSges(5,1) * t6 + mrSges(6,1) * t1 - mrSges(5,2) * t5 - mrSges(6,2) * t3 + t72 * t165 + (m(6) * t1 + t20) * pkin(4) - t240 * t100 * qJD(4) + t159; -mrSges(6,2) * t70 + t173 * t71 - t152 * t175 + ((-mrSges(5,1) * t122 - t213) * t137 + (mrSges(5,2) * t122 - t236) * t134) * qJD(4) + t182; -mrSges(6,2) * t88 + t173 * t89 + ((-mrSges(5,1) * pkin(8) - t213) * t137 + (mrSges(5,2) * pkin(8) - t236) * t134) * qJD(4) + t182; 0; m(6) * t16 + t17; m(6) * t108 + t102; m(6) * t174 + t102; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
