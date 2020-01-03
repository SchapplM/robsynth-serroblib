% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:21
% EndTime: 2019-12-31 18:21:26
% DurationCPUTime: 2.22s
% Computational Cost: add. (5124->294), mult. (10911->426), div. (0->0), fcn. (10577->8), ass. (0->155)
t248 = -m(6) / 0.2e1;
t165 = cos(qJ(3));
t153 = sin(pkin(8)) * pkin(1) + pkin(6);
t161 = sin(pkin(9));
t184 = pkin(4) * t161 + t153;
t126 = t184 * t165;
t162 = cos(pkin(9));
t230 = sin(qJ(5));
t231 = cos(qJ(5));
t133 = -t230 * t161 + t231 * t162;
t119 = t133 * t165;
t214 = t119 * mrSges(6,2);
t134 = -t161 * t231 - t162 * t230;
t117 = t134 * t165;
t216 = t117 * mrSges(6,1);
t207 = t214 / 0.2e1 - t216 / 0.2e1;
t247 = t126 * t248 - t207;
t160 = t162 ^ 2;
t246 = t160 / 0.2e1;
t164 = sin(qJ(3));
t245 = t164 * t165;
t211 = t162 * mrSges(5,2);
t213 = t161 * mrSges(5,1);
t174 = t213 / 0.2e1 + t211 / 0.2e1;
t227 = t164 * pkin(3);
t147 = -qJ(4) * t165 + t227;
t201 = t161 * t164;
t90 = t162 * t147 + t153 * t201;
t198 = t162 * t164;
t91 = t161 * t147 - t153 * t198;
t178 = -t90 * t161 + t91 * t162;
t244 = m(5) / 0.2e1;
t243 = m(6) / 0.2e1;
t226 = pkin(7) + qJ(4);
t144 = t226 * t161;
t146 = t226 * t162;
t89 = -t144 * t230 + t146 * t231;
t242 = -t89 / 0.2e1;
t116 = t134 * t164;
t241 = t116 / 0.2e1;
t240 = t117 / 0.2e1;
t118 = t133 * t164;
t239 = t118 / 0.2e1;
t238 = t119 / 0.2e1;
t237 = t133 / 0.2e1;
t236 = -t134 / 0.2e1;
t235 = -t161 / 0.2e1;
t234 = t162 / 0.2e1;
t233 = t164 / 0.2e1;
t232 = -t165 / 0.2e1;
t229 = m(5) * t153;
t225 = mrSges(6,1) * t116;
t224 = mrSges(6,2) * t118;
t223 = Ifges(5,4) * t161;
t222 = Ifges(5,4) * t162;
t221 = Ifges(6,4) * t118;
t220 = Ifges(6,4) * t134;
t219 = Ifges(5,5) * t162;
t218 = Ifges(5,6) * t161;
t92 = mrSges(6,2) * t165 + t116 * mrSges(6,3);
t217 = t116 * t92;
t94 = -mrSges(6,1) * t165 - t118 * mrSges(6,3);
t215 = t118 * t94;
t212 = t161 * Ifges(5,2);
t65 = t118 * mrSges(6,1) + mrSges(6,2) * t116;
t9 = -t215 / 0.2e1 + t217 / 0.2e1 + t65 * t232 + (-t116 ^ 2 / 0.2e1 - t118 ^ 2 / 0.2e1) * mrSges(6,3);
t210 = t9 * qJD(1);
t205 = -mrSges(5,1) * t162 + mrSges(5,2) * t161 - mrSges(4,1);
t137 = t165 * mrSges(5,2) - mrSges(5,3) * t201;
t139 = -t165 * mrSges(5,1) - mrSges(5,3) * t198;
t187 = -cos(pkin(8)) * pkin(1) - pkin(2);
t132 = -pkin(3) * t165 - t164 * qJ(4) + t187;
t121 = t162 * t132;
t64 = -pkin(7) * t198 + t121 + (-t153 * t161 - pkin(4)) * t165;
t196 = t165 * t153;
t81 = t161 * t132 + t162 * t196;
t72 = -pkin(7) * t201 + t81;
t26 = -t230 * t72 + t231 * t64;
t27 = t230 * t64 + t231 * t72;
t80 = -t161 * t196 + t121;
t13 = t217 - t215 + m(6) * (t116 * t27 - t118 * t26) + (m(5) * (-t161 * t81 - t162 * t80) - t161 * t137 - t162 * t139) * t164;
t204 = t13 * qJD(1);
t203 = t153 * t164;
t197 = t162 * t165;
t140 = t164 * mrSges(5,1) - mrSges(5,3) * t197;
t202 = t161 * t140;
t200 = t161 * t165;
t138 = -t164 * mrSges(5,2) - mrSges(5,3) * t200;
t199 = t162 * t138;
t195 = t65 * qJD(5);
t194 = Ifges(6,5) * t116 - Ifges(6,6) * t118;
t193 = Ifges(6,5) * t133 + Ifges(6,6) * t134;
t159 = t161 ^ 2;
t192 = t159 + t160;
t191 = qJD(3) * t165;
t190 = m(5) * t233;
t183 = t192 * mrSges(5,3);
t182 = t192 * qJ(4);
t114 = Ifges(5,6) * t164 + (-t212 + t222) * t165;
t115 = Ifges(5,5) * t164 + (t162 * Ifges(5,1) - t223) * t165;
t125 = t184 * t164;
t127 = (t211 + t213) * t165;
t175 = -Ifges(6,5) * t119 / 0.2e1 - Ifges(6,6) * t117 / 0.2e1;
t73 = t164 * pkin(4) - pkin(7) * t197 + t90;
t79 = -pkin(7) * t200 + t91;
t30 = -t230 * t79 + t231 * t73;
t31 = t230 * t73 + t231 * t79;
t60 = Ifges(6,2) * t116 - t165 * Ifges(6,6) + t221;
t61 = Ifges(6,4) * t119 + Ifges(6,2) * t117 + Ifges(6,6) * t164;
t107 = Ifges(6,4) * t116;
t62 = Ifges(6,1) * t118 - t165 * Ifges(6,5) + t107;
t63 = Ifges(6,1) * t119 + Ifges(6,4) * t117 + Ifges(6,5) * t164;
t66 = t214 - t216;
t93 = -t164 * mrSges(6,2) + mrSges(6,3) * t117;
t95 = t164 * mrSges(6,1) - mrSges(6,3) * t119;
t1 = t91 * t137 + t81 * t138 + t90 * t139 + t80 * t140 + t126 * (t224 - t225) + t63 * t239 + t62 * t238 + t125 * t66 + t61 * t241 + t60 * t240 + t31 * t92 + t27 * t93 + t30 * t94 + t26 * t95 + m(6) * (t125 * t126 + t26 * t30 + t27 * t31) + m(5) * (t80 * t90 + t81 * t91) + (Ifges(6,5) * t239 + Ifges(6,6) * t241 + t187 * mrSges(4,1) + t153 * t127 + t114 * t235 + t115 * t234 + (-Ifges(4,4) + t219 / 0.2e1 - t218 / 0.2e1) * t164) * t164 + (t187 * mrSges(4,2) + (-Ifges(4,2) - Ifges(6,3) + Ifges(4,1) - Ifges(5,3) + Ifges(5,1) * t246 + (t211 + t229) * t153 + (t153 * mrSges(5,1) - t222 + t212 / 0.2e1) * t161) * t164 + t175 + (Ifges(4,4) + t218 - t219) * t165) * t165;
t173 = t137 * t234 + t139 * t235;
t179 = -t161 * t80 + t162 * t81;
t6 = t94 * t240 + t95 * t241 + t92 * t238 + t93 * t239 + (t199 / 0.2e1 - t202 / 0.2e1 - t225 / 0.2e1 + t224 / 0.2e1 + t174 * t164) * t164 + (-t66 / 0.2e1 - t127 / 0.2e1 + t173) * t165 + (t30 * t116 + t26 * t117 + t31 * t118 + t27 * t119 + t125 * t164 - t126 * t165) * t243 + ((t179 - t196) * t165 + (t178 + t203) * t164) * t244;
t181 = t1 * qJD(1) + t6 * qJD(2);
t67 = -Ifges(6,2) * t118 + t107;
t68 = Ifges(6,1) * t116 - t221;
t5 = t26 * t92 - t27 * t94 + t194 * t232 + t125 * t65 + (-t27 * mrSges(6,3) + t68 / 0.2e1 - t60 / 0.2e1) * t118 + (-t26 * mrSges(6,3) + t62 / 0.2e1 + t67 / 0.2e1) * t116;
t180 = t5 * qJD(1) + t9 * qJD(2);
t21 = m(5) * (-0.1e1 + t192) * t245 + m(6) * (t116 * t117 + t118 * t119 - t245);
t177 = t6 * qJD(1) + t21 * qJD(2);
t82 = -t134 * mrSges(6,1) + mrSges(6,2) * t133;
t176 = qJD(1) * t65 + qJD(3) * t82;
t172 = (t116 * t134 + t118 * t133) * t243;
t154 = -pkin(4) * t162 - pkin(3);
t131 = Ifges(6,4) * t133;
t84 = Ifges(6,2) * t134 + t131;
t85 = Ifges(6,2) * t133 - t220;
t86 = Ifges(6,1) * t133 + t220;
t87 = -Ifges(6,1) * t134 + t131;
t14 = t154 * t82 + (-t86 / 0.2e1 + t85 / 0.2e1) * t134 + (t87 / 0.2e1 + t84 / 0.2e1) * t133;
t169 = t82 * t232;
t16 = t169 + t207;
t88 = -t144 * t231 - t146 * t230;
t166 = (t62 / 0.4e1 + t67 / 0.4e1) * t133 + (-t68 / 0.4e1 + t60 / 0.4e1) * t134 + (t86 / 0.4e1 - t85 / 0.4e1 + mrSges(6,3) * t242) * t118 + (t87 / 0.4e1 + t84 / 0.4e1 - t88 * mrSges(6,3) / 0.2e1) * t116 + t125 * t82 / 0.2e1 + t154 * t65 / 0.2e1 - t165 * t193 / 0.4e1 + t88 * t92 / 0.2e1 + t94 * t242;
t168 = Ifges(6,3) * t233 + t30 * mrSges(6,1) / 0.2e1 - t31 * mrSges(6,2) / 0.2e1 - t175;
t3 = t166 - t168;
t171 = t3 * qJD(1) + t16 * qJD(2) + t14 * qJD(3);
t20 = (t133 ^ 2 + t134 ^ 2) * mrSges(6,3) + t183 + m(6) * (t133 * t89 + t134 * t88) + m(5) * t182;
t25 = t172 + (t248 + (t159 / 0.2e1 + t246 - 0.1e1 / 0.2e1) * m(5)) * t164;
t167 = (t116 * t237 + t118 * t236) * mrSges(6,3) + t179 * t244 + (t89 * t116 - t88 * t118 + t133 * t27 + t134 * t26) * t243 + t92 * t237 + t134 * t94 / 0.2e1 + t173;
t8 = (-t229 / 0.2e1 - t174) * t165 + t167 + t247;
t170 = qJD(1) * t8 + qJD(2) * t25 + qJD(3) * t20;
t83 = -mrSges(6,1) * t133 - mrSges(6,2) * t134;
t24 = m(6) * t233 + t190 * t192 + t172 + t190;
t15 = t169 - t207;
t7 = t174 * t165 + t196 * t244 + t167 - t247;
t4 = t166 + t168;
t2 = t6 * qJD(3) + t9 * qJD(5);
t10 = [qJD(3) * t1 + qJD(4) * t13 + qJD(5) * t5, t2, t7 * qJD(4) + t4 * qJD(5) + (Ifges(4,5) + (Ifges(5,1) * t161 + t222) * t234 + (Ifges(5,2) * t162 + t223) * t235 + (-m(5) * pkin(3) + t205) * t153) * t191 + t181 + (-Ifges(4,6) * t164 + t161 * t115 / 0.2e1 + t114 * t234 + t154 * t66 + t61 * t237 + t63 * t236 + t126 * t83 - pkin(3) * t127 + t87 * t238 + t85 * t240 + t89 * t93 + t88 * t95 + m(6) * (t154 * t126 + t30 * t88 + t89 * t31) + mrSges(4,2) * t203 + (m(5) * t178 + t199 - t202) * qJ(4) + (Ifges(5,5) * t161 - Ifges(6,5) * t134 + Ifges(5,6) * t162 + Ifges(6,6) * t133) * t233 + (t31 * t133 + t30 * t134) * mrSges(6,3) + t178 * mrSges(5,3)) * qJD(3), t7 * qJD(3) + t204, t4 * qJD(3) + (-mrSges(6,1) * t27 - mrSges(6,2) * t26 + t194) * qJD(5) + t180; t2, t21 * qJD(3), t24 * qJD(4) + t15 * qJD(5) + (-mrSges(4,2) + t183) * t191 + t177 + ((t117 * t134 + t119 * t133) * mrSges(6,3) + (t83 + t205) * t164 + 0.2e1 * (t117 * t88 + t119 * t89 + t154 * t164) * t243 + 0.2e1 * (t165 * t182 - t227) * t244) * qJD(3), t24 * qJD(3), t15 * qJD(3) - t195 + t210; qJD(4) * t8 + qJD(5) * t3 - t181, qJD(4) * t25 + qJD(5) * t16 - t177, qJD(4) * t20 + qJD(5) * t14, t170, (-mrSges(6,1) * t89 - mrSges(6,2) * t88 + t193) * qJD(5) + t171; -t8 * qJD(3) + t195 - t204, -t25 * qJD(3), qJD(5) * t82 - t170, 0, t176; -qJD(3) * t3 - qJD(4) * t65 - t180, -t16 * qJD(3) - t210, -qJD(4) * t82 - t171, -t176, 0;];
Cq = t10;
