% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR5
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:32
% EndTime: 2020-01-03 11:41:42
% DurationCPUTime: 2.62s
% Computational Cost: add. (8093->268), mult. (17639->376), div. (0->0), fcn. (18764->8), ass. (0->146)
t156 = sin(pkin(9));
t160 = sin(qJ(3));
t162 = cos(qJ(3));
t202 = cos(pkin(9));
t136 = -t156 * t162 - t160 * t202;
t157 = sin(pkin(8));
t126 = t136 * t157;
t199 = t157 * t160;
t139 = pkin(3) * t199 + t157 * qJ(2);
t105 = -pkin(4) * t126 + t139;
t158 = cos(pkin(8));
t177 = t202 * t162;
t125 = t156 * t199 - t157 * t177;
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t176 = t125 * t159 + t161 * t126;
t83 = t125 * t161 - t126 * t159;
t254 = t83 * mrSges(6,1);
t179 = t176 * mrSges(6,2) - t254;
t228 = Ifges(6,4) * t83;
t238 = t83 / 0.2e1;
t226 = pkin(7) * t125;
t140 = -pkin(2) * t158 - pkin(6) * t157 - pkin(1);
t134 = t162 * t140;
t194 = t162 * t157;
t174 = -qJ(4) * t194 + t134;
t100 = (-qJ(2) * t160 - pkin(3)) * t158 + t174;
t197 = t158 * t162;
t120 = qJ(2) * t197 + t160 * t140;
t112 = -qJ(4) * t199 + t120;
t107 = t156 * t112;
t55 = t202 * t100 - t107;
t40 = -pkin(4) * t158 + t226 + t55;
t225 = pkin(7) * t126;
t178 = t202 * t112;
t56 = t156 * t100 + t178;
t42 = t56 + t225;
t25 = -t159 * t42 + t161 * t40;
t250 = t176 / 0.2e1;
t26 = t159 * t40 + t161 * t42;
t262 = (-t176 * t25 + t26 * t83) * mrSges(6,3) + (Ifges(6,2) * t176 - Ifges(6,6) * t158 - t228) * t238 - (Ifges(6,1) * t176 + t228) * t83 / 0.2e1 + (-Ifges(6,5) * t158 + 0.2e1 * Ifges(6,4) * t176 + (-Ifges(6,1) + Ifges(6,2)) * t83) * t250 + t105 * t179;
t187 = -t254 / 0.2e1;
t135 = -t156 * t160 + t177;
t104 = t135 * t159 - t136 * t161;
t232 = -t104 / 0.2e1;
t66 = -mrSges(6,1) * t158 + t83 * mrSges(6,3);
t261 = t66 * t232;
t175 = t161 * t135 + t136 * t159;
t16 = -t104 * mrSges(6,1) - t175 * mrSges(6,2);
t259 = t16 * qJD(5);
t127 = t136 * t158;
t128 = t135 * t158;
t182 = t202 * pkin(3);
t149 = t182 + pkin(4);
t227 = pkin(3) * t156;
t131 = t149 * t161 - t159 * t227;
t132 = t149 * t159 + t161 * t227;
t240 = m(5) * pkin(3);
t189 = t240 / 0.2e1;
t198 = t158 * t160;
t85 = t127 * t161 - t128 * t159;
t88 = t127 * t159 + t128 * t161;
t215 = t85 * mrSges(6,1) / 0.2e1 - t88 * mrSges(6,2) / 0.2e1;
t241 = m(6) / 0.2e1;
t258 = -(t131 * t85 + t132 * t88) * t241 - t127 * mrSges(5,1) / 0.2e1 + t128 * mrSges(5,2) / 0.2e1 - (t127 * t202 + t128 * t156) * t189 + mrSges(4,1) * t198 / 0.2e1 + mrSges(4,2) * t197 / 0.2e1 - t215;
t248 = Ifges(6,5) * t176;
t255 = Ifges(6,6) * t83;
t257 = t248 + t255;
t186 = -t255 / 0.2e1 - t248 / 0.2e1;
t244 = t131 * t176;
t234 = -t175 / 0.2e1;
t184 = t176 * t234;
t190 = Ifges(5,5) * t126 + Ifges(5,6) * t125;
t206 = t126 * mrSges(5,3);
t114 = mrSges(5,2) * t158 + t206;
t207 = t125 * mrSges(5,3);
t115 = -mrSges(5,1) * t158 + t207;
t137 = t158 * mrSges(4,2) - mrSges(4,3) * t199;
t195 = t162 * t137;
t138 = -t158 * mrSges(4,1) - mrSges(4,3) * t194;
t196 = t160 * t138;
t233 = t175 / 0.2e1;
t242 = m(5) / 0.2e1;
t185 = qJ(2) * t198;
t111 = t174 - t185;
t58 = -t111 * t156 - t178;
t47 = t58 - t225;
t59 = t202 * t111 - t107;
t48 = t59 + t226;
t31 = -t159 * t48 + t161 * t47;
t32 = t159 * t47 + t161 * t48;
t65 = mrSges(6,2) * t158 + mrSges(6,3) * t176;
t243 = ((t55 - t59) * t136 + (t56 + t58) * t135) * t242 + ((t26 + t31) * t175 + (-t25 + t32) * t104) * t241 + t261 + t65 * t233 + t135 * t114 / 0.2e1 + t136 * t115 / 0.2e1 - t196 / 0.2e1 + t195 / 0.2e1 + (-t136 * t125 / 0.2e1 - t135 * t126 / 0.2e1) * mrSges(5,3) + (t104 * t238 + t184) * mrSges(6,3);
t239 = t65 / 0.2e1;
t231 = -t132 / 0.2e1;
t230 = t157 / 0.2e1;
t229 = -t158 / 0.2e1;
t224 = t25 * mrSges(6,2);
t223 = t26 * mrSges(6,1);
t222 = t31 * mrSges(6,1);
t221 = t32 * mrSges(6,2);
t113 = pkin(3) * t194 - t125 * pkin(4);
t119 = t134 - t185;
t121 = t125 * mrSges(5,1);
t41 = -mrSges(6,1) * t176 - mrSges(6,2) * t83;
t89 = -mrSges(5,1) * t126 - mrSges(5,2) * t125;
t1 = t119 * t137 - t120 * t138 - t139 * t121 + t59 * t114 + t58 * t115 + t113 * t41 + t32 * t65 + t31 * t66 + m(6) * (t105 * t113 + t25 * t31 + t26 * t32) + m(5) * (t55 * t58 + t56 * t59) + (t139 * mrSges(5,2) - t55 * mrSges(5,3) + Ifges(5,4) * t126) * t126 + (t56 * mrSges(5,3) - Ifges(5,4) * t125 + (-Ifges(5,1) + Ifges(5,2)) * t126) * t125 + ((t119 * mrSges(4,3) + Ifges(4,5) * t158 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t160) * t157) * t160 + (-t120 * mrSges(4,3) + Ifges(4,6) * t158 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t162 + (-Ifges(4,1) + Ifges(4,2)) * t160) * t157 + (m(5) * t139 + t89) * pkin(3)) * t162) * t157 + (0.2e1 * t190 + t257) * t229 + t262;
t212 = t1 * qJD(1);
t2 = t229 * t257 + t25 * t65 - t26 * t66 + t262;
t205 = t2 * qJD(1);
t168 = (-t232 * t83 + t184) * mrSges(6,3) + t175 * t239 + t261;
t7 = t168 - t215;
t204 = t7 * qJD(1);
t154 = t157 ^ 2;
t152 = t154 * qJ(2);
t155 = t158 ^ 2;
t173 = t160 * mrSges(4,1) + t162 * mrSges(4,2);
t9 = t128 * t114 + t127 * t115 + t88 * t65 + t85 * t66 + (t195 - t196) * t158 + (t155 + t154) * mrSges(3,3) + (t157 * t173 + t41 + t89) * t157 + m(6) * (t105 * t157 + t25 * t85 + t26 * t88) + m(5) * (t127 * t55 + t128 * t56 + t139 * t157) + m(4) * (t152 + (-t119 * t160 + t120 * t162) * t158) + m(3) * (qJ(2) * t155 + t152);
t203 = t9 * qJD(1);
t10 = t126 * t114 + t125 * t115 + t176 * t65 + t83 * t66 + m(6) * (t176 * t26 + t25 * t83) + m(5) * (t125 * t55 + t126 * t56);
t201 = qJD(1) * t10;
t165 = (t131 * t83 + t132 * t176) * t241 + (t125 * t202 + t126 * t156) * t189;
t188 = m(5) * t230;
t170 = pkin(3) * t162 * t188 + t113 * t241;
t14 = t126 * mrSges(5,2) - t121 - t165 + t170 + t179;
t200 = t14 * qJD(1);
t167 = (t104 * t176 + t175 * t83) * t241 + (t125 * t135 - t126 * t136) * t242;
t19 = (-m(6) / 0.2e1 - m(5) / 0.2e1) * t157 + t167;
t193 = t19 * qJD(1);
t27 = 0.2e1 * t250 * mrSges(6,2) + 0.2e1 * t187;
t192 = t27 * qJD(1);
t37 = -mrSges(6,1) * t132 - mrSges(6,2) * t131;
t191 = t37 * qJD(5);
t4 = t243 + t258;
t172 = t4 * qJD(1);
t15 = (t234 + t233) * mrSges(6,2) + (t232 + t104 / 0.2e1) * mrSges(6,1);
t166 = (-t83 * t231 - t244 / 0.2e1) * mrSges(6,3) + t131 * t239 + t66 * t231 - t186;
t6 = (-t25 / 0.2e1 + t32 / 0.2e1) * mrSges(6,2) + (-t26 / 0.2e1 - t31 / 0.2e1) * mrSges(6,1) + t166 + t186;
t169 = t6 * qJD(1) + t15 * qJD(2) + t37 * qJD(3);
t28 = t187 + t254 / 0.2e1;
t20 = t165 + t170;
t18 = m(6) * t230 + t167 + t188;
t8 = t168 + t215;
t5 = -t224 / 0.2e1 - t223 / 0.2e1 - t221 / 0.2e1 + t222 / 0.2e1 + t166 - t186;
t3 = t243 - t258;
t11 = [qJD(2) * t9 + qJD(3) * t1 + qJD(4) * t10 + qJD(5) * t2, t203 + 0.2e1 * ((t104 * t88 + t175 * t85) * t241 + (t127 * t135 - t128 * t136) * t242) * qJD(2) + t3 * qJD(3) + t18 * qJD(4) + t8 * qJD(5), t212 + t3 * qJD(2) + (-Ifges(4,5) * t199 - Ifges(4,6) * t194 + t207 * t227 - t182 * t206 + m(6) * (t131 * t31 + t132 * t32) - t59 * mrSges(5,2) + t58 * mrSges(5,1) - t221 + t222 + (t156 * t59 + t202 * t58) * t240 - t119 * mrSges(4,2) - t120 * mrSges(4,1) + t190 + t257 + (t132 * t83 - t244) * mrSges(6,3)) * qJD(3) + t20 * qJD(4) + t5 * qJD(5), qJD(2) * t18 + qJD(3) * t20 + qJD(5) * t28 + t201, t205 + t8 * qJD(2) + t5 * qJD(3) + t28 * qJD(4) + (t257 - t223 - t224) * qJD(5); qJD(3) * t4 + qJD(4) * t19 + qJD(5) * t7 - t203, 0, (-t135 * mrSges(5,2) + t136 * mrSges(5,1) + (t135 * t156 + t136 * t202) * t240 + m(6) * (-t104 * t131 + t132 * t175) - t173 + t16) * qJD(3) + t259 + t172, t193, t16 * qJD(3) + t204 + t259; -qJD(2) * t4 - qJD(4) * t14 + qJD(5) * t6 - t212, qJD(5) * t15 - t172, t191, -t200, t169 + t191; -qJD(2) * t19 + qJD(3) * t14 + qJD(5) * t27 - t201, -t193, t200, 0, t192; -qJD(2) * t7 - qJD(3) * t6 - qJD(4) * t27 - t205, -t15 * qJD(3) - t204, -t169, -t192, 0;];
Cq = t11;
