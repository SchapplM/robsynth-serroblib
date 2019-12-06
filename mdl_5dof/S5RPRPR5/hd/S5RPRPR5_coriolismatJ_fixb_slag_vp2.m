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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:55:57
% EndTime: 2019-12-05 17:56:06
% DurationCPUTime: 2.56s
% Computational Cost: add. (8093->268), mult. (17639->377), div. (0->0), fcn. (18764->8), ass. (0->147)
t202 = cos(pkin(9));
t182 = t202 * pkin(3);
t149 = t182 + pkin(4);
t159 = sin(qJ(5));
t161 = cos(qJ(5));
t156 = sin(pkin(9));
t229 = pkin(3) * t156;
t132 = t149 * t159 + t161 * t229;
t157 = sin(pkin(8));
t162 = cos(qJ(3));
t177 = t202 * t162;
t160 = sin(qJ(3));
t199 = t157 * t160;
t125 = t156 * t199 - t157 * t177;
t136 = -t156 * t162 - t160 * t202;
t126 = t136 * t157;
t83 = t125 * t161 - t126 * t159;
t206 = t132 * t83;
t131 = t149 * t161 - t159 * t229;
t176 = t125 * t159 + t161 * t126;
t207 = t131 * t176;
t265 = (t206 - t207) * mrSges(6,3);
t139 = pkin(3) * t199 + t157 * qJ(2);
t105 = -pkin(4) * t126 + t139;
t158 = cos(pkin(8));
t256 = t83 * mrSges(6,1);
t179 = t176 * mrSges(6,2) - t256;
t230 = Ifges(6,4) * t83;
t239 = t83 / 0.2e1;
t228 = pkin(7) * t125;
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
t40 = -pkin(4) * t158 + t228 + t55;
t227 = pkin(7) * t126;
t178 = t202 * t112;
t56 = t156 * t100 + t178;
t42 = t56 + t227;
t25 = -t159 * t42 + t161 * t40;
t251 = t176 / 0.2e1;
t26 = t159 * t40 + t161 * t42;
t264 = (-t176 * t25 + t26 * t83) * mrSges(6,3) + (Ifges(6,2) * t176 - Ifges(6,6) * t158 - t230) * t239 - (Ifges(6,1) * t176 + t230) * t83 / 0.2e1 + (-Ifges(6,5) * t158 + 0.2e1 * Ifges(6,4) * t176 + (-Ifges(6,1) + Ifges(6,2)) * t83) * t251 + t105 * t179;
t187 = -t256 / 0.2e1;
t135 = -t156 * t160 + t177;
t104 = t135 * t159 - t136 * t161;
t175 = t161 * t135 + t136 * t159;
t16 = -t104 * mrSges(6,1) - t175 * mrSges(6,2);
t262 = t16 * qJD(5);
t127 = t136 * t158;
t128 = t135 * t158;
t242 = m(5) * pkin(3);
t189 = t242 / 0.2e1;
t198 = t158 * t160;
t85 = t127 * t161 - t128 * t159;
t88 = t127 * t159 + t128 * t161;
t217 = t85 * mrSges(6,1) / 0.2e1 - t88 * mrSges(6,2) / 0.2e1;
t243 = m(6) / 0.2e1;
t261 = -(t131 * t85 + t132 * t88) * t243 - t127 * mrSges(5,1) / 0.2e1 + t128 * mrSges(5,2) / 0.2e1 - (t127 * t202 + t128 * t156) * t189 + mrSges(4,1) * t198 / 0.2e1 + mrSges(4,2) * t197 / 0.2e1 - t217;
t249 = Ifges(6,5) * t176;
t257 = Ifges(6,6) * t83;
t259 = t249 + t257;
t186 = t257 / 0.2e1 + t249 / 0.2e1;
t234 = t175 / 0.2e1;
t65 = mrSges(6,2) * t158 + mrSges(6,3) * t176;
t254 = t65 * t234;
t235 = -t175 / 0.2e1;
t184 = t176 * t235;
t190 = Ifges(5,5) * t126 + Ifges(5,6) * t125;
t208 = t126 * mrSges(5,3);
t114 = mrSges(5,2) * t158 + t208;
t209 = t125 * mrSges(5,3);
t115 = -mrSges(5,1) * t158 + t209;
t137 = t158 * mrSges(4,2) - mrSges(4,3) * t199;
t195 = t162 * t137;
t138 = -t158 * mrSges(4,1) - mrSges(4,3) * t194;
t196 = t160 * t138;
t66 = -mrSges(6,1) * t158 + t83 * mrSges(6,3);
t240 = t66 / 0.2e1;
t244 = m(5) / 0.2e1;
t185 = qJ(2) * t198;
t111 = t174 - t185;
t58 = -t111 * t156 - t178;
t47 = t58 - t227;
t59 = t202 * t111 - t107;
t48 = t59 + t228;
t31 = -t159 * t48 + t161 * t47;
t32 = t159 * t47 + t161 * t48;
t245 = ((t55 - t59) * t136 + (t56 + t58) * t135) * t244 + ((t26 + t31) * t175 + (-t25 + t32) * t104) * t243 - t104 * t240 + t254 + t135 * t114 / 0.2e1 + t136 * t115 / 0.2e1 - t196 / 0.2e1 + t195 / 0.2e1 + (-t136 * t125 / 0.2e1 - t135 * t126 / 0.2e1) * mrSges(5,3) + (t104 * t239 + t184) * mrSges(6,3);
t233 = -t104 / 0.2e1;
t232 = t157 / 0.2e1;
t231 = -t158 / 0.2e1;
t226 = t25 * mrSges(6,2);
t225 = t26 * mrSges(6,1);
t224 = t31 * mrSges(6,1);
t223 = t32 * mrSges(6,2);
t113 = pkin(3) * t194 - t125 * pkin(4);
t119 = t134 - t185;
t121 = t125 * mrSges(5,1);
t41 = -mrSges(6,1) * t176 - mrSges(6,2) * t83;
t89 = -mrSges(5,1) * t126 - mrSges(5,2) * t125;
t1 = t119 * t137 - t120 * t138 - t139 * t121 + t113 * t41 + t59 * t114 + t58 * t115 + t32 * t65 + t31 * t66 + m(5) * (t55 * t58 + t56 * t59) + m(6) * (t105 * t113 + t25 * t31 + t26 * t32) + (t139 * mrSges(5,2) - t55 * mrSges(5,3) + Ifges(5,4) * t126) * t126 + (t56 * mrSges(5,3) - Ifges(5,4) * t125 + (-Ifges(5,1) + Ifges(5,2)) * t126) * t125 + ((t119 * mrSges(4,3) + Ifges(4,5) * t158 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t160) * t157) * t160 + (-t120 * mrSges(4,3) + Ifges(4,6) * t158 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t162 + (-Ifges(4,1) + Ifges(4,2)) * t160) * t157 + (m(5) * t139 + t89) * pkin(3)) * t162) * t157 + (0.2e1 * t190 + t259) * t231 + t264;
t214 = t1 * qJD(1);
t2 = t231 * t259 + t25 * t65 - t26 * t66 + t264;
t205 = t2 * qJD(1);
t167 = (-t233 * t83 + t184) * mrSges(6,3) + t254 + t66 * t233;
t7 = t167 - t217;
t204 = t7 * qJD(1);
t154 = t157 ^ 2;
t152 = t154 * qJ(2);
t155 = t158 ^ 2;
t173 = t160 * mrSges(4,1) + t162 * mrSges(4,2);
t9 = t128 * t114 + t127 * t115 + t88 * t65 + t85 * t66 + (t195 - t196) * t158 + (t155 + t154) * mrSges(3,3) + (t157 * t173 + t41 + t89) * t157 + m(6) * (t105 * t157 + t25 * t85 + t26 * t88) + m(5) * (t127 * t55 + t128 * t56 + t139 * t157) + m(4) * (t152 + (-t119 * t160 + t120 * t162) * t158) + m(3) * (qJ(2) * t155 + t152);
t203 = t9 * qJD(1);
t10 = t126 * t114 + t125 * t115 + t176 * t65 + t83 * t66 + m(6) * (t176 * t26 + t25 * t83) + m(5) * (t125 * t55 + t126 * t56);
t201 = qJD(1) * t10;
t165 = (t131 * t83 + t132 * t176) * t243 + (t125 * t202 + t126 * t156) * t189;
t188 = m(5) * t232;
t170 = pkin(3) * t162 * t188 + t113 * t243;
t14 = t126 * mrSges(5,2) - t121 - t165 + t170 + t179;
t200 = t14 * qJD(1);
t166 = (t104 * t176 + t175 * t83) * t243 + (t125 * t135 - t126 * t136) * t244;
t19 = (-m(6) / 0.2e1 - m(5) / 0.2e1) * t157 + t166;
t193 = t19 * qJD(1);
t27 = 0.2e1 * t251 * mrSges(6,2) + 0.2e1 * t187;
t192 = t27 * qJD(1);
t37 = -mrSges(6,1) * t132 - mrSges(6,2) * t131;
t191 = t37 * qJD(5);
t4 = t245 + t261;
t172 = t4 * qJD(1);
t15 = (t235 + t234) * mrSges(6,2) + (t233 + t104 / 0.2e1) * mrSges(6,1);
t168 = -t131 * t65 / 0.2e1 + t132 * t240 - t186;
t5 = (-t206 / 0.2e1 + t207 / 0.2e1) * mrSges(6,3) + (-t32 / 0.2e1 + t25 / 0.2e1) * mrSges(6,2) + (t31 / 0.2e1 + t26 / 0.2e1) * mrSges(6,1) + t168 + t186;
t169 = t5 * qJD(1) - t15 * qJD(2) - t37 * qJD(3);
t28 = t187 + t256 / 0.2e1;
t20 = t165 + t170;
t18 = m(6) * t232 + t166 + t188;
t8 = t167 + t217;
t6 = -t226 / 0.2e1 - t225 / 0.2e1 - t223 / 0.2e1 + t224 / 0.2e1 - t168 + t186 + t265 / 0.2e1;
t3 = t245 - t261;
t11 = [qJD(2) * t9 + qJD(3) * t1 + qJD(4) * t10 + qJD(5) * t2, t203 + 0.2e1 * ((t104 * t88 + t175 * t85) * t243 + (t127 * t135 - t128 * t136) * t244) * qJD(2) + t3 * qJD(3) + t18 * qJD(4) + t8 * qJD(5), t214 + t3 * qJD(2) + (-Ifges(4,5) * t199 - Ifges(4,6) * t194 + m(6) * (t131 * t31 + t132 * t32) + t209 * t229 - t182 * t208 - t223 + t224 - t59 * mrSges(5,2) + t58 * mrSges(5,1) + (t156 * t59 + t202 * t58) * t242 - t119 * mrSges(4,2) - t120 * mrSges(4,1) + t190 + t259 + t265) * qJD(3) + t20 * qJD(4) + t6 * qJD(5), qJD(2) * t18 + qJD(3) * t20 + qJD(5) * t28 + t201, t205 + t8 * qJD(2) + t6 * qJD(3) + t28 * qJD(4) + (t259 - t225 - t226) * qJD(5); qJD(3) * t4 + qJD(4) * t19 + qJD(5) * t7 - t203, 0, (-t135 * mrSges(5,2) + t136 * mrSges(5,1) + (t135 * t156 + t136 * t202) * t242 + m(6) * (-t104 * t131 + t132 * t175) - t173 + t16) * qJD(3) + t262 + t172, t193, t16 * qJD(3) + t204 + t262; -qJD(2) * t4 - qJD(4) * t14 - qJD(5) * t5 - t214, qJD(5) * t15 - t172, t191, -t200, -t169 + t191; -qJD(2) * t19 + qJD(3) * t14 + qJD(5) * t27 - t201, -t193, t200, 0, t192; -qJD(2) * t7 + qJD(3) * t5 - qJD(4) * t27 - t205, -t15 * qJD(3) - t204, t169, -t192, 0;];
Cq = t11;
