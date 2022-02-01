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
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:24:59
% EndTime: 2022-01-23 09:25:04
% DurationCPUTime: 2.42s
% Computational Cost: add. (8093->273), mult. (17639->389), div. (0->0), fcn. (18764->8), ass. (0->148)
t158 = cos(pkin(9));
t149 = pkin(3) * t158 + pkin(4);
t160 = sin(qJ(5));
t162 = cos(qJ(5));
t156 = sin(pkin(9));
t226 = pkin(3) * t156;
t132 = t149 * t160 + t162 * t226;
t157 = sin(pkin(8));
t161 = sin(qJ(3));
t189 = t161 * t157;
t163 = cos(qJ(3));
t193 = t158 * t163;
t125 = t156 * t189 - t157 * t193;
t136 = -t156 * t163 - t158 * t161;
t126 = t136 * t157;
t83 = t125 * t162 - t126 * t160;
t203 = t132 * t83;
t131 = t149 * t162 - t160 * t226;
t176 = t125 * t160 + t126 * t162;
t205 = t131 * t176;
t260 = (t203 - t205) * mrSges(6,3);
t139 = pkin(3) * t189 + qJ(2) * t157;
t105 = -pkin(4) * t126 + t139;
t159 = cos(pkin(8));
t244 = t176 * mrSges(6,2);
t178 = -t83 * mrSges(6,1) + t244;
t227 = Ifges(6,4) * t83;
t233 = -t83 / 0.2e1;
t235 = t83 / 0.2e1;
t225 = pkin(7) * t125;
t140 = -pkin(2) * t159 - pkin(6) * t157 - pkin(1);
t134 = t163 * t140;
t195 = t157 * t163;
t174 = -qJ(4) * t195 + t134;
t100 = (-qJ(2) * t161 - pkin(3)) * t159 + t174;
t191 = t159 * t163;
t120 = qJ(2) * t191 + t140 * t161;
t112 = -qJ(4) * t189 + t120;
t107 = t156 * t112;
t55 = t100 * t158 - t107;
t40 = -pkin(4) * t159 + t225 + t55;
t224 = pkin(7) * t126;
t194 = t158 * t112;
t56 = t100 * t156 + t194;
t42 = t56 + t224;
t25 = -t160 * t42 + t162 * t40;
t26 = t160 * t40 + t162 * t42;
t259 = (-t176 * t25 + t26 * t83) * mrSges(6,3) + (0.2e1 * Ifges(6,4) * t176 - Ifges(6,5) * t159 + (-Ifges(6,1) + Ifges(6,2)) * t83) * t176 / 0.2e1 + (Ifges(6,2) * t176 - Ifges(6,6) * t159 - t227) * t235 + (Ifges(6,1) * t176 + t227) * t233 + t105 * t178;
t135 = -t156 * t161 + t193;
t104 = t135 * t160 - t136 * t162;
t175 = t135 * t162 + t136 * t160;
t16 = -t104 * mrSges(6,1) - t175 * mrSges(6,2);
t257 = t16 * qJD(5);
t127 = t136 * t159;
t128 = t135 * t159;
t85 = t127 * t162 - t128 * t160;
t88 = t127 * t160 + t128 * t162;
t214 = t85 * mrSges(6,1) / 0.2e1 - t88 * mrSges(6,2) / 0.2e1;
t238 = m(5) * pkin(3);
t240 = -m(6) / 0.2e1;
t256 = (t131 * t85 + t132 * t88) * t240 - t127 * mrSges(5,1) / 0.2e1 + t128 * mrSges(5,2) / 0.2e1 - (t127 * t158 + t128 * t156) * t238 / 0.2e1 - t214;
t246 = Ifges(6,5) * t176;
t252 = Ifges(6,6) * t83;
t254 = t246 + t252;
t180 = t252 / 0.2e1 + t246 / 0.2e1;
t181 = t244 / 0.2e1;
t231 = t175 / 0.2e1;
t65 = mrSges(6,2) * t159 + mrSges(6,3) * t176;
t249 = t65 * t231;
t232 = -t175 / 0.2e1;
t242 = t176 * t232;
t184 = Ifges(5,5) * t126 + Ifges(5,6) * t125;
t241 = m(5) / 0.2e1;
t239 = m(6) / 0.2e1;
t66 = -mrSges(6,1) * t159 + t83 * mrSges(6,3);
t236 = t66 / 0.2e1;
t230 = -t104 / 0.2e1;
t228 = -t159 / 0.2e1;
t223 = t25 * mrSges(6,2);
t222 = t26 * mrSges(6,1);
t192 = t159 * t161;
t179 = qJ(2) * t192;
t111 = t174 - t179;
t58 = -t111 * t156 - t194;
t47 = t58 - t224;
t59 = t111 * t158 - t107;
t48 = t59 + t225;
t31 = -t160 * t48 + t162 * t47;
t221 = t31 * mrSges(6,1);
t32 = t160 * t47 + t162 * t48;
t220 = t32 * mrSges(6,2);
t182 = pkin(3) * t195;
t113 = -pkin(4) * t125 + t182;
t114 = mrSges(5,2) * t159 + t126 * mrSges(5,3);
t115 = -mrSges(5,1) * t159 + t125 * mrSges(5,3);
t119 = t134 - t179;
t121 = t126 * mrSges(5,2);
t137 = t159 * mrSges(4,2) - mrSges(4,3) * t189;
t138 = -t159 * mrSges(4,1) - mrSges(4,3) * t195;
t41 = -mrSges(6,1) * t176 - mrSges(6,2) * t83;
t89 = -mrSges(5,1) * t126 - mrSges(5,2) * t125;
t1 = t119 * t137 - t120 * t138 + t139 * t121 + t113 * t41 + t59 * t114 + t58 * t115 + t32 * t65 + t31 * t66 + m(5) * (t55 * t58 + t56 * t59) + m(6) * (t105 * t113 + t25 * t31 + t26 * t32) + (-t55 * mrSges(5,3) + Ifges(5,4) * t126) * t126 + (-t139 * mrSges(5,1) + t56 * mrSges(5,3) - Ifges(5,4) * t125 + (Ifges(5,2) - Ifges(5,1)) * t126) * t125 + ((t119 * mrSges(4,3) + Ifges(4,5) * t159 + (-mrSges(4,2) * qJ(2) + Ifges(4,4) * t161) * t157) * t161 + (-t120 * mrSges(4,3) + Ifges(4,6) * t159 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t163 + (-Ifges(4,1) + Ifges(4,2)) * t161) * t157 + (m(5) * t139 + t89) * pkin(3)) * t163) * t157 + (0.2e1 * t184 + t254) * t228 + t259;
t211 = t1 * qJD(1);
t206 = t131 * t83;
t204 = t132 * t176;
t202 = t161 * mrSges(4,1);
t201 = t163 * mrSges(4,2);
t2 = t228 * t254 + t25 * t65 - t26 * t66 + t259;
t200 = t2 * qJD(1);
t167 = (-t230 * t83 + t242) * mrSges(6,3) + t249 + t66 * t230;
t7 = t167 - t214;
t199 = t7 * qJD(1);
t154 = t157 ^ 2;
t152 = t154 * qJ(2);
t155 = t159 ^ 2;
t173 = t201 + t202;
t188 = t163 * t137;
t190 = t161 * t138;
t9 = t128 * t114 + t127 * t115 + t88 * t65 + t85 * t66 + (t188 - t190) * t159 + (t155 + t154) * mrSges(3,3) + (t157 * t173 + t41 + t89) * t157 + m(6) * (t105 * t157 + t25 * t85 + t26 * t88) + m(5) * (t127 * t55 + t128 * t56 + t139 * t157) + m(4) * (t152 + (-t119 * t161 + t120 * t163) * t159) + m(3) * (qJ(2) * t155 + t152);
t198 = t9 * qJD(1);
t10 = t126 * t114 + t125 * t115 + t176 * t65 + t83 * t66 + m(6) * (t176 * t26 + t25 * t83) + m(5) * (t125 * t55 + t126 * t56);
t197 = qJD(1) * t10;
t183 = t238 / 0.2e1;
t169 = (t125 * t158 + t126 * t156) * t183;
t177 = m(5) * t182;
t14 = t125 * mrSges(5,1) - t121 - t177 / 0.2e1 + 0.2e1 * (t206 / 0.4e1 + t204 / 0.4e1 - t113 / 0.4e1) * m(6) + t169 - t178;
t196 = t14 * qJD(1);
t166 = (t104 * t176 + t175 * t83) * t239 + (t125 * t135 - t126 * t136) * t241;
t19 = (t240 - m(5) / 0.2e1) * t157 + t166;
t187 = t19 * qJD(1);
t27 = 0.2e1 * mrSges(6,1) * t233 + 0.2e1 * t181;
t186 = t27 * qJD(1);
t37 = -mrSges(6,1) * t132 - mrSges(6,2) * t131;
t185 = t37 * qJD(5);
t164 = (t104 * t235 + t242) * mrSges(6,3) + (-t136 * t125 / 0.2e1 - t135 * t126 / 0.2e1) * mrSges(5,3) + ((t55 - t59) * t136 + (t56 + t58) * t135) * t241 + ((t26 + t31) * t175 + (-t25 + t32) * t104) * t239 - t104 * t236 + t249 + t135 * t114 / 0.2e1 + t136 * t115 / 0.2e1 - t190 / 0.2e1 + t188 / 0.2e1;
t4 = t164 + (t201 / 0.2e1 + t202 / 0.2e1) * t159 + t256;
t172 = t4 * qJD(1);
t15 = (t232 + t231) * mrSges(6,2) + (t230 + t104 / 0.2e1) * mrSges(6,1);
t168 = -t131 * t65 / 0.2e1 + t132 * t236 - t180;
t5 = (-t203 / 0.2e1 + t205 / 0.2e1) * mrSges(6,3) + (-t32 / 0.2e1 + t25 / 0.2e1) * mrSges(6,2) + (t31 / 0.2e1 + t26 / 0.2e1) * mrSges(6,1) + t168 + t180;
t170 = qJD(1) * t5 - qJD(2) * t15 - qJD(3) * t37;
t28 = t181 - t244 / 0.2e1;
t20 = t177 / 0.2e1 + t169 + (t113 + t204 + t206) * t239;
t18 = t166 + (m(5) + m(6)) * t157 / 0.2e1;
t8 = t167 + t214;
t6 = -t223 / 0.2e1 - t222 / 0.2e1 - t220 / 0.2e1 + t221 / 0.2e1 - t168 + t180 + t260 / 0.2e1;
t3 = t164 - mrSges(4,1) * t192 / 0.2e1 - mrSges(4,2) * t191 / 0.2e1 - t256;
t11 = [qJD(2) * t9 + qJD(3) * t1 + qJD(4) * t10 + qJD(5) * t2, t198 + 0.2e1 * ((t104 * t88 + t175 * t85) * t239 + (t127 * t135 - t128 * t136) * t241) * qJD(2) + t3 * qJD(3) + t18 * qJD(4) + t8 * qJD(5), t3 * qJD(2) + t20 * qJD(4) + t6 * qJD(5) + t211 + (-Ifges(4,5) * t189 - Ifges(4,6) * t195 + m(6) * (t131 * t31 + t132 * t32) - t220 + t221 + t58 * mrSges(5,1) - t59 * mrSges(5,2) - t119 * mrSges(4,2) - t120 * mrSges(4,1) + t184 + t254 + (m(5) * (t156 * t59 + t158 * t58) + (t125 * t156 - t126 * t158) * mrSges(5,3)) * pkin(3) + t260) * qJD(3), qJD(2) * t18 + qJD(3) * t20 + qJD(5) * t28 + t197, t200 + t8 * qJD(2) + t6 * qJD(3) + t28 * qJD(4) + (t254 - t222 - t223) * qJD(5); qJD(3) * t4 + qJD(4) * t19 + qJD(5) * t7 - t198, 0, (t136 * mrSges(5,1) - t135 * mrSges(5,2) + t16 - t173) * qJD(3) + t257 + 0.2e1 * ((-t104 * t131 + t132 * t175) * t239 + (t135 * t156 + t136 * t158) * t183) * qJD(3) + t172, t187, t16 * qJD(3) + t199 + t257; -qJD(2) * t4 + qJD(4) * t14 - qJD(5) * t5 - t211, qJD(5) * t15 - t172, t185, t196, -t170 + t185; -qJD(2) * t19 - qJD(3) * t14 + qJD(5) * t27 - t197, -t187, -t196, 0, t186; -qJD(2) * t7 + qJD(3) * t5 - qJD(4) * t27 - t200, -qJD(3) * t15 - t199, t170, -t186, 0;];
Cq = t11;
