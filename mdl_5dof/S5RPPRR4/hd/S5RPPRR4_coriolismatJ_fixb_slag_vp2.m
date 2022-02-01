% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:19
% EndTime: 2022-01-23 09:16:23
% DurationCPUTime: 1.84s
% Computational Cost: add. (7285->225), mult. (15891->318), div. (0->0), fcn. (17125->8), ass. (0->127)
t154 = sin(pkin(9));
t156 = cos(pkin(9));
t159 = sin(qJ(4));
t161 = cos(qJ(4));
t138 = t154 * t161 + t156 * t159;
t155 = sin(pkin(8));
t129 = t138 * t155;
t186 = t154 * t155;
t179 = pkin(3) * t186 + qJ(2) * t155;
t171 = -pkin(4) * t129 - t179;
t137 = -t154 * t159 + t156 * t161;
t131 = t137 * t155;
t158 = sin(qJ(5));
t160 = cos(qJ(5));
t103 = -t129 * t158 + t131 * t160;
t173 = -t129 * t160 - t131 * t158;
t237 = t173 * mrSges(6,2);
t175 = t103 * mrSges(6,1) + t237;
t250 = t171 * t175;
t118 = t137 * t158 + t138 * t160;
t172 = t137 * t160 - t138 * t158;
t28 = -t118 * mrSges(6,1) - t172 * mrSges(6,2);
t249 = t28 * qJD(5);
t157 = cos(pkin(8));
t130 = t138 * t157;
t132 = t137 * t157;
t101 = -t130 * t160 - t132 * t158;
t104 = -t130 * t158 + t132 * t160;
t229 = m(6) * pkin(4);
t178 = t229 / 0.2e1;
t212 = t101 * mrSges(6,1) / 0.2e1 - t104 * mrSges(6,2) / 0.2e1;
t248 = t130 * mrSges(5,1) / 0.2e1 + t132 * mrSges(5,2) / 0.2e1 - (t101 * t160 + t104 * t158) * t178 - t212;
t247 = Ifges(6,5) * t173 - Ifges(6,6) * t103;
t224 = -t118 / 0.2e1;
t226 = -t172 / 0.2e1;
t204 = t103 * mrSges(6,3);
t84 = -mrSges(6,1) * t157 - t204;
t246 = t84 * t224 + (t103 * t224 + t173 * t226) * mrSges(6,3);
t227 = t103 / 0.2e1;
t176 = t237 / 0.2e1;
t241 = t137 / 0.2e1;
t240 = Ifges(6,4) * t103;
t236 = t160 * t173;
t143 = -pkin(2) * t157 - qJ(3) * t155 - pkin(1);
t136 = t156 * t143;
t185 = t155 * t156;
t114 = -pkin(6) * t185 + t136 + (-qJ(2) * t154 - pkin(3)) * t157;
t190 = qJ(2) * t157;
t123 = t143 * t154 + t156 * t190;
t119 = -pkin(6) * t186 + t123;
t74 = t114 * t159 + t119 * t161;
t63 = -pkin(7) * t129 + t74;
t196 = t158 * t63;
t73 = t114 * t161 - t119 * t159;
t62 = -pkin(7) * t131 + t73;
t59 = -pkin(4) * t157 + t62;
t32 = t160 * t59 - t196;
t234 = t32 * mrSges(6,3) - Ifges(6,4) * t173;
t233 = t154 ^ 2;
t232 = t156 ^ 2;
t231 = m(5) / 0.2e1;
t230 = m(6) / 0.2e1;
t83 = mrSges(6,2) * t157 + mrSges(6,3) * t173;
t228 = t83 / 0.2e1;
t225 = t172 / 0.2e1;
t223 = -t138 / 0.2e1;
t222 = t155 / 0.2e1;
t221 = -t158 / 0.2e1;
t219 = t131 * pkin(4);
t218 = t32 * mrSges(6,2);
t195 = t160 * t63;
t33 = t158 * t59 + t195;
t216 = t33 * mrSges(6,1);
t36 = -t158 * t62 - t195;
t215 = t36 * mrSges(6,1);
t37 = t160 * t62 - t196;
t214 = t37 * mrSges(6,2);
t198 = t129 * mrSges(5,3);
t120 = mrSges(5,2) * t157 - t198;
t197 = t131 * mrSges(5,3);
t121 = -mrSges(5,1) * t157 - t197;
t167 = -Ifges(5,5) * t129 - Ifges(5,6) * t131 + t247;
t174 = t131 * mrSges(5,1) - mrSges(5,2) * t129;
t60 = -mrSges(6,1) * t173 + mrSges(6,2) * t103;
t1 = -t37 * t83 - t36 * t84 - t60 * t219 + t250 - m(6) * (-t171 * t219 + t32 * t36 + t33 * t37) - t179 * t174 + (-t129 ^ 2 + t131 ^ 2) * Ifges(5,4) + t167 * t157 + (-Ifges(6,1) * t103 + t234) * t173 - (-mrSges(6,3) * t33 - Ifges(6,2) * t173 - t240) * t103 + (t121 + t197) * t74 + (-t198 - t120) * t73 - (-Ifges(5,1) + Ifges(5,2)) * t129 * t131;
t208 = t1 * qJD(1);
t194 = t160 * t103;
t2 = -t32 * t83 + t250 + t103 * (-Ifges(6,6) * t157 + t240) + (Ifges(6,5) * t157 + (-Ifges(6,1) + Ifges(6,2)) * t103 + t234) * t173 + (t204 + t84) * t33;
t193 = t2 * qJD(1);
t165 = t172 * t228 + t246;
t7 = t165 - t212;
t192 = t7 * qJD(1);
t122 = -t154 * t190 + t136;
t140 = mrSges(4,2) * t157 - mrSges(4,3) * t186;
t141 = -mrSges(4,1) * t157 - mrSges(4,3) * t185;
t152 = t155 ^ 2;
t150 = t152 * qJ(2);
t153 = t157 ^ 2;
t9 = t101 * t84 + t104 * t83 + t132 * t120 - t130 * t121 + (t156 * t140 - t154 * t141) * t157 + (t153 + t152) * mrSges(3,3) + (t129 * mrSges(5,1) + t131 * mrSges(5,2) + t60 + (mrSges(4,1) * t154 + mrSges(4,2) * t156) * t155) * t155 + m(6) * (t32 * t101 + t33 * t104 - t155 * t171) + m(5) * (-t73 * t130 + t74 * t132 + t155 * t179) + m(4) * (t150 + (-t122 * t154 + t123 * t156) * t157) + m(3) * (qJ(2) * t153 + t150);
t191 = t9 * qJD(1);
t10 = t173 * t83 - t129 * t120 - t131 * t121 - t103 * t84 + m(6) * (-t103 * t32 + t173 * t33) + m(5) * (-t129 * t74 - t131 * t73) + (m(4) * (-t122 * t156 - t123 * t154) - t154 * t140 - t156 * t141) * t155;
t189 = t10 * qJD(1);
t188 = t173 * t158;
t164 = (-t103 * t172 + t118 * t173) * t230 + (-t129 * t138 - t131 * t137) * t231;
t19 = (-m(6) / 0.2e1 - m(5) / 0.2e1 + (-t233 / 0.2e1 - t232 / 0.2e1 - 0.1e1 / 0.2e1) * m(4)) * t155 + t164;
t182 = t19 * qJD(1);
t20 = (t188 / 0.2e1 - t194 / 0.2e1 - t131 / 0.2e1) * t229 - t174 - t175;
t181 = t20 * qJD(1);
t38 = 0.2e1 * mrSges(6,1) * t227 + 0.2e1 * t176;
t180 = t38 * qJD(1);
t177 = m(4) * t222;
t162 = (t129 * t241 + t131 * t223) * mrSges(5,3) + ((t33 + t36) * t172 + (-t32 + t37) * t118) * t230 + t83 * t225 + t120 * t241 + t121 * t223 + t246;
t3 = t162 + t248;
t170 = t3 * qJD(1);
t142 = (mrSges(6,1) * t158 + mrSges(6,2) * t160) * pkin(4);
t27 = (t226 + t225) * mrSges(6,2) + (t224 + t118 / 0.2e1) * mrSges(6,1);
t163 = (t160 * t228 + t84 * t221 + (t103 * t221 - t236 / 0.2e1) * mrSges(6,3)) * pkin(4);
t6 = (-t103 / 0.2e1 + t227) * Ifges(6,6) + (-t32 / 0.2e1 + t37 / 0.2e1) * mrSges(6,2) + (-t33 / 0.2e1 - t36 / 0.2e1) * mrSges(6,1) + t163;
t168 = -qJD(1) * t6 - qJD(2) * t27 + qJD(4) * t142;
t139 = t142 * qJD(5);
t52 = (t131 + t188 - t194) * t178;
t39 = t176 - t237 / 0.2e1;
t18 = t177 + (-t232 - t233) * t177 + t164 + (m(6) + m(5)) * t222;
t8 = t165 + t212;
t5 = -t218 / 0.2e1 - t216 / 0.2e1 - t214 / 0.2e1 + t215 / 0.2e1 + t163 + t247;
t4 = t162 - t248;
t11 = [qJD(2) * t9 + qJD(3) * t10 - qJD(4) * t1 - qJD(5) * t2, t191 + 0.2e1 * ((t101 * t172 + t104 * t118) * t230 + (-t130 * t137 + t132 * t138) * t231) * qJD(2) + t18 * qJD(3) + t4 * qJD(4) + t8 * qJD(5), qJD(2) * t18 + qJD(4) * t52 + qJD(5) * t39 + t189, t4 * qJD(2) + t52 * qJD(3) + t5 * qJD(5) - t208 + (-t74 * mrSges(5,1) - t73 * mrSges(5,2) + t167 - t214 + t215 + (m(6) * (t158 * t37 + t160 * t36) + (-t158 * t103 - t236) * mrSges(6,3)) * pkin(4)) * qJD(4), -t193 + t8 * qJD(2) + t39 * qJD(3) + t5 * qJD(4) + (-t216 - t218 + t247) * qJD(5); qJD(3) * t19 + qJD(4) * t3 + qJD(5) * t7 - t191, 0, t182, (-t138 * mrSges(5,1) - t137 * mrSges(5,2) + (-t118 * t160 + t158 * t172) * t229 + t28) * qJD(4) + t249 + t170, t28 * qJD(4) + t192 + t249; -qJD(2) * t19 - qJD(4) * t20 + qJD(5) * t38 - t189, -t182, 0, -t181, t180; -qJD(2) * t3 + qJD(3) * t20 + qJD(5) * t6 + t208, qJD(5) * t27 - t170, t181, -t139, -t139 - t168; -qJD(2) * t7 - qJD(3) * t38 - qJD(4) * t6 + t193, -qJD(4) * t27 - t192, -t180, t168, 0;];
Cq = t11;
