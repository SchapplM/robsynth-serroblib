% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:11
% EndTime: 2019-12-31 17:20:15
% DurationCPUTime: 2.00s
% Computational Cost: add. (2102->286), mult. (4823->386), div. (0->0), fcn. (3651->4), ass. (0->137)
t252 = Ifges(5,4) + Ifges(4,5);
t254 = Ifges(4,6) - Ifges(5,6);
t159 = cos(qJ(2));
t156 = sin(qJ(3));
t157 = sin(qJ(2));
t198 = t156 * t157;
t190 = mrSges(4,3) * t198;
t115 = mrSges(4,2) * t159 - t190;
t147 = t159 * mrSges(5,3);
t122 = -mrSges(5,2) * t198 - t147;
t253 = t122 + t115;
t251 = Ifges(5,2) + Ifges(4,3);
t149 = Ifges(5,5) * t156;
t158 = cos(qJ(3));
t174 = Ifges(5,1) * t158 + t149;
t76 = -t159 * Ifges(5,4) + t157 * t174;
t215 = Ifges(4,4) * t156;
t175 = Ifges(4,1) * t158 - t215;
t78 = -t159 * Ifges(4,5) + t157 * t175;
t250 = t76 + t78;
t196 = t157 * t158;
t142 = Ifges(5,5) * t196;
t72 = -t159 * Ifges(5,6) + Ifges(5,3) * t198 + t142;
t249 = t72 + t142;
t248 = t252 * t156 + t254 * t158;
t187 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t225 = t157 / 0.2e1;
t218 = pkin(5) * t158;
t181 = -qJ(4) + t218;
t125 = -pkin(2) * t159 - t157 * pkin(6) - pkin(1);
t199 = t156 * t125;
t44 = t159 * t181 + t199;
t193 = t158 * t159;
t60 = pkin(5) * t193 + t199;
t246 = -t44 + t60;
t184 = pkin(5) * t156 + pkin(3);
t195 = t158 * t125;
t45 = t159 * t184 - t195;
t197 = t156 * t159;
t59 = -pkin(5) * t197 + t195;
t245 = t45 + t59;
t129 = -t158 * Ifges(5,3) + t149;
t152 = Ifges(4,4) * t158;
t173 = -Ifges(4,2) * t156 + t152;
t169 = pkin(3) * t158 + qJ(4) * t156;
t217 = pkin(6) * t159;
t137 = t157 * pkin(2) - t217;
t200 = t137 * t158;
t62 = pkin(5) * t198 + t200;
t114 = t156 * t137;
t63 = -pkin(5) * t196 + t114;
t244 = -t156 * t62 + t158 * t63;
t46 = -t157 * t181 + t114;
t47 = -t157 * t184 - t200;
t243 = t156 * t47 + t158 * t46;
t242 = -t156 * Ifges(4,1) - t152;
t123 = -pkin(2) - t169;
t177 = t158 * mrSges(5,1) + t156 * mrSges(5,3);
t241 = -m(5) * t123 + t177;
t146 = m(5) * qJ(4) + mrSges(5,3);
t240 = m(5) / 0.2e1;
t239 = -mrSges(4,1) / 0.2e1;
t238 = mrSges(4,2) / 0.2e1;
t237 = -mrSges(5,3) / 0.2e1;
t234 = t47 / 0.2e1;
t233 = t60 / 0.2e1;
t74 = -t159 * Ifges(4,6) + t157 * t173;
t232 = -t74 / 0.4e1;
t127 = pkin(3) * t156 - qJ(4) * t158;
t165 = pkin(5) + t127;
t80 = t165 * t157;
t231 = t80 / 0.2e1;
t230 = -qJ(4) / 0.2e1;
t229 = t123 / 0.2e1;
t228 = -t156 / 0.2e1;
t227 = t156 / 0.2e1;
t224 = -t158 / 0.2e1;
t223 = t158 / 0.2e1;
t221 = -t159 / 0.4e1;
t214 = Ifges(5,5) * t158;
t211 = t156 * Ifges(5,1);
t208 = t157 * mrSges(5,1);
t176 = mrSges(5,1) * t156 - mrSges(5,3) * t158;
t101 = t176 * t157;
t102 = t176 * t159;
t128 = mrSges(4,1) * t156 + mrSges(4,2) * t158;
t103 = t159 * t128;
t116 = -t157 * mrSges(4,2) - mrSges(4,3) * t197;
t117 = -mrSges(4,1) * t159 - mrSges(4,3) * t196;
t118 = mrSges(5,1) * t159 + mrSges(5,2) * t196;
t119 = t157 * mrSges(4,1) - mrSges(4,3) * t193;
t189 = mrSges(5,2) * t193;
t120 = t189 - t208;
t121 = -mrSges(5,2) * t197 + t157 * mrSges(5,3);
t188 = Ifges(4,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t179 = t188 * t159;
t170 = Ifges(5,3) * t156 + t214;
t73 = Ifges(5,6) * t157 + t159 * t170;
t75 = Ifges(4,6) * t157 + t159 * t173;
t77 = Ifges(5,4) * t157 + t159 * t174;
t79 = Ifges(4,5) * t157 + t159 * t175;
t81 = t165 * t159;
t3 = t81 * t101 + t80 * t102 + t63 * t115 + t60 * t116 + t62 * t117 + t47 * t118 + t59 * t119 + t45 * t120 + t44 * t121 + t46 * t122 + m(4) * (t59 * t62 + t60 * t63) + m(5) * (t44 * t46 + t45 * t47 + t80 * t81) + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t159 + (t76 / 0.2e1 + t78 / 0.2e1 - t179) * t158 + (t72 / 0.2e1 - t74 / 0.2e1 - t187 * t159) * t156) * t159 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t157 + pkin(5) * t103 + (t77 / 0.2e1 + t79 / 0.2e1 + t188 * t157) * t158 + (t73 / 0.2e1 - t75 / 0.2e1 + t187 * t157) * t156 + (-Ifges(3,2) + Ifges(3,1) + (m(4) * pkin(5) + t128) * pkin(5) - t251) * t159) * t157;
t204 = t3 * qJD(1);
t100 = t169 * t157;
t130 = t158 * Ifges(4,2) + t215;
t178 = t158 * mrSges(4,1) - t156 * mrSges(4,2);
t4 = (-m(5) * t80 - t101) * t100 + (-m(5) * t45 + t117 - t118) * t60 + (-m(5) * t44 - t190 - t253) * t59 + (-t80 * t177 + t74 * t223 + t60 * t158 * mrSges(4,3) + (t156 * t45 + t158 * t44) * mrSges(5,2) - t248 * t159 / 0.2e1 + t250 * t227 + t249 * t224 + (-pkin(5) * t178 + t211 * t223 + t242 * t224 + (t130 - t129) * t228) * t157) * t157;
t203 = t4 * qJD(1);
t10 = m(5) * (-t159 * t44 - t196 * t80) - t159 * t122 - t101 * t196;
t201 = qJD(1) * t10;
t194 = t158 * t177;
t192 = Ifges(5,4) * t158 + Ifges(5,6) * t156;
t191 = m(5) * t234;
t186 = t45 / 0.2e1 + t59 / 0.2e1;
t185 = t233 - t44 / 0.2e1;
t150 = Ifges(4,5) * t158;
t180 = -Ifges(4,6) * t156 + t150;
t131 = t211 - t214;
t164 = -Ifges(4,2) / 0.4e1 - Ifges(5,3) / 0.4e1 + Ifges(4,1) / 0.4e1 + Ifges(5,1) / 0.4e1;
t160 = pkin(5) * t128 / 0.2e1 + (pkin(2) * t239 + mrSges(5,1) * t229 + t149 / 0.4e1 - t130 / 0.4e1 + t129 / 0.4e1 + t164 * t158) * t158 + (pkin(2) * t238 + mrSges(5,3) * t229 + t242 / 0.4e1 - t131 / 0.4e1 - t152 / 0.4e1 - t164 * t156 + (-0.3e1 / 0.4e1 * Ifges(4,4) + Ifges(5,5) / 0.2e1) * t158) * t156 + (mrSges(5,2) + mrSges(4,3)) * pkin(6) * (-t158 ^ 2 / 0.2e1 - t156 ^ 2 / 0.2e1);
t161 = -m(5) * (-pkin(3) * t47 + qJ(4) * t46) / 0.2e1 + pkin(3) * t120 / 0.2e1 + t121 * t230 + t46 * t237 + mrSges(5,1) * t234 + t62 * t239 + t63 * t238;
t162 = (t100 * t123 + t127 * t80) * t240 - t100 * t177 / 0.2e1 + t127 * t101 / 0.2e1 + t192 * t221;
t2 = t150 * t221 + (t78 / 0.4e1 + t76 / 0.4e1 + t80 * t237 - t179 + t186 * mrSges(5,2) + (t245 * t240 + t118 / 0.2e1 - t117 / 0.2e1) * pkin(6)) * t158 + (mrSges(5,1) * t231 + t232 + t72 / 0.4e1 + t142 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(4,6) - Ifges(5,6) / 0.2e1) * t159 + t185 * mrSges(5,2) + (t246 * t240 - t122 / 0.2e1 - t115 / 0.2e1) * pkin(6)) * t156 + (-Ifges(5,2) / 0.2e1 - Ifges(4,3) / 0.2e1 + t160) * t157 + t161 + t162;
t5 = pkin(2) * t128 - t123 * t176 + t130 * t227 + t170 * t223 + t241 * t127 + (t129 + t175 + t174) * t228 + (t131 - t242 + t173) * t224;
t168 = -t2 * qJD(1) + t5 * qJD(2);
t28 = t241 * t156;
t163 = (-t156 * t80 + (-t123 * t157 - t217) * t158) * t240 + t101 * t228;
t8 = t189 + (-mrSges(5,1) / 0.2e1 - t194 / 0.2e1) * t157 + t191 - t163;
t167 = qJD(1) * t8 - qJD(2) * t28;
t13 = -t147 + 0.2e1 * (t199 / 0.4e1 - t60 / 0.4e1 + (t218 / 0.4e1 + t230) * t159) * m(5);
t166 = qJD(1) * t13 + qJD(3) * t146;
t124 = (m(5) * pkin(6) + mrSges(5,2)) * t158;
t12 = (t199 + (-0.2e1 * qJ(4) + t218) * t159) * t240 + m(5) * t233 + t122;
t9 = t194 * t225 - t208 / 0.2e1 + t191 + t163;
t1 = t180 * t221 + t176 * t231 + t162 - t161 + t156 * t232 + ((t156 * t246 + t158 * t245) * t240 + t118 * t223 + t117 * t224 + t253 * t228) * pkin(6) + (t156 * t185 + t158 * t186) * mrSges(5,2) + t160 * t157 + t249 * t156 / 0.4e1 + t251 * t225 + t250 * t158 / 0.4e1 + t187 * t197 + t252 * t193 / 0.2e1;
t6 = [qJD(2) * t3 - qJD(3) * t4 + qJD(4) * t10, t1 * qJD(3) + t9 * qJD(4) + t204 + (-pkin(2) * t103 + t123 * t102 + t75 * t223 + t73 * t224 + ((t116 + t121) * t158 + (-t119 + t120) * t156 + m(4) * t244 + m(5) * t243) * pkin(6) + (Ifges(3,5) + (t131 / 0.2e1 - t242 / 0.2e1) * t158 + (t129 / 0.2e1 - t130 / 0.2e1) * t156 + (-m(4) * pkin(2) - mrSges(3,1) - t178) * pkin(5)) * t159 - t241 * t81 + (t77 + t79) * t227 + t248 * t225 + (pkin(5) * mrSges(3,2) - Ifges(3,6)) * t157 + t244 * mrSges(4,3) + t243 * mrSges(5,2)) * qJD(2), t1 * qJD(2) + t12 * qJD(4) - t203 + (((-qJ(4) * mrSges(5,2) - t254) * t158 + (pkin(3) * mrSges(5,2) - t252) * t156) * t157 + (-m(5) * pkin(3) - mrSges(4,1) - mrSges(5,1)) * t60 + (-mrSges(4,2) + t146) * t59) * qJD(3), qJD(2) * t9 + qJD(3) * t12 + t201; qJD(3) * t2 - qJD(4) * t8 - t204, -qJD(3) * t5 + qJD(4) * t28, t124 * qJD(4) - t168 + (t180 + t192 + (-m(5) * t169 - t177 - t178) * pkin(6) - t169 * mrSges(5,2)) * qJD(3), qJD(3) * t124 - t167; -qJD(2) * t2 + qJD(4) * t13 + t203, t168, t146 * qJD(4), t166; qJD(2) * t8 - qJD(3) * t13 - t201, t167, -t166, 0;];
Cq = t6;
