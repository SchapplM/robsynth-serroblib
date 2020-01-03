% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:55
% EndTime: 2019-12-31 18:51:04
% DurationCPUTime: 3.91s
% Computational Cost: add. (3351->337), mult. (9063->441), div. (0->0), fcn. (6367->6), ass. (0->157)
t238 = Ifges(5,4) + Ifges(6,4);
t239 = Ifges(5,1) + Ifges(6,1);
t229 = Ifges(5,5) + Ifges(6,5);
t237 = Ifges(5,2) + Ifges(6,2);
t228 = Ifges(6,6) + Ifges(5,6);
t240 = Ifges(4,2) / 0.2e1;
t122 = sin(pkin(8));
t123 = cos(pkin(8));
t125 = sin(qJ(3));
t127 = cos(qJ(3));
t110 = t122 * t127 + t123 * t125;
t105 = t110 * qJD(1);
t124 = sin(qJ(4));
t126 = cos(qJ(4));
t90 = qJD(3) * t126 - t105 * t124;
t236 = t238 * t90;
t235 = t238 * t126;
t234 = t238 * t124;
t91 = qJD(3) * t124 + t105 * t126;
t233 = t238 * t91;
t109 = t122 * t125 - t127 * t123;
t104 = t109 * qJD(1);
t163 = -pkin(2) * t123 - pkin(1);
t113 = qJD(1) * t163 + qJD(2);
t61 = pkin(3) * t104 - pkin(7) * t105 + t113;
t190 = pkin(6) + qJ(2);
t114 = t190 * t122;
t111 = qJD(1) * t114;
t115 = t190 * t123;
t112 = qJD(1) * t115;
t84 = -t111 * t125 + t112 * t127;
t79 = qJD(3) * pkin(7) + t84;
t22 = -t124 * t79 + t126 * t61;
t23 = t124 * t61 + t126 * t79;
t137 = t124 * t23 + t126 * t22;
t15 = qJ(5) * t90 + t23;
t150 = mrSges(6,1) * t124 + mrSges(6,2) * t126;
t152 = mrSges(5,1) * t124 + mrSges(5,2) * t126;
t197 = t126 / 0.2e1;
t200 = -t124 / 0.2e1;
t223 = qJD(4) + t104;
t202 = -t223 / 0.2e1;
t204 = t91 / 0.2e1;
t207 = -t90 / 0.2e1;
t216 = t239 * t126 - t234;
t217 = -t237 * t124 + t235;
t218 = -t228 * t124 + t229 * t126;
t220 = t229 * t223 + t239 * t91 + t236;
t227 = t228 * t223 + t237 * t90 + t233;
t83 = -t111 * t127 - t125 * t112;
t78 = -qJD(3) * pkin(3) - t83;
t40 = -pkin(4) * t90 + qJD(5) + t78;
t14 = -qJ(5) * t91 + t22;
t9 = pkin(4) * t223 + t14;
t212 = t137 * mrSges(5,3) + (t15 * t124 + t9 * t126) * mrSges(6,3) - t150 * t40 - t152 * t78 + t217 * t207 - t216 * t204 + t218 * t202 - t227 * t200 - t220 * t197;
t231 = t105 * Ifges(4,1) / 0.2e1;
t232 = t113 * mrSges(4,2) - t83 * mrSges(4,3) - Ifges(4,4) * t104 + Ifges(4,5) * qJD(3) - t212 + t231;
t230 = t104 * t240;
t226 = t229 * t124 + t228 * t126;
t225 = t237 * t126 + t234;
t224 = t239 * t124 + t235;
t205 = -t91 / 0.2e1;
t106 = t109 * qJD(3);
t97 = qJD(1) * t106;
t53 = qJD(4) * t90 - t126 * t97;
t54 = -qJD(4) * t91 + t124 * t97;
t107 = t110 * qJD(3);
t98 = qJD(1) * t107;
t222 = t228 * t98 + t237 * t54 + t238 * t53;
t221 = t229 * t98 + t238 * t54 + t239 * t53;
t219 = -t127 * t114 - t115 * t125;
t168 = qJD(4) * t126;
t169 = qJD(4) * t124;
t133 = t109 * qJD(2);
t42 = -qJD(1) * t133 + t83 * qJD(3);
t67 = pkin(3) * t98 + pkin(7) * t97;
t5 = t124 * t67 + t126 * t42 + t61 * t168 - t169 * t79;
t6 = -qJD(4) * t23 - t124 * t42 + t126 * t67;
t215 = -t124 * t6 + t126 * t5;
t213 = (m(3) * qJ(2) + mrSges(3,3)) * (t122 ^ 2 + t123 ^ 2);
t164 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t165 = Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1;
t166 = Ifges(6,5) / 0.2e1 + Ifges(5,5) / 0.2e1;
t210 = t164 * t223 + t165 * t90 + t166 * t91 - t15 * mrSges(6,2) - t23 * mrSges(5,2) - t84 * mrSges(4,3) - Ifges(4,6) * qJD(3) - t105 * Ifges(4,4) + t230 + t113 * mrSges(4,1) + t22 * mrSges(5,1) + t9 * mrSges(6,1) - t228 * t207 - t229 * t205 - (Ifges(6,3) + Ifges(5,3)) * t202;
t209 = t53 / 0.2e1;
t208 = t54 / 0.2e1;
t203 = t98 / 0.2e1;
t134 = t110 * qJD(2);
t43 = qJD(1) * t134 + qJD(3) * t84;
t194 = t43 * t219;
t189 = -qJ(5) - pkin(7);
t80 = pkin(3) * t105 + pkin(7) * t104;
t27 = t124 * t80 + t126 * t83;
t82 = pkin(3) * t109 - pkin(7) * t110 + t163;
t87 = -t114 * t125 + t115 * t127;
t85 = t126 * t87;
t39 = t124 * t82 + t85;
t186 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t90 - mrSges(5,2) * t91 - t105 * mrSges(4,3);
t175 = qJ(5) * t126;
t174 = t104 * t124;
t173 = t110 * t124;
t62 = t219 * qJD(3) - t133;
t81 = pkin(3) * t107 + pkin(7) * t106;
t167 = t124 * t81 + t126 * t62 + t82 * t168;
t162 = t110 * t168;
t161 = t98 * mrSges(4,1) - t97 * mrSges(4,2);
t18 = -t54 * mrSges(6,1) + t53 * mrSges(6,2);
t160 = -t124 * t62 + t126 * t81;
t26 = -t124 * t83 + t126 * t80;
t38 = -t124 * t87 + t126 * t82;
t159 = qJD(4) * t189;
t1 = pkin(4) * t98 - qJ(5) * t53 - qJD(5) * t91 + t6;
t2 = qJ(5) * t54 + qJD(5) * t90 + t5;
t157 = -t1 * t126 - t124 * t2;
t156 = -t124 * t5 - t126 * t6;
t155 = t124 * t9 - t126 * t15;
t153 = mrSges(5,1) * t126 - mrSges(5,2) * t124;
t151 = mrSges(6,1) * t126 - mrSges(6,2) * t124;
t136 = t124 * t22 - t126 * t23;
t135 = qJ(5) * t106 - qJD(5) * t110;
t132 = t6 * mrSges(5,1) + t1 * mrSges(6,1) - t5 * mrSges(5,2) - t2 * mrSges(6,2);
t63 = qJD(3) * t87 + t134;
t119 = -pkin(4) * t126 - pkin(3);
t117 = t189 * t126;
t116 = t189 * t124;
t102 = -qJD(5) * t124 + t126 * t159;
t101 = qJD(5) * t126 + t124 * t159;
t96 = Ifges(5,3) * t98;
t95 = Ifges(6,3) * t98;
t92 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t104;
t64 = pkin(4) * t173 - t219;
t60 = mrSges(5,1) * t223 - mrSges(5,3) * t91;
t59 = mrSges(6,1) * t223 - mrSges(6,3) * t91;
t58 = -mrSges(5,2) * t223 + mrSges(5,3) * t90;
t57 = -mrSges(6,2) * t223 + mrSges(6,3) * t90;
t51 = Ifges(5,5) * t53;
t50 = Ifges(6,5) * t53;
t49 = Ifges(5,6) * t54;
t48 = Ifges(6,6) * t54;
t46 = -pkin(4) * t174 + t84;
t44 = -mrSges(6,1) * t90 + mrSges(6,2) * t91;
t37 = -mrSges(5,2) * t98 + mrSges(5,3) * t54;
t36 = -mrSges(6,2) * t98 + mrSges(6,3) * t54;
t35 = mrSges(5,1) * t98 - mrSges(5,3) * t53;
t34 = mrSges(6,1) * t98 - mrSges(6,3) * t53;
t25 = (-t106 * t124 + t162) * pkin(4) + t63;
t24 = -qJ(5) * t173 + t39;
t21 = qJ(5) * t174 + t27;
t20 = pkin(4) * t109 - t110 * t175 + t38;
t19 = -mrSges(5,1) * t54 + mrSges(5,2) * t53;
t17 = -pkin(4) * t54 + t43;
t16 = pkin(4) * t105 + t104 * t175 + t26;
t8 = -qJD(4) * t39 + t160;
t7 = -t169 * t87 + t167;
t4 = -qJ(5) * t162 + (-qJD(4) * t87 + t135) * t124 + t167;
t3 = pkin(4) * t107 + t135 * t126 + (-t85 + (qJ(5) * t110 - t82) * t124) * qJD(4) + t160;
t10 = [t62 * t92 - t219 * t19 + t4 * t57 + t7 * t58 + t3 * t59 + t8 * t60 + t64 * t18 + t25 * t44 + t20 * t34 + t24 * t36 + t38 * t35 + t39 * t37 + t163 * t161 - t186 * t63 + (t219 * t97 - t87 * t98) * mrSges(4,3) + m(5) * (t22 * t8 + t23 * t7 + t38 * t6 + t39 * t5 + t63 * t78 - t194) + m(4) * (t42 * t87 + t62 * t84 - t63 * t83 - t194) + m(6) * (t1 * t20 + t15 * t4 + t17 * t64 + t2 * t24 + t25 * t40 + t3 * t9) + 0.2e1 * t213 * qJD(2) * qJD(1) + (t230 + t210) * t107 - (t231 + t232) * t106 + (-t42 * mrSges(4,3) + t50 / 0.2e1 + t48 / 0.2e1 + t95 / 0.2e1 + t51 / 0.2e1 + t49 / 0.2e1 + t96 / 0.2e1 + Ifges(4,4) * t97 + t165 * t54 + t166 * t53 + (Ifges(4,2) + t164) * t98 + t132) * t109 + (t17 * t150 - Ifges(4,1) * t97 - Ifges(4,4) * t98 + (mrSges(4,3) + t152) * t43 + t157 * mrSges(6,3) + t156 * mrSges(5,3) + (mrSges(5,3) * t136 + mrSges(6,3) * t155 + t151 * t40 + t153 * t78 + t225 * t207 + t224 * t205 + t226 * t202 - t227 * t126 / 0.2e1) * qJD(4) + t216 * t209 + t217 * t208 + t218 * t203 + t221 * t197 + (t220 * qJD(4) + t222) * t200) * t110; t104 * t92 + (-t44 + t186) * t105 + (t34 + t35 + t223 * (t57 + t58)) * t126 + (t36 + t37 - t223 * (t59 + t60)) * t124 - m(4) * (-t104 * t84 - t105 * t83) + t161 - t213 * qJD(1) ^ 2 + (-t105 * t40 - t155 * t223 - t157) * m(6) + (-t105 * t78 - t136 * t223 - t156) * m(5); -t210 * t105 + m(6) * (t1 * t116 + t101 * t15 + t102 * t9 - t117 * t2 + t119 * t17) - Ifges(4,6) * t98 - t83 * t92 - Ifges(4,5) * t97 - t27 * t58 - t26 * t60 - t42 * mrSges(4,2) - t46 * t44 - pkin(3) * t19 + (-t212 + (m(6) * t40 + t44) * t124 * pkin(4)) * qJD(4) + t224 * t209 + t225 * t208 + t226 * t203 + (-t1 * t124 + t126 * t2) * mrSges(6,3) - ((t240 - Ifges(4,1) / 0.2e1) * t105 - t232) * t104 + t186 * t84 + (-pkin(3) * t43 - t22 * t26 - t23 * t27 - t78 * t84) * m(5) + t221 * t124 / 0.2e1 + t222 * t197 + (-t124 * t35 + t126 * t37 + m(5) * t215 + (-m(5) * t137 - t124 * t58 - t126 * t60) * qJD(4)) * pkin(7) + t215 * mrSges(5,3) - t17 * t151 + (-mrSges(4,1) - t153) * t43 + (-t16 + t102) * t59 + (-t21 + t101) * t57 - m(6) * (t15 * t21 + t16 * t9 + t40 * t46) + t116 * t34 - t117 * t36 + t119 * t18; t96 + t95 + t51 + t50 + t48 + t49 - t78 * (mrSges(5,1) * t91 + mrSges(5,2) * t90) - t40 * (mrSges(6,1) * t91 + mrSges(6,2) * t90) - t14 * t57 - t22 * t58 + t15 * t59 + t23 * t60 + (-t44 * t91 + t34) * pkin(4) + t132 + (-(t14 - t9) * t15 + (-t40 * t91 + t1) * pkin(4)) * m(6) + (t15 * t91 + t9 * t90) * mrSges(6,3) + (t22 * t90 + t23 * t91) * mrSges(5,3) + (t239 * t90 - t233) * t205 + t227 * t204 + (-t228 * t91 + t229 * t90) * t202 + (-t237 * t91 + t220 + t236) * t207; -t90 * t57 + t91 * t59 + 0.2e1 * (t17 / 0.2e1 + t15 * t207 + t9 * t204) * m(6) + t18;];
tauc = t10(:);
