% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP11
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
% Datum: 2019-12-31 18:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:53:27
% EndTime: 2019-12-31 18:53:36
% DurationCPUTime: 3.72s
% Computational Cost: add. (3349->345), mult. (9018->446), div. (0->0), fcn. (6291->6), ass. (0->163)
t227 = Ifges(5,1) + Ifges(6,1);
t223 = Ifges(6,4) + Ifges(5,5);
t228 = Ifges(5,6) - Ifges(6,6);
t115 = sin(pkin(8));
t116 = cos(pkin(8));
t118 = sin(qJ(3));
t120 = cos(qJ(3));
t103 = t115 * t118 - t120 * t116;
t98 = t103 * qJD(1);
t95 = qJD(4) + t98;
t158 = -pkin(2) * t116 - pkin(1);
t107 = qJD(1) * t158 + qJD(2);
t117 = sin(qJ(4));
t119 = cos(qJ(4));
t104 = t115 * t120 + t116 * t118;
t99 = t104 * qJD(1);
t61 = pkin(3) * t98 - pkin(7) * t99 + t107;
t180 = pkin(6) + qJ(2);
t109 = t180 * t115;
t105 = qJD(1) * t109;
t110 = t180 * t116;
t106 = qJD(1) * t110;
t79 = -t105 * t118 + t106 * t120;
t74 = qJD(3) * pkin(7) + t79;
t21 = -t117 * t74 + t119 * t61;
t22 = t117 * t61 + t119 * t74;
t133 = t117 * t22 + t119 * t21;
t217 = qJD(5) - t21;
t15 = -pkin(4) * t95 + t217;
t16 = qJ(5) * t95 + t22;
t135 = t117 * t16 - t119 * t15;
t170 = Ifges(6,5) * t119;
t139 = Ifges(6,3) * t117 + t170;
t172 = Ifges(5,4) * t119;
t145 = -Ifges(5,2) * t117 + t172;
t150 = mrSges(6,1) * t117 - mrSges(6,3) * t119;
t152 = mrSges(5,1) * t117 + mrSges(5,2) * t119;
t194 = t119 / 0.2e1;
t196 = t117 / 0.2e1;
t197 = -t117 / 0.2e1;
t199 = -t95 / 0.2e1;
t86 = qJD(3) * t117 + t119 * t99;
t201 = t86 / 0.2e1;
t128 = t119 * qJD(3) - t117 * t99;
t203 = -t128 / 0.2e1;
t204 = t128 / 0.2e1;
t171 = Ifges(6,5) * t117;
t173 = Ifges(5,4) * t117;
t214 = t227 * t119 + t171 - t173;
t215 = -t117 * t228 + t223 * t119;
t190 = Ifges(6,5) * t128;
t84 = Ifges(5,4) * t128;
t219 = t223 * t95 + t227 * t86 - t190 + t84;
t78 = -t105 * t120 - t118 * t106;
t73 = -qJD(3) * pkin(3) - t78;
t23 = -pkin(4) * t128 - qJ(5) * t86 + t73;
t83 = Ifges(6,5) * t86;
t29 = t95 * Ifges(6,6) - Ifges(6,3) * t128 + t83;
t191 = Ifges(5,4) * t86;
t32 = Ifges(5,2) * t128 + t95 * Ifges(5,6) + t191;
t210 = t135 * mrSges(6,2) + t133 * mrSges(5,3) - t139 * t203 - t145 * t204 - t150 * t23 - t152 * t73 - t219 * t194 - t196 * t29 - t197 * t32 + t215 * t199 - t214 * t201;
t224 = t99 * Ifges(4,1) / 0.2e1;
t225 = t107 * mrSges(4,2) - t78 * mrSges(4,3) - Ifges(4,4) * t98 + Ifges(4,5) * qJD(3) - t210 + t224;
t222 = t223 * t117 + t228 * t119;
t221 = t227 * t117 - t170 + t172;
t202 = -t86 / 0.2e1;
t100 = t103 * qJD(3);
t92 = qJD(1) * t100;
t54 = qJD(4) * t128 - t119 * t92;
t55 = qJD(4) * t86 - t117 * t92;
t101 = t104 * qJD(3);
t93 = qJD(1) * t101;
t220 = t223 * t93 + (-Ifges(5,4) + Ifges(6,5)) * t55 + t227 * t54;
t136 = pkin(4) * t117 - qJ(5) * t119;
t218 = -qJD(5) * t117 + t95 * t136 - t79;
t216 = -t120 * t109 - t110 * t118;
t162 = qJD(4) * t119;
t163 = qJD(4) * t117;
t126 = t103 * qJD(2);
t43 = -qJD(1) * t126 + t78 * qJD(3);
t65 = pkin(3) * t93 + pkin(7) * t92;
t4 = t117 * t65 + t119 * t43 + t61 * t162 - t163 * t74;
t5 = -qJD(4) * t22 - t117 * t43 + t119 * t65;
t213 = -t117 * t5 + t119 * t4;
t1 = qJ(5) * t93 + qJD(5) * t95 + t4;
t2 = -pkin(4) * t93 - t5;
t212 = t1 * t119 + t117 * t2;
t211 = (m(3) * qJ(2) + mrSges(3,3)) * (t115 ^ 2 + t116 ^ 2);
t77 = pkin(3) * t103 - pkin(7) * t104 + t158;
t82 = -t109 * t118 + t110 * t120;
t175 = t117 * t77 + t119 * t82;
t62 = t216 * qJD(3) - t126;
t76 = pkin(3) * t101 + pkin(7) * t100;
t9 = -qJD(4) * t175 - t117 * t62 + t119 * t76;
t159 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t160 = Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t161 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t183 = t98 * Ifges(4,2);
t208 = t159 * t128 - t160 * t95 - t161 * t86 + t15 * mrSges(6,1) + t22 * mrSges(5,2) + t79 * mrSges(4,3) + Ifges(5,6) * t203 + Ifges(6,6) * t204 + Ifges(4,6) * qJD(3) + t99 * Ifges(4,4) - t183 / 0.2e1 - t107 * mrSges(4,1) - t16 * mrSges(6,3) - t21 * mrSges(5,1) + t223 * t202 + (Ifges(5,3) + Ifges(6,2)) * t199;
t207 = t54 / 0.2e1;
t206 = -t55 / 0.2e1;
t205 = t55 / 0.2e1;
t200 = t93 / 0.2e1;
t195 = -t119 / 0.2e1;
t193 = mrSges(5,3) * t128;
t192 = mrSges(5,3) * t86;
t127 = t104 * qJD(2);
t44 = qJD(1) * t127 + qJD(3) * t79;
t185 = t44 * t216;
t35 = mrSges(5,1) * t93 - mrSges(5,3) * t54;
t36 = -t93 * mrSges(6,1) + t54 * mrSges(6,2);
t179 = -t35 + t36;
t37 = -mrSges(5,2) * t93 - mrSges(5,3) * t55;
t38 = -mrSges(6,2) * t55 + mrSges(6,3) * t93;
t178 = t37 + t38;
t57 = mrSges(6,2) * t128 + mrSges(6,3) * t95;
t58 = -mrSges(5,2) * t95 + t193;
t177 = t57 + t58;
t59 = mrSges(5,1) * t95 - t192;
t60 = -mrSges(6,1) * t95 + mrSges(6,2) * t86;
t176 = t59 - t60;
t75 = pkin(3) * t99 + pkin(7) * t98;
t27 = t117 * t75 + t119 * t78;
t174 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t128 - mrSges(5,2) * t86 - t99 * mrSges(4,3);
t157 = t93 * mrSges(4,1) - t92 * mrSges(4,2);
t155 = -t1 * t117 + t119 * t2;
t154 = -t117 * t4 - t119 * t5;
t153 = mrSges(5,1) * t119 - mrSges(5,2) * t117;
t151 = mrSges(6,1) * t119 + mrSges(6,3) * t117;
t144 = Ifges(5,2) * t119 + t173;
t138 = -Ifges(6,3) * t119 + t171;
t137 = pkin(4) * t119 + qJ(5) * t117;
t134 = t117 * t15 + t119 * t16;
t132 = t117 * t21 - t119 * t22;
t26 = -t117 * t78 + t119 * t75;
t39 = -t117 * t82 + t119 * t77;
t8 = t117 * t76 + t119 * t62 + t77 * t162 - t163 * t82;
t125 = t5 * mrSges(5,1) - t2 * mrSges(6,1) - t4 * mrSges(5,2) + t1 * mrSges(6,3);
t63 = qJD(3) * t82 + t127;
t108 = -pkin(3) - t137;
t91 = Ifges(6,2) * t93;
t90 = Ifges(5,3) * t93;
t87 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t98;
t52 = Ifges(6,4) * t54;
t51 = Ifges(5,5) * t54;
t50 = Ifges(5,6) * t55;
t49 = Ifges(6,6) * t55;
t46 = -mrSges(6,1) * t128 - mrSges(6,3) * t86;
t45 = pkin(4) * t86 - qJ(5) * t128;
t41 = t104 * t136 - t216;
t25 = -pkin(4) * t103 - t39;
t24 = qJ(5) * t103 + t175;
t20 = -pkin(4) * t99 - t26;
t19 = qJ(5) * t99 + t27;
t18 = mrSges(5,1) * t55 + mrSges(5,2) * t54;
t17 = mrSges(6,1) * t55 - mrSges(6,3) * t54;
t12 = t54 * Ifges(5,4) - t55 * Ifges(5,2) + t93 * Ifges(5,6);
t11 = t54 * Ifges(6,5) + t93 * Ifges(6,6) + t55 * Ifges(6,3);
t10 = -t136 * t100 + (qJD(4) * t137 - qJD(5) * t119) * t104 + t63;
t7 = -pkin(4) * t101 - t9;
t6 = pkin(4) * t55 - qJ(5) * t54 - qJD(5) * t86 + t44;
t3 = qJ(5) * t101 + qJD(5) * t103 + t8;
t13 = [t62 * t87 - t216 * t18 + t3 * t57 + t8 * t58 + t9 * t59 + t7 * t60 + t10 * t46 + t25 * t36 + t24 * t38 + t39 * t35 + t175 * t37 + t41 * t17 + t158 * t157 - t174 * t63 + (t216 * t92 - t82 * t93) * mrSges(4,3) + m(5) * (t175 * t4 + t21 * t9 + t22 * t8 + t39 * t5 + t63 * t73 - t185) + m(4) * (t43 * t82 + t62 * t79 - t63 * t78 - t185) + m(6) * (t1 * t24 + t10 * t23 + t15 * t7 + t16 * t3 + t2 * t25 + t41 * t6) + 0.2e1 * t211 * qJD(2) * qJD(1) + (t183 / 0.2e1 - t208) * t101 - (t224 + t225) * t100 + (t51 / 0.2e1 - t50 / 0.2e1 + t90 / 0.2e1 + t52 / 0.2e1 + t91 / 0.2e1 + t49 / 0.2e1 + Ifges(4,4) * t92 - t43 * mrSges(4,3) + t159 * t55 + t161 * t54 + (Ifges(4,2) + t160) * t93 + t125) * t103 + (t139 * t205 + t145 * t206 + t6 * t150 + t11 * t196 - Ifges(4,1) * t92 - Ifges(4,4) * t93 + (mrSges(4,3) + t152) * t44 + t154 * mrSges(5,3) + t155 * mrSges(6,2) + (-mrSges(6,2) * t134 + mrSges(5,3) * t132 + t138 * t204 + t144 * t203 + t151 * t23 + t153 * t73 + t195 * t32 + t222 * t199 + t221 * t202) * qJD(4) + t214 * t207 + t215 * t200 + (t219 * qJD(4) + t12) * t197 + (qJD(4) * t29 + t220) * t194) * t104; t98 * t87 + (-t46 + t174) * t99 + (t95 * t177 - t179) * t119 + (-t95 * t176 + t178) * t117 - m(4) * (-t78 * t99 - t79 * t98) + t157 - t211 * qJD(1) ^ 2 + (t95 * t134 - t23 * t99 - t155) * m(6) + (-t95 * t132 - t73 * t99 - t154) * m(5); t108 * t17 - Ifges(4,5) * t92 - Ifges(4,6) * t93 - t78 * t87 - t19 * t57 - t27 * t58 - t26 * t59 - t20 * t60 - t43 * mrSges(4,2) - pkin(3) * t18 + t174 * t79 + t218 * t46 + (t108 * t6 - t15 * t20 - t16 * t19 + t218 * t23) * m(6) + t220 * t196 + t212 * mrSges(6,2) + t213 * mrSges(5,3) + (m(5) * t213 + m(6) * t212 + t117 * t179 + t119 * t178 + (-m(5) * t133 - m(6) * t135 - t117 * t177 - t119 * t176) * qJD(4)) * pkin(7) - t210 * qJD(4) + t138 * t205 + t144 * t206 + t12 * t194 + t11 * t195 + (-pkin(3) * t44 - t21 * t26 - t22 * t27 - t73 * t79) * m(5) - t6 * t151 + (-mrSges(4,1) - t153) * t44 + t208 * t99 - ((-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t99 - t225) * t98 + t221 * t207 + t222 * t200; t91 + t90 + t52 - t50 + t51 + t49 - t73 * (mrSges(5,1) * t86 + mrSges(5,2) * t128) - t23 * (mrSges(6,1) * t86 - mrSges(6,3) * t128) + qJD(5) * t57 - t45 * t46 + qJ(5) * t38 - pkin(4) * t36 + (t176 + t192) * t22 + (-t177 + t193) * t21 + t125 + (Ifges(6,3) * t86 + t190) * t204 + t32 * t201 + (-t128 * t15 + t16 * t86) * mrSges(6,2) + (t223 * t128 - t228 * t86) * t199 + (-pkin(4) * t2 + qJ(5) * t1 - t15 * t22 + t217 * t16 - t23 * t45) * m(6) + (-Ifges(5,2) * t86 + t219 + t84) * t203 + (t227 * t128 - t191 + t29 + t83) * t202; t86 * t46 - t95 * t57 + 0.2e1 * (t2 / 0.2e1 + t16 * t199 + t23 * t201) * m(6) + t36;];
tauc = t13(:);
