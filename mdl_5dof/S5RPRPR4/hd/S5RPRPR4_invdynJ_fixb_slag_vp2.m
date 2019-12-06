% Calculate vector of inverse dynamics joint torques for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:10
% EndTime: 2019-12-05 17:53:26
% DurationCPUTime: 6.63s
% Computational Cost: add. (3408->390), mult. (7405->537), div. (0->0), fcn. (5001->16), ass. (0->175)
t251 = -qJD(3) / 0.2e1;
t151 = qJDD(3) + qJDD(5);
t153 = qJD(3) + qJD(5);
t156 = sin(pkin(9));
t158 = cos(pkin(9));
t162 = sin(qJ(3));
t165 = cos(qJ(3));
t112 = -t156 * t162 + t158 * t165;
t104 = t112 * qJD(1);
t189 = t165 * qJD(1);
t190 = t162 * qJD(1);
t105 = -t156 * t189 - t158 * t190;
t161 = sin(qJ(5));
t164 = cos(qJ(5));
t179 = t164 * t104 + t105 * t161;
t188 = qJD(1) * qJD(3);
t118 = qJDD(1) * t165 - t162 * t188;
t119 = qJDD(1) * t162 + t165 * t188;
t68 = t118 * t158 - t119 * t156;
t69 = t118 * t156 + t119 * t158;
t19 = qJD(5) * t179 + t161 * t68 + t164 * t69;
t157 = sin(pkin(8));
t135 = pkin(1) * t157 + pkin(6);
t124 = t135 * qJD(1);
t147 = t165 * qJDD(2);
t187 = qJD(1) * qJD(4);
t191 = qJD(3) * t165;
t122 = t135 * qJDD(1);
t247 = qJD(2) * qJD(3) + t122;
t37 = -t124 * t191 + qJDD(3) * pkin(3) - qJ(4) * t119 + t147 + (-t187 - t247) * t162;
t192 = qJD(3) * t162;
t57 = t162 * qJDD(2) - t124 * t192 + t247 * t165;
t40 = qJ(4) * t118 + t165 * t187 + t57;
t12 = -t156 * t40 + t158 * t37;
t6 = qJDD(3) * pkin(4) - pkin(7) * t69 + t12;
t214 = pkin(7) * t105;
t177 = qJ(4) * qJD(1) + t124;
t193 = qJD(2) * t162;
t83 = t165 * t177 + t193;
t72 = t156 * t83;
t149 = t165 * qJD(2);
t82 = -t162 * t177 + t149;
t76 = qJD(3) * pkin(3) + t82;
t38 = t158 * t76 - t72;
t26 = qJD(3) * pkin(4) + t214 + t38;
t215 = pkin(7) * t104;
t203 = t158 * t83;
t39 = t156 * t76 + t203;
t27 = t39 + t215;
t7 = -t161 * t27 + t164 * t26;
t13 = t156 * t37 + t158 * t40;
t9 = pkin(7) * t68 + t13;
t2 = qJD(5) * t7 + t161 * t6 + t164 * t9;
t56 = t104 * t161 - t105 * t164;
t20 = -qJD(5) * t56 - t161 * t69 + t164 * t68;
t219 = Ifges(6,4) * t56;
t50 = Ifges(6,4) * t179;
t24 = Ifges(6,1) * t56 + Ifges(6,5) * t153 + t50;
t154 = qJ(3) + pkin(9);
t148 = qJ(5) + t154;
t133 = sin(t148);
t134 = cos(t148);
t246 = mrSges(6,1) * t133 + mrSges(6,2) * t134;
t155 = qJ(1) + pkin(8);
t143 = sin(t155);
t249 = g(2) * t143;
t8 = t161 * t26 + t164 * t27;
t3 = -qJD(5) * t8 - t161 * t9 + t164 * t6;
t159 = cos(pkin(8));
t137 = -pkin(1) * t159 - pkin(2);
t150 = t165 * pkin(3);
t121 = t137 - t150;
t103 = qJD(1) * t121 + qJD(4);
t64 = -pkin(4) * t104 + t103;
t250 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * t151 - t246 * t249 - (Ifges(6,5) * t179 - Ifges(6,6) * t56) * t153 / 0.2e1 + (t179 * t7 + t56 * t8) * mrSges(6,3) - (-Ifges(6,2) * t56 + t24 + t50) * t179 / 0.2e1 - t64 * (mrSges(6,1) * t56 + mrSges(6,2) * t179) - (Ifges(6,1) * t179 - t219) * t56 / 0.2e1;
t231 = m(3) + m(6) + m(5) + m(4);
t248 = pkin(1) * t231 + mrSges(2,1);
t128 = -t165 * mrSges(4,1) + t162 * mrSges(4,2);
t142 = sin(t154);
t144 = cos(t154);
t180 = t134 * mrSges(6,1) - t133 * mrSges(6,2);
t245 = mrSges(5,1) * t144 - mrSges(5,2) * t142 - t128 + t180;
t244 = -Ifges(4,4) * t189 / 0.2e1 + Ifges(4,5) * t251;
t23 = Ifges(6,2) * t179 + Ifges(6,6) * t153 + t219;
t242 = t23 / 0.2e1;
t136 = pkin(3) * t158 + pkin(4);
t216 = pkin(3) * t156;
t101 = t136 * t164 - t161 * t216;
t42 = -t156 * t82 - t203;
t30 = t42 - t215;
t43 = t158 * t82 - t72;
t31 = t43 + t214;
t240 = t101 * qJD(5) - t161 * t30 - t164 * t31;
t102 = t136 * t161 + t164 * t216;
t239 = -t102 * qJD(5) + t161 * t31 - t164 * t30;
t200 = qJDD(3) / 0.2e1;
t94 = t124 * t165 + t193;
t58 = -t94 * qJD(3) - t122 * t162 + t147;
t235 = -t58 * t162 + t57 * t165;
t207 = Ifges(4,4) * t162;
t173 = t165 * Ifges(4,2) + t207;
t198 = Ifges(4,6) * qJD(3);
t234 = -t94 * mrSges(4,3) - qJD(1) * t173 / 0.2e1 - t198 / 0.2e1;
t93 = -t124 * t162 + t149;
t233 = -t93 * mrSges(4,3) + Ifges(4,1) * t190 / 0.2e1 - t244;
t232 = 0.2e1 * t200;
t209 = mrSges(4,2) * t165;
t226 = m(5) * pkin(3);
t230 = -t162 * (mrSges(4,1) + t226) + m(6) * (-pkin(3) * t162 - pkin(4) * t142) - mrSges(5,1) * t142 - mrSges(5,2) * t144 - t209;
t194 = pkin(4) * t144 + t150;
t229 = mrSges(3,1) + m(6) * (pkin(2) + t194) + m(5) * (t150 + pkin(2)) + m(4) * pkin(2) + t245;
t160 = -qJ(4) - pkin(6);
t227 = m(4) * pkin(6) - m(5) * t160 - m(6) * (-pkin(7) + t160) - mrSges(3,2) + mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t223 = t56 / 0.2e1;
t221 = -t105 / 0.2e1;
t145 = cos(t155);
t213 = g(3) * t145;
t195 = qJ(4) + t135;
t178 = qJD(3) * t195;
t86 = qJD(4) * t165 - t162 * t178;
t87 = -qJD(4) * t162 - t165 * t178;
t45 = t156 * t87 + t158 * t86;
t110 = t195 * t162;
t111 = t195 * t165;
t60 = -t156 * t110 + t158 * t111;
t206 = Ifges(4,4) * t165;
t205 = Ifges(5,4) * t105;
t185 = pkin(3) * t192;
t183 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t182 = -t20 * mrSges(6,1) + t19 * mrSges(6,2);
t44 = -t156 * t86 + t158 * t87;
t59 = -t158 * t110 - t111 * t156;
t176 = -g(2) * t145 - g(3) * t143;
t123 = t137 * qJDD(1);
t113 = t156 * t165 + t158 * t162;
t46 = -pkin(7) * t113 + t59;
t47 = pkin(7) * t112 + t60;
t21 = -t161 * t47 + t164 * t46;
t22 = t161 * t46 + t164 * t47;
t65 = t112 * t164 - t113 * t161;
t66 = t112 * t161 + t113 * t164;
t81 = -pkin(3) * t118 + qJDD(4) + t123;
t166 = cos(qJ(1));
t163 = sin(qJ(1));
t127 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t189;
t126 = t137 * qJD(1);
t125 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t190;
t107 = t112 * qJD(3);
t106 = t113 * qJD(3);
t100 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t99 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t118;
t98 = Ifges(5,4) * t104;
t89 = qJD(3) * mrSges(5,1) + mrSges(5,3) * t105;
t88 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t104;
t85 = pkin(4) * t106 + t185;
t84 = pkin(3) * t190 - pkin(4) * t105;
t80 = -pkin(4) * t112 + t121;
t63 = -mrSges(5,1) * t104 - mrSges(5,2) * t105;
t62 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t69;
t61 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t68;
t52 = -t105 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t98;
t51 = t104 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t205;
t49 = mrSges(6,1) * t153 - mrSges(6,3) * t56;
t48 = -mrSges(6,2) * t153 + mrSges(6,3) * t179;
t41 = -pkin(4) * t68 + t81;
t33 = -pkin(7) * t106 + t45;
t32 = -pkin(7) * t107 + t44;
t29 = -qJD(5) * t66 - t106 * t164 - t107 * t161;
t28 = qJD(5) * t65 - t106 * t161 + t107 * t164;
t25 = -mrSges(6,1) * t179 + mrSges(6,2) * t56;
t15 = -mrSges(6,2) * t151 + mrSges(6,3) * t20;
t14 = mrSges(6,1) * t151 - mrSges(6,3) * t19;
t5 = -qJD(5) * t22 - t161 * t33 + t164 * t32;
t4 = qJD(5) * t21 + t161 * t32 + t164 * t33;
t1 = [qJD(3) ^ 2 * (Ifges(4,5) * t165 - Ifges(4,6) * t162) / 0.2e1 + t119 * t162 * Ifges(4,1) + t126 * (mrSges(4,1) * t162 + t209) * qJD(3) + (-mrSges(2,2) * t163 + t143 * t227 + t145 * t229 + t248 * t166) * g(2) + (mrSges(2,2) * t166 + t143 * t229 - t145 * t227 + t248 * t163) * g(3) + t165 * (Ifges(4,4) * t119 + Ifges(4,2) * t118 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t153 * (Ifges(6,5) * t28 + Ifges(6,6) * t29) / 0.2e1 + t137 * (-mrSges(4,1) * t118 + mrSges(4,2) * t119) - t106 * t51 / 0.2e1 + t107 * t52 / 0.2e1 + t44 * t89 + t118 * t173 / 0.2e1 + t85 * t25 + t45 * t88 + t60 * t61 + t59 * t62 + t64 * (-mrSges(6,1) * t29 + mrSges(6,2) * t28) + t4 * t48 + t5 * t49 + (mrSges(6,2) * t41 - mrSges(6,3) * t3 + Ifges(6,1) * t19 + Ifges(6,4) * t20 + Ifges(6,5) * t151) * t66 + t28 * t24 / 0.2e1 + t21 * t14 + t22 * t15 + t63 * t185 + (-t81 * mrSges(5,1) + t13 * mrSges(5,3) + Ifges(5,4) * t69 + t68 * Ifges(5,2) + Ifges(5,6) * t232) * t112 + (t81 * mrSges(5,2) - t12 * mrSges(5,3) + t69 * Ifges(5,1) + Ifges(5,4) * t68 + Ifges(5,5) * t232) * t113 + t233 * t191 + t234 * t192 + t162 * (Ifges(4,4) * t118 + Ifges(4,5) * qJDD(3)) / 0.2e1 + m(5) * (t103 * t185 + t12 * t59 + t121 * t81 + t13 * t60 + t38 * t44 + t39 * t45) + t80 * t182 + t121 * t183 + (t165 * (-Ifges(4,2) * t162 + t206) + t162 * (Ifges(4,1) * t165 - t207)) * t188 / 0.2e1 + (m(4) * t137 + t128) * t123 + t179 * (Ifges(6,4) * t28 + Ifges(6,2) * t29) / 0.2e1 + t103 * (mrSges(5,1) * t106 + mrSges(5,2) * t107) + t104 * (Ifges(5,4) * t107 - Ifges(5,2) * t106) / 0.2e1 + qJD(3) * (Ifges(5,5) * t107 - Ifges(5,6) * t106) / 0.2e1 + (-t106 * t39 - t107 * t38) * mrSges(5,3) + (Ifges(5,1) * t107 - Ifges(5,4) * t106) * t221 + (-t28 * t7 + t29 * t8) * mrSges(6,3) + t29 * t242 + (-t125 * t191 - t127 * t192 + t165 * t99 - t162 * t100 + m(4) * ((-t94 * t162 - t93 * t165) * qJD(3) + t235)) * t135 + t235 * mrSges(4,3) + t119 * t206 / 0.2e1 + (Ifges(4,5) * t162 + Ifges(4,6) * t165) * t200 + (Ifges(6,1) * t28 + Ifges(6,4) * t29) * t223 + m(6) * (t2 * t22 + t21 * t3 + t4 * t8 + t41 * t80 + t5 * t7 + t64 * t85) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t159 * mrSges(3,1) - 0.2e1 * t157 * mrSges(3,2) + m(3) * (t157 ^ 2 + t159 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (-mrSges(6,1) * t41 + mrSges(6,3) * t2 + Ifges(6,4) * t19 + Ifges(6,2) * t20 + Ifges(6,6) * t151) * t65; m(3) * qJDD(2) + t165 * t100 - t106 * t89 + t107 * t88 + t112 * t62 + t113 * t61 + t65 * t14 + t66 * t15 + t162 * t99 + t28 * t48 + t29 * t49 + (-t125 * t162 + t127 * t165) * qJD(3) + m(4) * (t162 * t57 + t165 * t58 + (-t162 * t93 + t165 * t94) * qJD(3)) + m(5) * (-t106 * t38 + t107 * t39 + t112 * t12 + t113 * t13) + m(6) * (t2 * t66 + t28 * t8 + t29 * t7 + t3 * t65) - t231 * g(1); (t104 * t38 - t105 * t39) * mrSges(5,3) + Ifges(4,5) * t119 + t94 * t125 - t93 * t127 + Ifges(4,6) * t118 + t101 * t14 + t102 * t15 - t103 * (-mrSges(5,1) * t105 + mrSges(5,2) * t104) - t42 * t89 - m(5) * (t38 * t42 + t39 * t43) + ((-t126 * mrSges(4,2) - t233 + t244) * t165 + (-t126 * mrSges(4,1) + t198 / 0.2e1 + (t207 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t165) * qJD(1) + (-m(5) * t103 - t63) * pkin(3) - t234) * t162) * qJD(1) - t84 * t25 - t43 * t88 + Ifges(5,6) * t68 + Ifges(5,5) * t69 - t57 * mrSges(4,2) + t58 * mrSges(4,1) + t12 * mrSges(5,1) - t13 * mrSges(5,2) + (t246 - t230) * t213 + t105 * (Ifges(5,1) * t104 + t205) / 0.2e1 - (Ifges(5,2) * t105 + t52 + t98) * t104 / 0.2e1 + (-m(5) * t150 - t245) * g(1) + t250 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t239 * t49 + t240 * t48 + (-t194 * g(1) + t101 * t3 + t102 * t2 + t239 * t7 + t240 * t8 - t64 * t84) * m(6) + t230 * t249 + (t156 * t61 + t158 * t62) * pkin(3) + (Ifges(5,5) * t104 + Ifges(5,6) * t105) * t251 + t51 * t221 + (t12 * t158 + t13 * t156) * t226 + t56 * t242; -t104 * t88 - t105 * t89 - t179 * t48 + t56 * t49 + t182 + t183 + (-t179 * t8 + t56 * t7 + t176 + t41) * m(6) + (-t104 * t39 - t105 * t38 + t176 + t81) * m(5); -g(1) * t180 + t213 * t246 + t23 * t223 - t7 * t48 + t8 * t49 + t250;];
tau = t1;
