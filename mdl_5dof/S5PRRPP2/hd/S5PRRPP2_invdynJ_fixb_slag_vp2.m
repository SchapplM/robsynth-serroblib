% Calculate vector of inverse dynamics joint torques for
% S5PRRPP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:44
% EndTime: 2019-12-05 16:09:05
% DurationCPUTime: 8.36s
% Computational Cost: add. (1773->362), mult. (3986->479), div. (0->0), fcn. (2553->10), ass. (0->166)
t243 = mrSges(5,1) + mrSges(6,1);
t245 = -mrSges(5,2) + mrSges(6,3);
t207 = m(5) + m(6);
t242 = -mrSges(6,2) - mrSges(5,3);
t234 = Ifges(5,1) + Ifges(6,1);
t232 = Ifges(6,4) + Ifges(5,5);
t116 = sin(pkin(8));
t120 = sin(qJ(3));
t122 = cos(qJ(3));
t119 = -qJ(4) - pkin(6);
t154 = qJD(3) * t119;
t129 = -qJD(4) * t120 + t122 * t154;
t123 = cos(qJ(2));
t171 = qJD(1) * t123;
t179 = cos(pkin(8));
t71 = qJD(4) * t122 + t120 * t154;
t153 = t179 * t120;
t83 = t116 * t122 + t153;
t230 = t116 * t71 - t129 * t179 - t83 * t171;
t178 = t116 * t120;
t133 = t179 * t122 - t178;
t128 = t123 * t133;
t229 = -qJD(1) * t128 + t116 * t129 + t179 * t71;
t191 = Ifges(4,4) * t120;
t237 = pkin(3) * t207;
t241 = mrSges(4,1) + t237;
t117 = sin(pkin(7));
t118 = cos(pkin(7));
t240 = g(1) * t118 + g(2) * t117;
t239 = -m(4) * pkin(6) + mrSges(3,2) - mrSges(4,3) + t242;
t113 = qJ(3) + pkin(8);
t110 = sin(t113);
t111 = cos(t113);
t145 = -mrSges(4,1) * t122 + mrSges(4,2) * t120;
t238 = m(4) * pkin(2) + t245 * t110 + t243 * t111 + mrSges(3,1) - t145;
t233 = -Ifges(5,4) + Ifges(6,5);
t208 = t83 / 0.2e1;
t163 = qJD(2) * qJD(3);
t86 = qJDD(2) * t122 - t120 * t163;
t87 = qJDD(2) * t120 + t122 * t163;
t41 = t116 * t87 - t179 * t86;
t42 = t116 * t86 + t179 * t87;
t6 = t41 * mrSges(6,1) - t42 * mrSges(6,3);
t7 = t41 * mrSges(5,1) + t42 * mrSges(5,2);
t235 = -t6 - t7;
t194 = qJD(3) / 0.2e1;
t231 = Ifges(6,6) - Ifges(5,6);
t149 = qJD(2) * t179;
t170 = qJD(2) * t120;
t72 = t116 * t170 - t122 * t149;
t205 = Ifges(6,5) * t72;
t68 = Ifges(5,4) * t72;
t74 = t83 * qJD(2);
t228 = t232 * qJD(3) + t234 * t74 + t205 - t68;
t30 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t42;
t31 = -qJDD(3) * mrSges(6,1) + t42 * mrSges(6,2);
t227 = t31 - t30;
t29 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t41;
t32 = -mrSges(6,2) * t41 + qJDD(3) * mrSges(6,3);
t226 = t32 + t29;
t53 = -mrSges(6,2) * t72 + qJD(3) * mrSges(6,3);
t193 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t72 + t53;
t192 = t243 * qJD(3) + t242 * t74;
t121 = sin(qJ(2));
t225 = t121 * t240;
t93 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t170;
t168 = qJD(2) * t122;
t94 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t168;
t224 = t120 * t94 + t122 * t93;
t223 = t122 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t86) - t120 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t87);
t167 = qJD(3) * t120;
t164 = qJD(1) * qJD(2);
t107 = t123 * t164;
t89 = t121 * qJDD(1) + t107;
t81 = qJDD(2) * pkin(6) + t89;
t172 = qJD(1) * t121;
t99 = qJD(2) * pkin(6) + t172;
t39 = t122 * t81 - t167 * t99;
t165 = qJD(3) * t122;
t158 = t99 * t165;
t40 = -t120 * t81 - t158;
t138 = -t120 * t40 + t122 * t39;
t222 = t120 * t93 - t122 * t94;
t33 = mrSges(6,1) * t72 - mrSges(6,3) * t74;
t34 = mrSges(5,1) * t72 + mrSges(5,2) * t74;
t220 = t145 * qJD(2) + t33 + t34;
t219 = m(6) * pkin(4) + t243;
t217 = -m(6) * qJ(5) - t245;
t108 = pkin(3) * t122 + pkin(2);
t76 = -qJD(2) * t108 + qJD(4) - t171;
t14 = pkin(4) * t72 - qJ(5) * t74 + t76;
t148 = qJ(4) * qJD(2) + t99;
t65 = t148 * t122;
t48 = t179 * t65;
t64 = t148 * t120;
t54 = qJD(3) * pkin(3) - t64;
t17 = t116 * t54 + t48;
t15 = qJD(3) * qJ(5) + t17;
t216 = t76 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) - t17 * mrSges(5,3);
t189 = t116 * t65;
t16 = t179 * t54 - t189;
t12 = -qJD(3) * pkin(4) + qJD(5) - t16;
t215 = mrSges(5,2) * t76 + t12 * mrSges(6,2) - t16 * mrSges(5,3) - mrSges(6,3) * t14;
t124 = qJD(2) ^ 2;
t213 = -t72 / 0.2e1;
t212 = t72 / 0.2e1;
t210 = t74 / 0.2e1;
t162 = qJD(2) * qJD(4);
t11 = -t158 + qJDD(3) * pkin(3) - qJ(4) * t87 + (-t81 - t162) * t120;
t13 = qJ(4) * t86 + t122 * t162 + t39;
t4 = t116 * t11 + t179 * t13;
t206 = Ifges(5,4) * t74;
t204 = pkin(3) * t116;
t200 = g(3) * t121;
t190 = Ifges(4,4) * t122;
t176 = t117 * t123;
t175 = t118 * t123;
t174 = t120 * t123;
t173 = t122 * t123;
t169 = qJD(2) * t121;
t166 = qJD(3) * t121;
t161 = pkin(3) * t170;
t160 = pkin(3) * t167;
t156 = t179 * pkin(3);
t106 = t121 * t164;
t88 = qJDD(1) * t123 - t106;
t144 = mrSges(4,1) * t120 + mrSges(4,2) * t122;
t141 = t122 * Ifges(4,2) + t191;
t140 = Ifges(4,5) * t122 - Ifges(4,6) * t120;
t139 = pkin(4) * t111 + qJ(5) * t110;
t3 = t11 * t179 - t116 * t13;
t100 = -qJD(2) * pkin(2) - t171;
t135 = t100 * t144;
t134 = t120 * (Ifges(4,1) * t122 - t191);
t80 = -qJDD(2) * pkin(2) - t88;
t130 = t100 * t121 + (t120 ^ 2 + t122 ^ 2) * t99 * t123;
t43 = -pkin(3) * t86 + qJDD(4) + t80;
t109 = Ifges(4,4) * t168;
t105 = -t156 - pkin(4);
t103 = qJ(5) + t204;
t95 = t119 * t122;
t78 = Ifges(4,1) * t170 + Ifges(4,5) * qJD(3) + t109;
t77 = Ifges(4,6) * qJD(3) + qJD(2) * t141;
t75 = t133 * qJD(3);
t73 = t83 * qJD(3);
t67 = Ifges(6,5) * t74;
t63 = t133 * t121;
t62 = t83 * t121;
t60 = t110 * t175 - t117 * t111;
t58 = t110 * t176 + t111 * t118;
t46 = t119 * t178 - t179 * t95;
t45 = -t116 * t95 - t119 * t153;
t44 = -mrSges(4,1) * t86 + mrSges(4,2) * t87;
t35 = -pkin(4) * t133 - qJ(5) * t83 - t108;
t26 = -t72 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t206;
t25 = Ifges(6,6) * qJD(3) + t72 * Ifges(6,3) + t67;
t22 = pkin(4) * t74 + qJ(5) * t72 + t161;
t21 = qJD(2) * t128 - t166 * t83;
t20 = -t116 * t123 * t168 - t133 * t166 - t149 * t174;
t19 = -t179 * t64 - t189;
t18 = -t116 * t64 + t48;
t9 = pkin(4) * t73 - qJ(5) * t75 - qJD(5) * t83 + t160;
t5 = pkin(4) * t41 - qJ(5) * t42 - qJD(5) * t74 + t43;
t2 = -qJDD(3) * pkin(4) + qJDD(5) - t3;
t1 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t4;
t8 = [m(2) * qJDD(1) + t226 * t63 + t227 * t62 + t193 * t21 + t192 * t20 + (-m(2) - m(3) - m(4) - t207) * g(3) + (qJDD(2) * mrSges(3,1) - t124 * mrSges(3,2) - qJD(2) * t222 + t235 - t44) * t123 + (-t124 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + qJD(2) * t220 - qJD(3) * t224 + t223) * t121 + m(4) * (qJD(2) * t130 + t121 * t138 - t123 * t80) + m(3) * (t121 * t89 + t123 * t88) + m(5) * (-t123 * t43 + t16 * t20 + t169 * t76 + t17 * t21 - t3 * t62 + t4 * t63) + m(6) * (t1 * t63 - t12 * t20 - t123 * t5 + t14 * t169 + t15 * t21 + t2 * t62); (t191 + t141) * t86 / 0.2e1 + t122 * (Ifges(4,4) * t87 + Ifges(4,2) * t86) / 0.2e1 + (-(Ifges(5,2) + Ifges(6,3)) * t133 + 0.2e1 * t233 * t208) * t41 + (-t133 * t231 + t232 * t83) * qJDD(3) / 0.2e1 + (t25 / 0.2e1 - t26 / 0.2e1 + Ifges(6,3) * t212 - Ifges(5,2) * t213 + t233 * t210 + t231 * t194 + t216) * t73 + (t232 * qJDD(3) + t234 * t42) * t208 + (t140 * t194 + t135) * qJD(3) + t226 * t46 + t227 * t45 + (t1 * t46 + t12 * t230 + t14 * t9 + t15 * t229 + t2 * t45 + t35 * t5) * m(6) + (-t108 * t43 - t16 * t230 + t160 * t76 + t17 * t229 - t3 * t45 + t4 * t46) * m(5) + t122 * qJDD(3) * Ifges(4,6) + t240 * ((t207 * t119 + t239) * t123 + (m(5) * t108 - m(6) * (-t108 - t139) + t238) * t121) + (t122 * (-Ifges(4,2) * t120 + t190) + t134) * t163 / 0.2e1 - t108 * t7 + (t88 + t106) * mrSges(3,1) + t87 * t120 * Ifges(4,1) + (-pkin(2) * t80 - qJD(1) * t130) * m(4) + (-t89 + t107) * mrSges(3,2) + (Ifges(5,4) * t213 + Ifges(6,5) * t212 + t232 * t194 + t234 * t210 + t215 + t228 / 0.2e1) * t75 + t222 * t171 + t138 * mrSges(4,3) + (m(4) * t138 - t165 * t93 - t167 * t94 + t223) * pkin(6) + (-m(5) * t76 - m(6) * t14 - t220) * t172 + t87 * t190 / 0.2e1 - pkin(2) * t44 + t9 * t33 + t35 * t6 + t80 * t145 + (-t207 * (t123 * t108 - t119 * t121) + (-m(6) * t139 - t238) * t123 + t239 * t121) * g(3) + t43 * (-mrSges(5,1) * t133 + mrSges(5,2) * t83) + t5 * (-mrSges(6,1) * t133 - mrSges(6,3) * t83) - t133 * (Ifges(6,5) * t42 + Ifges(6,6) * qJDD(3)) / 0.2e1 + t133 * (Ifges(5,4) * t42 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t78 * t165 / 0.2e1 + t34 * t160 + t120 * Ifges(4,5) * qJDD(3) + (-t133 * t233 + t234 * t83) * t42 / 0.2e1 + (t1 * t133 + t2 * t83) * mrSges(6,2) + (t133 * t4 - t3 * t83) * mrSges(5,3) + t229 * t193 - t230 * t192 - t77 * t167 / 0.2e1 + Ifges(3,3) * qJDD(2); -qJD(2) * t135 - t205 * t213 + t232 * t42 - (t231 * t74 - t232 * t72) * qJD(3) / 0.2e1 - (-t234 * t72 - t206 + t25 + t67) * t74 / 0.2e1 + t192 * t18 - t193 * t19 + (-Ifges(5,2) * t74 + t228 - t68) * t212 + t231 * t41 + t215 * t72 + (Ifges(6,3) * t213 - t216) * t74 + (t110 * t219 + t111 * t217 + t120 * t237 + t144) * t200 + (-(-t117 * t120 - t118 * t173) * mrSges(4,2) + t217 * (t110 * t117 + t111 * t175) + t219 * t60 + t241 * (-t117 * t122 + t118 * t174)) * g(1) + (-(-t117 * t173 + t118 * t120) * mrSges(4,2) + t217 * (-t118 * t110 + t111 * t176) + t219 * t58 - t241 * (-t117 * t174 - t118 * t122)) * g(2) - (-Ifges(4,2) * t170 + t109 + t78) * t168 / 0.2e1 + t103 * t32 + t105 * t31 + Ifges(4,6) * t86 + Ifges(4,5) * t87 + (t1 * t103 + t105 * t2 - t12 * t18 - t14 * t22 + (-t19 + qJD(5)) * t15) * m(6) + t224 * t99 + qJD(5) * t53 - t22 * t33 - t39 * mrSges(4,2) + t40 * mrSges(4,1) + t1 * mrSges(6,3) - t2 * mrSges(6,1) + t3 * mrSges(5,1) - t4 * mrSges(5,2) + t29 * t204 + t26 * t210 + t30 * t156 + (t16 * t18 - t76 * t161 - t17 * t19 + (t116 * t4 + t179 * t3) * pkin(3)) * m(5) + (Ifges(6,2) + Ifges(5,3) + Ifges(4,3)) * qJDD(3) - t140 * t163 / 0.2e1 - t34 * t161 + t77 * t170 / 0.2e1 - t124 * t134 / 0.2e1; t207 * t123 * g(3) + t192 * t74 + t193 * t72 + (-t12 * t74 + t15 * t72 - t225 + t5) * m(6) + (t16 * t74 + t17 * t72 - t225 + t43) * m(5) - t235; -qJD(3) * t53 + t74 * t33 + (-g(1) * t60 - g(2) * t58 - t15 * qJD(3) - t110 * t200 + t14 * t74 + t2) * m(6) + t31;];
tau = t8;
