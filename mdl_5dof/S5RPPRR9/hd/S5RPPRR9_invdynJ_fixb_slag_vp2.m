% Calculate vector of inverse dynamics joint torques for
% S5RPPRR9
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:19
% EndTime: 2019-12-31 18:02:29
% DurationCPUTime: 5.34s
% Computational Cost: add. (2351->381), mult. (4338->531), div. (0->0), fcn. (2373->8), ass. (0->178)
t110 = sin(qJ(4));
t185 = qJD(1) * t110;
t160 = mrSges(5,3) * t185;
t109 = sin(qJ(5));
t111 = cos(qJ(5));
t180 = qJD(4) * t111;
t68 = t109 * t185 + t180;
t172 = t109 * qJD(4);
t69 = t111 * t185 - t172;
t207 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t68 - mrSges(6,2) * t69 - t160;
t112 = cos(qJ(4));
t106 = sin(pkin(8));
t107 = cos(pkin(8));
t186 = qJ(2) * qJD(1);
t113 = -pkin(1) - pkin(2);
t86 = qJD(1) * t113 + qJD(2);
t52 = t106 * t86 + t107 * t186;
t47 = -qJD(1) * pkin(6) + t52;
t34 = qJD(3) * t112 - t110 * t47;
t30 = -qJD(4) * pkin(4) - t34;
t244 = m(6) * t30;
t226 = t207 + t244;
t251 = t68 / 0.2e1;
t219 = -t69 / 0.2e1;
t184 = qJD(1) * t112;
t88 = qJD(5) + t184;
t250 = t88 / 0.2e1;
t249 = -qJD(1) / 0.2e1;
t179 = qJD(4) * t112;
t181 = qJD(4) * t110;
t85 = qJDD(1) * t113 + qJDD(2);
t171 = qJD(1) * qJD(2);
t87 = qJDD(1) * qJ(2) + t171;
t43 = t106 * t85 + t107 * t87;
t41 = -qJDD(1) * pkin(6) + t43;
t12 = qJD(3) * t179 + t110 * qJDD(3) + t112 * t41 - t181 * t47;
t10 = qJDD(4) * pkin(7) + t12;
t42 = -t106 * t87 + t107 * t85;
t40 = qJDD(1) * pkin(3) - t42;
t170 = qJD(1) * qJD(4);
t75 = -qJDD(1) * t112 + t110 * t170;
t76 = -qJDD(1) * t110 - t112 * t170;
t14 = -pkin(4) * t75 - pkin(7) * t76 + t40;
t35 = qJD(3) * t110 + t112 * t47;
t31 = qJD(4) * pkin(7) + t35;
t146 = pkin(4) * t112 + pkin(7) * t110;
t51 = -t106 * t186 + t107 * t86;
t46 = qJD(1) * pkin(3) - t51;
t32 = qJD(1) * t146 + t46;
t8 = -t109 * t31 + t111 * t32;
t1 = qJD(5) * t8 + t10 * t111 + t109 * t14;
t9 = t109 * t32 + t111 * t31;
t2 = -qJD(5) * t9 - t10 * t109 + t111 * t14;
t144 = t1 * t111 - t109 * t2;
t176 = qJD(5) * t111;
t178 = qJD(5) * t109;
t248 = -t8 * t176 - t9 * t178 + t144;
t247 = -t75 / 0.2e1;
t246 = -t76 / 0.2e1;
t245 = -m(5) - m(6);
t28 = qJD(5) * t68 + qJDD(4) * t109 + t111 * t76;
t29 = qJD(5) * t69 + qJDD(4) * t111 - t109 * t76;
t7 = -mrSges(6,1) * t29 + mrSges(6,2) * t28;
t215 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t76 + t7;
t243 = mrSges(3,1) + mrSges(2,1);
t242 = -mrSges(3,3) + mrSges(2,2);
t241 = -mrSges(5,3) + mrSges(4,2);
t155 = t106 * t181;
t191 = t109 * t112;
t189 = t111 * t112;
t57 = t106 * t189 - t107 * t109;
t240 = -qJD(5) * t57 + t109 * t155 - (t106 * t111 - t107 * t191) * qJD(1);
t56 = -t106 * t191 - t107 * t111;
t239 = qJD(5) * t56 - t111 * t155 - (t106 * t109 + t107 * t189) * qJD(1);
t205 = Ifges(5,4) * t112;
t140 = -t110 * Ifges(5,1) - t205;
t64 = Ifges(6,4) * t68;
t23 = -Ifges(6,1) * t69 + Ifges(6,5) * t88 + t64;
t238 = Ifges(5,5) * qJD(4) + qJD(1) * t140 + t111 * t23;
t206 = Ifges(5,4) * t110;
t237 = t110 * (-Ifges(5,1) * t112 + t206) + t112 * (Ifges(5,2) * t110 - t205);
t188 = t35 * qJD(4);
t13 = qJDD(3) * t112 - t110 * t41 - t188;
t182 = qJD(2) * t107;
t78 = t107 * qJ(2) + t106 * t113;
t71 = -pkin(6) + t78;
t236 = t112 * t182 - t71 * t181;
t234 = -t110 * t13 + t112 * t12;
t65 = qJDD(5) - t75;
t15 = mrSges(6,1) * t65 - mrSges(6,3) * t28;
t16 = -mrSges(6,2) * t65 + mrSges(6,3) * t29;
t233 = -t109 * t15 + t111 * t16;
t231 = m(4) - t245;
t137 = -t112 * Ifges(5,2) - t206;
t230 = -Ifges(5,6) * qJD(4) / 0.2e1 + t137 * t249 + Ifges(6,5) * t219 + Ifges(6,6) * t251 + Ifges(6,3) * t250;
t141 = t109 * mrSges(6,1) + t111 * mrSges(6,2);
t212 = Ifges(6,4) * t69;
t22 = Ifges(6,2) * t68 + Ifges(6,6) * t88 - t212;
t228 = -t109 * t22 / 0.2e1 + t141 * t30;
t227 = t107 * t112 * t35 + t106 * t46;
t225 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t142 = -mrSges(6,1) * t111 + mrSges(6,2) * t109;
t125 = m(6) * pkin(4) - t142;
t164 = m(6) * pkin(7) + mrSges(6,3);
t84 = mrSges(5,1) * t112 - mrSges(5,2) * t110;
t224 = t110 * t164 + t112 * t125 + t84;
t114 = qJD(1) ^ 2;
t223 = t28 / 0.2e1;
t222 = t29 / 0.2e1;
t221 = t65 / 0.2e1;
t218 = t69 / 0.2e1;
t214 = cos(qJ(1));
t213 = sin(qJ(1));
t204 = Ifges(6,4) * t109;
t203 = Ifges(6,4) * t111;
t194 = t112 * t71;
t193 = t107 * t110;
t192 = t109 * t110;
t190 = t110 * t111;
t187 = t214 * pkin(1) + t213 * qJ(2);
t183 = qJD(2) * t106;
t177 = qJD(5) * t110;
t175 = qJDD(1) * mrSges(3,1);
t174 = qJDD(1) * mrSges(4,1);
t173 = qJDD(1) * mrSges(4,2);
t167 = Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * t65;
t163 = mrSges(5,3) * t184;
t158 = t214 * pkin(2) + t187;
t154 = t112 * t172;
t150 = t176 / 0.2e1;
t149 = -t170 / 0.2e1;
t77 = -t106 * qJ(2) + t107 * t113;
t70 = pkin(3) - t77;
t148 = -pkin(1) * t213 + t214 * qJ(2);
t145 = -pkin(4) * t110 + pkin(7) * t112;
t147 = qJD(4) * t145 - qJD(5) * t194 + t183;
t143 = mrSges(5,1) * t110 + mrSges(5,2) * t112;
t139 = Ifges(6,1) * t111 - t204;
t138 = Ifges(6,1) * t109 + t203;
t136 = -Ifges(6,2) * t109 + t203;
t135 = Ifges(6,2) * t111 + t204;
t134 = -Ifges(5,5) * t112 + Ifges(5,6) * t110;
t133 = Ifges(6,5) * t111 - Ifges(6,6) * t109;
t132 = Ifges(6,5) * t109 + Ifges(6,6) * t111;
t131 = -t106 * t51 + t107 * t52;
t128 = t46 * t143;
t124 = -t109 * t177 + t111 * t179;
t123 = t110 * t176 + t154;
t122 = -pkin(2) * t213 + t148;
t120 = -Ifges(6,5) * t110 - t112 * t139;
t119 = -Ifges(6,6) * t110 - t112 * t136;
t118 = -Ifges(6,3) * t110 - t112 * t133;
t48 = t146 + t70;
t117 = qJD(5) * t48 + t236;
t116 = (-t110 * t35 - t112 * t34) * qJD(4) + t234;
t96 = -qJDD(1) * pkin(1) + qJDD(2);
t83 = -qJD(4) * mrSges(5,2) - t163;
t74 = t145 * qJD(1);
t72 = t84 * qJD(1);
t67 = t106 * t214 - t107 * t213;
t66 = -t106 * t213 - t107 * t214;
t63 = t141 * t110;
t53 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t75;
t45 = mrSges(6,1) * t88 + mrSges(6,3) * t69;
t44 = -mrSges(6,2) * t88 + mrSges(6,3) * t68;
t39 = -mrSges(5,1) * t75 + mrSges(5,2) * t76;
t27 = t109 * t67 - t189 * t66;
t26 = t111 * t67 + t191 * t66;
t20 = t109 * t48 + t189 * t71;
t19 = t111 * t48 - t191 * t71;
t18 = t109 * t74 + t111 * t34;
t17 = -t109 * t34 + t111 * t74;
t11 = -qJDD(4) * pkin(4) - t13;
t6 = t28 * Ifges(6,1) + t29 * Ifges(6,4) + t65 * Ifges(6,5);
t5 = t28 * Ifges(6,4) + t29 * Ifges(6,2) + t65 * Ifges(6,6);
t4 = -t109 * t117 + t111 * t147;
t3 = t109 * t147 + t111 * t117;
t21 = [(-m(3) * t187 - m(4) * t158 - t27 * mrSges(6,1) - t26 * mrSges(6,2) + t241 * t67 - t243 * t214 + t242 * t213 + t245 * (-t66 * pkin(3) + pkin(6) * t67 + t158) + (m(6) * t146 + mrSges(4,1) + t84) * t66) * g(2) + t22 * t154 / 0.2e1 + (qJD(5) * t23 + t5) * t192 / 0.2e1 + (t107 * t171 + t43) * mrSges(4,2) + (-Ifges(5,6) * qJDD(4) + t167 / 0.2e1 + Ifges(5,4) * t246 + Ifges(5,2) * t247 + Ifges(6,3) * t221 + Ifges(6,6) * t222 + Ifges(6,5) * t223 + t225) * t112 + t76 * t140 / 0.2e1 + m(4) * (qJD(2) * t131 + t42 * t77 + t43 * t78) + t75 * t137 / 0.2e1 + (t118 * t250 + t119 * t251 + t120 * t219 - t128 + t134 * qJD(4) / 0.2e1) * qJD(4) + (t1 * t192 + t123 * t9 + t124 * t8 + t190 * t2) * mrSges(6,3) - t6 * t190 / 0.2e1 + (t132 * t250 + t135 * t251 + t138 * t219) * t177 + (Ifges(2,3) + Ifges(4,3) + Ifges(3,2)) * qJDD(1) - t96 * mrSges(3,1) + m(3) * (-pkin(1) * t96 + (t87 + t171) * qJ(2)) - t77 * t174 + t40 * t84 + 0.2e1 * t87 * mrSges(3,3) + t70 * t39 - t11 * t63 + t3 * t44 + t4 * t45 + t19 * t15 + t20 * t16 + (mrSges(6,3) * t66 * g(2) + Ifges(5,1) * t246 + Ifges(5,4) * t247 - Ifges(5,5) * qJDD(4) - t133 * t221 - t136 * t222 - t139 * t223 + t22 * t150 + (m(6) * t11 + t215) * t71 + (-m(5) * t34 + t226) * t182) * t110 + (t106 * t171 - t42) * mrSges(4,1) + t30 * (-mrSges(6,1) * t123 - mrSges(6,2) * t124) + m(6) * (t1 * t20 + t19 * t2 + t3 * t9 + t4 * t8) + t236 * t83 + t237 * t149 + (-t8 * mrSges(6,1) + t9 * mrSges(6,2) - t230) * t181 + (-m(3) * t148 - m(4) * t122 + t242 * t214 + t243 * t213 + t245 * (t67 * pkin(3) + t122) + (-mrSges(4,1) - t224) * t67 + (pkin(6) * t245 - t141 + t241) * t66) * g(1) + m(5) * (qJD(2) * t227 + t116 * t71 + t40 * t70) + t78 * t173 + pkin(1) * t175 + t72 * t183 + t53 * t194 + (t181 * t35 - t234) * mrSges(5,3) + (t226 * t71 - t238 / 0.2e1 + t34 * mrSges(5,3)) * t179; -t175 + t56 * t15 + t57 * t16 + t240 * t45 + t239 * t44 + (-m(3) * qJ(2) - mrSges(3,3)) * t114 + (-t114 * mrSges(4,2) - t174 - t39) * t107 + (-t114 * mrSges(4,1) + t173 + t112 * t53 + t215 * t110 + (-t110 * t83 + t112 * t207) * qJD(4)) * t106 + m(5) * (t106 * t116 - t107 * t40) + m(4) * (t106 * t43 + t107 * t42) + m(3) * t96 + (t1 * t57 + t2 * t56 + (t11 * t110 + t179 * t30) * t106 + t239 * t9 + t240 * t8) * m(6) + ((-t110 * t207 - t112 * t83) * t107 - t106 * t72 - m(5) * (-t193 * t34 + t227) - t193 * t244 - m(4) * t131) * qJD(1) + (-t213 * g(1) + t214 * g(2)) * (m(3) + t231); m(4) * qJDD(3) + t231 * g(3) + ((-t109 * t45 + t111 * t44 + t83) * qJD(4) + m(5) * (t13 + t188) + m(6) * (-t172 * t8 + t180 * t9 - t11) - t215) * t112 + (t53 + (-t109 * t44 - t111 * t45) * qJD(5) + t207 * qJD(4) + m(5) * (-t34 * qJD(4) + t12) + m(6) * (qJD(4) * t30 + t248) + t233) * t110; (t133 * t88 + t136 * t68) * qJD(5) / 0.2e1 + (-t8 * (-mrSges(6,1) * t110 + mrSges(6,3) * t189) - t9 * (mrSges(6,2) * t110 + mrSges(6,3) * t191) + t128 + t120 * t218) * qJD(1) + t248 * mrSges(6,3) + (-t163 - t83) * t34 + t11 * t142 + (-pkin(4) * t11 - t17 * t8 - t18 * t9) * m(6) + t111 * t5 / 0.2e1 + t109 * t6 / 0.2e1 + (t238 / 0.2e1 + t228) * t184 + Ifges(5,6) * t75 + Ifges(5,5) * t76 - t18 * t44 - t17 * t45 - t12 * mrSges(5,2) + t13 * mrSges(5,1) - pkin(4) * t7 + (t110 * t125 - t112 * t164 + t143) * (-g(1) * t66 - g(2) * t67) + t237 * t114 / 0.2e1 + (-t45 * t176 - t44 * t178 + m(6) * ((-t9 * t109 - t8 * t111) * qJD(5) + t144) + t233) * pkin(7) + t230 * t185 + (t139 * t219 + t228) * qJD(5) + (-t160 - t226) * t35 + t224 * g(3) + Ifges(5,3) * qJDD(4) + t134 * t149 + t23 * t150 + (t88 * t118 + t68 * t119) * t249 + t132 * t221 + t135 * t222 + t138 * t223; -t30 * (-mrSges(6,1) * t69 + mrSges(6,2) * t68) + (Ifges(6,1) * t68 + t212) * t218 + t22 * t219 - t88 * (Ifges(6,5) * t68 + Ifges(6,6) * t69) / 0.2e1 - t8 * t44 + t9 * t45 - g(1) * (mrSges(6,1) * t26 - mrSges(6,2) * t27) - g(2) * ((-t111 * t66 + t191 * t67) * mrSges(6,1) + (t109 * t66 + t189 * t67) * mrSges(6,2)) - g(3) * t63 + (t68 * t8 - t69 * t9) * mrSges(6,3) + t167 - (Ifges(6,2) * t69 + t23 + t64) * t68 / 0.2e1 + t225;];
tau = t21;
