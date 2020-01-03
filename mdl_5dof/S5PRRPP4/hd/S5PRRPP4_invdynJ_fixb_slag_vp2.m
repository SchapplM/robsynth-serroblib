% Calculate vector of inverse dynamics joint torques for
% S5PRRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:40:56
% DurationCPUTime: 4.90s
% Computational Cost: add. (915->310), mult. (1869->363), div. (0->0), fcn. (840->4), ass. (0->139)
t205 = Ifges(5,1) + Ifges(6,1);
t203 = Ifges(6,4) + Ifges(5,5);
t204 = Ifges(5,4) + Ifges(4,5);
t186 = -Ifges(6,5) + t204;
t82 = pkin(7) + qJ(2);
t72 = cos(t82);
t170 = g(1) * t72;
t212 = Ifges(6,2) + Ifges(5,3);
t211 = Ifges(4,6) - Ifges(5,6);
t201 = Ifges(5,6) - Ifges(6,6);
t86 = sin(qJ(3));
t77 = t86 * qJ(4);
t118 = pkin(2) + t77;
t173 = pkin(3) + pkin(4);
t87 = cos(qJ(3));
t11 = qJD(5) + (t173 * t87 + t118) * qJD(2);
t156 = t87 * mrSges(5,3);
t157 = t87 * mrSges(6,2);
t81 = t87 * pkin(3);
t100 = -t118 - t81;
t32 = t100 * qJD(2);
t210 = t32 * (mrSges(5,1) * t86 - t156) + t11 * (-t86 * mrSges(6,1) + t157);
t159 = Ifges(5,5) * t87;
t161 = Ifges(6,4) * t87;
t209 = t205 * t86 - t159 - t161;
t71 = sin(t82);
t188 = g(2) * t71 + t170;
t174 = m(5) + m(6);
t208 = -m(4) - t174;
t143 = t86 * qJD(2);
t207 = t203 * t143;
t151 = t81 + t77;
t108 = t87 * mrSges(5,1) + t86 * mrSges(5,3);
t110 = mrSges(4,1) * t87 - mrSges(4,2) * t86;
t206 = -t108 - t110;
t137 = qJD(2) * qJD(3);
t41 = qJDD(2) * t86 + t87 * t137;
t36 = t41 * mrSges(5,2);
t18 = -qJDD(3) * mrSges(5,1) + t36;
t200 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t41 + t18;
t125 = t86 * t137;
t139 = t87 * qJDD(2);
t40 = t125 - t139;
t19 = -mrSges(5,2) * t40 + qJDD(3) * mrSges(5,3);
t199 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t40 + t19;
t142 = t87 * qJD(2);
t198 = qJD(3) * t201 - t142 * t212 + t207;
t127 = mrSges(4,3) * t143;
t131 = mrSges(5,2) * t143;
t197 = t127 + t131 + (-mrSges(4,1) - mrSges(5,1)) * qJD(3);
t128 = mrSges(6,3) * t142;
t48 = qJD(3) * mrSges(6,2) - t128;
t126 = mrSges(5,2) * t142;
t50 = qJD(3) * mrSges(5,3) + t126;
t196 = t48 + t50;
t130 = mrSges(4,3) * t142;
t49 = -qJD(3) * mrSges(4,2) + t130;
t195 = -t49 - t50;
t78 = t86 * mrSges(6,2);
t152 = t87 * mrSges(6,1) + t78;
t194 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t192 = t204 * t87 - t211 * t86;
t191 = t203 * t86;
t138 = qJD(1) * qJD(3);
t135 = pkin(6) * t139 + qJDD(1) * t86 + t138 * t87;
t7 = -pkin(6) * t125 + t135;
t8 = -pkin(6) * t41 + qJDD(1) * t87 - t138 * t86;
t190 = t7 * t87 - t8 * t86;
t5 = t7 + t194;
t92 = qJDD(4) - t8;
t6 = -qJDD(3) * pkin(3) + t92;
t189 = t5 * t87 + t6 * t86;
t141 = qJ(5) * qJD(2);
t42 = -pkin(6) * t143 + qJD(1) * t87;
t20 = t141 * t86 + t42;
t187 = -t20 + qJD(4);
t68 = Ifges(4,4) * t142;
t185 = Ifges(4,1) * t143 + t209 * qJD(2) + qJD(3) * t186 + t68;
t184 = -mrSges(3,1) - t78 + t206;
t150 = qJ(4) * t87;
t183 = t188 * t150;
t175 = m(2) + m(3);
t164 = Ifges(4,4) * t86;
t163 = Ifges(4,4) * t87;
t158 = t41 * mrSges(6,3);
t155 = pkin(6) - qJ(5);
t153 = pkin(2) * t72 + pkin(6) * t71;
t43 = pkin(6) * t142 + qJD(1) * t86;
t148 = qJD(3) * t86;
t147 = qJD(3) * t87;
t146 = qJD(4) * t86;
t145 = qJD(5) * t87;
t144 = qJDD(2) * pkin(2);
t134 = pkin(4) * t87 + t151;
t133 = pkin(6) * t148;
t129 = mrSges(6,3) * t143;
t53 = t155 * t87;
t119 = -t40 * mrSges(6,1) + mrSges(6,2) * t41;
t117 = -t137 / 0.2e1;
t116 = t137 / 0.2e1;
t115 = t151 * t72 + t153;
t114 = -m(6) * t173 - mrSges(6,1);
t21 = -t141 * t87 + t43;
t85 = qJD(3) * qJ(4);
t12 = t21 + t85;
t9 = -qJD(3) * t173 + t187;
t111 = t12 * t87 + t86 * t9;
t109 = mrSges(4,1) * t86 + mrSges(4,2) * t87;
t105 = t87 * Ifges(4,2) + t164;
t102 = Ifges(6,5) * t87 + Ifges(6,6) * t86;
t101 = pkin(3) * t86 - t150;
t99 = pkin(2) * t109;
t98 = -t173 * t86 + t150;
t95 = t86 * (Ifges(4,1) * t87 - t164);
t94 = t87 * (Ifges(6,2) * t86 + t161);
t93 = t87 * (Ifges(5,3) * t86 + t159);
t89 = qJ(4) * t41 + qJD(4) * t143 + t144;
t51 = t155 * t86;
t45 = -qJD(3) * mrSges(6,1) - t129;
t44 = -pkin(2) - t151;
t39 = t101 * qJD(2);
t38 = t152 * qJD(2);
t37 = t108 * qJD(2);
t34 = t85 + t43;
t33 = pkin(2) + t134;
t28 = Ifges(4,6) * qJD(3) + qJD(2) * t105;
t25 = -qJD(3) * pkin(3) + qJD(4) - t42;
t24 = qJD(3) * t53 - qJD(5) * t86;
t23 = qJD(3) * t101 - t146;
t22 = -t148 * t155 - t145;
t16 = -qJDD(3) * mrSges(6,1) - t158;
t14 = qJDD(3) * mrSges(6,2) + mrSges(6,3) * t40;
t13 = t98 * qJD(2);
t10 = qJD(3) * t98 + t146;
t4 = pkin(3) * t40 - t89;
t3 = qJ(5) * t40 + (-t133 - t145) * qJD(2) + t135 + t194;
t2 = -t173 * t40 + qJDD(5) + t89;
t1 = -qJ(5) * t41 - qJD(5) * t143 - qJDD(3) * t173 + t92;
t15 = [t175 * qJDD(1) + (-t16 - t200) * t87 + (t14 + t199) * t86 + ((t49 + t196) * t87 + (t45 + t197) * t86) * qJD(3) + m(4) * (t7 * t86 + t8 * t87 + (-t42 * t86 + t43 * t87) * qJD(3)) + m(5) * (t5 * t86 - t6 * t87 + (t25 * t86 + t34 * t87) * qJD(3)) + m(6) * (qJD(3) * t111 - t1 * t87 + t3 * t86) + (-t175 + t208) * g(3); (m(4) * pkin(2) ^ 2 + Ifges(3,3)) * qJDD(2) + (-t34 * mrSges(5,2) - t43 * mrSges(4,3) + t12 * mrSges(6,3) - t28 / 0.2e1 + t198 / 0.2e1) * t148 + m(6) * (t1 * t51 + t10 * t11 + t12 * t22 + t2 * t33 + t24 * t9 + t3 * t53) - t4 * t108 + (t147 * t25 - t188 + t189) * mrSges(5,2) + (-t147 * t42 - t188 + t190) * mrSges(4,3) + (-t1 * t86 - t147 * t9 - t3 * t87 + t188) * mrSges(6,3) + (t86 * Ifges(4,1) + t163 + t209) * t41 / 0.2e1 - t40 * t105 / 0.2e1 + t195 * t133 + (t192 / 0.2e1 - t102 / 0.2e1) * qJD(3) ^ 2 + t185 * t147 / 0.2e1 + (-m(6) * (-qJ(5) * t71 + t115) - m(4) * t153 + mrSges(3,2) * t71 - m(5) * t115 + (-(m(6) * pkin(4) + mrSges(6,1)) * t87 + t184) * t72) * g(2) + (t94 + t93) * t117 + t191 * t40 / 0.2e1 - (t201 * qJDD(3) + t203 * t41) * t87 / 0.2e1 + ((m(6) * qJ(5) + mrSges(3,2)) * t72 + (m(4) * pkin(2) - m(5) * t100 + m(6) * t118 - t114 * t87 - t184) * t71) * g(1) + t2 * t152 - t99 * t137 + m(5) * (t23 * t32 + t4 * t44) - qJDD(3) * (Ifges(6,5) * t86 - Ifges(6,6) * t87) / 0.2e1 + t87 * (Ifges(4,4) * t41 - Ifges(4,2) * t40 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t44 * (mrSges(5,1) * t40 - mrSges(5,3) * t41) + t24 * t45 + t22 * t48 + t51 * t16 + t53 * t14 - pkin(2) * (mrSges(4,1) * t40 + mrSges(4,2) * t41) - t212 * t87 * t40 - t23 * t37 + t10 * t38 + t33 * t119 + (m(5) * ((t25 * t87 - t34 * t86) * qJD(3) + t189) + m(4) * ((-t42 * t87 - t43 * t86) * qJD(3) + t190) + t199 * t87 + t200 * t86 + t197 * t147 + t208 * t170) * pkin(6) + t110 * t144 + ((Ifges(4,1) + t205) * t41 + (-Ifges(4,4) + t203) * t40 + t186 * qJDD(3)) * t86 / 0.2e1 + (t87 * (-Ifges(4,2) * t86 + t163) + t95 + (t205 * t87 + t191) * t86) * t116 + t210 * qJD(3) + (t204 * t86 + t211 * t87) * qJDD(3) / 0.2e1; t34 * t131 + (-m(5) * t151 - m(6) * t134 - t152 + t206) * g(3) + t9 * t128 + (-pkin(3) * t6 + qJ(4) * t5 + qJD(4) * t34 - t32 * t39 - t183) * m(5) + (qJ(4) * t3 - t1 * t173 - t11 * t13 + t12 * t187 - t21 * t9 - t183) * m(6) - (t205 * t142 + t198 + t207) * t143 / 0.2e1 + (-m(5) * t34 + t130 + t195) * t42 + t196 * qJD(4) + (-m(5) * t25 + t127 - t197) * t43 + t192 * t117 - (-Ifges(4,2) * t143 + t185 + t68) * t142 / 0.2e1 + t186 * t41 - t12 * t129 - t25 * t126 + (Ifges(6,3) + Ifges(5,2) + Ifges(4,3)) * qJDD(3) + t28 * t143 / 0.2e1 + (-t157 - t156 + t109 + (m(5) * pkin(3) + mrSges(5,1) - t114) * t86) * t188 - t21 * t45 - t20 * t48 + (t14 + t19) * qJ(4) - t13 * t38 + t39 * t37 - pkin(3) * t18 - t7 * mrSges(4,2) + t8 * mrSges(4,1) - t1 * mrSges(6,1) + t3 * mrSges(6,2) + t5 * mrSges(5,3) - t6 * mrSges(5,1) - t173 * t16 + (-Ifges(4,6) + t201) * t40 + t102 * t116 + ((-t95 / 0.2e1 + t94 / 0.2e1 + t93 / 0.2e1 + t99) * qJD(2) - t210) * qJD(2); -t158 + t36 + (-mrSges(5,1) - mrSges(6,1)) * qJDD(3) - t196 * qJD(3) + t174 * t87 * g(3) + ((-t37 - t38) * qJD(2) - t174 * t188) * t86 + (-qJD(3) * t12 - t11 * t143 + t1) * m(6) + (-qJD(3) * t34 + t143 * t32 + t6) * m(5); (t45 * t86 + t48 * t87) * qJD(2) + (g(1) * t71 - g(2) * t72 + qJD(2) * t111 + t2) * m(6) + t119;];
tau = t15;
