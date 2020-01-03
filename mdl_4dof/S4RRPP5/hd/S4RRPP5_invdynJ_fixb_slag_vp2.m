% Calculate vector of inverse dynamics joint torques for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:18
% DurationCPUTime: 4.36s
% Computational Cost: add. (652->273), mult. (1501->325), div. (0->0), fcn. (640->4), ass. (0->124)
t182 = Ifges(3,5) + Ifges(5,5);
t170 = -Ifges(4,4) + t182;
t183 = Ifges(5,4) + Ifges(4,5);
t169 = -Ifges(3,6) + t183;
t74 = cos(qJ(2));
t140 = t74 * mrSges(4,3);
t141 = t74 * mrSges(5,2);
t72 = sin(qJ(2));
t64 = t72 * qJ(3);
t109 = -pkin(1) - t64;
t66 = t74 * pkin(2);
t85 = t109 - t66;
t23 = t85 * qJD(1);
t134 = qJD(1) * t72;
t58 = pkin(5) * t134;
t37 = -qJD(2) * pkin(2) + qJD(3) + t58;
t128 = qJD(2) * qJ(3);
t133 = qJD(1) * t74;
t60 = pkin(5) * t133;
t40 = -t60 - t128;
t139 = pkin(2) + qJ(4);
t6 = (-t139 * t74 + t109) * qJD(1);
t191 = -m(4) * (t37 * t74 + t40 * t72) * pkin(5) - t23 * (-t72 * mrSges(4,2) - t140) - t6 * (t72 * mrSges(5,3) - t141);
t189 = Ifges(5,2) + Ifges(4,3);
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t173 = g(1) * t75 + g(2) * t73;
t184 = -m(5) - m(4);
t187 = -m(3) + t184;
t138 = t66 + t64;
t100 = mrSges(3,1) * t74 - mrSges(3,2) * t72;
t98 = t74 * mrSges(4,2) - t72 * mrSges(4,3);
t186 = t98 - t100;
t185 = pkin(3) * t133 + qJD(4);
t125 = qJD(1) * qJD(2);
t36 = qJDD(1) * t72 + t74 * t125;
t57 = Ifges(3,4) * t133;
t145 = Ifges(5,6) * t74;
t90 = t72 * Ifges(5,3) - t145;
t181 = Ifges(3,1) * t134 + qJD(1) * t90 + qJD(2) * t182 + t57;
t56 = Ifges(5,6) * t134;
t148 = Ifges(4,6) * t72;
t91 = -t74 * Ifges(4,3) - t148;
t180 = -Ifges(5,2) * t133 + qJD(1) * t91 + qJD(2) * t183 + t56;
t118 = mrSges(4,1) * t133;
t43 = -qJD(2) * mrSges(4,3) - t118;
t179 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t133 - t43;
t117 = mrSges(5,1) * t133;
t44 = qJD(2) * mrSges(5,2) + t117;
t178 = t43 - t44;
t120 = mrSges(4,1) * t134;
t177 = -mrSges(3,3) * t134 - t120 + (mrSges(3,1) - mrSges(4,2)) * qJD(2);
t115 = t72 * t125;
t126 = t74 * qJDD(1);
t55 = pkin(5) * t126;
t25 = -pkin(5) * t115 + t55;
t26 = t36 * pkin(5);
t175 = t25 * t74 + t26 * t72;
t7 = qJD(2) * (-qJD(3) + t58) - qJDD(2) * qJ(3) - t55;
t116 = qJDD(3) + t26;
t9 = -qJDD(2) * pkin(2) + t116;
t174 = -t7 * t74 + t72 * t9;
t31 = -pkin(3) * t134 - t58;
t172 = -t31 + qJD(3);
t171 = -t60 - t185;
t147 = Ifges(4,6) * t74;
t168 = t72 * t148 + (t145 - t147 + (-Ifges(4,2) + t189) * t72) * t74;
t167 = t169 * t72 + t170 * t74;
t142 = t72 * mrSges(5,2);
t166 = -mrSges(2,1) - t142 + t186;
t165 = -m(5) * pkin(3) + mrSges(2,2) - mrSges(3,3);
t136 = qJ(3) * t74;
t164 = t173 * t136;
t157 = pkin(3) + pkin(5);
t156 = pkin(5) * t72;
t155 = pkin(5) * t74;
t150 = Ifges(3,4) * t72;
t149 = Ifges(3,4) * t74;
t146 = Ifges(5,6) * t72;
t132 = qJD(2) * t72;
t130 = qJD(3) * t72;
t129 = qJDD(1) * pkin(1);
t48 = t157 * t74;
t119 = mrSges(5,1) * t134;
t14 = t36 * mrSges(4,1) + qJDD(2) * mrSges(4,2);
t35 = t115 - t126;
t13 = -t35 * mrSges(5,1) + qJDD(2) * mrSges(5,2);
t108 = -t125 / 0.2e1;
t11 = t36 * mrSges(5,1) - qJDD(2) * mrSges(5,3);
t104 = qJ(4) * t74 + t138;
t102 = -m(5) * t139 - mrSges(5,3);
t99 = mrSges(3,1) * t72 + mrSges(3,2) * t74;
t97 = -t74 * mrSges(5,3) - t142;
t96 = t74 * Ifges(3,2) + t150;
t92 = -t72 * Ifges(4,2) - t147;
t88 = qJ(4) * t72 - t136;
t87 = -qJ(3) * t36 - t129;
t86 = -qJD(4) * t74 - t130;
t84 = pkin(1) * t99;
t81 = t72 * (Ifges(3,1) * t74 - t150);
t62 = pkin(2) * t132;
t59 = pkin(2) * t134;
t47 = t157 * t72;
t42 = -qJD(2) * mrSges(5,3) + t119;
t38 = -pkin(1) - t138;
t34 = qJD(2) * t48;
t33 = t157 * t132;
t30 = t97 * qJD(1);
t29 = -qJ(3) * t133 + t59;
t28 = t98 * qJD(1);
t24 = -pkin(1) - t104;
t22 = Ifges(4,4) * qJD(2) + qJD(1) * t92;
t17 = Ifges(3,6) * qJD(2) + qJD(1) * t96;
t16 = -t40 + t185;
t15 = -t128 * t74 - t130 + t62;
t12 = mrSges(4,1) * t35 - qJDD(2) * mrSges(4,3);
t10 = qJD(1) * t88 + t59;
t8 = -qJD(2) * t139 + t172;
t5 = qJD(2) * t88 + t62 + t86;
t4 = -pkin(3) * t35 + qJDD(4) - t7;
t3 = pkin(2) * t35 - qJD(1) * t130 + t87;
t2 = pkin(3) * t36 - qJD(2) * qJD(4) - qJDD(2) * t139 + t116;
t1 = qJD(1) * t86 + t139 * t35 + t87;
t18 = [(t165 * t73 + t187 * (t75 * pkin(1) + t73 * pkin(5)) + (t184 * t138 - (m(5) * qJ(4) + mrSges(5,3)) * t74 + t166) * t75) * g(2) + ((m(3) * pkin(1) - m(4) * t85 - m(5) * t109 - t102 * t74 - t166) * t73 + (pkin(5) * t187 + t165) * t75) * g(1) + (-t16 * mrSges(5,1) - t179 * pkin(5) + t180 / 0.2e1 + t40 * mrSges(4,1) - t17 / 0.2e1) * t132 + m(4) * (t174 * pkin(5) + t15 * t23 + t3 * t38) + t100 * t129 + (-qJDD(2) * mrSges(3,2) - t12) * t155 + ((t8 * mrSges(5,1) - t177 * pkin(5) + t181 / 0.2e1 + t37 * mrSges(4,1) - t22 / 0.2e1) * t74 + t167 * qJD(2) / 0.2e1 - t191) * qJD(2) + (t72 * (Ifges(5,3) * t74 + t146) + t74 * (-Ifges(3,2) * t72 + t149) + t81) * t125 / 0.2e1 + (t72 * Ifges(3,1) + t149 + t90) * t36 / 0.2e1 + (-t74 * Ifges(5,2) + t146 + t91) * t35 / 0.2e1 + m(5) * (t1 * t24 - t16 * t33 + t2 * t47 + t34 * t8 + t4 * t48 + t5 * t6) - ((-Ifges(4,6) + Ifges(5,6)) * t36 + t189 * t35 + t183 * qJDD(2)) * t74 / 0.2e1 + ((Ifges(3,1) + Ifges(5,3)) * t36 + (-Ifges(3,4) + Ifges(5,6)) * t35 + t182 * qJDD(2)) * t72 / 0.2e1 + (-t169 * t74 + t170 * t72) * qJDD(2) / 0.2e1 + (-t155 * t35 + t156 * t36 + t175) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + t175 * pkin(5)) + t168 * t108 - t84 * t125 + t74 * (Ifges(3,4) * t36 - Ifges(3,2) * t35 + Ifges(3,6) * qJDD(2)) / 0.2e1 - t72 * (Ifges(4,4) * qJDD(2) - Ifges(4,2) * t36 + Ifges(4,6) * t35) / 0.2e1 + (-qJDD(2) * mrSges(3,1) + t14) * t156 + t47 * t11 + t48 * t13 + t38 * (-mrSges(4,2) * t35 - mrSges(4,3) * t36) + t34 * t42 - t33 * t44 + t24 * (-mrSges(5,2) * t36 + mrSges(5,3) * t35) - pkin(1) * (mrSges(3,1) * t35 + mrSges(3,2) * t36) + t15 * t28 + t5 * t30 + (t2 * t72 + t4 * t74 - t173) * mrSges(5,1) + (-t173 + t174) * mrSges(4,1) - t36 * t92 / 0.2e1 + Ifges(2,3) * qJDD(1) - t35 * t96 / 0.2e1 + t1 * t97 + t3 * t98; (-m(4) * t138 - m(5) * t104 + t186 + t97) * g(3) - t139 * t11 + (Ifges(5,1) + Ifges(4,1) + Ifges(3,3)) * qJDD(2) + ((t84 - t81 / 0.2e1 + t168 / 0.2e1) * qJD(1) + t191) * qJD(1) + t177 * t60 - t178 * qJD(3) + t179 * t58 - (Ifges(5,3) * t133 + t180 + t56) * t134 / 0.2e1 - (-Ifges(3,2) * t134 + t181 + t57) * t133 / 0.2e1 + t170 * t36 + t171 * t42 + t167 * t108 + t169 * t35 + t17 * t134 / 0.2e1 + t22 * t133 / 0.2e1 - t40 * t120 - t8 * t117 - t37 * t118 - t31 * t44 - t25 * mrSges(3,2) - t26 * mrSges(3,1) - t29 * t28 - t10 * t30 + t9 * mrSges(4,2) - pkin(2) * t14 - t7 * mrSges(4,3) - t2 * mrSges(5,3) + t4 * mrSges(5,2) + (-t141 - t140 + t99 + (m(4) * pkin(2) - mrSges(4,2) - t102) * t72) * t173 + t16 * t119 + (qJ(3) * t4 - t10 * t6 - t139 * t2 + t16 * t172 + t171 * t8 - t164) * m(5) + (-pkin(2) * t9 - qJ(3) * t7 - qJD(3) * t40 - t23 * t29 - t164) * m(4) + (-t12 + t13) * qJ(3); t178 * qJD(2) - t184 * t74 * g(3) + ((t28 + t30) * qJD(1) + t184 * t173) * t72 + t11 + t14 + (-qJD(2) * t16 + t134 * t6 + t2) * m(5) + (qJD(2) * t40 + t134 * t23 + t9) * m(4); t30 * t133 + qJD(2) * t42 + (-g(3) * t72 + t8 * qJD(2) + t133 * t6 - t173 * t74 + t4) * m(5) + t13;];
tau = t18;
