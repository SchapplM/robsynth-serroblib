% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:37
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:37:41
% EndTime: 2018-11-23 15:37:43
% DurationCPUTime: 2.07s
% Computational Cost: add. (3183->264), mult. (6834->362), div. (0->0), fcn. (4478->8), ass. (0->128)
t186 = -Ifges(6,1) / 0.2e1;
t84 = sin(pkin(9)) * pkin(1) + qJ(3);
t184 = qJD(1) * t84;
t185 = m(4) * t184;
t98 = sin(qJ(6));
t99 = cos(qJ(6));
t115 = Ifges(7,5) * t99 - Ifges(7,6) * t98;
t155 = Ifges(7,4) * t99;
t117 = -Ifges(7,2) * t98 + t155;
t156 = Ifges(7,4) * t98;
t119 = Ifges(7,1) * t99 - t156;
t120 = mrSges(7,1) * t98 + mrSges(7,2) * t99;
t159 = sin(qJ(5));
t160 = cos(qJ(5));
t83 = -cos(pkin(9)) * pkin(1) - pkin(2) - qJ(4);
t180 = qJD(1) * t83;
t69 = qJD(3) + t180;
t94 = sin(pkin(10));
t96 = cos(pkin(10));
t49 = -qJD(2) * t94 + t96 * t69;
t44 = -pkin(7) * qJD(1) * t96 + t49;
t139 = qJD(1) * t94;
t50 = t96 * qJD(2) + t94 * t69;
t45 = -pkin(7) * t139 + t50;
t24 = t159 * t44 + t160 * t45;
t22 = qJD(5) * pkin(8) + t24;
t73 = qJD(4) + t184;
t62 = pkin(4) * t139 + t73;
t71 = t159 * t96 + t160 * t94;
t65 = t71 * qJD(1);
t70 = t159 * t94 - t160 * t96;
t66 = t70 * qJD(1);
t28 = pkin(5) * t65 + pkin(8) * t66 + t62;
t5 = -t22 * t98 + t28 * t99;
t6 = t22 * t99 + t28 * t98;
t124 = t5 * t99 + t6 * t98;
t168 = t99 / 0.2e1;
t48 = qJD(5) * t98 - t66 * t99;
t157 = Ifges(7,4) * t48;
t47 = qJD(5) * t99 + t66 * t98;
t61 = qJD(6) + t65;
t17 = Ifges(7,2) * t47 + Ifges(7,6) * t61 + t157;
t170 = t61 / 0.2e1;
t171 = t48 / 0.2e1;
t172 = t47 / 0.2e1;
t46 = Ifges(7,4) * t47;
t18 = Ifges(7,1) * t48 + Ifges(7,5) * t61 + t46;
t23 = -t159 * t45 + t160 * t44;
t21 = -qJD(5) * pkin(5) - t23;
t183 = -t124 * mrSges(7,3) + t18 * t168 - t98 * t17 / 0.2e1 + t115 * t170 + t117 * t172 + t119 * t171 + t21 * t120;
t60 = Ifges(6,4) * t65;
t177 = t66 * t186 - t60 / 0.2e1;
t182 = t62 * mrSges(6,2) + Ifges(6,5) * qJD(5) + t177 + t183;
t144 = t94 ^ 2 + t96 ^ 2;
t181 = mrSges(5,3) * t144;
t179 = t70 * qJD(4);
t123 = -t5 * t98 + t6 * t99;
t106 = t71 * qJD(4);
t14 = -qJD(1) * t106 + qJD(5) * t23;
t138 = qJD(1) * qJD(3);
t56 = qJD(5) * t65;
t57 = qJD(5) * t66;
t35 = -pkin(5) * t57 + pkin(8) * t56 + t138;
t1 = qJD(6) * t5 + t14 * t99 + t35 * t98;
t2 = -qJD(6) * t6 - t14 * t98 + t35 * t99;
t126 = t1 * t99 - t2 * t98;
t31 = qJD(6) * t47 - t56 * t99;
t32 = -qJD(6) * t48 + t56 * t98;
t178 = t2 * mrSges(7,1) - t1 * mrSges(7,2) + Ifges(7,5) * t31 + Ifges(7,6) * t32;
t33 = -mrSges(7,2) * t61 + mrSges(7,3) * t47;
t34 = mrSges(7,1) * t61 - mrSges(7,3) * t48;
t110 = t98 * t33 + t99 * t34;
t176 = -m(7) * t124 - t110;
t111 = t99 * t33 - t98 * t34;
t51 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t65;
t108 = t111 + t51;
t175 = m(7) * t123 + t108;
t19 = -mrSges(7,1) * t57 - mrSges(7,3) * t31;
t20 = mrSges(7,2) * t57 + mrSges(7,3) * t32;
t113 = -t98 * t19 + t99 * t20;
t174 = m(7) * (qJD(6) * t124 - t126) + t110 * qJD(6) - m(6) * t14 - t57 * mrSges(6,3) - t113;
t169 = t98 / 0.2e1;
t161 = -pkin(7) + t83;
t15 = -qJD(1) * t179 + t24 * qJD(5);
t63 = t161 * t94;
t64 = t161 * t96;
t37 = t159 * t63 - t160 * t64;
t154 = t15 * t37;
t153 = t15 * t70;
t152 = t15 * t71;
t148 = t65 * Ifges(6,2);
t146 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t47 + mrSges(7,2) * t48 - t66 * mrSges(6,3);
t122 = mrSges(5,1) * t94 + mrSges(5,2) * t96;
t145 = mrSges(6,1) * t65 - mrSges(6,2) * t66 + qJD(1) * t122;
t143 = m(5) * qJD(4);
t74 = t94 * pkin(4) + t84;
t135 = -t57 * mrSges(6,1) - t56 * mrSges(6,2);
t11 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t134 = -t56 * mrSges(6,3) + t11;
t133 = qJD(5) * t160;
t132 = qJD(5) * t159;
t125 = t1 * t98 + t2 * t99;
t121 = -mrSges(7,1) * t99 + mrSges(7,2) * t98;
t118 = Ifges(7,1) * t98 + t155;
t116 = Ifges(7,2) * t99 + t156;
t114 = Ifges(7,5) * t98 + Ifges(7,6) * t99;
t67 = -t132 * t96 - t133 * t94;
t68 = -t132 * t94 + t133 * t96;
t112 = t23 * t67 + t24 * t68;
t36 = pkin(5) * t71 + pkin(8) * t70 + t74;
t38 = t159 * t64 + t160 * t63;
t12 = t36 * t99 - t38 * t98;
t13 = t36 * t98 + t38 * t99;
t109 = -t49 * t96 - t50 * t94;
t103 = t6 * mrSges(7,2) - t61 * Ifges(7,3) - t48 * Ifges(7,5) - t47 * Ifges(7,6) + Ifges(6,6) * qJD(5) - t148 / 0.2e1 - Ifges(6,4) * t66 - t5 * mrSges(7,1) - t62 * mrSges(6,1);
t100 = qJD(1) ^ 2;
t54 = Ifges(7,3) * t57;
t43 = -pkin(5) * t66 + pkin(8) * t65;
t41 = pkin(5) * t68 - pkin(8) * t67 + qJD(3);
t26 = qJD(5) * t38 - t179;
t25 = -qJD(5) * t37 - t106;
t10 = t23 * t99 + t43 * t98;
t9 = -t23 * t98 + t43 * t99;
t8 = t31 * Ifges(7,1) + t32 * Ifges(7,4) - t57 * Ifges(7,5);
t7 = t31 * Ifges(7,4) + t32 * Ifges(7,2) - t57 * Ifges(7,6);
t4 = -qJD(6) * t13 - t25 * t98 + t41 * t99;
t3 = qJD(6) * t12 + t25 * t99 + t41 * t98;
t16 = [t74 * t135 + t25 * t51 + t3 * t33 + t4 * t34 + t37 * t11 + t12 * t19 + t13 * t20 + (t148 / 0.2e1 - t103) * t68 + (t177 + t182) * t67 + t146 * t26 + 0.2e1 * qJD(4) * t181 * qJD(1) + (-t37 * t56 + t38 * t57 - t112) * mrSges(6,3) + m(6) * (t14 * t38 - t23 * t26 + t24 * t25 + t154) + m(7) * (t1 * t13 + t12 * t2 + t21 * t26 + t3 * t6 + t4 * t5 + t154) + (-t144 * t180 + t109) * t143 + (Ifges(6,4) * t56 - t54 / 0.2e1 + mrSges(6,1) * t138 - t14 * mrSges(6,3) - (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t57 + t178) * t71 + (t7 * t169 - t99 * t8 / 0.2e1 - t31 * t119 / 0.2e1 - t32 * t117 / 0.2e1 + Ifges(6,1) * t56 - mrSges(6,2) * t138 + (-mrSges(6,3) - t120) * t15 + t125 * mrSges(7,3) + (mrSges(7,3) * t123 + t114 * t170 + t116 * t172 + t118 * t171 + t121 * t21 + t168 * t17 + t169 * t18) * qJD(6) - (-t115 / 0.2e1 + Ifges(6,4)) * t57) * t70 + (t145 + ((2 * mrSges(4,3)) + t122) * qJD(1) + m(5) * (t73 + t184) + 0.2e1 * t185 + m(6) * (qJD(1) * t74 + t62)) * qJD(3); t134 * t71 + t146 * t68 + m(7) * (t21 * t68 + t152) + m(6) * (-t23 * t68 + t152) + t174 * t70 + (m(6) * t24 + t175) * t67; -t100 * mrSges(4,3) + t134 * t70 - t146 * t67 + m(7) * (-t21 * t67 + t153) + m(6) * (t112 + t153) - t174 * t71 + (-m(5) * t73 - m(6) * t62 - t143 * t144 - t145 + t176 - t185) * qJD(1) + t175 * t68; t99 * t19 + t98 * t20 + t146 * t66 + t111 * qJD(6) + (m(5) + m(6)) * t138 - t100 * t181 + t108 * t65 - m(6) * (t23 * t66 - t24 * t65) - m(5) * t109 * qJD(1) + t135 + (t61 * t123 + t21 * t66 + t125) * m(7); t31 * t118 / 0.2e1 + t32 * t116 / 0.2e1 + t7 * t168 + t8 * t169 - Ifges(6,5) * t56 - t23 * t51 - t10 * t33 - t9 * t34 - t14 * mrSges(6,2) - pkin(5) * t11 - t146 * t24 + (-mrSges(6,1) + t121) * t15 + t126 * mrSges(7,3) - (t24 * mrSges(6,3) + t103) * t66 - (t60 / 0.2e1 + t23 * mrSges(6,3) - (t186 + Ifges(6,2) / 0.2e1) * t66 - t182) * t65 + t183 * qJD(6) - (t114 / 0.2e1 - Ifges(6,6)) * t57 + (-pkin(5) * t15 - t10 * t6 - t21 * t24 - t5 * t9) * m(7) + (m(7) * t126 + qJD(6) * t176 + t113) * pkin(8); -t54 - t21 * (mrSges(7,1) * t48 + mrSges(7,2) * t47) - t48 * (Ifges(7,1) * t47 - t157) / 0.2e1 + t17 * t171 - t61 * (Ifges(7,5) * t47 - Ifges(7,6) * t48) / 0.2e1 - t5 * t33 + t6 * t34 + (t47 * t5 + t48 * t6) * mrSges(7,3) - (-Ifges(7,2) * t48 + t18 + t46) * t47 / 0.2e1 + t178;];
tauc  = t16(:);
