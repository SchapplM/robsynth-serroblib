% Calculate time derivative of joint inertia matrix for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:18
% EndTime: 2019-03-09 02:53:22
% DurationCPUTime: 1.83s
% Computational Cost: add. (3000->295), mult. (6298->461), div. (0->0), fcn. (5898->8), ass. (0->130)
t168 = -Ifges(4,1) + Ifges(4,2);
t101 = cos(pkin(10));
t99 = sin(pkin(10));
t134 = t101 ^ 2 + t99 ^ 2;
t167 = -Ifges(6,5) * t101 + Ifges(6,6) * t99 + (2 * Ifges(5,4));
t102 = sin(qJ(6));
t104 = cos(qJ(6));
t84 = t101 * t104 - t102 * t99;
t77 = t84 * qJD(6);
t105 = -pkin(1) - pkin(7);
t124 = qJ(4) - t105;
t151 = cos(qJ(3));
t115 = qJD(3) * t151;
t103 = sin(qJ(3));
t123 = qJD(3) * t103;
t166 = -mrSges(4,1) * t123 - mrSges(4,2) * t115;
t138 = t99 * mrSges(6,3);
t100 = sin(pkin(9));
t125 = cos(pkin(9));
t112 = t125 * t151;
t75 = -qJD(3) * t112 + t100 * t123;
t83 = t100 * t151 + t103 * t125;
t76 = t83 * qJD(3);
t49 = mrSges(6,2) * t75 + t138 * t76;
t128 = t101 * t76;
t50 = -mrSges(6,1) * t75 + mrSges(6,3) * t128;
t165 = t101 * t49 - t99 * t50;
t164 = t100 * t75 + t125 * t76;
t82 = t100 * t103 - t112;
t58 = -mrSges(6,2) * t83 + t138 * t82;
t127 = t101 * t82;
t59 = mrSges(6,1) * t83 + mrSges(6,3) * t127;
t163 = -t101 * t58 + t99 * t59;
t162 = qJD(5) * t134;
t93 = -pkin(3) * t125 - pkin(4);
t161 = m(6) * t93 - mrSges(6,1) * t101 + mrSges(6,2) * t99 - mrSges(5,1);
t160 = 0.2e1 * m(6);
t159 = 2 * m(7);
t158 = -2 * mrSges(5,3);
t157 = 2 * qJD(2);
t156 = m(5) * pkin(3);
t91 = pkin(3) * t100 + qJ(5);
t153 = pkin(8) + t91;
t145 = t76 * t99;
t42 = -mrSges(6,1) * t145 - mrSges(6,2) * t128;
t85 = t101 * t102 + t104 * t99;
t44 = t85 * t82;
t22 = qJD(6) * t44 - t76 * t84;
t24 = t76 * t85 + t77 * t82;
t7 = -t24 * mrSges(7,1) + t22 * mrSges(7,2);
t152 = t42 + t7;
t150 = Ifges(6,4) * t99;
t106 = -qJD(4) * t151 + t123 * t124;
t107 = t124 * t151;
t70 = -qJD(3) * t107 - t103 * qJD(4);
t40 = t100 * t70 - t125 * t106;
t86 = t124 * t103;
t64 = -t100 * t86 + t125 * t107;
t149 = t40 * t64;
t148 = t64 * t76;
t65 = -t100 * t107 - t125 * t86;
t147 = t65 * t75;
t144 = t77 * mrSges(7,3);
t78 = t85 * qJD(6);
t143 = t78 * mrSges(7,3);
t66 = t82 * t76;
t142 = t82 * t99;
t141 = t83 * mrSges(5,3);
t140 = t84 * mrSges(7,3);
t139 = t85 * mrSges(7,3);
t89 = pkin(3) * t115 + qJD(2);
t33 = -pkin(4) * t75 + qJ(5) * t76 + qJD(5) * t82 + t89;
t41 = t100 * t106 + t125 * t70;
t11 = t101 * t41 + t99 * t33;
t94 = t103 * pkin(3) + qJ(2);
t57 = pkin(4) * t83 + qJ(5) * t82 + t94;
t27 = t101 * t65 + t99 * t57;
t135 = Ifges(7,5) * t77 - Ifges(7,6) * t78;
t132 = Ifges(6,4) * t101;
t126 = t11 * t101;
t122 = t82 * t158;
t121 = Ifges(7,5) * t22 + Ifges(7,6) * t24 - Ifges(7,3) * t75;
t117 = t75 * t134;
t116 = -t75 * mrSges(5,1) - t76 * mrSges(5,2);
t10 = t101 * t33 - t41 * t99;
t26 = t101 * t57 - t65 * t99;
t111 = t40 * t82 + t148;
t110 = -t101 * Ifges(6,1) + t150;
t109 = t99 * Ifges(6,2) - t132;
t14 = pkin(5) * t83 + pkin(8) * t127 + t26;
t17 = pkin(8) * t142 + t27;
t3 = -t102 * t17 + t104 * t14;
t4 = t102 * t14 + t104 * t17;
t79 = t153 * t99;
t80 = t153 * t101;
t51 = -t102 * t80 - t104 * t79;
t52 = -t102 * t79 + t104 * t80;
t87 = -t101 * pkin(5) + t93;
t63 = Ifges(7,1) * t85 + Ifges(7,4) * t84;
t62 = Ifges(7,4) * t85 + Ifges(7,2) * t84;
t61 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t56 = Ifges(7,1) * t77 - Ifges(7,4) * t78;
t55 = Ifges(7,4) * t77 - Ifges(7,2) * t78;
t54 = mrSges(7,1) * t78 + mrSges(7,2) * t77;
t53 = (-mrSges(6,1) * t99 - mrSges(6,2) * t101) * t82;
t46 = t84 * t82;
t45 = t84 * t83;
t43 = t85 * t83;
t39 = -pkin(5) * t142 + t64;
t37 = -qJD(5) * t85 - qJD(6) * t52;
t36 = qJD(5) * t84 + qJD(6) * t51;
t35 = mrSges(7,1) * t83 + mrSges(7,3) * t46;
t34 = -mrSges(7,2) * t83 + mrSges(7,3) * t44;
t32 = -t75 * Ifges(6,5) + t110 * t76;
t31 = -t75 * Ifges(6,6) + t109 * t76;
t28 = -pkin(5) * t145 + t40;
t25 = -mrSges(7,1) * t44 - mrSges(7,2) * t46;
t23 = t75 * t85 - t77 * t83;
t21 = -t75 * t84 - t78 * t83;
t16 = -Ifges(7,1) * t46 + Ifges(7,4) * t44 + Ifges(7,5) * t83;
t15 = -Ifges(7,4) * t46 + Ifges(7,2) * t44 + Ifges(7,6) * t83;
t13 = mrSges(7,2) * t75 + mrSges(7,3) * t24;
t12 = -mrSges(7,1) * t75 - mrSges(7,3) * t22;
t9 = pkin(8) * t145 + t11;
t8 = -pkin(5) * t75 + pkin(8) * t128 + t10;
t6 = Ifges(7,1) * t22 + Ifges(7,4) * t24 - t75 * Ifges(7,5);
t5 = Ifges(7,4) * t22 + Ifges(7,2) * t24 - t75 * Ifges(7,6);
t2 = -qJD(6) * t4 - t102 * t9 + t104 * t8;
t1 = qJD(6) * t3 + t102 * t8 + t104 * t9;
t18 = [t148 * t158 + (t1 * t4 + t2 * t3 + t28 * t39) * t159 + (t10 * t26 + t11 * t27 + t149) * t160 + t31 * t142 + (Ifges(6,6) * t83 + t109 * t82) * t145 + (((m(4) + m(3)) * t157) + 0.2e1 * (mrSges(4,1) * t151 - mrSges(4,2) * t103) * qJD(3)) * qJ(2) + t83 * t121 + 0.2e1 * t94 * t116 + (-0.2e1 * Ifges(4,4) * t151 + t103 * t168) * t115 + (0.2e1 * Ifges(4,4) * t103 + t151 * t168) * t123 + (0.2e1 * Ifges(5,1) * t82 + t167 * t83) * t76 + (Ifges(7,5) * t46 - Ifges(7,6) * t44 - t167 * t82 + (-(2 * Ifges(5,2)) - (2 * Ifges(6,3)) - Ifges(7,3)) * t83) * t75 + (0.2e1 * t53 + t122) * t40 + (t103 * mrSges(4,1) + mrSges(4,2) * t151 + mrSges(3,3)) * t157 + 0.2e1 * t3 * t12 + 0.2e1 * t4 * t13 + t22 * t16 - t32 * t127 + t24 * t15 - (Ifges(6,5) * t83 + t110 * t82) * t128 + 0.2e1 * t28 * t25 + 0.2e1 * t1 * t34 + 0.2e1 * t2 * t35 + 0.2e1 * t39 * t7 - 0.2e1 * t41 * t141 + t44 * t5 - t46 * t6 + 0.2e1 * t27 * t49 + 0.2e1 * t26 * t50 + 0.2e1 * t11 * t58 + 0.2e1 * t10 * t59 + 0.2e1 * t64 * t42 + 0.2e1 * mrSges(5,3) * t147 + 0.2e1 * m(5) * (t41 * t65 + t89 * t94 + t149) + 0.2e1 * t89 * (t83 * mrSges(5,1) - t82 * mrSges(5,2)); -t43 * t12 + t45 * t13 + t21 * t34 + t23 * t35 + t165 * t83 + t152 * t82 + (t25 + t53 + t122) * t76 + (0.2e1 * t141 + t163) * t75 + m(7) * (t1 * t45 - t2 * t43 + t21 * t4 + t23 * t3 + t28 * t82 + t39 * t76) + m(6) * ((-t10 * t83 + t26 * t75) * t99 + (t11 * t83 - t27 * t75) * t101 + t111) + m(5) * (t41 * t83 + t111 - t147); 0.2e1 * m(6) * (-t117 * t83 + t66) + 0.2e1 * m(7) * (t21 * t45 - t23 * t43 + t66) + 0.2e1 * m(5) * (-t75 * t83 + t66); m(7) * (t1 * t52 + t2 * t51 + t28 * t87 + t3 * t37 + t36 * t4) + mrSges(6,3) * t126 + t1 * t140 - t4 * t143 - Ifges(4,6) * t115 + t164 * mrSges(5,3) * pkin(3) + (t100 * t156 - mrSges(5,2)) * t41 + (m(6) * (-t10 * t99 + t126) + t165) * t91 + t166 * t105 + (m(6) * (t101 * t27 - t26 * t99) - t163) * qJD(5) + (-t125 * t156 + t161) * t40 - (Ifges(6,5) * t99 + Ifges(7,5) * t85 + Ifges(6,6) * t101 + Ifges(7,6) * t84) * t75 / 0.2e1 - Ifges(4,5) * t123 - (Ifges(6,1) * t99 + t132) * t128 / 0.2e1 + t83 * t135 / 0.2e1 + t36 * t34 + t37 * t35 - t10 * t138 - t2 * t139 + t51 * t12 - t3 * t144 + t52 * t13 + t39 * t54 + t44 * t55 / 0.2e1 - t46 * t56 / 0.2e1 + t28 * t61 + t24 * t62 / 0.2e1 + t22 * t63 / 0.2e1 + Ifges(5,6) * t75 + (Ifges(6,2) * t101 + t150) * t145 / 0.2e1 - Ifges(5,5) * t76 + t77 * t16 / 0.2e1 - t78 * t15 / 0.2e1 + t84 * t5 / 0.2e1 + t85 * t6 / 0.2e1 + t87 * t7 + t93 * t42 + t99 * t32 / 0.2e1 + t101 * t31 / 0.2e1; m(6) * (-t117 * t91 + t162 * t83) + m(7) * (t21 * t52 + t23 * t51 + t36 * t45 - t37 * t43) + t21 * t140 - t45 * t143 - t23 * t139 + t43 * t144 + t82 * t54 - t164 * t156 + (m(7) * t87 + t161 + t61) * t76 + t166 + (-mrSges(6,3) * t134 + mrSges(5,2)) * t75; 0.2e1 * t87 * t54 + t77 * t63 + t85 * t56 - t78 * t62 + t84 * t55 + (t36 * t52 + t37 * t51) * t159 + 0.2e1 * (t36 * t84 - t37 * t85 - t51 * t77 - t52 * t78) * mrSges(7,3) + (t160 * t91 + 0.2e1 * mrSges(6,3)) * t162; t101 * t50 + t84 * t12 + t85 * t13 + t77 * t34 - t78 * t35 + t99 * t49 + m(7) * (t1 * t85 + t2 * t84 - t3 * t78 + t4 * t77) + m(6) * (t10 * t101 + t11 * t99) + m(5) * t89 + t116; m(7) * (t21 * t85 + t23 * t84 + t43 * t78 + t45 * t77); m(7) * (t36 * t85 + t37 * t84 - t51 * t78 + t52 * t77); (t77 * t85 - t78 * t84) * t159; m(6) * t40 + m(7) * t28 + t152; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t76; t54; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t121; mrSges(7,1) * t23 - mrSges(7,2) * t21; mrSges(7,1) * t37 - t36 * mrSges(7,2) + t135; -t54; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
