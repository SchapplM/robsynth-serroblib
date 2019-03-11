% Calculate time derivative of joint inertia matrix for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:01
% EndTime: 2019-03-08 19:43:03
% DurationCPUTime: 1.43s
% Computational Cost: add. (1823->250), mult. (4512->370), div. (0->0), fcn. (4368->10), ass. (0->118)
t155 = m(6) + m(5);
t77 = cos(qJ(6));
t113 = qJD(6) * t77;
t143 = cos(qJ(4));
t71 = sin(pkin(11));
t73 = cos(pkin(11));
t75 = sin(qJ(4));
t53 = t143 * t71 + t75 * t73;
t49 = t53 * qJD(4);
t107 = t143 * t73;
t133 = t71 * t75;
t52 = -t107 + t133;
t74 = sin(qJ(6));
t87 = t52 * t113 + t74 * t49;
t154 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t118 = t71 ^ 2 + t73 ^ 2;
t117 = cos(pkin(6));
t72 = sin(pkin(6));
t76 = sin(qJ(2));
t132 = t72 * t76;
t47 = t117 * t71 + t73 * t132;
t84 = -t117 * t73 + t71 * t132;
t152 = t73 * t47 + t71 * t84;
t146 = pkin(4) + pkin(9);
t106 = qJD(4) * t143;
t48 = qJD(4) * t133 - t73 * t106;
t91 = qJ(5) * t48 - qJD(5) * t53;
t13 = t146 * t49 + t91;
t120 = pkin(8) + qJ(3);
t58 = t120 * t71;
t59 = t120 * t73;
t38 = t143 * t59 - t75 * t58;
t25 = t53 * qJD(3) + t38 * qJD(4);
t15 = -t48 * pkin(5) + t25;
t66 = -pkin(3) * t73 - pkin(2);
t89 = -qJ(5) * t53 + t66;
t21 = t146 * t52 + t89;
t37 = t143 * t58 + t59 * t75;
t26 = pkin(5) * t53 + t37;
t6 = -t21 * t74 + t26 * t77;
t1 = t6 * qJD(6) + t13 * t77 + t15 * t74;
t7 = t21 * t77 + t26 * t74;
t2 = -t7 * qJD(6) - t13 * t74 + t15 * t77;
t151 = t1 * t74 + t2 * t77;
t78 = cos(qJ(2));
t115 = qJD(2) * t78;
t110 = t72 * t115;
t81 = t143 * t84;
t11 = qJD(4) * t81 - t107 * t110 + (qJD(4) * t47 + t71 * t110) * t75;
t29 = t143 * t47 - t75 * t84;
t12 = t29 * qJD(4) + t53 * t110;
t150 = t11 * t52 + t12 * t53 - t29 * t49;
t24 = (qJD(3) * t71 + qJD(4) * t59) * t75 - qJD(3) * t107 + t58 * t106;
t131 = t72 * t78;
t28 = t47 * t75 + t81;
t18 = t74 * t131 + t28 * t77;
t116 = qJD(2) * t76;
t111 = t72 * t116;
t88 = t77 * t131 - t28 * t74;
t3 = t88 * qJD(6) - t74 * t111 + t12 * t77;
t4 = t18 * qJD(6) + t77 * t111 + t12 * t74;
t149 = (t18 * t74 + t77 * t88) * qJD(6) - t3 * t77 - t4 * t74;
t148 = t74 / 0.2e1;
t147 = t77 / 0.2e1;
t142 = mrSges(7,3) * t52;
t141 = Ifges(7,4) * t74;
t140 = Ifges(7,4) * t77;
t5 = t29 * t11;
t136 = t48 * mrSges(6,1);
t135 = t53 * Ifges(7,6);
t134 = t72 ^ 2 * t76;
t98 = Ifges(7,1) * t74 + t140;
t23 = t53 * Ifges(7,5) + t98 * t52;
t129 = t74 * t23;
t62 = Ifges(7,1) * t77 - t141;
t127 = t74 * t62;
t97 = Ifges(7,2) * t77 + t141;
t22 = t97 * t52 + t135;
t126 = t77 * t22;
t125 = t77 * t49;
t61 = -Ifges(7,2) * t74 + t140;
t124 = t77 * t61;
t122 = -mrSges(6,2) + mrSges(5,1);
t121 = -mrSges(6,3) + mrSges(5,2);
t114 = qJD(6) * t74;
t112 = 0.2e1 * t49;
t109 = t52 * t114;
t104 = t87 * Ifges(7,5) + Ifges(7,6) * t125 - Ifges(7,3) * t48;
t99 = t6 * t74 - t7 * t77;
t60 = mrSges(7,1) * t74 + mrSges(7,2) * t77;
t96 = -Ifges(7,5) * t74 - Ifges(7,6) * t77;
t94 = -t24 * t38 + t25 * t37;
t34 = mrSges(7,1) * t53 - t74 * t142;
t35 = -mrSges(7,2) * t53 + t77 * t142;
t93 = -t74 * t34 + t77 * t35;
t92 = -qJ(5) * t11 + t29 * qJD(5);
t86 = t109 - t125;
t82 = -t11 * t38 + t12 * t37 - t24 * t29 + t25 * t28;
t80 = m(7) * t149;
t57 = t98 * qJD(6);
t56 = t97 * qJD(6);
t55 = -mrSges(7,1) * t113 + mrSges(7,2) * t114;
t45 = t48 * mrSges(5,2);
t44 = t48 * mrSges(6,3);
t36 = -mrSges(6,2) * t52 - mrSges(6,3) * t53;
t33 = (-mrSges(7,1) * t77 + mrSges(7,2) * t74) * t52;
t32 = pkin(4) * t52 + t89;
t31 = t49 * mrSges(5,1) - t45;
t30 = -t49 * mrSges(6,2) + t44;
t27 = -t52 * pkin(5) + t38;
t20 = pkin(4) * t49 + t91;
t17 = -mrSges(7,1) * t48 - t87 * mrSges(7,3);
t16 = mrSges(7,2) * t48 - t86 * mrSges(7,3);
t14 = -pkin(5) * t49 - t24;
t10 = t86 * mrSges(7,1) + t87 * mrSges(7,2);
t9 = t87 * Ifges(7,1) - t86 * Ifges(7,4) - t48 * Ifges(7,5);
t8 = t87 * Ifges(7,4) - t86 * Ifges(7,2) - t48 * Ifges(7,6);
t19 = [0.2e1 * m(7) * (t18 * t3 - t4 * t88 - t5) + 0.2e1 * m(4) * (t152 * t72 - t134) * t115 + 0.2e1 * t155 * (-t115 * t134 + t12 * t28 - t5); -t11 * t33 + t3 * t34 + t4 * t35 + t29 * t10 + t18 * t17 - t88 * t16 + m(7) * (-t1 * t88 - t11 * t27 + t14 * t29 + t18 * t2 + t3 * t6 + t4 * t7) + m(5) * t82 + m(6) * ((t32 * t116 - t20 * t78) * t72 + t82) + m(4) * (t152 * qJD(3) + (t118 * t78 * qJ(3) - pkin(2) * t76) * t72 * qJD(2)) - t28 * t136 + (-t31 - t30) * t131 + t150 * mrSges(6,1) + (m(5) * t66 - mrSges(4,1) * t73 + mrSges(5,1) * t52 + mrSges(4,2) * t71 + mrSges(5,2) * t53 - mrSges(3,1) + t36) * t111 + (-t28 * t48 + t150) * mrSges(5,3) + (t118 * mrSges(4,3) - mrSges(3,2)) * t110; 0.2e1 * t1 * t35 + 0.2e1 * t27 * t10 + 0.2e1 * t14 * t33 + 0.2e1 * t7 * t16 + 0.2e1 * t6 * t17 + 0.2e1 * t2 * t34 + 0.2e1 * t20 * t36 + 0.2e1 * t32 * t30 + 0.2e1 * t66 * t31 - t48 * t37 * t154 + (-t154 * t38 + t126 + t129) * t49 + 0.2e1 * m(5) * t94 + 0.2e1 * m(6) * (t20 * t32 + t94) + 0.2e1 * m(7) * (t1 * t7 + t14 * t27 + t2 * t6) + ((-Ifges(5,4) - Ifges(6,6)) * t112 + t25 * t154 + (-(2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t48 + t104) * t53 + (t74 * t9 + t77 * t8 + (Ifges(5,2) + Ifges(6,3)) * t112 + t24 * t154 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) + t96) * t48 + (t77 * t23 + (-t22 - t135) * t74) * qJD(6)) * t52 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * t118 * qJD(3); m(7) * (-t74 * t3 + t77 * t4 + (-t18 * t77 + t74 * t88) * qJD(6)) + (m(4) + t155) * t111; t77 * t16 - t74 * t17 + t44 - t45 + t122 * t49 + (-t77 * t34 - t74 * t35) * qJD(6) + m(7) * (t1 * t77 - t2 * t74 + (-t6 * t77 - t7 * t74) * qJD(6)) + m(6) * t20; 0; -t29 * t55 - t122 * t12 + (-t60 + t121) * t11 + m(7) * t92 + m(6) * (-pkin(4) * t12 + t92) + t146 * t80 + t149 * mrSges(7,3); qJ(5) * t10 + qJD(5) * t33 + t14 * t60 - t27 * t55 - t122 * t25 + t121 * t24 + (t9 / 0.2e1 - t2 * mrSges(7,3) - t146 * t17) * t77 + (-t8 / 0.2e1 - t146 * t16 - t1 * mrSges(7,3)) * t74 + m(6) * (-pkin(4) * t25 - qJ(5) * t24 + qJD(5) * t38) + m(7) * (qJ(5) * t14 + qJD(5) * t27 - t146 * t151) + (-qJD(5) * mrSges(6,1) - t56 * t147 - t57 * t148) * t52 + (-Ifges(7,5) * t77 / 0.2e1 + Ifges(7,6) * t148 + Ifges(6,4) - Ifges(5,5) + pkin(4) * mrSges(6,1)) * t48 + (-Ifges(5,6) + Ifges(6,5) + t124 / 0.2e1 + t127 / 0.2e1 - qJ(5) * mrSges(6,1)) * t49 + (t53 * t96 / 0.2e1 - t126 / 0.2e1 - t129 / 0.2e1 + (t62 * t147 - t74 * t61 / 0.2e1) * t52 + t99 * mrSges(7,3) - (-m(7) * t99 + t93) * t146) * qJD(6); 0; -0.2e1 * qJ(5) * t55 + t56 * t74 - t57 * t77 + (-t124 - t127) * qJD(6) + 0.2e1 * (mrSges(6,3) + t60 + (m(6) + m(7)) * qJ(5)) * qJD(5); m(6) * t12 - t80; -t136 + t74 * t16 + t77 * t17 + t93 * qJD(6) + m(7) * (-qJD(6) * t99 + t151) + m(6) * t25; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t109 + t104; t55; ((mrSges(7,2) * t146 - Ifges(7,6)) * t77 + (mrSges(7,1) * t146 - Ifges(7,5)) * t74) * qJD(6); -t60 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t19(1) t19(2) t19(4) t19(7) t19(11) t19(16); t19(2) t19(3) t19(5) t19(8) t19(12) t19(17); t19(4) t19(5) t19(6) t19(9) t19(13) t19(18); t19(7) t19(8) t19(9) t19(10) t19(14) t19(19); t19(11) t19(12) t19(13) t19(14) t19(15) t19(20); t19(16) t19(17) t19(18) t19(19) t19(20) t19(21);];
Mq  = res;
