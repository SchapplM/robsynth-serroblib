% Calculate time derivative of joint inertia matrix for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:50
% EndTime: 2022-01-20 11:16:54
% DurationCPUTime: 1.26s
% Computational Cost: add. (2008->222), mult. (4467->318), div. (0->0), fcn. (3678->8), ass. (0->130)
t95 = sin(qJ(5));
t96 = sin(qJ(4));
t98 = cos(qJ(5));
t99 = cos(qJ(4));
t75 = t95 * t99 + t96 * t98;
t93 = sin(pkin(9));
t53 = t75 * t93;
t170 = -0.2e1 * t53;
t106 = t95 * t96 - t98 * t99;
t54 = t106 * t93;
t169 = -0.2e1 * t54;
t100 = cos(qJ(2));
t132 = pkin(1) * qJD(2);
t85 = t100 * t132 + qJD(3);
t97 = sin(qJ(2));
t90 = pkin(1) * t97 + qJ(3);
t168 = qJ(3) * t85 + qJD(3) * t90;
t115 = t97 * t132;
t125 = qJD(4) * t99;
t94 = cos(pkin(9));
t138 = t94 * t99;
t145 = pkin(1) * t100;
t79 = -t94 * pkin(3) - t93 * pkin(7) - pkin(2);
t65 = t79 - t145;
t116 = t96 * t115 + t65 * t125 + t85 * t138;
t139 = t94 * t96;
t117 = t90 * t139;
t140 = t93 * t99;
t118 = pkin(8) * t140;
t17 = (-t117 - t118) * qJD(4) + t116;
t39 = t90 * t138 + t96 * t65;
t22 = -t39 * qJD(4) + t99 * t115 - t85 * t139;
t141 = t93 * t96;
t119 = pkin(8) * t141;
t81 = qJD(4) * t119;
t18 = t81 + t22;
t59 = t99 * t65;
t27 = -t118 + t59 + (-t90 * t96 - pkin(4)) * t94;
t33 = -t119 + t39;
t9 = t27 * t98 - t33 * t95;
t2 = t9 * qJD(5) + t17 * t98 + t18 * t95;
t10 = t27 * t95 + t33 * t98;
t3 = -t10 * qJD(5) - t17 * t95 + t18 * t98;
t167 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t130 = qJ(3) * t96;
t70 = t99 * t79;
t34 = -t118 + t70 + (-pkin(4) - t130) * t94;
t48 = qJ(3) * t138 + t96 * t79;
t43 = -t119 + t48;
t15 = t34 * t98 - t43 * t95;
t114 = t94 * t130;
t128 = qJD(3) * t94;
t136 = t79 * t125 + t99 * t128;
t30 = (-t114 - t118) * qJD(4) + t136;
t36 = -t48 * qJD(4) - t96 * t128;
t31 = t81 + t36;
t5 = t15 * qJD(5) + t30 * t98 + t31 * t95;
t16 = t34 * t95 + t43 * t98;
t6 = -t16 * qJD(5) - t30 * t95 + t31 * t98;
t166 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t127 = qJD(4) * t93;
t162 = qJD(4) + qJD(5);
t42 = t162 * t75;
t28 = t42 * t93;
t29 = t162 * t54;
t14 = -Ifges(6,5) * t28 + Ifges(6,6) * t29;
t165 = -(-Ifges(5,5) * t96 - Ifges(5,6) * t99) * t127 - t14;
t91 = t93 ^ 2;
t92 = t94 ^ 2;
t164 = (t91 + t92) * mrSges(4,3);
t163 = qJD(3) + t85;
t161 = ((-mrSges(4,1) * t94 + mrSges(4,2) * t93 - mrSges(3,1)) * t97 - mrSges(3,2) * t100) * t132;
t160 = 2 * m(4);
t159 = 2 * m(5);
t158 = 2 * m(6);
t13 = -mrSges(6,1) * t29 - mrSges(6,2) * t28;
t157 = 0.2e1 * t13;
t32 = mrSges(6,1) * t53 - mrSges(6,2) * t54;
t156 = 0.2e1 * t32;
t44 = mrSges(6,2) * t94 - mrSges(6,3) * t53;
t155 = 0.2e1 * t44;
t45 = -mrSges(6,1) * t94 + mrSges(6,3) * t54;
t154 = 0.2e1 * t45;
t55 = (mrSges(5,1) * t99 - mrSges(5,2) * t96) * t127;
t153 = 0.2e1 * t55;
t110 = mrSges(5,1) * t96 + mrSges(5,2) * t99;
t63 = t110 * t93;
t152 = 0.2e1 * t63;
t71 = mrSges(5,2) * t94 - mrSges(5,3) * t141;
t151 = 0.2e1 * t71;
t72 = -mrSges(5,1) * t94 - mrSges(5,3) * t140;
t150 = 0.2e1 * t72;
t147 = Ifges(5,4) * t96;
t146 = Ifges(5,4) * t99;
t144 = t28 * mrSges(6,3);
t143 = t29 * mrSges(6,3);
t142 = t90 * t85;
t137 = t96 * (-Ifges(5,2) * t99 - t147) * t127;
t135 = t168 * t91;
t126 = qJD(4) * t96;
t124 = qJD(5) * t95;
t123 = qJD(5) * t98;
t122 = qJ(3) * qJD(3);
t121 = 0.2e1 * mrSges(5,3);
t120 = 0.2e1 * mrSges(6,3);
t113 = t94 * t126;
t41 = t162 * t106;
t111 = -t42 * mrSges(6,1) + t41 * mrSges(6,2);
t38 = t59 - t117;
t109 = t38 * t96 - t39 * t99;
t47 = t70 - t114;
t108 = t47 * t96 - t48 * t99;
t107 = -(-Ifges(5,6) * t94 + (-Ifges(5,2) * t96 + t146) * t93) * t99 - (-Ifges(5,5) * t94 + (Ifges(5,1) * t99 - t147) * t93) * t96;
t105 = 0.2e1 * t164;
t104 = -t45 * t124 + t98 * t144;
t103 = t91 * (-Ifges(5,1) * t96 - t146) * t125 + t165 * t94 + (Ifges(6,4) * t169 + Ifges(6,2) * t170 - Ifges(6,6) * t94) * t29 - (Ifges(6,1) * t169 + Ifges(6,4) * t170 - Ifges(6,5) * t94) * t28;
t102 = -t106 * t144 + t71 * t125 - t72 * t126 + t75 * t143 - t41 * t44 - t42 * t45;
t101 = -t165 + (t123 * t44 + t143 * t95) * pkin(4);
t89 = t91 * t122;
t88 = pkin(4) * t141;
t82 = t93 * pkin(4) * t125;
t73 = t93 * qJ(3) + t88;
t67 = (-mrSges(6,1) * t95 - mrSges(6,2) * t98) * qJD(5) * pkin(4);
t66 = t93 * qJD(3) + t82;
t61 = t91 * t142;
t60 = t93 * t90 + t88;
t49 = t93 * t85 + t82;
t35 = -qJ(3) * t113 + t136;
t21 = -t90 * t113 + t116;
t1 = [t21 * t151 + t22 * t150 + t49 * t156 + t60 * t157 + t2 * t155 + t3 * t154 + t103 + t85 * t105 + (t90 * t153 - t137 + t85 * t152 + (t109 * t121 + t107) * qJD(4)) * t93 + 0.2e1 * t161 + (t10 * t29 + t9 * t28) * t120 + (t92 * t142 + t61 + (-pkin(2) - t145) * t115) * t160 + (t10 * t2 + t3 * t9 + t49 * t60) * t158 + (t21 * t39 + t22 * t38 + t61) * t159; (-t137 + t163 * t63 + (qJ(3) + t90) * t55 + (((-t39 - t48) * t99 + (t38 + t47) * t96) * mrSges(5,3) + t107) * qJD(4)) * t93 + ((t10 + t16) * t29 - (-t15 - t9) * t28) * mrSges(6,3) + t103 + m(4) * (-pkin(2) * t115 + t168 * t92 + t135) + m(6) * (t10 * t5 + t15 * t3 + t16 * t2 + t49 * t73 + t6 * t9 + t60 * t66) + m(5) * (t21 * t48 + t22 * t47 + t35 * t39 + t36 * t38 + t135) + (t36 + t22) * t72 + (t35 + t21) * t71 + (t6 + t3) * t45 + (t5 + t2) * t44 + (t49 + t66) * t32 + (t73 + t60) * t13 + t161 + t163 * t164; (t35 * t48 + t36 * t47 + t89) * t159 + (t92 * t122 + t89) * t160 + (t15 * t6 + t16 * t5 + t66 * t73) * t158 + t35 * t151 + t36 * t150 + t73 * t157 + t66 * t156 + t5 * t155 + t6 * t154 + t103 + qJD(3) * t105 + (qJ(3) * t153 + qJD(3) * t152 - t137 + (t108 * t121 + t107) * qJD(4)) * t93 + (t15 * t28 + t16 * t29) * t120; m(4) * t115 + m(6) * (-t10 * t41 - t106 * t3 + t2 * t75 - t42 * t9) + m(5) * (-t109 * qJD(4) + t21 * t96 + t22 * t99) + t102; m(6) * (-t106 * t6 - t15 * t42 - t16 * t41 + t5 * t75) + m(5) * (-t108 * qJD(4) + t35 * t96 + t36 * t99) + t102; (t106 * t42 - t41 * t75) * t158; t22 * mrSges(5,1) - t21 * mrSges(5,2) + (m(6) * (t10 * t123 - t9 * t124 + t2 * t95 + t3 * t98) + t104) * pkin(4) + t101 + t167; t36 * mrSges(5,1) - t35 * mrSges(5,2) + (m(6) * (t16 * t123 - t15 * t124 + t5 * t95 + t6 * t98) + t104) * pkin(4) + t101 + t166; -t110 * qJD(4) + m(6) * (-t41 * t95 - t42 * t98 + (t106 * t95 + t75 * t98) * qJD(5)) * pkin(4) + t111; 0.2e1 * t67; t14 + t167; t14 + t166; t111; t67; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
