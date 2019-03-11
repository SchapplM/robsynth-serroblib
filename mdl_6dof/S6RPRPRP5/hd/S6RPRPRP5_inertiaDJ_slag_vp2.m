% Calculate time derivative of joint inertia matrix for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:34
% EndTime: 2019-03-09 03:14:41
% DurationCPUTime: 2.64s
% Computational Cost: add. (4011->331), mult. (9255->479), div. (0->0), fcn. (9035->8), ass. (0->133)
t193 = Ifges(7,4) + Ifges(6,5);
t191 = Ifges(7,6) - Ifges(6,6);
t192 = Ifges(7,2) + Ifges(6,3);
t189 = m(6) + m(7);
t131 = sin(pkin(9));
t175 = pkin(7) + qJ(2);
t120 = t175 * t131;
t133 = cos(pkin(9));
t122 = t175 * t133;
t135 = sin(qJ(3));
t137 = cos(qJ(3));
t188 = -t137 * t120 - t122 * t135;
t130 = sin(pkin(10));
t132 = cos(pkin(10));
t134 = sin(qJ(5));
t136 = cos(qJ(5));
t116 = t130 * t134 - t136 * t132;
t110 = t116 * qJD(5);
t118 = t130 * t136 + t132 * t134;
t111 = t118 * qJD(5);
t59 = pkin(5) * t111 + qJ(6) * t110 - qJD(6) * t118;
t75 = t111 * mrSges(7,1) + t110 * mrSges(7,3);
t76 = t111 * mrSges(6,1) - t110 * mrSges(6,2);
t187 = -m(7) * t59 - t75 - t76;
t186 = -t110 * t193 + t191 * t111;
t185 = m(7) * qJ(6) + mrSges(7,3);
t119 = t131 * t137 + t135 * t133;
t113 = t119 * qJD(3);
t117 = t131 * t135 - t133 * t137;
t112 = t117 * qJD(3);
t162 = t112 * t132;
t60 = pkin(3) * t113 + qJ(4) * t112 - qJD(4) * t119;
t64 = -t117 * qJD(2) + qJD(3) * t188;
t27 = -t130 * t64 + t132 * t60;
t17 = pkin(4) * t113 + pkin(8) * t162 + t27;
t160 = t119 * t132;
t151 = -pkin(2) * t133 - pkin(1);
t82 = pkin(3) * t117 - qJ(4) * t119 + t151;
t95 = -t135 * t120 + t122 * t137;
t47 = -t130 * t95 + t132 * t82;
t29 = pkin(4) * t117 - pkin(8) * t160 + t47;
t161 = t119 * t130;
t48 = t130 * t82 + t132 * t95;
t35 = -pkin(8) * t161 + t48;
t171 = t134 * t29 + t136 * t35;
t163 = t112 * t130;
t28 = t130 * t60 + t132 * t64;
t20 = pkin(8) * t163 + t28;
t4 = -qJD(5) * t171 - t134 * t20 + t136 * t17;
t184 = 2 * m(5);
t183 = 2 * m(6);
t182 = 0.2e1 * m(7);
t128 = t132 ^ 2;
t179 = t113 / 0.2e1;
t65 = qJD(2) * t119 + qJD(3) * t95;
t176 = t65 * t188;
t174 = pkin(8) + qJ(4);
t43 = -t111 * t119 + t112 * t116;
t23 = mrSges(6,1) * t113 - mrSges(6,3) * t43;
t24 = -t113 * mrSges(7,1) + t43 * mrSges(7,2);
t173 = -t23 + t24;
t154 = qJD(5) * t136;
t155 = qJD(5) * t134;
t44 = -t112 * t118 + t154 * t160 - t155 * t161;
t25 = -mrSges(6,2) * t113 - mrSges(6,3) * t44;
t26 = -mrSges(7,2) * t44 + mrSges(7,3) * t113;
t172 = t25 + t26;
t69 = t118 * t119;
t53 = -mrSges(6,2) * t117 - mrSges(6,3) * t69;
t56 = -mrSges(7,2) * t69 + mrSges(7,3) * t117;
t170 = t53 + t56;
t70 = t116 * t119;
t54 = mrSges(6,1) * t117 + mrSges(6,3) * t70;
t55 = -mrSges(7,1) * t117 - mrSges(7,2) * t70;
t169 = -t54 + t55;
t68 = -mrSges(5,1) * t163 - mrSges(5,2) * t162;
t168 = Ifges(5,4) * t130;
t167 = Ifges(5,4) * t132;
t166 = t113 * Ifges(5,5);
t165 = t113 * Ifges(5,6);
t164 = t130 * Ifges(5,2);
t153 = t134 * qJD(4);
t152 = t136 * qJD(4);
t125 = -pkin(4) * t132 - pkin(3);
t15 = t44 * mrSges(6,1) + t43 * mrSges(6,2);
t14 = t44 * mrSges(7,1) - t43 * mrSges(7,3);
t121 = t174 * t132;
t148 = qJD(5) * t174;
t62 = t132 * t152 - t121 * t155 + (-t136 * t148 - t153) * t130;
t63 = t132 * t153 + t121 * t154 + (-t134 * t148 + t152) * t130;
t149 = t174 * t130;
t92 = t134 * t121 + t136 * t149;
t94 = t136 * t121 - t134 * t149;
t150 = t94 * t62 + t63 * t92;
t147 = t113 * mrSges(4,1) - t112 * mrSges(4,2);
t66 = pkin(4) * t161 - t188;
t144 = Ifges(5,5) * t132 - Ifges(5,6) * t130;
t8 = -t134 * t35 + t136 * t29;
t141 = t192 * t113 + t191 * t44 + t193 * t43;
t3 = t134 * t17 + t136 * t20 + t29 * t154 - t155 * t35;
t49 = -pkin(4) * t163 + t65;
t90 = Ifges(6,1) * t118 - Ifges(6,4) * t116;
t89 = Ifges(7,1) * t118 + Ifges(7,5) * t116;
t88 = Ifges(6,4) * t118 - Ifges(6,2) * t116;
t87 = Ifges(7,5) * t118 + Ifges(7,3) * t116;
t86 = mrSges(7,1) * t116 - mrSges(7,3) * t118;
t84 = mrSges(5,1) * t117 - mrSges(5,3) * t160;
t83 = -mrSges(5,2) * t117 - mrSges(5,3) * t161;
t81 = pkin(5) * t116 - qJ(6) * t118 + t125;
t80 = -Ifges(6,1) * t110 - Ifges(6,4) * t111;
t79 = -Ifges(7,1) * t110 + Ifges(7,5) * t111;
t78 = -Ifges(6,4) * t110 - Ifges(6,2) * t111;
t77 = -Ifges(7,5) * t110 + Ifges(7,3) * t111;
t74 = mrSges(5,1) * t113 + mrSges(5,3) * t162;
t73 = -mrSges(5,2) * t113 + mrSges(5,3) * t163;
t51 = t166 - (t132 * Ifges(5,1) - t168) * t112;
t50 = t165 - (-t164 + t167) * t112;
t45 = mrSges(7,1) * t69 + mrSges(7,3) * t70;
t33 = -Ifges(6,1) * t70 - Ifges(6,4) * t69 + Ifges(6,5) * t117;
t32 = -Ifges(7,1) * t70 + Ifges(7,4) * t117 + Ifges(7,5) * t69;
t31 = -Ifges(6,4) * t70 - Ifges(6,2) * t69 + Ifges(6,6) * t117;
t30 = -Ifges(7,5) * t70 + Ifges(7,6) * t117 + Ifges(7,3) * t69;
t19 = pkin(5) * t69 + qJ(6) * t70 + t66;
t13 = Ifges(6,1) * t43 - Ifges(6,4) * t44 + t113 * Ifges(6,5);
t12 = Ifges(7,1) * t43 + t113 * Ifges(7,4) + Ifges(7,5) * t44;
t11 = Ifges(6,4) * t43 - Ifges(6,2) * t44 + t113 * Ifges(6,6);
t10 = Ifges(7,5) * t43 + t113 * Ifges(7,6) + Ifges(7,3) * t44;
t7 = -pkin(5) * t117 - t8;
t6 = qJ(6) * t117 + t171;
t5 = t44 * pkin(5) - t43 * qJ(6) + t70 * qJD(6) + t49;
t2 = -pkin(5) * t113 - t4;
t1 = qJ(6) * t113 + qJD(6) * t117 + t3;
t9 = [(-t130 * t50 + t132 * t51 + (-(2 * Ifges(4,4)) + t144) * t113 - (Ifges(5,1) * t128 + (2 * Ifges(4,1)) + (t164 - 0.2e1 * t167) * t130) * t112 + 0.2e1 * (mrSges(5,1) * t130 + mrSges(5,2) * t132 + mrSges(4,3)) * t65) * t119 + (t171 * t3 + t4 * t8 + t49 * t66) * t183 + 0.2e1 * t171 * t25 + t113 * (-Ifges(7,4) * t70 + Ifges(7,6) * t69) + t113 * (-Ifges(6,5) * t70 - Ifges(6,6) * t69) + 0.2e1 * t49 * (mrSges(6,1) * t69 - mrSges(6,2) * t70) + 0.2e1 * (t112 * t188 - t113 * t95) * mrSges(4,3) - 0.2e1 * t188 * t68 + 0.2e1 * m(4) * (t64 * t95 - t176) + (t27 * t47 + t28 * t48 - t176) * t184 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t131 ^ 2 + t133 ^ 2) * qJD(2) + (t30 - t31) * t44 + 0.2e1 * t151 * t147 + (t32 + t33) * t43 + (-0.2e1 * t64 * mrSges(4,3) - 0.2e1 * (-Ifges(4,4) + t144) * t112 + ((2 * Ifges(4,2)) + (2 * Ifges(5,3)) + t192) * t113 + t141) * t117 + (t1 * t6 + t19 * t5 + t2 * t7) * t182 + 0.2e1 * t19 * t14 + 0.2e1 * t8 * t23 + 0.2e1 * t7 * t24 + 0.2e1 * t6 * t26 + 0.2e1 * t5 * t45 + 0.2e1 * t3 * t53 + 0.2e1 * t4 * t54 + 0.2e1 * t2 * t55 + 0.2e1 * t1 * t56 + 0.2e1 * t66 * t15 + t69 * t10 - t69 * t11 - t70 * t12 - t70 * t13 + 0.2e1 * t48 * t73 + 0.2e1 * t47 * t74 + 0.2e1 * t28 * t83 + 0.2e1 * t27 * t84; t130 * t73 + t132 * t74 + t172 * t118 + t173 * t116 + t169 * t111 - t170 * t110 + m(6) * (-t110 * t171 - t111 * t8 - t116 * t4 + t118 * t3) + m(7) * (t1 * t118 - t110 * t6 + t111 * t7 + t116 * t2) + m(5) * (t130 * t28 + t132 * t27) + t147; 0.2e1 * t189 * (-t118 * t110 + t111 * t116); (t110 * t8 - t111 * t171) * mrSges(6,3) + (-t110 * t7 - t111 * t6) * mrSges(7,2) + (t13 + t12) * t118 / 0.2e1 + m(6) * (t125 * t49 + t171 * t62 + t3 * t94 - t4 * t92 - t63 * t8) + t186 * t117 / 0.2e1 - (t79 / 0.2e1 + t80 / 0.2e1) * t70 + (t77 / 0.2e1 - t78 / 0.2e1) * t69 + t169 * t63 + t170 * t62 + t172 * t94 + t173 * t92 - (t132 * (Ifges(5,1) * t130 + t167) / 0.2e1 - t130 * (Ifges(5,2) * t132 + t168) / 0.2e1 + Ifges(4,5)) * t112 + (t28 * mrSges(5,3) + qJ(4) * t73 + qJD(4) * t83 + t50 / 0.2e1 - t65 * mrSges(5,1) + t165 / 0.2e1) * t132 + (-t27 * mrSges(5,3) - qJ(4) * t74 - qJD(4) * t84 + t51 / 0.2e1 + t65 * mrSges(5,2) + t166 / 0.2e1) * t130 + (t87 / 0.2e1 - t88 / 0.2e1) * t44 + (t89 / 0.2e1 + t90 / 0.2e1) * t43 + (-t3 * mrSges(6,3) - t1 * mrSges(7,2) + t10 / 0.2e1 - t11 / 0.2e1 + t49 * mrSges(6,1) + t191 * t179) * t116 + (t49 * mrSges(6,2) + t2 * mrSges(7,2) - t4 * mrSges(6,3) + t193 * t179) * t118 + m(5) * (-pkin(3) * t65 + (-t130 * t47 + t132 * t48) * qJD(4) + (-t130 * t27 + t132 * t28) * qJ(4)) + m(7) * (t1 * t94 + t19 * t59 + t2 * t92 + t5 * t81 + t6 * t62 + t63 * t7) + (t30 / 0.2e1 - t31 / 0.2e1) * t111 - (t32 / 0.2e1 + t33 / 0.2e1) * t110 + t59 * t45 - t64 * mrSges(4,2) - t65 * mrSges(4,1) - pkin(3) * t68 + t19 * t75 + t66 * t76 + t81 * t14 + t5 * t86 - Ifges(4,6) * t113 + t125 * t15; t189 * (-t94 * t110 + t111 * t92 + t116 * t63 + t62 * t118); 0.2e1 * t125 * t76 + 0.2e1 * t59 * t86 + 0.2e1 * t81 * t75 + (t79 + t80) * t118 + (t77 - t78) * t116 + (t87 - t88) * t111 - (t89 + t90) * t110 + (t59 * t81 + t150) * t182 + t150 * t183 + (qJ(4) * t184 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t130 ^ 2 + t128) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t110 * t92 - t111 * t94 - t116 * t62 + t118 * t63); m(5) * t65 + m(6) * t49 + m(7) * t5 + t14 + t15 + t68; 0; -t187; 0; -pkin(5) * t24 + qJD(6) * t56 + qJ(6) * t26 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) + t1 * mrSges(7,3) - t3 * mrSges(6,2) - t2 * mrSges(7,1) + t4 * mrSges(6,1) + t141; t187; m(7) * qJD(6) * t94 + (pkin(5) * t110 - qJ(6) * t111 - qJD(6) * t116) * mrSges(7,2) + (-m(7) * pkin(5) - mrSges(6,1) - mrSges(7,1)) * t63 + (-mrSges(6,2) + t185) * t62 + t186; 0; 0.2e1 * t185 * qJD(6); m(7) * t2 + t24; m(7) * t111; m(7) * t63 - t110 * mrSges(7,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
