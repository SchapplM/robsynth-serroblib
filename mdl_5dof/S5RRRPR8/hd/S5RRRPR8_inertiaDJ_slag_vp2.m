% Calculate time derivative of joint inertia matrix for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:09
% EndTime: 2019-12-31 21:19:14
% DurationCPUTime: 1.79s
% Computational Cost: add. (1993->221), mult. (4455->331), div. (0->0), fcn. (3786->6), ass. (0->108)
t172 = mrSges(5,2) - mrSges(4,1);
t83 = cos(qJ(5));
t124 = qJD(5) * t83;
t162 = qJD(2) + qJD(3);
t147 = cos(qJ(3));
t82 = sin(qJ(2));
t109 = t147 * t82;
t148 = cos(qJ(2));
t81 = sin(qJ(3));
t60 = t148 * t81 + t109;
t45 = t162 * t60;
t102 = t147 * t148;
t133 = t81 * t82;
t59 = -t102 + t133;
t80 = sin(qJ(5));
t94 = t59 * t124 + t45 * t80;
t144 = Ifges(6,4) * t83;
t66 = -Ifges(6,2) * t80 + t144;
t145 = Ifges(6,4) * t80;
t67 = Ifges(6,1) * t83 - t145;
t171 = t66 * t83 + t67 * t80;
t125 = qJD(5) * t80;
t165 = -t59 * t125 + t45 * t83;
t130 = t80 ^ 2 + t83 ^ 2;
t170 = (-mrSges(4,2) * t147 + (-t130 * mrSges(6,3) + t172) * t81) * pkin(2) * qJD(3);
t160 = t162 * t133;
t44 = -t162 * t102 + t160;
t169 = -0.2e1 * t44;
t168 = 0.2e1 * t59;
t151 = -pkin(7) - pkin(6);
t65 = mrSges(6,1) * t80 + mrSges(6,2) * t83;
t167 = mrSges(5,3) + t65;
t166 = Ifges(3,1) - Ifges(3,2);
t135 = t80 * mrSges(6,3);
t37 = mrSges(6,1) * t60 - t135 * t59;
t132 = t83 * mrSges(6,3);
t38 = -mrSges(6,2) * t60 + t132 * t59;
t152 = pkin(3) + pkin(8);
t76 = -pkin(2) * t148 - pkin(1);
t91 = -t60 * qJ(4) + t76;
t27 = t152 * t59 + t91;
t61 = t151 * t109;
t68 = t151 * t148;
t51 = -t68 * t81 - t61;
t31 = pkin(4) * t60 + t51;
t11 = -t27 * t80 + t31 * t83;
t12 = t27 * t83 + t31 * t80;
t105 = qJD(3) * t147;
t89 = qJD(2) * t68;
t21 = -t68 * t105 - t147 * t89 + t151 * t160;
t14 = -t44 * pkin(4) + t21;
t128 = qJD(2) * t82;
t77 = pkin(2) * t128;
t92 = qJ(4) * t44 - qJD(4) * t60 + t77;
t7 = t152 * t45 + t92;
t1 = qJD(5) * t11 + t14 * t80 + t7 * t83;
t149 = t1 * t80;
t15 = mrSges(6,2) * t44 + mrSges(6,3) * t165;
t16 = -mrSges(6,1) * t44 - mrSges(6,3) * t94;
t126 = qJD(5) * t12;
t2 = t14 * t83 - t7 * t80 - t126;
t86 = m(6) * (t149 + t2 * t83 + (-t11 * t80 + t12 * t83) * qJD(5)) + t83 * t16 + t80 * t15;
t161 = t38 * t124 - t37 * t125 + t86;
t156 = 2 * m(5);
t155 = 0.2e1 * m(6);
t101 = mrSges(6,1) * t83 - mrSges(6,2) * t80;
t62 = t101 * qJD(5);
t154 = 0.2e1 * t62;
t150 = pkin(2) * t81;
t143 = t44 * mrSges(5,1);
t140 = t59 * t80;
t139 = t59 * t83;
t103 = pkin(2) * t105;
t69 = t103 + qJD(4);
t73 = qJ(4) + t150;
t136 = t69 * t73;
t127 = qJD(3) * t81;
t118 = pkin(2) * t127;
t117 = t147 * pkin(2);
t112 = t139 / 0.2e1;
t106 = qJD(2) * t148;
t104 = t60 * t118;
t75 = -t117 - pkin(3);
t100 = Ifges(6,1) * t80 + t144;
t99 = Ifges(6,2) * t83 + t145;
t98 = Ifges(6,5) * t80 + Ifges(6,6) * t83;
t20 = -t127 * t68 - t162 * t61 - t81 * t89;
t52 = t133 * t151 - t147 * t68;
t97 = -t20 * t52 + t21 * t51;
t96 = t130 * t118;
t95 = qJ(4) * t69 + qJD(4) * t73;
t90 = t94 * Ifges(6,5) + t165 * Ifges(6,6) - Ifges(6,3) * t44;
t63 = t99 * qJD(5);
t64 = t100 * qJD(5);
t88 = -t171 * qJD(5) + t80 * t63 - t83 * t64;
t13 = -pkin(4) * t45 - t20;
t25 = Ifges(6,6) * t60 + t59 * t99;
t26 = Ifges(6,5) * t60 + t100 * t59;
t32 = -t59 * pkin(4) + t52;
t5 = Ifges(6,4) * t94 + Ifges(6,2) * t165 - Ifges(6,6) * t44;
t6 = Ifges(6,1) * t94 + Ifges(6,4) * t165 - Ifges(6,5) * t44;
t85 = t13 * t65 - t25 * t124 / 0.2e1 + t32 * t62 - t80 * t5 / 0.2e1 + t83 * t6 / 0.2e1 - t64 * t140 / 0.2e1 - t63 * t112 + t11 * mrSges(6,3) * t125 + (-Ifges(6,5) * t83 / 0.2e1 + Ifges(6,6) * t80 / 0.2e1 + Ifges(5,4) - Ifges(4,5)) * t44 + t172 * t21 - (t59 * t66 + t26) * t125 / 0.2e1 + (t67 * t112 - t60 * t98 / 0.2e1) * qJD(5) + (Ifges(5,5) - Ifges(4,6) + t171 / 0.2e1) * t45;
t72 = -pkin(8) + t75;
t34 = t59 * pkin(3) + t91;
t33 = t101 * t59;
t17 = pkin(3) * t45 + t92;
t8 = -mrSges(6,1) * t165 + mrSges(6,2) * t94;
t3 = [(Ifges(5,3) + Ifges(4,2)) * t45 * t168 + (Ifges(5,2) + Ifges(4,1)) * t60 * t169 + (mrSges(4,3) + mrSges(5,1)) * (t20 * t168 + t51 * t169 + 0.2e1 * t21 * t60 - 0.2e1 * t45 * t52) + 0.2e1 * (Ifges(5,6) + Ifges(4,4)) * (t44 * t59 - t45 * t60) + 0.2e1 * t76 * (mrSges(4,1) * t45 - mrSges(4,2) * t44) + 0.2e1 * t17 * (-mrSges(5,2) * t59 - mrSges(5,3) * t60) + 0.2e1 * t34 * (-mrSges(5,2) * t45 + mrSges(5,3) * t44) + 0.2e1 * t32 * t8 - 0.2e1 * t13 * t33 + 0.2e1 * t2 * t37 + 0.2e1 * t1 * t38 + 0.2e1 * t12 * t15 + 0.2e1 * t11 * t16 + (0.2e1 * Ifges(3,4) * t148 + t166 * t82) * t106 + (-0.2e1 * Ifges(3,4) * t82 + t166 * t148) * t128 + t165 * t25 + t94 * t26 + (t1 * t12 + t11 * t2 + t13 * t32) * t155 + (t17 * t34 + t97) * t156 - 0.2e1 * pkin(1) * (mrSges(3,1) * t82 + mrSges(3,2) * t148) * qJD(2) + t60 * t90 + t5 * t139 + t6 * t140 - t44 * (Ifges(6,3) * t60 + t59 * t98) + 0.2e1 * (mrSges(4,1) * t59 + mrSges(4,2) * t60) * t77 + 0.2e1 * m(4) * (t76 * t77 + t97); t85 - t69 * t33 + t73 * t8 - t20 * mrSges(5,3) + t20 * mrSges(4,2) - t12 * mrSges(6,3) * t124 + m(4) * (-t147 * t21 - t20 * t81 + (t147 * t52 + t51 * t81) * qJD(3)) * pkin(2) + Ifges(3,5) * t106 + m(6) * (t13 * t73 + t32 * t69) + m(5) * (-t20 * t73 + t21 * t75 + t52 * t69) - Ifges(3,6) * t128 - t2 * t132 - t1 * t135 - t75 * t143 + (t83 * t37 + t80 * t38 + m(6) * (t11 * t83 + t12 * t80) + m(5) * t51) * t118 + (-mrSges(3,1) * t106 + mrSges(3,2) * t128) * pkin(6) + (-t45 * t73 - t59 * t69 + t104) * mrSges(5,1) + t161 * t72 + (-t103 * t59 + t117 * t44 - t150 * t45 + t104) * mrSges(4,3); t73 * t154 + 0.2e1 * t167 * t69 + 0.2e1 * t170 + (t72 * t96 + t136) * t155 + (t118 * t75 + t136) * t156 + t88; t85 - qJD(4) * t33 + qJ(4) * t8 - t161 * t152 + (-t149 + (-t2 - t126) * t83) * mrSges(6,3) + (pkin(3) * t44 - qJ(4) * t45 - qJD(4) * t59) * mrSges(5,1) + m(5) * (-pkin(3) * t21 - qJ(4) * t20 + qJD(4) * t52) + m(6) * (qJ(4) * t13 + qJD(4) * t32) + (mrSges(4,2) - mrSges(5,3)) * t20; (t73 + qJ(4)) * t62 + t170 + m(6) * (-t152 * t96 + t95) + m(5) * (-pkin(3) * t118 + t95) + t88 + t167 * (t69 + qJD(4)); qJ(4) * t154 + 0.2e1 * ((m(5) + m(6)) * qJ(4) + t167) * qJD(4) + t88; -t143 + (-t80 * t37 + t83 * t38) * qJD(5) + m(5) * t21 + t86; (m(6) * t130 + m(5)) * t118; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t90; t101 * t118 + (-t65 * t72 - t98) * qJD(5); ((mrSges(6,2) * t152 - Ifges(6,6)) * t83 + (mrSges(6,1) * t152 - Ifges(6,5)) * t80) * qJD(5); -t65 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
