% Calculate time derivative of joint inertia matrix for
% S5RRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:50:56
% EndTime: 2019-12-31 19:50:58
% DurationCPUTime: 0.62s
% Computational Cost: add. (927->131), mult. (2102->174), div. (0->0), fcn. (1671->6), ass. (0->73)
t106 = mrSges(5,3) + mrSges(6,2);
t107 = 2 * mrSges(4,3);
t57 = cos(qJ(2));
t87 = pkin(1) * qJD(2);
t78 = t57 * t87;
t46 = qJD(3) + t78;
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t88 = t53 ^ 2 + t54 ^ 2;
t72 = t88 * t46;
t55 = sin(qJ(4));
t90 = t55 * t54;
t91 = cos(qJ(4));
t42 = t91 * t53 + t90;
t75 = t91 * t54;
t64 = -t55 * t53 + t75;
t105 = -mrSges(4,1) * t54 - mrSges(5,1) * t64 + mrSges(4,2) * t53 + mrSges(5,2) * t42;
t104 = m(6) * qJD(5);
t103 = m(6) * qJ(5) + mrSges(6,3);
t101 = -m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1);
t100 = -mrSges(5,2) + t103;
t99 = 2 * m(4);
t98 = 2 * m(5);
t97 = 0.2e1 * m(6);
t37 = t64 * qJD(4);
t38 = t42 * qJD(4);
t22 = t38 * mrSges(6,1) - t37 * mrSges(6,3);
t96 = 0.2e1 * t22;
t23 = t38 * mrSges(5,1) + t37 * mrSges(5,2);
t95 = 0.2e1 * t23;
t25 = -mrSges(6,1) * t64 - mrSges(6,3) * t42;
t94 = 0.2e1 * t25;
t93 = pkin(1) * t57;
t56 = sin(qJ(2));
t48 = pkin(1) * t56 + qJ(3);
t92 = -pkin(7) - t48;
t89 = -pkin(7) - qJ(3);
t86 = qJD(4) * t55;
t85 = t55 * qJD(3);
t79 = t56 * t87;
t49 = -t54 * pkin(3) - pkin(2);
t50 = t54 * pkin(7);
t39 = t48 * t54 + t50;
t66 = t92 * t91;
t20 = t55 * t39 - t53 * t66;
t76 = t55 * t92;
t21 = t91 * t39 + t53 * t76;
t5 = t46 * t75 - t39 * t86 + (qJD(4) * t66 - t55 * t46) * t53;
t71 = qJD(4) * t91;
t6 = t39 * t71 + t46 * t90 + (qJD(4) * t76 + t91 * t46) * t53;
t77 = t20 * t6 + t21 * t5;
t45 = qJ(3) * t54 + t50;
t65 = t89 * t91;
t69 = t91 * qJD(3);
t17 = t54 * t69 - t45 * t86 + (qJD(4) * t65 - t85) * t53;
t73 = t55 * t89;
t18 = t45 * t71 + t54 * t85 + (qJD(4) * t73 + t69) * t53;
t27 = t55 * t45 - t53 * t65;
t28 = t91 * t45 + t53 * t73;
t74 = t28 * t17 + t18 * t27;
t70 = t88 * qJ(3);
t67 = t88 * qJD(3);
t63 = t17 * t21 + t18 * t20 + t27 * t6 + t28 * t5;
t62 = t22 + t23;
t61 = 0.2e1 * (-Ifges(5,2) - Ifges(6,3)) * t64 * t38 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t42 * t37 + 0.2e1 * (-Ifges(6,5) + Ifges(5,4)) * (t37 * t64 - t42 * t38);
t8 = pkin(4) * t38 - qJ(5) * t37 - qJD(5) * t42;
t24 = -pkin(4) * t64 - t42 * qJ(5) + t49;
t59 = qJD(5) * t64 * mrSges(6,2) + (-qJ(5) * mrSges(6,2) - Ifges(5,6) + Ifges(6,6)) * t38 + (-pkin(4) * mrSges(6,2) + Ifges(6,4) + Ifges(5,5)) * t37;
t43 = t49 - t93;
t30 = t37 * mrSges(6,2);
t19 = t24 - t93;
t7 = t8 + t79;
t1 = [t77 * t98 + (t19 * t7 + t77) * t97 + t43 * t95 + t7 * t94 + t19 * t96 - 0.2e1 * mrSges(3,2) * t78 + t61 + (t48 * t99 + t107) * t72 + (t43 * t98 - (2 * mrSges(3,1)) + (-pkin(2) - t93) * t99 + 0.2e1 * t105) * t79 + 0.2e1 * t106 * (t20 * t37 - t21 * t38 + t42 * t6 + t5 * t64); m(5) * (t49 * t79 + t63) + m(6) * (t19 * t8 + t24 * t7 + t63) + (t8 + t7) * t25 + (t43 + t49) * t23 + (t19 + t24) * t22 + t61 + m(4) * (-pkin(2) * t79 + t46 * t70 + t48 * t67) + (t72 + t67) * mrSges(4,3) + (-mrSges(3,2) * t57 + (-mrSges(3,1) + t105) * t56) * t87 + t106 * ((t18 + t6) * t42 - (-t17 - t5) * t64 + (-t21 - t28) * t38 + (t20 + t27) * t37); (t24 * t8 + t74) * t97 + t74 * t98 + t49 * t95 + t8 * t94 + t24 * t96 + t61 + (t88 * t107 + t70 * t99) * qJD(3) + 0.2e1 * t106 * (t17 * t64 + t18 * t42 + t27 * t37 - t28 * t38); m(6) * t7 + (m(4) + m(5)) * t79 + t62; m(6) * t8 + t62; 0; t100 * t5 + t101 * t6 + t21 * t104 + t59; t100 * t17 + t101 * t18 + t28 * t104 + t59; 0; 0.2e1 * t103 * qJD(5); m(6) * t6 + t30; m(6) * t18 + t30; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
