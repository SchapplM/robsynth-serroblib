% Calculate time derivative of joint inertia matrix for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:30
% EndTime: 2019-12-31 20:51:32
% DurationCPUTime: 0.73s
% Computational Cost: add. (460->145), mult. (1106->182), div. (0->0), fcn. (568->4), ass. (0->66)
t106 = Ifges(6,4) + Ifges(5,5);
t57 = sin(qJ(3));
t59 = cos(qJ(3));
t105 = t57 ^ 2 + t59 ^ 2;
t95 = -0.2e1 * t57;
t94 = -0.2e1 * t59;
t103 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1) - Ifges(4,2) - Ifges(6,2) - Ifges(5,3);
t60 = cos(qJ(2));
t83 = pkin(1) * qJD(2);
t73 = t60 * t83;
t102 = t105 * t73;
t80 = qJD(3) * t59;
t101 = t59 * pkin(3) + t57 * qJ(4);
t100 = (-mrSges(3,2) + (mrSges(5,2) + mrSges(4,3)) * t105) * t73;
t99 = 2 * m(6);
t81 = qJD(3) * t57;
t26 = -mrSges(6,1) * t81 + mrSges(6,2) * t80;
t98 = 0.2e1 * t26;
t32 = -t59 * mrSges(5,1) - t57 * mrSges(5,3);
t97 = 0.2e1 * t32;
t33 = mrSges(6,1) * t59 + mrSges(6,2) * t57;
t96 = 0.2e1 * t33;
t87 = pkin(7) - qJ(5);
t58 = sin(qJ(2));
t42 = pkin(1) * t58 + pkin(7);
t86 = t102 * t42;
t85 = t102 * pkin(7);
t82 = -qJ(5) + t42;
t79 = qJD(4) * t59;
t78 = qJ(4) * qJD(3);
t77 = pkin(2) + t101;
t8 = pkin(3) * t81 - t57 * qJD(4) - t59 * t78;
t76 = 0.2e1 * qJD(3);
t75 = m(5) * t79;
t74 = t58 * t83;
t43 = -pkin(1) * t60 - pkin(2);
t35 = t87 * t59;
t17 = t82 * t59;
t69 = -qJD(5) + t73;
t34 = -t59 * mrSges(4,1) + t57 * mrSges(4,2);
t68 = t57 * mrSges(4,1) + t59 * mrSges(4,2);
t67 = t57 * mrSges(5,1) - t59 * mrSges(5,3);
t14 = t43 - t101;
t66 = (-mrSges(3,1) + t34) * t74;
t4 = -pkin(4) * t81 - t8;
t65 = mrSges(5,2) * t79 + Ifges(5,6) * t81 + (Ifges(5,4) + Ifges(4,5)) * t80 + (t57 * t78 - t79) * mrSges(6,3);
t61 = -pkin(3) - pkin(4);
t64 = (-qJ(4) * mrSges(5,2) - Ifges(4,6) - Ifges(6,6)) * t57 + (-pkin(3) * mrSges(5,2) - t61 * mrSges(6,3) - Ifges(6,5)) * t59;
t63 = (0.2e1 * Ifges(4,4) * t59 + t103 * t57 + t106 * t94) * t80 + (Ifges(4,4) * t95 + t103 * t59 + 0.2e1 * t106 * t57) * t81;
t62 = -m(5) * t101 + t32 + t34;
t53 = t59 * pkin(4);
t46 = mrSges(5,2) * t80;
t40 = qJ(5) * t81;
t31 = t87 * t57;
t27 = t68 * qJD(3);
t25 = t67 * qJD(3);
t16 = t82 * t57;
t15 = t53 + t77;
t9 = qJD(3) * t35 - qJD(5) * t57;
t7 = -pkin(7) * t81 - qJD(5) * t59 + t40;
t6 = -t14 + t53;
t5 = t74 + t8;
t3 = t4 - t74;
t2 = qJD(3) * t17 + t69 * t57;
t1 = -t42 * t81 + t69 * t59 + t40;
t10 = [t63 + 0.2e1 * m(4) * (t43 * t74 + t86) + (t1 * t17 + t16 * t2 + t3 * t6) * t99 + 0.2e1 * m(5) * (t14 * t5 + t86) + t5 * t97 + t3 * t96 + 0.2e1 * t43 * t27 + 0.2e1 * t14 * t25 + t6 * t98 + 0.2e1 * t66 + (t1 * t94 + t2 * t95 + (-t16 * t59 + t17 * t57) * t76) * mrSges(6,3) + 0.2e1 * t100; t63 + m(4) * (-pkin(2) * t74 + t85) + m(6) * (t1 * t35 + t15 * t3 + t16 * t9 + t17 * t7 + t2 * t31 + t4 * t6) + m(5) * (t14 * t8 - t5 * t77 + t85) + t100 + t66 + (t4 + t3) * t33 + (t8 + t5) * t32 + (t43 - pkin(2)) * t27 + (t15 + t6) * t26 + (-t77 + t14) * t25 + ((-t1 - t7) * t59 + (-t2 - t9) * t57 + ((-t16 - t31) * t59 + (t17 + t35) * t57) * qJD(3)) * mrSges(6,3); t63 + t8 * t97 + t4 * t96 + t15 * t98 - 0.2e1 * pkin(2) * t27 + (t9 * t95 + t7 * t94 + (-t31 * t59 + t35 * t57) * t76) * mrSges(6,3) - 0.2e1 * (m(5) * t8 + t25) * t77 + (t15 * t4 + t31 * t9 + t35 * t7) * t99; -t2 * mrSges(6,1) + t1 * mrSges(6,2) + t42 * t75 + m(6) * (qJ(4) * t1 + qJD(4) * t17 + t2 * t61) + (m(5) * (-pkin(3) * t57 + qJ(4) * t59) - t67 - t68) * t73 + (t62 * t42 + t64) * qJD(3) + t65; -t9 * mrSges(6,1) + t7 * mrSges(6,2) + pkin(7) * t75 + m(6) * (qJ(4) * t7 + qJD(4) * t35 + t61 * t9) + (t62 * pkin(7) + t64) * qJD(3) + t65; 0.2e1 * (mrSges(6,2) + mrSges(5,3) + (m(5) + m(6)) * qJ(4)) * qJD(4); t46 - mrSges(6,3) * t80 + m(6) * t2 + (t42 * t80 + t57 * t73) * m(5); t46 + m(6) * t9 + (m(5) * pkin(7) - mrSges(6,3)) * t80; 0; 0; m(6) * t3 + t26; m(6) * t4 + t26; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
