% Calculate time derivative of joint inertia matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR16_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:45
% EndTime: 2019-12-31 18:38:47
% DurationCPUTime: 0.81s
% Computational Cost: add. (570->166), mult. (1181->248), div. (0->0), fcn. (747->4), ass. (0->84)
t39 = sin(qJ(5));
t40 = sin(qJ(3));
t41 = cos(qJ(5));
t65 = qJD(5) * t41;
t42 = cos(qJ(3));
t68 = qJD(3) * t42;
t46 = t39 * t68 + t40 * t65;
t99 = 2 * qJ(2);
t98 = 2 * qJ(4);
t97 = m(5) + m(6);
t96 = -t40 * pkin(3) + qJ(4) * t42;
t95 = t40 * mrSges(4,1) + t42 * mrSges(4,2);
t94 = 2 * m(6);
t93 = -2 * mrSges(5,2);
t92 = 2 * mrSges(5,3);
t91 = -t39 / 0.2e1;
t90 = t39 / 0.2e1;
t89 = t41 / 0.2e1;
t44 = -pkin(1) - pkin(6);
t88 = pkin(4) - t44;
t87 = Ifges(6,4) * t39;
t86 = Ifges(6,4) * t41;
t85 = Ifges(6,6) * t42;
t50 = Ifges(6,1) * t39 + t86;
t12 = Ifges(6,5) * t42 + t50 * t40;
t84 = t39 * t12;
t28 = Ifges(6,1) * t41 - t87;
t83 = t39 * t28;
t43 = -pkin(3) - pkin(7);
t82 = t39 * t43;
t80 = t40 * mrSges(6,3);
t49 = Ifges(6,2) * t41 + t87;
t11 = t49 * t40 + t85;
t79 = t41 * t11;
t27 = -Ifges(6,2) * t39 + t86;
t78 = t41 * t27;
t77 = t41 * t43;
t75 = t42 * mrSges(5,3);
t73 = t39 ^ 2 + t41 ^ 2;
t71 = qJD(3) * t39;
t70 = qJD(3) * t40;
t69 = qJD(3) * t41;
t67 = qJD(5) * t39;
t66 = qJD(5) * t40;
t64 = qJD(5) * t42;
t63 = qJ(4) * qJD(3);
t62 = 0.2e1 * t42;
t58 = t41 * t68;
t61 = t46 * Ifges(6,5) + Ifges(6,6) * t58;
t57 = pkin(3) * t68 + t40 * t63 + qJD(2);
t56 = m(6) * t73;
t55 = m(5) * t44 - mrSges(5,1);
t54 = qJD(3) * t88;
t16 = t40 * t54;
t23 = qJ(2) - t96;
t14 = pkin(7) * t40 + t23;
t25 = t88 * t42;
t6 = -t14 * t39 + t25 * t41;
t8 = (pkin(7) * qJD(3) - qJD(4)) * t42 + t57;
t1 = qJD(5) * t6 - t39 * t16 + t41 * t8;
t7 = t14 * t41 + t25 * t39;
t2 = -qJD(5) * t7 - t41 * t16 - t39 * t8;
t53 = t39 * t1 + t41 * t2;
t52 = t39 * t6 - t41 * t7;
t51 = mrSges(6,1) * t41 - mrSges(6,2) * t39;
t26 = mrSges(6,1) * t39 + mrSges(6,2) * t41;
t48 = -Ifges(6,5) * t39 - Ifges(6,6) * t41;
t47 = -t39 * t66 + t58;
t10 = -mrSges(6,1) * t70 - t46 * mrSges(6,3);
t21 = mrSges(6,1) * t42 - t39 * t80;
t22 = -mrSges(6,2) * t42 + t41 * t80;
t9 = mrSges(6,2) * t70 + t47 * mrSges(6,3);
t45 = t41 * t10 - t21 * t67 + t22 * t65 + t39 * t9;
t24 = t88 * t40;
t20 = t50 * qJD(5);
t19 = t49 * qJD(5);
t18 = t51 * qJD(5);
t17 = t42 * t54;
t15 = t51 * t40;
t13 = -qJD(4) * t42 + t57;
t5 = -t47 * mrSges(6,1) + t46 * mrSges(6,2);
t4 = t28 * t66 + (-Ifges(6,5) * t40 + t50 * t42) * qJD(3);
t3 = t27 * t66 + (-Ifges(6,6) * t40 + t49 * t42) * qJD(3);
t29 = [t42 * t61 + 0.2e1 * t17 * t15 + 0.2e1 * t2 * t21 + 0.2e1 * t1 * t22 - 0.2e1 * t24 * t5 + 0.2e1 * t7 * t9 + 0.2e1 * t6 * t10 + (t1 * t7 + t17 * t24 + t2 * t6) * t94 + (t13 * t93 + t41 * t3 + t39 * t4 + (t12 * t41 + (-t11 - t85) * t39) * qJD(5)) * t40 + ((mrSges(4,1) * t99 + t23 * t93 + t79 + t84 + (-Ifges(4,4) - Ifges(5,6)) * t62) * t42 + (-0.2e1 * qJ(2) * mrSges(4,2) + t23 * t92 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + t48) * t40 + (-Ifges(4,1) + Ifges(4,2) - Ifges(5,2) + Ifges(5,3) - Ifges(6,3)) * t62) * t40) * qJD(3) + 0.2e1 * (m(5) * t23 - t75) * t13 + ((2 * mrSges(3,3)) + 0.2e1 * t95 + ((m(3) + m(4)) * t99)) * qJD(2); (t5 + m(6) * (t6 * t69 + t7 * t71 - t17) + t22 * t71 + t21 * t69) * t40 + (-qJD(3) * t15 + m(6) * (-qJD(3) * t24 + t6 * t67 - t7 * t65 - t53) - t45) * t42; (0.1e1 - t73) * t40 * t68 * t94; t3 * t91 + t4 * t89 - qJD(4) * t15 - t24 * t18 - t17 * t26 + qJ(4) * t5 + m(6) * (-qJ(4) * t17 - qJD(4) * t24 + t1 * t82 + t2 * t77) + t10 * t77 + t9 * t82 - t53 * mrSges(6,3) + (t55 * qJD(4) - t19 * t89 - t20 * t90) * t40 + (t42 * t48 / 0.2e1 - t79 / 0.2e1 - t84 / 0.2e1 + (t27 * t91 + t28 * t89) * t40 + t52 * mrSges(6,3) + (-m(6) * t52 - t39 * t21 + t41 * t22) * t43) * qJD(5) + ((-Ifges(6,5) * t41 / 0.2e1 + Ifges(6,6) * t90 + Ifges(5,4) - Ifges(4,5) + pkin(3) * mrSges(5,1)) * t40 + (Ifges(5,5) - Ifges(4,6) + t78 / 0.2e1 + t83 / 0.2e1 - qJ(4) * mrSges(5,1)) * t42 + (m(5) * t96 + t40 * mrSges(5,2) + t75 - t95) * t44) * qJD(3); t40 * t18 + ((-mrSges(4,2) + mrSges(5,3) + t26) * t42 + (-m(5) * pkin(3) - t73 * mrSges(6,3) + t43 * t56 - mrSges(4,1) + mrSges(5,2)) * t40) * qJD(3) + t97 * (qJD(4) * t40 + t42 * t63); t18 * t98 + t19 * t39 - t20 * t41 + (-t78 - t83) * qJD(5) + (t97 * t98 + 0.2e1 * t26 + t92) * qJD(4); m(6) * (-t52 * qJD(5) + t53) + t55 * t70 + t45; (m(5) + t56) * t70; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + (-Ifges(6,6) * t67 - Ifges(6,3) * qJD(3)) * t40 + t61; (-t39 * t70 + t41 * t64) * mrSges(6,2) + (t39 * t64 + t40 * t69) * mrSges(6,1); ((-mrSges(6,2) * t43 - Ifges(6,6)) * t41 + (-mrSges(6,1) * t43 - Ifges(6,5)) * t39) * qJD(5); -t26 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t29(1), t29(2), t29(4), t29(7), t29(11); t29(2), t29(3), t29(5), t29(8), t29(12); t29(4), t29(5), t29(6), t29(9), t29(13); t29(7), t29(8), t29(9), t29(10), t29(14); t29(11), t29(12), t29(13), t29(14), t29(15);];
Mq = res;
