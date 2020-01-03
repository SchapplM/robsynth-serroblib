% Calculate time derivative of joint inertia matrix for
% S5RRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5]';
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
% Datum: 2019-12-31 19:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:40:42
% EndTime: 2019-12-31 19:40:44
% DurationCPUTime: 0.76s
% Computational Cost: add. (627->169), mult. (1312->245), div. (0->0), fcn. (828->4), ass. (0->77)
t47 = sin(qJ(5));
t49 = cos(qJ(5));
t50 = cos(qJ(2));
t69 = qJD(5) * t50;
t48 = sin(qJ(2));
t74 = qJD(2) * t48;
t54 = t47 * t69 + t49 * t74;
t86 = -t49 / 0.2e1;
t77 = pkin(6) - qJ(4);
t73 = qJD(2) * t50;
t17 = -qJD(4) * t48 + t77 * t73;
t25 = -t50 * pkin(2) - t48 * qJ(3) - pkin(1);
t19 = t50 * pkin(3) - t25;
t11 = pkin(4) * t48 + pkin(7) * t50 + t19;
t27 = t77 * t48;
t5 = t11 * t49 - t27 * t47;
t51 = -pkin(2) - pkin(3);
t45 = -pkin(7) + t51;
t76 = qJ(3) * t73 + t48 * qJD(3);
t8 = (pkin(4) * t50 + t45 * t48) * qJD(2) + t76;
t1 = qJD(5) * t5 + t17 * t49 + t47 * t8;
t6 = t11 * t47 + t27 * t49;
t2 = -qJD(5) * t6 - t17 * t47 + t49 * t8;
t95 = -t1 * t49 + t2 * t47 + (t47 * t6 + t49 * t5) * qJD(5);
t94 = -0.2e1 * pkin(1);
t93 = 2 * mrSges(5,1);
t92 = -2 * mrSges(5,3);
t16 = pkin(2) * t74 - t76;
t91 = -0.2e1 * t16;
t90 = 0.2e1 * t25;
t28 = t77 * t50;
t89 = 0.2e1 * t28;
t46 = qJ(3) + pkin(4);
t88 = 0.2e1 * t46;
t87 = t47 / 0.2e1;
t85 = mrSges(6,3) * t50;
t84 = Ifges(6,4) * t47;
t83 = Ifges(6,4) * t49;
t82 = Ifges(6,5) * t49;
t15 = qJD(4) * t50 + t77 * t74;
t81 = t15 * t28;
t56 = Ifges(6,2) * t49 + t84;
t80 = t47 * t56;
t58 = Ifges(6,1) * t47 + t83;
t79 = t49 * t58;
t78 = -mrSges(4,2) + mrSges(5,3);
t75 = mrSges(5,1) * t73 + mrSges(5,2) * t74;
t72 = qJD(3) * t28;
t71 = qJD(5) * t47;
t70 = qJD(5) * t49;
t68 = t47 * t74;
t66 = t49 * t69;
t64 = t54 * Ifges(6,5) + Ifges(6,6) * t66 + Ifges(6,3) * t73;
t63 = m(4) * pkin(6) - t78;
t26 = mrSges(6,1) * t49 - mrSges(6,2) * t47;
t60 = -mrSges(6,1) * t47 - mrSges(6,2) * t49;
t59 = Ifges(6,1) * t49 - t84;
t57 = -Ifges(6,2) * t47 + t83;
t55 = -Ifges(6,6) * t47 - (2 * Ifges(3,4)) - (2 * Ifges(5,4)) + (2 * Ifges(4,5));
t53 = -t66 + t68;
t10 = -mrSges(6,2) * t73 - t53 * mrSges(6,3);
t9 = mrSges(6,1) * t73 - t54 * mrSges(6,3);
t52 = -m(6) * t95 + t49 * t10 - t47 * t9;
t39 = Ifges(6,6) * t71;
t24 = mrSges(6,1) * t48 + t49 * t85;
t23 = -mrSges(6,2) * t48 + t47 * t85;
t22 = t59 * qJD(5);
t21 = t57 * qJD(5);
t20 = t60 * qJD(5);
t18 = t60 * t50;
t14 = Ifges(6,5) * t48 - t59 * t50;
t13 = Ifges(6,6) * t48 - t57 * t50;
t12 = t51 * t74 + t76;
t7 = t53 * mrSges(6,1) + t54 * mrSges(6,2);
t4 = t58 * t69 + (t50 * Ifges(6,5) + t59 * t48) * qJD(2);
t3 = t56 * t69 + (t50 * Ifges(6,6) + t57 * t48) * qJD(2);
t29 = [m(4) * t16 * t90 + 0.2e1 * t19 * t75 - 0.2e1 * t15 * t18 + 0.2e1 * t1 * t23 + 0.2e1 * t2 * t24 + t7 * t89 + 0.2e1 * t5 * t9 + 0.2e1 * t6 * t10 + 0.2e1 * m(6) * (t1 * t6 + t2 * t5 - t81) + 0.2e1 * m(5) * (t12 * t19 + t17 * t27 - t81) + (mrSges(4,3) * t91 + t12 * t93 + t17 * t92 + t64) * t48 + (mrSges(4,1) * t91 - 0.2e1 * t12 * mrSges(5,2) + 0.2e1 * t15 * mrSges(5,3) + t47 * t3 - t49 * t4 + (t13 * t49 + t14 * t47) * qJD(5)) * t50 + ((mrSges(3,2) * t94 - 0.2e1 * t25 * mrSges(4,3) + t27 * t92 + (-t55 - t82) * t50) * t50 + (mrSges(3,1) * t94 + mrSges(4,1) * t90 + mrSges(5,3) * t89 - t47 * t13 + t49 * t14 + t55 * t48 + ((2 * Ifges(3,1)) + (2 * Ifges(4,1)) - (2 * Ifges(5,1)) - (2 * Ifges(3,2)) + (2 * Ifges(5,2)) - (2 * Ifges(4,3)) + Ifges(6,3)) * t50) * t48) * qJD(2); -t14 * t70 / 0.2e1 + t48 * (-Ifges(6,5) * t70 + t39) / 0.2e1 + t3 * t86 + t46 * t7 + t17 * mrSges(5,2) + qJD(3) * t18 + t28 * t20 + (qJD(5) * t13 / 0.2e1 - t4 / 0.2e1) * t47 + (-t26 - mrSges(5,1)) * t15 + m(5) * (-qJ(3) * t15 + t17 * t51 + t72) + m(6) * (-t15 * t46 + t72) + t95 * mrSges(6,3) + (-t23 * t71 - t24 * t70 + t52) * t45 + (-t22 * t86 - t21 * t87 + (t56 * t86 - t58 * t87) * qJD(5) + t63 * qJD(3)) * t50 + ((-Ifges(6,5) * t47 / 0.2e1 + Ifges(6,6) * t86 + Ifges(5,6) + Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,2) - t51 * mrSges(5,3) + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * pkin(6)) * t50 + (-Ifges(5,5) + Ifges(4,6) - Ifges(3,6) - t79 / 0.2e1 + t80 / 0.2e1 + t78 * qJ(3) + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(6)) * t48) * qJD(2); t20 * t88 + t21 * t49 + t22 * t47 + (t79 - t80) * qJD(5) + (m(6) * t88 + t93 + 0.2e1 * mrSges(4,3) + 0.2e1 * t26 + 0.2e1 * (m(4) + m(5)) * qJ(3)) * qJD(3); (-t47 * t23 - t49 * t24) * qJD(5) + m(5) * t17 + t63 * t73 + t52; 0; 0; t47 * t10 + t49 * t9 + (t49 * t23 - t47 * t24) * qJD(5) + m(6) * (t1 * t47 + t2 * t49 + (-t47 * t5 + t49 * t6) * qJD(5)) + m(5) * t12 + t75; 0; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,6) * t68 + t64; t39 + (-t26 * t45 - t82) * qJD(5); -t26 * qJD(5); t20; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t29(1), t29(2), t29(4), t29(7), t29(11); t29(2), t29(3), t29(5), t29(8), t29(12); t29(4), t29(5), t29(6), t29(9), t29(13); t29(7), t29(8), t29(9), t29(10), t29(14); t29(11), t29(12), t29(13), t29(14), t29(15);];
Mq = res;
