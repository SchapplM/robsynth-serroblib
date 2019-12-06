% Calculate time derivative of joint inertia matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:47
% EndTime: 2019-12-05 16:47:49
% DurationCPUTime: 0.74s
% Computational Cost: add. (810->137), mult. (2142->208), div. (0->0), fcn. (1760->6), ass. (0->65)
t89 = mrSges(5,1) + mrSges(6,1);
t88 = 2 * mrSges(5,3);
t87 = 2 * mrSges(6,3);
t54 = cos(qJ(4));
t86 = pkin(3) * t54;
t53 = sin(qJ(2));
t51 = sin(qJ(4));
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t61 = t51 * t52 - t54 * t55;
t32 = t61 * t53;
t81 = -pkin(7) - pkin(6);
t45 = t81 * t52;
t46 = t81 * t55;
t28 = t51 * t45 - t54 * t46;
t76 = t52 ^ 2 + t55 ^ 2;
t84 = pkin(3) * qJD(4);
t83 = qJD(3) + qJD(4);
t82 = m(6) * pkin(4);
t80 = m(6) * t53;
t38 = t51 * t55 + t52 * t54;
t24 = t83 * t38;
t56 = cos(qJ(2));
t75 = qJD(2) * t56;
t10 = -t24 * t53 - t61 * t75;
t71 = qJD(4) * t54;
t78 = (t10 * t51 - t32 * t71) * pkin(3);
t77 = -mrSges(5,2) - mrSges(6,2);
t74 = qJD(3) * t52;
t73 = qJD(3) * t55;
t72 = qJD(4) * t51;
t70 = 0.2e1 * t52;
t69 = pkin(3) * t74;
t68 = t53 * t75;
t31 = t38 * t53;
t67 = t31 * t72;
t48 = -pkin(3) * t55 - pkin(2);
t66 = qJD(3) * t81;
t65 = t77 * t54;
t23 = t83 * t61;
t6 = mrSges(6,1) * t24 - mrSges(6,2) * t23;
t27 = t45 * t54 + t46 * t51;
t64 = 2 * Ifges(5,4) + 2 * Ifges(6,4);
t43 = t52 * t66;
t44 = t55 * t66;
t14 = -qJD(4) * t28 - t43 * t51 + t54 * t44;
t4 = qJ(5) * t23 - qJD(5) * t38 + t14;
t63 = m(6) * t4 + t23 * mrSges(6,3);
t62 = -mrSges(4,1) * t55 + mrSges(4,2) * t52;
t11 = t32 * t83 - t38 * t75;
t60 = t77 * t10 + t11 * t89;
t13 = t43 * t54 + t44 * t51 + t45 * t71 + t46 * t72;
t59 = -t51 * t24 + (t38 * t51 - t54 * t61) * qJD(4);
t3 = -qJ(5) * t24 - qJD(5) * t61 + t13;
t57 = t14 * mrSges(5,1) + t4 * mrSges(6,1) - t13 * mrSges(5,2) - t3 * mrSges(6,2) - (Ifges(5,6) + Ifges(6,6)) * t24 - (Ifges(5,5) + Ifges(6,5)) * t23;
t47 = pkin(4) + t86;
t42 = (mrSges(4,1) * t52 + mrSges(4,2) * t55) * qJD(3);
t30 = pkin(4) * t61 + t48;
t26 = mrSges(5,1) * t61 + mrSges(5,2) * t38;
t25 = mrSges(6,1) * t61 + mrSges(6,2) * t38;
t17 = pkin(4) * t24 + t69;
t16 = -qJ(5) * t61 + t28;
t15 = -qJ(5) * t38 + t27;
t7 = mrSges(5,1) * t24 - mrSges(5,2) * t23;
t1 = [0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t10 * t32 - t11 * t31 - t68) + 0.2e1 * m(4) * (-0.1e1 + t76) * t68; (-t42 - t6 - t7) * t56 + m(5) * (t10 * t28 + t11 * t27 - t13 * t32 - t14 * t31 - t56 * t69) + m(6) * (t10 * t16 + t11 * t15 - t17 * t56 - t3 * t32 - t31 * t4) + (t30 * t80 + (-m(4) * pkin(2) + m(5) * t48 - mrSges(3,1) + t25 + t26 + t62) * t53 + (-mrSges(3,2) + (m(4) * pkin(6) + mrSges(4,3)) * t76) * t56) * qJD(2) + (mrSges(6,3) + mrSges(5,3)) * (-t10 * t61 - t11 * t38 - t23 * t31 + t24 * t32); -0.2e1 * pkin(2) * t42 + 0.2e1 * t17 * t25 + 0.2e1 * t30 * t6 + 0.2e1 * t48 * t7 + (-Ifges(4,4) * t52 + pkin(3) * t26) * qJD(3) * t70 + 0.2e1 * m(5) * (t13 * t28 + t14 * t27 + t48 * t69) + 0.2e1 * m(6) * (t15 * t4 + t16 * t3 + t17 * t30) + (0.2e1 * Ifges(4,4) * t55 + (Ifges(4,1) - Ifges(4,2)) * t70) * t73 + (-t3 * t61 - t38 * t4) * t87 + (-t13 * t61 - t14 * t38) * t88 - (t16 * t87 + t28 * t88 + t38 * t64 - 0.2e1 * (Ifges(5,2) + Ifges(6,2)) * t61) * t24 - (-0.2e1 * mrSges(5,3) * t27 - 0.2e1 * mrSges(6,3) * t15 - t61 * t64 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t38) * t23; (t53 * t74 - t55 * t75) * mrSges(4,2) + (-t52 * t75 - t53 * t73) * mrSges(4,1) + m(5) * ((t11 * t54 + t67) * pkin(3) + t78) + m(6) * (pkin(3) * t67 + t11 * t47 + t78) + t60; t63 * t47 + (Ifges(4,5) * t55 - Ifges(4,6) * t52 + pkin(6) * t62) * qJD(3) + (m(6) * (-t15 * t72 + t16 * t71 + t3 * t51) + m(5) * (t13 * t51 + t14 * t54 - t27 * t72 + t28 * t71) + t59 * mrSges(6,3) + (t54 * t23 + t59) * mrSges(5,3)) * pkin(3) + t57; 0.2e1 * (t65 + ((-t47 + t86) * m(6) - t89) * t51) * t84; t11 * t82 + t60; pkin(4) * t63 + t57; (t65 + (-t82 - t89) * t51) * t84; 0; qJD(2) * t80; m(6) * t17 + t6; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
