% Calculate time derivative of joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:08
% EndTime: 2020-01-03 11:45:09
% DurationCPUTime: 0.54s
% Computational Cost: add. (431->113), mult. (979->152), div. (0->0), fcn. (602->6), ass. (0->61)
t85 = Ifges(6,4) + Ifges(5,4);
t84 = -Ifges(5,2) - Ifges(6,2);
t43 = cos(qJ(4));
t76 = 0.2e1 * t43;
t83 = Ifges(5,1) + Ifges(6,1);
t41 = sin(qJ(4));
t56 = qJD(4) * t43;
t57 = qJD(4) * t41;
t81 = t85 * t41;
t45 = (t83 * t43 - t81) * t57 + ((t83 + t84) * t41 + t85 * t76) * t56;
t62 = t84 * t43 - t81;
t82 = t62 * t57 + t45;
t30 = cos(pkin(8)) * pkin(1) + pkin(2);
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t73 = pkin(1) * sin(pkin(8));
t49 = t44 * t30 - t42 * t73;
t61 = t42 * t30 + t44 * t73;
t79 = 2 * m(6);
t19 = mrSges(6,1) * t57 + mrSges(6,2) * t56;
t78 = 0.2e1 * t19;
t77 = -0.2e1 * t41;
t75 = m(6) * pkin(4);
t74 = mrSges(6,3) * pkin(4);
t72 = t43 * pkin(4);
t22 = -t43 * mrSges(6,1) + t41 * mrSges(6,2);
t54 = pkin(4) * t57;
t8 = t61 * qJD(3);
t5 = t8 + t54;
t71 = t5 * t22;
t70 = mrSges(5,2) * t43;
t64 = -Ifges(5,6) - Ifges(6,6);
t63 = -qJ(5) - pkin(7);
t60 = (Ifges(5,5) + Ifges(6,5)) * t56;
t59 = t41 ^ 2 + t43 ^ 2;
t10 = pkin(7) + t61;
t58 = -qJ(5) - t10;
t55 = 0.2e1 * qJD(4);
t53 = mrSges(6,1) + t75;
t52 = (-t43 * mrSges(5,1) + t41 * mrSges(5,2) - mrSges(4,1)) * t8;
t7 = t49 * qJD(3);
t51 = t59 * t7;
t50 = t59 * mrSges(5,3);
t48 = qJD(4) * t63;
t47 = qJD(4) * t58;
t9 = -pkin(3) - t49;
t46 = mrSges(5,1) * t41 + t70;
t37 = t43 * qJ(5);
t36 = t43 * qJD(5);
t31 = -pkin(3) - t72;
t24 = t43 * pkin(7) + t37;
t21 = t63 * t41;
t20 = t46 * qJD(4);
t12 = -t41 * qJD(5) + t43 * t48;
t11 = t41 * t48 + t36;
t6 = t9 - t72;
t4 = t43 * t10 + t37;
t3 = t58 * t41;
t2 = (-qJD(5) - t7) * t41 + t43 * t47;
t1 = t41 * t47 + t43 * t7 + t36;
t13 = [t6 * t78 + 0.2e1 * t9 * t20 + 0.2e1 * t71 - 0.2e1 * t7 * mrSges(4,2) + 0.2e1 * m(4) * (-t49 * t8 + t61 * t7) + (t4 * t1 + t3 * t2 + t6 * t5) * t79 + 0.2e1 * m(5) * (t10 * t51 + t9 * t8) + (t1 * t76 + t2 * t77 + (-t3 * t43 - t4 * t41) * t55) * mrSges(6,3) + 0.2e1 * t50 * t7 + 0.2e1 * t52 + t82; m(6) * (t1 * t41 + t2 * t43 + (-t3 * t41 + t4 * t43) * qJD(4)); 0; t71 + t52 + (-pkin(3) + t9) * t20 + (t31 + t6) * t19 + (-mrSges(4,2) + t50) * t7 + (pkin(4) * t22 + t62) * t57 + m(6) * (t24 * t1 + t11 * t4 + t12 * t3 + t21 * t2 + t31 * t5 + t54 * t6) + m(5) * (-pkin(3) * t8 + pkin(7) * t51) + ((t1 + t11) * t43 + (-t12 - t2) * t41 + ((-t21 - t3) * t43 + (-t24 - t4) * t41) * qJD(4)) * mrSges(6,3) + t45; m(6) * (t11 * t41 + t12 * t43 + (-t21 * t41 + t24 * t43) * qJD(4)); -0.2e1 * pkin(3) * t20 + 0.2e1 * t22 * t54 + t31 * t78 + (t24 * t11 + t21 * t12 + t31 * t54) * t79 + (t11 * t76 + t12 * t77 + (-t21 * t43 - t24 * t41) * t55) * mrSges(6,3) + t82; -t1 * mrSges(6,2) - t46 * t7 + t53 * t2 + ((-mrSges(5,1) * t10 - t74) * t43 + (mrSges(5,2) * t10 + t64) * t41) * qJD(4) + t60; (-t70 + (-mrSges(5,1) - t75) * t41) * qJD(4) - t19; -t11 * mrSges(6,2) + t53 * t12 + ((-mrSges(5,1) * pkin(7) - t74) * t43 + (mrSges(5,2) * pkin(7) + t64) * t41) * qJD(4) + t60; 0; m(6) * t5 + t19; 0; m(6) * t54 + t19; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t13(1), t13(2), t13(4), t13(7), t13(11); t13(2), t13(3), t13(5), t13(8), t13(12); t13(4), t13(5), t13(6), t13(9), t13(13); t13(7), t13(8), t13(9), t13(10), t13(14); t13(11), t13(12), t13(13), t13(14), t13(15);];
Mq = res;
