% Calculate time derivative of joint inertia matrix for
% S5PRRRP6
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:48
% EndTime: 2019-12-05 16:50:51
% DurationCPUTime: 0.77s
% Computational Cost: add. (786->141), mult. (2121->217), div. (0->0), fcn. (1689->6), ass. (0->64)
t81 = (mrSges(6,2) + mrSges(5,3));
t88 = 2 * t81;
t87 = -mrSges(5,1) - mrSges(6,1);
t86 = -mrSges(5,2) + mrSges(6,3);
t46 = sin(qJ(4));
t47 = sin(qJ(3));
t49 = cos(qJ(4));
t50 = cos(qJ(3));
t31 = t46 * t50 + t47 * t49;
t48 = sin(qJ(2));
t27 = t31 * t48;
t80 = t47 ^ 2 + t50 ^ 2;
t78 = qJD(3) * t47;
t51 = cos(qJ(2));
t79 = qJD(2) * t51;
t85 = -t48 * t78 + t50 * t79;
t84 = qJD(3) + qJD(4);
t83 = -pkin(7) - pkin(6);
t77 = qJD(3) * t50;
t76 = qJD(4) * t46;
t75 = qJD(4) * t49;
t74 = 0.2e1 * t47;
t73 = pkin(3) * t78;
t72 = pkin(3) * t76;
t71 = pkin(3) * t75;
t70 = t83 * t47;
t69 = t47 * t79;
t68 = t48 * t79;
t37 = t83 * t50;
t25 = -t46 * t37 - t49 * t70;
t65 = t25 * t76;
t64 = t27 * t76;
t42 = -pkin(3) * t50 - pkin(2);
t62 = qJD(3) * t83;
t34 = t47 * t62;
t60 = t50 * t62;
t14 = -qJD(4) * t25 + t49 * t34 + t46 * t60;
t26 = -t49 * t37 + t46 * t70;
t15 = qJD(4) * t26 + t46 * t34 - t49 * t60;
t63 = t14 * t26 + t15 * t25;
t61 = -2 * Ifges(5,4) + 2 * Ifges(6,5);
t59 = -t50 * mrSges(4,1) + t47 * mrSges(4,2);
t30 = t46 * t47 - t49 * t50;
t10 = -t27 * t84 - t30 * t79;
t11 = -t47 * t48 * t76 + (t48 * t50 * t84 + t69) * t49 + t85 * t46;
t58 = t10 * t86 + t11 * t87;
t28 = t30 * t48;
t57 = t10 * t26 + t11 * t25 - t14 * t28 + t15 * t27;
t21 = t84 * t30;
t22 = t84 * t31;
t53 = (-Ifges(5,6) + Ifges(6,6)) * t22 - (Ifges(6,4) + Ifges(5,5)) * t21 + t87 * t15 + t86 * t14;
t38 = qJD(5) + t71;
t52 = -mrSges(5,2) * t71 + t38 * mrSges(6,3) + t72 * t87;
t45 = qJD(5) * mrSges(6,3);
t41 = -pkin(3) * t49 - pkin(4);
t39 = pkin(3) * t46 + qJ(5);
t33 = (mrSges(4,1) * t47 + mrSges(4,2) * t50) * qJD(3);
t24 = mrSges(5,1) * t30 + mrSges(5,2) * t31;
t23 = mrSges(6,1) * t30 - mrSges(6,3) * t31;
t16 = pkin(4) * t30 - qJ(5) * t31 + t42;
t7 = mrSges(5,1) * t22 - mrSges(5,2) * t21;
t6 = mrSges(6,1) * t22 + mrSges(6,3) * t21;
t3 = pkin(4) * t22 + qJ(5) * t21 - qJD(5) * t31 + t73;
t1 = [0.2e1 * m(4) * (-0.1e1 + t80) * t68 + 0.2e1 * (m(5) + m(6)) * (-t10 * t28 + t11 * t27 - t68); (-t33 - t6 - t7) * t51 + m(5) * (-t51 * t73 + t57) + m(6) * (-t3 * t51 + t57) + ((-m(4) * pkin(2) + m(5) * t42 + m(6) * t16 - mrSges(3,1) + t23 + t24 + t59) * t48 + (-mrSges(3,2) + (m(4) * pkin(6) + mrSges(4,3)) * t80) * t51) * qJD(2) + t81 * (-t10 * t30 + t11 * t31 - t21 * t27 + t22 * t28); -0.2e1 * pkin(2) * t33 + 0.2e1 * t16 * t6 + 0.2e1 * t3 * t23 + 0.2e1 * t42 * t7 + (-Ifges(4,4) * t47 + pkin(3) * t24) * qJD(3) * t74 + 0.2e1 * m(5) * (t42 * t73 + t63) + 0.2e1 * m(6) * (t16 * t3 + t63) + (t31 * t61 + 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t30 - 0.2e1 * t81 * t26) * t22 - (t30 * t61 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t31 + t25 * t88) * t21 + (0.2e1 * Ifges(4,4) * t50 + (Ifges(4,1) - Ifges(4,2)) * t74) * t77 + (-t14 * t30 + t15 * t31) * t88; -t85 * mrSges(4,2) + (-t48 * t77 - t69) * mrSges(4,1) + m(5) * (t10 * t46 - t11 * t49 - t28 * t75 + t64) * pkin(3) + t58 + (pkin(3) * t64 + t10 * t39 + t11 * t41 - t28 * t38) * m(6); m(6) * (t14 * t39 + t15 * t41 + t26 * t38) + (-t21 * t41 - t22 * t39 - t30 * t38) * mrSges(6,2) + (Ifges(4,5) * t50 - Ifges(4,6) * t47 + pkin(6) * t59) * qJD(3) + (t31 * mrSges(6,2) * t76 + m(6) * t65 + m(5) * (t14 * t46 - t15 * t49 + t26 * t75 + t65) + (t49 * t21 - t46 * t22 + (-t30 * t49 + t31 * t46) * qJD(4)) * mrSges(5,3)) * pkin(3) + t53; 0.2e1 * m(6) * (t38 * t39 + t41 * t72) + 0.2e1 * t52; m(6) * (-pkin(4) * t11 + qJ(5) * t10 - qJD(5) * t28) + t58; m(6) * (-pkin(4) * t15 + qJ(5) * t14 + qJD(5) * t26) + (pkin(4) * t21 - qJ(5) * t22 - qJD(5) * t30) * mrSges(6,2) + t53; t45 + m(6) * (-pkin(4) * t72 + qJ(5) * t38 + qJD(5) * t39) + t52; 0.2e1 * m(6) * qJ(5) * qJD(5) + 0.2e1 * t45; m(6) * t11; m(6) * t15 - t21 * mrSges(6,2); m(6) * t72; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
