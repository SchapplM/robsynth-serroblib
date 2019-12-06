% Calculate time derivative of joint inertia matrix for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:11:18
% EndTime: 2019-12-05 18:11:21
% DurationCPUTime: 0.98s
% Computational Cost: add. (3105->152), mult. (6719->255), div. (0->0), fcn. (6819->8), ass. (0->71)
t68 = sin(pkin(9));
t86 = pkin(6) + qJ(2);
t62 = t86 * t68;
t75 = cos(qJ(3));
t58 = t75 * t62;
t69 = cos(pkin(9));
t63 = t86 * t69;
t72 = sin(qJ(3));
t47 = -t63 * t72 - t58;
t61 = t68 * t75 + t72 * t69;
t40 = -pkin(7) * t61 + t47;
t48 = -t72 * t62 + t75 * t63;
t89 = t69 * t75;
t60 = -t72 * t68 + t89;
t41 = pkin(7) * t60 + t48;
t71 = sin(qJ(4));
t74 = cos(qJ(4));
t19 = t74 * t40 - t41 * t71;
t20 = t71 * t40 + t74 * t41;
t93 = 2 * m(6);
t54 = t61 * qJD(3);
t92 = pkin(3) * t54;
t65 = pkin(3) * t74 + pkin(4);
t73 = cos(qJ(5));
t84 = qJD(5) * t73;
t70 = sin(qJ(5));
t85 = qJD(5) * t70;
t88 = t70 * t71;
t43 = t65 * t84 + (-t71 * t85 + (t73 * t74 - t88) * qJD(4)) * pkin(3);
t90 = t43 * mrSges(6,2);
t87 = t71 * t73;
t83 = -pkin(2) * t69 - pkin(1);
t45 = t60 * t74 - t61 * t71;
t53 = t60 * qJD(3);
t24 = t45 * qJD(4) + t74 * t53 - t71 * t54;
t46 = t60 * t71 + t61 * t74;
t25 = -t46 * qJD(4) - t71 * t53 - t74 * t54;
t26 = t45 * t73 - t46 * t70;
t11 = t26 * qJD(5) + t73 * t24 + t70 * t25;
t27 = t45 * t70 + t46 * t73;
t12 = -t27 * qJD(5) - t70 * t24 + t73 * t25;
t82 = -t12 * mrSges(6,1) + t11 * mrSges(6,2);
t81 = -t25 * mrSges(5,1) + t24 * mrSges(5,2);
t44 = -t65 * t85 + (-t71 * t84 + (-t70 * t74 - t87) * qJD(4)) * pkin(3);
t42 = t44 * mrSges(6,1);
t80 = t42 - t90;
t35 = -qJD(3) * t58 + qJD(2) * t89 + (-qJD(2) * t68 - qJD(3) * t63) * t72;
t31 = -t54 * pkin(7) + t35;
t36 = -t61 * qJD(2) - t48 * qJD(3);
t32 = -t53 * pkin(7) + t36;
t14 = t19 * qJD(4) + t74 * t31 + t71 * t32;
t4 = t25 * pkin(8) + t14;
t15 = -t20 * qJD(4) - t71 * t31 + t74 * t32;
t5 = -t24 * pkin(8) + t15;
t16 = -t46 * pkin(8) + t19;
t17 = pkin(8) * t45 + t20;
t6 = t16 * t73 - t17 * t70;
t2 = t6 * qJD(5) + t73 * t4 + t70 * t5;
t7 = t16 * t70 + t17 * t73;
t3 = -t7 * qJD(5) - t70 * t4 + t73 * t5;
t78 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t11 + Ifges(6,6) * t12;
t49 = -pkin(3) * t60 + t83;
t77 = (-mrSges(5,1) * t71 - mrSges(5,2) * t74) * qJD(4) * pkin(3);
t76 = t15 * mrSges(5,1) - t14 * mrSges(5,2) + Ifges(5,5) * t24 + Ifges(5,6) * t25 + t78;
t56 = (-mrSges(6,1) * t70 - mrSges(6,2) * t73) * qJD(5) * pkin(4);
t52 = pkin(3) * t87 + t65 * t70;
t51 = -pkin(3) * t88 + t65 * t73;
t50 = t53 * mrSges(4,2);
t33 = -t45 * pkin(4) + t49;
t18 = -pkin(4) * t25 + t92;
t1 = [0.2e1 * (-mrSges(5,1) * t45 + mrSges(5,2) * t46) * t92 + 0.2e1 * t26 * Ifges(6,2) * t12 + 0.2e1 * t27 * t11 * Ifges(6,1) + 0.2e1 * t18 * (-t26 * mrSges(6,1) + t27 * mrSges(6,2)) + 0.2e1 * t33 * t82 + 0.2e1 * t45 * Ifges(5,2) * t25 + 0.2e1 * t46 * t24 * Ifges(5,1) + 0.2e1 * t49 * t81 - 0.2e1 * t60 * Ifges(4,2) * t54 + 0.2e1 * t61 * t53 * Ifges(4,1) + 0.2e1 * t83 * (t54 * mrSges(4,1) + t50) + (t18 * t33 + t2 * t7 + t3 * t6) * t93 + 0.2e1 * m(4) * (t35 * t48 + t36 * t47) + 0.2e1 * m(5) * (t14 * t20 + t15 * t19 + t49 * t92) + 0.2e1 * (t11 * t26 + t12 * t27) * Ifges(6,4) + 0.2e1 * (t24 * t45 + t25 * t46) * Ifges(5,4) + 0.2e1 * (t53 * t60 - t54 * t61) * Ifges(4,4) + 0.2e1 * (-t11 * t6 + t12 * t7 + t2 * t26 - t27 * t3) * mrSges(6,3) + 0.2e1 * (t14 * t45 - t15 * t46 - t19 * t24 + t20 * t25) * mrSges(5,3) + 0.2e1 * (t35 * t60 - t36 * t61 - t47 * t53 - t48 * t54) * mrSges(4,3) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t68 ^ 2 + t69 ^ 2) * qJD(2); m(6) * t18 + t50 - (-m(5) * pkin(3) - mrSges(4,1)) * t54 + t81 + t82; 0; m(6) * (t2 * t52 + t3 * t51 + t43 * t7 + t44 * t6) + Ifges(4,5) * t53 - Ifges(4,6) * t54 - t35 * mrSges(4,2) + t36 * mrSges(4,1) + (-t11 * t51 + t12 * t52 + t26 * t43 - t27 * t44) * mrSges(6,3) + (m(5) * (t14 * t71 + t15 * t74 + (-t19 * t71 + t20 * t74) * qJD(4)) + (-t74 * t24 + t71 * t25 + (t45 * t74 + t46 * t71) * qJD(4)) * mrSges(5,3)) * pkin(3) + t76; 0; -0.2e1 * t90 + 0.2e1 * t42 + (t43 * t52 + t44 * t51) * t93 + 0.2e1 * t77; (m(6) * (t2 * t70 + t3 * t73 + (-t6 * t70 + t7 * t73) * qJD(5)) + (-t73 * t11 + t70 * t12 + (t26 * t73 + t27 * t70) * qJD(5)) * mrSges(6,3)) * pkin(4) + t76; 0; t77 + (m(6) * (t43 * t70 + t44 * t73 - t51 * t85 + t52 * t84) - mrSges(6,2) * t84 - mrSges(6,1) * t85) * pkin(4) + t80; 0.2e1 * t56; t78; 0; t80; t56; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
