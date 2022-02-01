% Calculate time derivative of joint inertia matrix for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% m [6x1]
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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:24:57
% EndTime: 2022-01-23 09:25:00
% DurationCPUTime: 1.20s
% Computational Cost: add. (1713->181), mult. (4125->304), div. (0->0), fcn. (3776->8), ass. (0->95)
t80 = sin(pkin(8));
t84 = sin(qJ(3));
t104 = t84 * t80;
t67 = pkin(3) * t104 + t80 * qJ(2);
t113 = -0.2e1 * t67;
t112 = 0.2e1 * t80;
t79 = sin(pkin(9));
t81 = cos(pkin(9));
t86 = cos(qJ(3));
t88 = t79 * t84 - t81 * t86;
t98 = qJD(3) * t80;
t49 = t88 * t98;
t64 = t79 * t86 + t81 * t84;
t57 = t64 * qJD(3);
t50 = t80 * t57;
t111 = Ifges(5,5) * t50 - Ifges(5,6) * t49;
t82 = cos(pkin(8));
t68 = -pkin(2) * t82 - pkin(6) * t80 - pkin(1);
t71 = t86 * t82 * qJ(2);
t48 = t84 * t68 + t71;
t53 = t64 * t80;
t54 = t88 * t80;
t83 = sin(qJ(5));
t85 = cos(qJ(5));
t29 = -t53 * t85 + t54 * t83;
t110 = 0.2e1 * t29;
t30 = -t53 * t83 - t54 * t85;
t109 = 0.2e1 * t30;
t108 = -0.2e1 * t53;
t107 = -0.2e1 * t54;
t106 = 0.2e1 * t82;
t105 = pkin(3) * t79;
t13 = qJD(5) * t29 + t49 * t83 - t50 * t85;
t14 = -qJD(5) * t30 + t49 * t85 + t50 * t83;
t103 = -Ifges(6,5) * t13 - Ifges(6,6) * t14;
t97 = qJD(3) * t86;
t99 = qJD(2) * t82;
t102 = t68 * t97 + t86 * t99;
t101 = qJ(2) * t84;
t93 = t82 * t101;
t100 = qJ(4) * t80;
t94 = t86 * t100;
t96 = qJD(4) * t80;
t27 = -t84 * t96 + (-t93 - t94) * qJD(3) + t102;
t92 = t84 * t99;
t28 = -t92 - t86 * t96 + (-t71 + (-t68 + t100) * t84) * qJD(3);
t9 = t81 * t27 + t79 * t28;
t62 = t86 * t68;
t34 = -t94 + t62 + (-pkin(3) - t101) * t82;
t41 = -t84 * t100 + t48;
t21 = t79 * t34 + t81 * t41;
t75 = t80 * qJD(2);
t60 = t80 * pkin(3) * t97 + t75;
t95 = qJ(2) * qJD(2);
t35 = -t64 * t83 - t85 * t88;
t58 = t88 * qJD(3);
t18 = qJD(5) * t35 - t57 * t83 - t58 * t85;
t36 = t64 * t85 - t83 * t88;
t19 = -qJD(5) * t36 - t57 * t85 + t58 * t83;
t91 = t19 * mrSges(6,1) - t18 * mrSges(6,2);
t8 = -t27 * t79 + t81 * t28;
t20 = t81 * t34 - t41 * t79;
t90 = mrSges(4,1) * t84 + mrSges(4,2) * t86;
t74 = pkin(3) * t81 + pkin(4);
t55 = -t83 * t105 + t74 * t85;
t51 = t55 * qJD(5);
t56 = t85 * t105 + t74 * t83;
t52 = t56 * qJD(5);
t89 = -t52 * mrSges(6,1) - t51 * mrSges(6,2);
t15 = -pkin(4) * t82 + t54 * pkin(7) + t20;
t16 = -pkin(7) * t53 + t21;
t4 = t15 * t85 - t16 * t83;
t5 = t15 * t83 + t16 * t85;
t6 = t50 * pkin(7) + t8;
t7 = pkin(7) * t49 + t9;
t2 = qJD(5) * t4 + t6 * t83 + t7 * t85;
t3 = -qJD(5) * t5 + t6 * t85 - t7 * t83;
t87 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t103;
t78 = t82 ^ 2;
t77 = t80 ^ 2;
t73 = t77 * t95;
t66 = -mrSges(4,3) * t80 * t86 - t82 * mrSges(4,1);
t65 = mrSges(4,2) * t82 - mrSges(4,3) * t104;
t47 = t62 - t93;
t44 = t50 * mrSges(5,2);
t43 = -mrSges(5,1) * t82 + t54 * mrSges(5,3);
t42 = mrSges(5,2) * t82 - t53 * mrSges(5,3);
t40 = -qJD(3) * t48 - t92;
t39 = -qJD(3) * t93 + t102;
t37 = t53 * pkin(4) + t67;
t31 = -t49 * pkin(4) + t60;
t23 = -mrSges(6,1) * t82 - t30 * mrSges(6,3);
t22 = mrSges(6,2) * t82 + t29 * mrSges(6,3);
t10 = t13 * mrSges(6,2);
t1 = [0.2e1 * t60 * (mrSges(5,1) * t53 - mrSges(5,2) * t54) + 0.2e1 * t39 * t65 + 0.2e1 * t40 * t66 + t44 * t113 + 0.2e1 * t8 * t43 + 0.2e1 * t31 * (-t29 * mrSges(6,1) + t30 * mrSges(6,2)) + 0.2e1 * t37 * t10 + 0.2e1 * t9 * t42 + 0.2e1 * t2 * t22 + 0.2e1 * t3 * t23 + (t103 + t111) * t82 + 0.2e1 * (t77 + t78) * qJD(2) * mrSges(3,3) + (mrSges(5,1) * t113 + 0.2e1 * t21 * mrSges(5,3) + Ifges(5,4) * t107 + Ifges(5,2) * t108 - Ifges(5,6) * t82) * t49 - (-0.2e1 * t20 * mrSges(5,3) + Ifges(5,1) * t107 + Ifges(5,4) * t108 - Ifges(5,5) * t82) * t50 + (-0.2e1 * t37 * mrSges(6,1) + 0.2e1 * t5 * mrSges(6,3) + Ifges(6,4) * t109 + Ifges(6,2) * t110 - Ifges(6,6) * t82) * t14 + (-0.2e1 * t4 * mrSges(6,3) + Ifges(6,1) * t109 + Ifges(6,4) * t110 - Ifges(6,5) * t82) * t13 + 0.2e1 * m(4) * (t48 * t39 + t47 * t40 + t73) + 0.2e1 * m(3) * (t78 * t95 + t73) + 0.2e1 * m(5) * (t20 * t8 + t21 * t9 + t60 * t67) + 0.2e1 * m(6) * (t2 * t5 + t3 * t4 + t31 * t37) + (0.2e1 * t90 * t75 + ((-0.2e1 * t48 * mrSges(4,3) + Ifges(4,6) * t106 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t86) * t112) * t86 + (0.2e1 * t47 * mrSges(4,3) + Ifges(4,5) * t106 + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t84 + (-Ifges(4,1) + Ifges(4,2)) * t86) * t112) * t84) * qJD(3)) * t80; t18 * t22 + t19 * t23 - t58 * t42 - t57 * t43 + (t86 * t65 - t84 * t66) * qJD(3) + (-t13 * t35 + t14 * t36) * mrSges(6,3) + (t49 * t64 - t50 * t88) * mrSges(5,3) + m(6) * (t18 * t5 + t19 * t4 + t2 * t36 + t3 * t35) + m(5) * (-t57 * t20 - t58 * t21 + t64 * t9 - t8 * t88) + m(4) * (t84 * t39 + t40 * t86 + (-t47 * t84 + t48 * t86) * qJD(3)); 0.2e1 * m(5) * (t57 * t88 - t58 * t64) + 0.2e1 * m(6) * (t18 * t36 + t19 * t35); m(6) * (t2 * t56 + t3 * t55 - t4 * t52 + t5 * t51) + t51 * t22 - t52 * t23 + t8 * mrSges(5,1) - t9 * mrSges(5,2) - t39 * mrSges(4,2) + t40 * mrSges(4,1) + (-Ifges(4,5) * t84 - Ifges(4,6) * t86) * t98 + (-t13 * t55 + t14 * t56) * mrSges(6,3) + (m(5) * (t79 * t9 + t8 * t81) + (t49 * t79 + t50 * t81) * mrSges(5,3)) * pkin(3) + t87 - t111; -t57 * mrSges(5,1) + t58 * mrSges(5,2) - t90 * qJD(3) + m(6) * (t18 * t56 + t19 * t55 - t35 * t52 + t36 * t51) + m(5) * (-t57 * t81 - t58 * t79) * pkin(3) + t91; 0.2e1 * m(6) * (t51 * t56 - t52 * t55) + 0.2e1 * t89; m(5) * t60 + m(6) * t31 - t49 * mrSges(5,1) - t14 * mrSges(6,1) + t10 - t44; 0; 0; 0; t87; t91; t89; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
