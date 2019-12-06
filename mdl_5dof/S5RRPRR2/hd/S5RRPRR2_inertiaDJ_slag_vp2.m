% Calculate time derivative of joint inertia matrix for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:27:06
% EndTime: 2019-12-05 18:27:10
% DurationCPUTime: 1.19s
% Computational Cost: add. (3322->168), mult. (7129->281), div. (0->0), fcn. (6983->8), ass. (0->76)
t82 = sin(qJ(2));
t98 = -qJ(3) - pkin(6);
t72 = t98 * t82;
t85 = cos(qJ(2));
t73 = t98 * t85;
t78 = sin(pkin(9));
t79 = cos(pkin(9));
t49 = t79 * t72 + t73 * t78;
t68 = t78 * t85 + t79 * t82;
t45 = -pkin(7) * t68 + t49;
t50 = t78 * t72 - t79 * t73;
t67 = -t78 * t82 + t79 * t85;
t46 = pkin(7) * t67 + t50;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t22 = t84 * t45 - t46 * t81;
t23 = t81 * t45 + t84 * t46;
t100 = pkin(2) * t78;
t75 = pkin(2) * t79 + pkin(3);
t60 = -t81 * t100 + t84 * t75;
t76 = -pkin(2) * t85 - pkin(1);
t101 = 0.2e1 * t76;
t89 = qJD(2) * t98;
t62 = qJD(3) * t85 + t82 * t89;
t63 = -t82 * qJD(3) + t85 * t89;
t44 = t79 * t62 + t78 * t63;
t80 = sin(qJ(5));
t97 = qJD(5) * t80;
t83 = cos(qJ(5));
t96 = qJD(5) * t83;
t95 = 0.2e1 * t85;
t77 = qJD(2) * t82 * pkin(2);
t47 = t67 * t84 - t68 * t81;
t64 = t68 * qJD(2);
t65 = t67 * qJD(2);
t27 = t47 * qJD(4) - t64 * t81 + t65 * t84;
t48 = t67 * t81 + t68 * t84;
t28 = -t48 * qJD(4) - t64 * t84 - t65 * t81;
t29 = t47 * t83 - t48 * t80;
t11 = t29 * qJD(5) + t27 * t83 + t28 * t80;
t30 = t47 * t80 + t48 * t83;
t12 = -t30 * qJD(5) - t27 * t80 + t28 * t83;
t93 = -t12 * mrSges(6,1) + t11 * mrSges(6,2);
t51 = pkin(3) * t64 + t77;
t92 = t64 * mrSges(4,1) + t65 * mrSges(4,2);
t91 = -t28 * mrSges(5,1) + t27 * mrSges(5,2);
t58 = pkin(4) + t60;
t61 = t84 * t100 + t75 * t81;
t41 = t58 * t83 - t61 * t80;
t56 = t60 * qJD(4);
t57 = t61 * qJD(4);
t20 = qJD(5) * t41 + t56 * t83 - t57 * t80;
t42 = t58 * t80 + t61 * t83;
t21 = -qJD(5) * t42 - t56 * t80 - t57 * t83;
t90 = t21 * mrSges(6,1) - t20 * mrSges(6,2);
t43 = -t62 * t78 + t79 * t63;
t34 = -pkin(7) * t65 + t43;
t35 = -pkin(7) * t64 + t44;
t14 = t22 * qJD(4) + t81 * t34 + t84 * t35;
t4 = t28 * pkin(8) + t14;
t15 = -t23 * qJD(4) + t84 * t34 - t35 * t81;
t5 = -t27 * pkin(8) + t15;
t16 = -t48 * pkin(8) + t22;
t17 = pkin(8) * t47 + t23;
t6 = t16 * t83 - t17 * t80;
t2 = t6 * qJD(5) + t4 * t83 + t5 * t80;
t7 = t16 * t80 + t17 * t83;
t3 = -t7 * qJD(5) - t4 * t80 + t5 * t83;
t88 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t11 + Ifges(6,6) * t12;
t52 = -t67 * pkin(3) + t76;
t87 = -t57 * mrSges(5,1) - t56 * mrSges(5,2) + t90;
t86 = t15 * mrSges(5,1) - t14 * mrSges(5,2) + Ifges(5,5) * t27 + Ifges(5,6) * t28 + t88;
t66 = (-t80 * mrSges(6,1) - t83 * mrSges(6,2)) * qJD(5) * pkin(4);
t36 = -t47 * pkin(4) + t52;
t18 = -t28 * pkin(4) + t51;
t1 = [0.2e1 * t29 * Ifges(6,2) * t12 + 0.2e1 * t30 * t11 * Ifges(6,1) + 0.2e1 * t18 * (-t29 * mrSges(6,1) + t30 * mrSges(6,2)) + 0.2e1 * t36 * t93 + 0.2e1 * t47 * Ifges(5,2) * t28 + 0.2e1 * t48 * t27 * Ifges(5,1) + 0.2e1 * t51 * (-t47 * mrSges(5,1) + t48 * mrSges(5,2)) + 0.2e1 * t52 * t91 - 0.2e1 * t67 * Ifges(4,2) * t64 + 0.2e1 * t68 * t65 * Ifges(4,1) + t92 * t101 + 0.2e1 * m(5) * (t14 * t23 + t15 * t22 + t51 * t52) + 0.2e1 * m(6) * (t18 * t36 + t2 * t7 + t3 * t6) + 0.2e1 * m(4) * (t49 * t43 + t50 * t44) + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t85) * t95 + (m(4) * pkin(2) * t101 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t67 + mrSges(4,2) * t68) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t82 + (Ifges(3,1) - Ifges(3,2)) * t95) * t82) * qJD(2) + 0.2e1 * (t29 * t11 + t30 * t12) * Ifges(6,4) + 0.2e1 * (t47 * t27 + t48 * t28) * Ifges(5,4) + 0.2e1 * (-t68 * t64 + t67 * t65) * Ifges(4,4) + 0.2e1 * (-t6 * t11 + t7 * t12 + t2 * t29 - t3 * t30) * mrSges(6,3) + 0.2e1 * (t14 * t47 - t15 * t48 - t22 * t27 + t23 * t28) * mrSges(5,3) + 0.2e1 * (-t43 * t68 + t44 * t67 - t49 * t65 - t50 * t64) * mrSges(4,3); t86 + (m(4) * (t43 * t79 + t44 * t78) + (-t64 * t78 - t65 * t79) * mrSges(4,3)) * pkin(2) + (Ifges(3,5) * t85 - Ifges(3,6) * t82 + (-mrSges(3,1) * t85 + mrSges(3,2) * t82) * pkin(6)) * qJD(2) + m(5) * (t14 * t61 + t15 * t60 - t57 * t22 + t56 * t23) + m(6) * (t2 * t42 + t20 * t7 + t21 * t6 + t3 * t41) + t43 * mrSges(4,1) - t44 * mrSges(4,2) - Ifges(4,6) * t64 + Ifges(4,5) * t65 + (-t27 * t60 + t28 * t61 + t47 * t56 + t48 * t57) * mrSges(5,3) + (-t11 * t41 + t12 * t42 + t20 * t29 - t21 * t30) * mrSges(6,3); 0.2e1 * m(6) * (t20 * t42 + t21 * t41) + 0.2e1 * m(5) * (t56 * t61 - t57 * t60) + 0.2e1 * t87; m(4) * t77 + m(5) * t51 + m(6) * t18 + t91 + t92 + t93; 0; 0; (m(6) * (t2 * t80 + t3 * t83 + (-t6 * t80 + t7 * t83) * qJD(5)) + (-t83 * t11 + t80 * t12 + (t29 * t83 + t30 * t80) * qJD(5)) * mrSges(6,3)) * pkin(4) + t86; (m(6) * (t20 * t80 + t21 * t83 - t41 * t97 + t42 * t96) - mrSges(6,2) * t96 - mrSges(6,1) * t97) * pkin(4) + t87; 0; 0.2e1 * t66; t88; t90; 0; t66; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
