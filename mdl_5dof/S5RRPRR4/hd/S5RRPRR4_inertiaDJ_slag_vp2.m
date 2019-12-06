% Calculate time derivative of joint inertia matrix for
% S5RRPRR4
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:50
% EndTime: 2019-12-05 18:31:52
% DurationCPUTime: 0.62s
% Computational Cost: add. (1196->147), mult. (2636->219), div. (0->0), fcn. (2057->8), ass. (0->86)
t105 = qJD(4) + qJD(5);
t61 = sin(qJ(5));
t62 = sin(qJ(4));
t64 = cos(qJ(5));
t65 = cos(qJ(4));
t43 = -t61 * t62 + t64 * t65;
t22 = t105 * t43;
t44 = t61 * t65 + t62 * t64;
t109 = t22 * t44;
t23 = t105 * t44;
t108 = t23 * t43;
t59 = sin(pkin(9));
t52 = pkin(2) * t59 + pkin(7);
t97 = -pkin(8) - t52;
t39 = t97 * t62;
t56 = t65 * pkin(8);
t40 = t52 * t65 + t56;
t16 = t39 * t64 - t40 * t61;
t72 = qJD(4) * t97;
t35 = t62 * t72;
t36 = t65 * t72;
t5 = qJD(5) * t16 + t35 * t64 + t36 * t61;
t17 = t39 * t61 + t40 * t64;
t6 = -qJD(5) * t17 - t35 * t61 + t36 * t64;
t107 = t6 * mrSges(6,1) - t5 * mrSges(6,2);
t66 = cos(qJ(2));
t54 = pkin(1) * t66 + pkin(2);
t60 = cos(pkin(9));
t63 = sin(qJ(2));
t91 = t60 * t63;
t89 = pkin(1) * t91 + t59 * t54;
t32 = pkin(7) + t89;
t98 = -pkin(8) - t32;
t25 = t98 * t62;
t26 = t32 * t65 + t56;
t10 = t25 * t64 - t26 * t61;
t87 = pkin(1) * qJD(2);
t92 = t59 * t63;
t34 = (t60 * t66 - t92) * t87;
t73 = qJD(4) * t98;
t14 = t34 * t65 + t62 * t73;
t15 = -t34 * t62 + t65 * t73;
t2 = qJD(5) * t10 + t14 * t64 + t15 * t61;
t11 = t25 * t61 + t26 * t64;
t3 = -qJD(5) * t11 - t14 * t61 + t15 * t64;
t106 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t71 = mrSges(5,1) * t62 + mrSges(5,2) * t65;
t45 = t71 * qJD(4);
t104 = 2 * m(6);
t9 = t23 * mrSges(6,1) + t22 * mrSges(6,2);
t103 = 0.2e1 * t9;
t102 = 0.2e1 * t45;
t101 = pkin(4) * t65;
t96 = Ifges(5,4) * t62;
t94 = Ifges(5,6) * t62;
t24 = -mrSges(6,1) * t43 + mrSges(6,2) * t44;
t33 = (t59 * t66 + t91) * t87;
t85 = qJD(4) * t62;
t78 = pkin(4) * t85;
t27 = t33 + t78;
t93 = t27 * t24;
t90 = Ifges(6,5) * t22 - Ifges(6,6) * t23;
t88 = t62 ^ 2 + t65 ^ 2;
t86 = pkin(4) * qJD(5);
t84 = qJD(4) * t65;
t83 = qJD(5) * t61;
t82 = qJD(5) * t64;
t81 = 2 * mrSges(6,3);
t80 = t64 * t22 * mrSges(6,3);
t79 = mrSges(6,3) * t86;
t53 = -pkin(2) * t60 - pkin(3);
t47 = -mrSges(5,1) * t65 + mrSges(5,2) * t62;
t77 = (t47 - mrSges(4,1)) * t33;
t76 = t88 * mrSges(5,3);
t75 = t88 * t34;
t74 = -pkin(1) * t92 + t54 * t60;
t31 = -pkin(3) - t74;
t70 = t64 * t43 * t79 + Ifges(5,5) * t84 + t90 + (-mrSges(6,3) * pkin(4) * t23 + t44 * t79) * t61;
t69 = 0.2e1 * Ifges(6,1) * t109 - 0.2e1 * Ifges(6,2) * t108 + (Ifges(5,1) * t65 - t96) * t85 + (0.2e1 * Ifges(5,4) * t65 + (Ifges(5,1) - Ifges(5,2)) * t62) * t84 + 0.2e1 * (t22 * t43 - t23 * t44) * Ifges(6,4);
t68 = (-mrSges(3,1) * t63 - mrSges(3,2) * t66) * t87;
t48 = Ifges(5,2) * t65 + t96;
t67 = -t48 * t85 + t69;
t46 = t53 - t101;
t42 = (-mrSges(6,1) * t61 - mrSges(6,2) * t64) * t86;
t28 = t31 - t101;
t1 = [0.2e1 * m(4) * (-t74 * t33 + t89 * t34) + (t10 * t3 + t11 * t2 + t27 * t28) * t104 + t31 * t102 + 0.2e1 * t93 + t28 * t103 + t67 - 0.2e1 * t34 * mrSges(4,2) + (-t10 * t22 - t11 * t23 + t2 * t43 - t3 * t44) * t81 + 0.2e1 * m(5) * (t31 * t33 + t32 * t75) + 0.2e1 * t76 * t34 + 0.2e1 * t68 + 0.2e1 * t77; t77 + m(4) * (-t33 * t60 + t34 * t59) * pkin(2) + (pkin(4) * t24 - t48) * t85 + m(6) * (t10 * t6 + t11 * t5 + t16 * t3 + t17 * t2 + t27 * t46 + t28 * t78) + (t46 + t28) * t9 + (t53 + t31) * t45 + t93 + t69 + t68 + (-mrSges(4,2) + t76) * t34 + ((-t3 - t6) * t44 + (t2 + t5) * t43 - (t11 + t17) * t23 + (-t10 - t16) * t22) * mrSges(6,3) + m(5) * (t33 * t53 + t52 * t75); 0.2e1 * t24 * t78 + t46 * t103 + t53 * t102 + (t16 * t6 + t17 * t5 + t46 * t78) * t104 + (-t16 * t22 - t17 * t23 + t43 * t5 - t44 * t6) * t81 + t67; m(6) * (-t10 * t23 + t11 * t22 + t2 * t44 + t3 * t43); m(6) * (-t16 * t23 + t17 * t22 + t43 * t6 + t44 * t5); (-t108 + t109) * t104; -t71 * t34 + (t47 * t32 - t94) * qJD(4) + (-t80 + m(6) * (-t10 * t83 + t11 * t82 + t2 * t61 + t3 * t64)) * pkin(4) + t70 + t106; (t47 * t52 - t94) * qJD(4) + (-t80 + m(6) * (-t16 * t83 + t17 * t82 + t5 * t61 + t6 * t64)) * pkin(4) + t70 + t107; -t45 + m(6) * (t22 * t61 - t23 * t64 + (-t43 * t61 + t44 * t64) * qJD(5)) * pkin(4) - t9; 0.2e1 * t42; t90 + t106; t90 + t107; -t9; t42; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
