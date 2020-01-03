% Calculate time derivative of joint inertia matrix for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:42
% EndTime: 2019-12-31 17:25:44
% DurationCPUTime: 0.94s
% Computational Cost: add. (1240->163), mult. (2924->262), div. (0->0), fcn. (2429->6), ass. (0->83)
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t70 = cos(qJ(3));
t117 = (t66 ^ 2 + t69 ^ 2) * t70;
t116 = qJD(2) + qJD(3);
t67 = sin(qJ(3));
t68 = sin(qJ(2));
t71 = cos(qJ(2));
t47 = t67 * t68 - t70 * t71;
t34 = t116 * t47;
t105 = t34 * t66;
t48 = t67 * t71 + t70 * t68;
t91 = qJD(4) * t69;
t76 = t48 * t91 - t105;
t35 = t116 * t48;
t11 = pkin(2) * qJD(2) * t68 + t35 * pkin(3) + t34 * pkin(7);
t62 = -pkin(2) * t71 - pkin(1);
t28 = t47 * pkin(3) - t48 * pkin(7) + t62;
t111 = -pkin(6) - pkin(5);
t57 = t111 * t68;
t58 = t111 * t71;
t42 = t57 * t67 - t58 * t70;
t13 = t28 * t69 - t42 * t66;
t86 = qJD(2) * t111;
t53 = t68 * t86;
t77 = t70 * t57 + t58 * t67;
t83 = t71 * t86;
t16 = qJD(3) * t77 + t70 * t53 + t67 * t83;
t2 = qJD(4) * t13 + t11 * t66 + t16 * t69;
t14 = t28 * t66 + t42 * t69;
t3 = -qJD(4) * t14 + t11 * t69 - t16 * t66;
t118 = t2 * t69 - t3 * t66;
t115 = 2 * m(5);
t17 = qJD(3) * t42 + t53 * t67 - t70 * t83;
t114 = 0.2e1 * t17;
t113 = 0.2e1 * t62;
t108 = Ifges(5,4) * t66;
t107 = Ifges(5,4) * t69;
t106 = Ifges(5,6) * t66;
t104 = t34 * t69;
t103 = t77 * t17;
t101 = t48 * t66;
t100 = t48 * t69;
t97 = t67 * mrSges(4,1);
t54 = -mrSges(5,1) * t69 + mrSges(5,2) * t66;
t96 = t67 * t54;
t95 = t70 * mrSges(4,2);
t94 = -Ifges(5,5) * t104 + Ifges(5,3) * t35;
t93 = pkin(2) * qJD(3);
t92 = qJD(4) * t66;
t90 = 0.2e1 * t71;
t89 = t48 * t92;
t75 = t89 + t104;
t8 = mrSges(5,1) * t76 - mrSges(5,2) * t75;
t87 = m(5) * t17 + t8;
t85 = -t92 / 0.2e1;
t84 = -(2 * Ifges(4,4)) - t106;
t82 = mrSges(5,3) * t117;
t81 = mrSges(5,1) * t66 + mrSges(5,2) * t69;
t80 = Ifges(5,1) * t69 - t108;
t79 = -Ifges(5,2) * t66 + t107;
t78 = Ifges(5,5) * t66 + Ifges(5,6) * t69;
t51 = t79 * qJD(4);
t52 = t80 * qJD(4);
t55 = Ifges(5,2) * t69 + t108;
t56 = Ifges(5,1) * t66 + t107;
t74 = t69 * t51 + t66 * t52 - t55 * t92 + t56 * t91;
t10 = -t35 * mrSges(5,2) - mrSges(5,3) * t76;
t29 = -t47 * mrSges(5,2) - mrSges(5,3) * t101;
t30 = t47 * mrSges(5,1) - mrSges(5,3) * t100;
t9 = t35 * mrSges(5,1) + mrSges(5,3) * t75;
t73 = -t30 * t91 - t29 * t92 + m(5) * (-t13 * t91 - t14 * t92 + t118) + t69 * t10 - t66 * t9;
t20 = Ifges(5,6) * t47 + t48 * t79;
t21 = Ifges(5,5) * t47 + t48 * t80;
t50 = t81 * qJD(4);
t6 = -Ifges(5,4) * t75 - Ifges(5,2) * t76 + Ifges(5,6) * t35;
t63 = Ifges(5,5) * t91;
t7 = -Ifges(5,1) * t75 - Ifges(5,4) * t76 + Ifges(5,5) * t35;
t72 = t20 * t85 + t21 * t91 / 0.2e1 - t77 * t50 - Ifges(4,5) * t34 - t51 * t101 / 0.2e1 + t52 * t100 / 0.2e1 + t66 * t7 / 0.2e1 + t47 * (-Ifges(5,6) * t92 + t63) / 0.2e1 + t69 * t6 / 0.2e1 - t16 * mrSges(4,2) + (-t104 / 0.2e1 + t48 * t85) * t56 + (t78 / 0.2e1 - Ifges(4,6)) * t35 + (t54 - mrSges(4,1)) * t17 - t76 * t55 / 0.2e1 + ((-t13 * t69 - t14 * t66) * qJD(4) + t118) * mrSges(5,3);
t61 = -pkin(2) * t70 - pkin(3);
t60 = pkin(2) * t67 + pkin(7);
t26 = t81 * t48;
t1 = [(t35 * mrSges(4,1) - t34 * mrSges(4,2)) * t113 - 0.2e1 * t77 * t8 + t26 * t114 + 0.2e1 * t2 * t29 + 0.2e1 * t3 * t30 + 0.2e1 * t13 * t9 + 0.2e1 * t14 * t10 - t21 * t104 + t20 * t105 + 0.2e1 * (t34 * t77 - t42 * t35) * mrSges(4,3) + 0.2e1 * m(4) * (t42 * t16 - t103) + (t13 * t3 + t14 * t2 - t103) * t115 + (-0.2e1 * mrSges(4,3) * t16 + ((2 * Ifges(4,2)) + Ifges(5,3)) * t35 - t84 * t34 + t94) * t47 + (mrSges(4,3) * t114 - 0.2e1 * Ifges(4,1) * t34 - t66 * t6 + t69 * t7 + (Ifges(5,5) * t69 + t84) * t35 + (-t69 * t20 - t66 * t21 - t47 * t78) * qJD(4)) * t48 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t71) * t90 + (0.2e1 * pkin(2) * (mrSges(4,1) * t47 + mrSges(4,2) * t48) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t113 - 0.2e1 * Ifges(3,4) * t68 + (Ifges(3,1) - Ifges(3,2)) * t90) * t68) * qJD(2); t72 + t87 * t61 + t73 * t60 + (m(4) * (t16 * t67 - t17 * t70) + (t70 * t34 - t67 * t35) * mrSges(4,3) + ((-t47 * mrSges(4,3) + t69 * t29 - t66 * t30 + m(5) * (-t13 * t66 + t14 * t69) + m(4) * t42) * t70 + (t48 * mrSges(4,3) + t26 - (m(5) + m(4)) * t77) * t67) * qJD(3)) * pkin(2) + (Ifges(3,5) * t71 - Ifges(3,6) * t68 + (-mrSges(3,1) * t71 + mrSges(3,2) * t68) * pkin(5)) * qJD(2); 0.2e1 * t61 * t50 + (-0.2e1 * t97 - 0.2e1 * t95 + (t117 * t60 + t61 * t67) * t115 + 0.2e1 * t96 + 0.2e1 * t82) * t93 + t74; -pkin(3) * t87 + pkin(7) * t73 + t72; (-pkin(3) + t61) * t50 + (-t97 - t95 + m(5) * (-pkin(3) * t67 + t117 * pkin(7)) + t96 + t82) * t93 + t74; -0.2e1 * pkin(3) * t50 + t74; t3 * mrSges(5,1) - t2 * mrSges(5,2) - Ifges(5,5) * t89 - Ifges(5,6) * t76 + t94; t63 - t81 * t70 * t93 + (t54 * t60 - t106) * qJD(4); t63 + (pkin(7) * t54 - t106) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
