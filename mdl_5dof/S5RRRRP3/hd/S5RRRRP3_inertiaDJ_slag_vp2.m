% Calculate time derivative of joint inertia matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:08
% EndTime: 2019-12-31 21:49:09
% DurationCPUTime: 0.54s
% Computational Cost: add. (710->128), mult. (1827->186), div. (0->0), fcn. (1042->6), ass. (0->89)
t120 = 2 * Ifges(5,4) - 2 * Ifges(6,5);
t73 = cos(qJ(4));
t69 = t73 ^ 2;
t74 = cos(qJ(3));
t96 = qJD(3) * t74;
t90 = pkin(2) * t96;
t86 = t69 * t90;
t70 = sin(qJ(4));
t68 = t70 ^ 2;
t87 = t68 * t90;
t119 = t86 + t87;
t71 = sin(qJ(3));
t72 = sin(qJ(2));
t102 = t71 * t72;
t75 = cos(qJ(2));
t61 = t75 * pkin(1) + pkin(2);
t97 = qJD(3) * t71;
t16 = t61 * t96 + (-t72 * t97 + (t74 * t75 - t102) * qJD(2)) * pkin(1);
t103 = t16 * t69;
t104 = t16 * t68;
t118 = t103 + t104;
t117 = Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,3);
t91 = pkin(2) * t97;
t94 = qJD(4) * t73;
t85 = t70 * mrSges(5,1) + t73 * mrSges(5,2);
t44 = t85 * qJD(4);
t111 = pkin(3) * t44;
t95 = qJD(4) * t70;
t30 = pkin(4) * t95 - qJ(5) * t94 - t70 * qJD(5);
t47 = -t73 * mrSges(6,1) - t70 * mrSges(6,3);
t21 = t30 * t47;
t84 = t70 * mrSges(6,1) - t73 * mrSges(6,3);
t43 = t84 * qJD(4);
t83 = -t73 * pkin(4) - t70 * qJ(5);
t46 = -pkin(3) + t83;
t25 = t46 * t43;
t116 = t21 + t25 - t111;
t26 = t30 + t91;
t18 = t26 * t47;
t110 = t74 * pkin(2);
t37 = t46 - t110;
t20 = t37 * t43;
t60 = -pkin(3) - t110;
t27 = t60 * t44;
t48 = -t73 * mrSges(5,1) + t70 * mrSges(5,2);
t36 = t48 * t91;
t50 = mrSges(5,3) * t87;
t51 = mrSges(6,2) * t87;
t52 = mrSges(5,3) * t86;
t53 = mrSges(6,2) * t86;
t115 = t18 + t20 + t27 + t36 + t50 + t51 + t52 + t53;
t114 = 2 * m(5);
t113 = 2 * m(6);
t101 = t72 * t74;
t32 = pkin(1) * t101 + t71 * t61;
t29 = pkin(8) + t32;
t112 = t118 * t29;
t105 = t16 * mrSges(4,2);
t100 = t118 * pkin(8);
t59 = t71 * pkin(2) + pkin(8);
t99 = t119 * t59;
t98 = t119 * pkin(8);
t93 = qJD(5) * t73;
t92 = m(6) * t93;
t89 = t118 * t59 + t119 * t29;
t31 = -pkin(1) * t102 + t74 * t61;
t88 = mrSges(6,2) * t93 + Ifges(6,6) * t95 + (Ifges(6,4) + Ifges(5,5)) * t94;
t82 = (-mrSges(3,1) * t72 - mrSges(3,2) * t75) * qJD(2) * pkin(1);
t81 = (-mrSges(4,1) * t71 - mrSges(4,2) * t74) * qJD(3) * pkin(2);
t80 = t83 * mrSges(6,2) - Ifges(5,6) * t70;
t79 = (t117 * t73 - t120 * t70) * t95 + (t117 * t70 + t120 * t73) * t94;
t17 = t61 * t97 + (t72 * t96 + (t71 * t75 + t101) * qJD(2)) * pkin(1);
t78 = m(6) * t83 + t47 + t48;
t77 = m(6) * (-pkin(4) * t70 + qJ(5) * t73) - t84 - t85;
t4 = t17 + t30;
t1 = t4 * t47;
t10 = mrSges(5,3) * t103;
t11 = mrSges(6,2) * t103;
t22 = -t31 + t46;
t14 = t22 * t43;
t15 = t17 * mrSges(4,1);
t28 = -pkin(3) - t31;
t19 = t28 * t44;
t5 = t17 * t48;
t8 = mrSges(5,3) * t104;
t9 = mrSges(6,2) * t104;
t76 = t1 + t10 + t11 + t14 - t15 + t19 + t5 + t79 + t8 + t9;
t63 = mrSges(6,2) * t94;
t2 = [(t28 * t17 + t112) * t114 + (t22 * t4 + t112) * t113 + 0.2e1 * m(4) * (t32 * t16 - t31 * t17) + 0.2e1 * t82 + 0.2e1 * t19 + 0.2e1 * t14 - 0.2e1 * t15 + 0.2e1 * t11 + 0.2e1 * t10 + 0.2e1 * t9 + 0.2e1 * t8 + 0.2e1 * t5 - 0.2e1 * t105 + 0.2e1 * t1 + t79; m(4) * (t16 * t71 - t17 * t74 - t31 * t97 + t32 * t96) * pkin(2) + m(6) * (t26 * t22 + t37 * t4 + t89) + (-t16 - t90) * mrSges(4,2) + t82 - mrSges(4,1) * t91 + t76 + (t60 * t17 + t28 * t91 + t89) * m(5) + t115; 0.2e1 * t81 + (t37 * t26 + t99) * t113 + (t60 * t91 + t99) * t114 + 0.2e1 * t50 + 0.2e1 * t51 + 0.2e1 * t52 + 0.2e1 * t53 + 0.2e1 * t36 + 0.2e1 * t27 + 0.2e1 * t20 + 0.2e1 * t18 + t79; m(5) * (-pkin(3) * t17 + t100) + m(6) * (t30 * t22 + t46 * t4 + t100) - t105 + t76 + t116; t81 + m(6) * (t46 * t26 + t30 * t37 + t98) + m(5) * (-pkin(3) * t91 + t98) + t79 + t115 + t116; t46 * t30 * t113 - 0.2e1 * t111 + 0.2e1 * t21 + 0.2e1 * t25 + t79; t29 * t92 + t77 * t16 + (t78 * t29 + t80) * qJD(4) + t88; t59 * t92 + t77 * t90 + (t78 * t59 + t80) * qJD(4) + t88; t80 * qJD(4) + (t78 * qJD(4) + t92) * pkin(8) + t88; 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t63 + (t70 * t16 + t29 * t94) * m(6); t63 + (t59 * t94 + t70 * t90) * m(6); m(6) * pkin(8) * t94 + t63; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t2(1), t2(2), t2(4), t2(7), t2(11); t2(2), t2(3), t2(5), t2(8), t2(12); t2(4), t2(5), t2(6), t2(9), t2(13); t2(7), t2(8), t2(9), t2(10), t2(14); t2(11), t2(12), t2(13), t2(14), t2(15);];
Mq = res;
