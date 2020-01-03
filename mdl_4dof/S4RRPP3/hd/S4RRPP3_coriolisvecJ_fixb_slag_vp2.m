% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:27
% DurationCPUTime: 1.84s
% Computational Cost: add. (762->205), mult. (2186->279), div. (0->0), fcn. (1289->4), ass. (0->101)
t124 = mrSges(5,1) + mrSges(4,1);
t122 = Ifges(4,1) + Ifges(5,1);
t123 = Ifges(5,4) + Ifges(4,5);
t94 = -Ifges(4,4) + Ifges(5,5);
t120 = -Ifges(4,6) + Ifges(5,6);
t66 = sin(pkin(6));
t68 = cos(qJ(2));
t90 = cos(pkin(6));
t75 = t90 * t68;
t67 = sin(qJ(2));
t86 = qJD(1) * t67;
t42 = -qJD(1) * t75 + t66 * t86;
t103 = Ifges(5,5) * t42;
t40 = Ifges(4,4) * t42;
t76 = t90 * t67;
t52 = t66 * t68 + t76;
t44 = t52 * qJD(1);
t119 = t123 * qJD(2) + t122 * t44 + t103 - t40;
t100 = t42 * mrSges(4,3);
t101 = t42 * mrSges(5,2);
t31 = qJD(2) * mrSges(5,3) - t101;
t92 = -qJD(2) * mrSges(4,2) - t100 + t31;
t98 = t44 * mrSges(4,3);
t99 = t44 * mrSges(5,2);
t91 = t124 * qJD(2) - t98 - t99;
t105 = Ifges(3,4) * t67;
t85 = qJD(1) * t68;
t88 = Ifges(3,6) * qJD(2);
t118 = t88 / 0.2e1 + (t68 * Ifges(3,2) + t105) * qJD(1) / 0.2e1 + pkin(5) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t85);
t117 = -2 * pkin(1);
t115 = -t42 / 0.2e1;
t114 = t42 / 0.2e1;
t112 = t44 / 0.2e1;
t110 = pkin(2) * t66;
t109 = pkin(5) * (qJD(2) * mrSges(3,1) - mrSges(3,3) * t86);
t93 = -qJ(3) - pkin(5);
t59 = t93 * t68;
t23 = -t59 * t66 - t76 * t93;
t74 = qJD(2) * t93;
t41 = qJD(3) * t68 + t67 * t74;
t32 = t41 * qJD(1);
t71 = -t67 * qJD(3) + t68 * t74;
t69 = t90 * t71;
t4 = -qJD(1) * t69 + t32 * t66;
t107 = t23 * t4;
t106 = t4 * t52;
t104 = Ifges(4,4) * t44;
t72 = -t66 * t67 + t75;
t45 = t72 * qJD(2);
t37 = qJD(1) * t45;
t102 = t37 * mrSges(5,2);
t55 = qJD(1) * t59;
t97 = t66 * t55;
t96 = -qJD(2) / 0.2e1;
t95 = mrSges(5,2) + mrSges(4,3);
t70 = t66 * t71;
t5 = qJD(1) * t70 + t90 * t32;
t46 = t90 * t55;
t81 = t93 * t67;
t54 = qJD(1) * t81;
t50 = qJD(2) * pkin(2) + t54;
t20 = t66 * t50 - t46;
t89 = Ifges(3,5) * qJD(2);
t64 = -pkin(2) * t68 - pkin(1);
t87 = qJD(1) * t64;
t84 = qJD(2) * t67;
t83 = qJD(1) * qJD(2);
t82 = pkin(2) * t86;
t80 = t90 * pkin(2);
t79 = t67 * t83;
t78 = t68 * t83;
t73 = pkin(2) * t79;
t19 = t50 * t90 + t97;
t56 = qJD(3) + t87;
t43 = t52 * qJD(2);
t65 = Ifges(3,4) * t85;
t63 = -t80 - pkin(3);
t61 = qJ(4) + t110;
t49 = Ifges(3,1) * t86 + t65 + t89;
t39 = Ifges(5,5) * t44;
t36 = qJD(1) * t43;
t35 = t37 * mrSges(4,2);
t34 = t36 * mrSges(5,1);
t24 = -t59 * t90 + t66 * t81;
t22 = t54 * t90 + t97;
t21 = t54 * t66 - t46;
t18 = -pkin(3) * t72 - t52 * qJ(4) + t64;
t17 = mrSges(4,1) * t42 + mrSges(4,2) * t44;
t16 = mrSges(5,1) * t42 - mrSges(5,3) * t44;
t15 = qJD(2) * qJ(4) + t20;
t12 = -Ifges(4,2) * t42 + Ifges(4,6) * qJD(2) + t104;
t11 = Ifges(5,6) * qJD(2) + Ifges(5,3) * t42 + t39;
t10 = -qJD(2) * pkin(3) + qJD(4) - t19;
t9 = t41 * t90 + t70;
t8 = t41 * t66 - t69;
t7 = pkin(3) * t44 + qJ(4) * t42 + t82;
t6 = t42 * pkin(3) - t44 * qJ(4) + t56;
t3 = qJD(2) * qJD(4) + t5;
t2 = pkin(2) * t84 + pkin(3) * t43 - qJ(4) * t45 - qJD(4) * t52;
t1 = pkin(3) * t36 - qJ(4) * t37 - qJD(4) * t44 + t73;
t13 = [t18 * t34 + t1 * (-mrSges(5,1) * t72 - mrSges(5,3) * t52) + t64 * t35 + t2 * t16 + t92 * t9 - t91 * t8 + m(4) * (-t19 * t8 + t20 * t9 + t24 * t5 + t107) + m(5) * (t1 * t18 + t10 * t8 + t15 * t9 + t2 * t6 + t24 * t3 + t107) + (t5 * t72 + t106) * mrSges(4,3) + (t3 * t72 + t106) * mrSges(5,2) + (-t109 + t49 / 0.2e1 + t89 / 0.2e1 + (mrSges(3,2) * t117 + 0.3e1 / 0.2e1 * Ifges(3,4) * t68) * qJD(1)) * t68 * qJD(2) + (-mrSges(5,3) * t18 + t122 * t52 + t95 * t23 - t72 * t94) * t37 + (mrSges(4,1) * t64 + t94 * t52 - (Ifges(4,2) + Ifges(5,3)) * t72 - t95 * t24) * t36 + (-t88 / 0.2e1 + (mrSges(3,1) * t117 - 0.3e1 / 0.2e1 * t105 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t68) * qJD(1) + (qJD(1) * (-mrSges(4,1) * t72 + mrSges(4,2) * t52) + t17 + m(4) * (t56 + t87)) * pkin(2) - t118) * t84 + (t56 * mrSges(4,1) + Ifges(5,3) * t114 + t6 * mrSges(5,1) - Ifges(4,2) * t115 + t11 / 0.2e1 - t12 / 0.2e1 - t20 * mrSges(4,3) - t15 * mrSges(5,2) + (Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1) * qJD(2) + t94 * t112) * t43 + (t56 * mrSges(4,2) + Ifges(5,5) * t114 - t6 * mrSges(5,3) + Ifges(4,4) * t115 - t19 * mrSges(4,3) + t10 * mrSges(5,2) + (Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1) * qJD(2) + t122 * t112 + t119 / 0.2e1) * t45; -(-Ifges(3,2) * t86 + t49 + t65) * t85 / 0.2e1 + (-t67 * (Ifges(3,1) * t68 - t105) / 0.2e1 + pkin(1) * (mrSges(3,1) * t67 + mrSges(3,2) * t68)) * qJD(1) ^ 2 - t124 * t4 + (-mrSges(3,1) * t78 + mrSges(3,2) * t79) * pkin(5) + (-mrSges(4,3) * t80 + t123) * t37 + (t120 * t44 - t123 * t42) * t96 - t19 * t100 + (-Ifges(4,2) * t44 + t119 - t40) * t114 + (-mrSges(5,2) * t61 - mrSges(4,3) * t110 + t120) * t36 - (-t122 * t42 - t104 + t11 + t39) * t44 / 0.2e1 + t118 * t86 + t91 * t21 - t92 * t22 + ((-t4 * t90 + t5 * t66) * pkin(2) + t19 * t21 - t20 * t22 - t56 * t82) * m(4) - (Ifges(3,5) * t68 - Ifges(3,6) * t67) * t83 / 0.2e1 - t17 * t82 - Ifges(3,6) * t79 - t56 * (t44 * mrSges(4,1) - t42 * mrSges(4,2)) - t6 * (t44 * mrSges(5,1) + t42 * mrSges(5,3)) + qJD(4) * t31 - t7 * t16 + t3 * mrSges(5,3) - t5 * mrSges(4,2) + (-t10 * t21 + t3 * t61 + t4 * t63 - t6 * t7 + (-t22 + qJD(4)) * t15) * m(5) + Ifges(3,5) * t78 + t20 * t98 + t15 * t99 + t10 * t101 + t63 * t102 + t85 * t109 + t12 * t112 + (Ifges(5,3) * t44 - t103) * t115; t36 * mrSges(4,1) - t37 * mrSges(5,3) + t92 * t42 + t91 * t44 + t34 + t35 + (-t10 * t44 + t15 * t42 + t1) * m(5) + (t19 * t44 + t20 * t42 + t73) * m(4); t102 - qJD(2) * t31 + t44 * t16 + 0.2e1 * (t4 / 0.2e1 + t15 * t96 + t6 * t112) * m(5);];
tauc = t13(:);
