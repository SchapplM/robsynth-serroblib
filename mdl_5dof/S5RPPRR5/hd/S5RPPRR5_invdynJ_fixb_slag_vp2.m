% Calculate vector of inverse dynamics joint torques for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:33
% DurationCPUTime: 1.82s
% Computational Cost: add. (1588->204), mult. (2403->263), div. (0->0), fcn. (1143->10), ass. (0->97)
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t47 = -mrSges(6,1) * t72 + mrSges(6,2) * t69;
t153 = mrSges(5,1) - t47;
t114 = t69 * mrSges(6,3);
t66 = -qJD(1) + qJD(4);
t41 = qJD(5) * mrSges(6,1) - t114 * t66;
t111 = t72 * mrSges(6,3);
t42 = -qJD(5) * mrSges(6,2) + t111 * t66;
t68 = cos(pkin(8));
t61 = -pkin(1) * t68 - pkin(2);
t55 = -pkin(3) + t61;
t35 = qJD(1) * t55 + qJD(3);
t67 = sin(pkin(8));
t58 = pkin(1) * t67 + qJ(3);
t45 = t58 * qJD(1);
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t15 = t35 * t70 + t45 * t73;
t12 = pkin(7) * t66 + t15;
t7 = -qJD(2) * t72 - t12 * t69;
t83 = t69 * qJD(2) - t12 * t72;
t89 = t69 * t7 + t72 * t83;
t132 = m(6) * t89 + t66 * mrSges(5,2) + t69 * t41 - t72 * t42;
t152 = -m(5) * t15 + t132;
t151 = m(6) * pkin(4) + t153;
t14 = t35 * t73 - t45 * t70;
t11 = -pkin(4) * t66 - t14;
t94 = t153 * t66;
t150 = -m(6) * t11 + t94;
t130 = t69 / 0.2e1;
t149 = m(4) + m(5);
t148 = m(5) * t14 + t150;
t147 = qJD(4) * t14;
t86 = mrSges(6,1) * t69 + mrSges(6,2) * t72;
t143 = t11 * t86;
t141 = mrSges(4,1) + mrSges(3,1);
t140 = -mrSges(4,3) + mrSges(3,2);
t23 = t70 * t55 + t73 * t58;
t105 = qJD(5) * t72;
t106 = qJD(5) * t69;
t139 = -t7 * t105 + t106 * t83;
t107 = qJD(4) * t15;
t34 = qJDD(1) * t55 + qJDD(3);
t102 = qJD(1) * qJD(3);
t36 = qJDD(1) * t58 + t102;
t6 = t34 * t73 - t36 * t70 - t107;
t65 = -qJDD(1) + qJDD(4);
t32 = -t106 * t66 + t65 * t72;
t18 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t32;
t33 = t105 * t66 + t65 * t69;
t19 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t33;
t138 = t72 * t18 - t69 * t19;
t4 = -pkin(4) * t65 - t6;
t137 = m(6) * t4 - mrSges(6,1) * t32 + mrSges(6,2) * t33;
t101 = qJ(1) + pkin(8);
t62 = cos(t101);
t92 = sin(t101);
t28 = -t62 * t73 - t70 * t92;
t29 = t62 * t70 - t73 * t92;
t136 = t29 * mrSges(5,2) + t151 * t28;
t135 = t28 * mrSges(5,2) - t151 * t29;
t90 = -t69 * t83 + t7 * t72;
t5 = t70 * t34 + t73 * t36 + t147;
t3 = pkin(7) * t65 + t5;
t1 = qJD(5) * t7 - qJDD(2) * t69 + t3 * t72;
t2 = qJD(5) * t83 - qJDD(2) * t72 - t3 * t69;
t91 = t1 * t72 - t2 * t69;
t131 = m(6) * (-qJD(5) * t90 + t91) + t138 - t41 * t105 - t42 * t106;
t71 = sin(qJ(1));
t127 = t71 * pkin(1);
t74 = cos(qJ(1));
t64 = t74 * pkin(1);
t125 = Ifges(6,4) * t69;
t124 = Ifges(6,4) * t72;
t123 = Ifges(6,2) * t72;
t120 = t65 * mrSges(5,1);
t119 = t65 * mrSges(5,2);
t116 = t66 * t69;
t115 = t66 * t72;
t104 = qJDD(1) * mrSges(4,1);
t103 = m(3) + t149;
t100 = t62 * pkin(2) + t92 * qJ(3) + t64;
t95 = -m(6) * pkin(7) - mrSges(6,3);
t93 = t62 * pkin(3) + t100;
t85 = Ifges(6,1) * t72 - t125;
t49 = t123 + t125;
t22 = t55 * t73 - t58 * t70;
t38 = qJD(5) * (Ifges(6,5) * t72 - Ifges(6,6) * t69);
t79 = -pkin(2) * t92 + t62 * qJ(3) - t127;
t77 = -pkin(3) * t92 + t79;
t24 = Ifges(6,6) * qJD(5) + t49 * t66;
t54 = Ifges(6,4) * t115;
t25 = Ifges(6,1) * t116 + Ifges(6,5) * qJD(5) + t54;
t76 = Ifges(5,3) * t65 + qJD(5) * t143 + t4 * t47 + t6 * mrSges(5,1) + t1 * t111 + t32 * t49 / 0.2e1 + t33 * (Ifges(6,1) * t69 + t124) / 0.2e1 - t5 * mrSges(5,2) + (Ifges(6,1) * t33 + Ifges(6,4) * t32) * t130 + t72 * (Ifges(6,4) * t33 + Ifges(6,2) * t32) / 0.2e1 - t24 * t106 / 0.2e1 + t25 * t105 / 0.2e1 - t2 * t114 + t139 * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t130 + Ifges(6,6) * t72) * qJDD(5) + (t38 + t85 * t116 + (-Ifges(6,2) * t69 + t124) * t115) * qJD(5) / 0.2e1;
t43 = qJDD(1) * t61 + qJDD(3);
t8 = [-t76 + m(4) * (qJD(3) * t45 + t36 * t58 + t43 * t61) + m(5) * (t22 * t6 + t23 * t5) + t22 * t120 - t43 * mrSges(4,1) - t61 * t104 - t23 * t119 + t137 * (pkin(4) - t22) - t148 * (qJD(3) * t70 + t23 * qJD(4)) + (t102 + t36) * mrSges(4,3) - t152 * (qJD(3) * t73 + qJD(4) * t22) + t131 * (-pkin(7) + t23) + (-mrSges(2,1) * t74 + mrSges(2,2) * t71 - m(6) * (pkin(7) * t29 + t93) - t29 * mrSges(6,3) - m(5) * t93 - m(4) * t100 - m(3) * t64 + t140 * t92 - t141 * t62 + t136) * g(2) + (mrSges(2,1) * t71 + mrSges(2,2) * t74 - m(6) * (t28 * pkin(7) + t77) - t28 * mrSges(6,3) - m(5) * t77 - m(4) * t79 + m(3) * t127 + t141 * t92 + t140 * t62 + t135) * g(1) + (t58 * mrSges(4,3) + Ifges(4,2) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t68 * mrSges(3,1) - 0.2e1 * t67 * mrSges(3,2) + m(3) * (t67 ^ 2 + t68 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1); -t42 * t105 + t41 * t106 - t69 * t18 - t72 * t19 + m(6) * (qJD(5) * t89 - t1 * t69 - t2 * t72) + t103 * qJDD(2) + (-m(6) - t103) * g(3); -t104 + m(4) * t43 + (t120 + m(5) * (t6 + t107) - t132 * qJD(4) - t137) * t73 + (-t119 + (-t72 * t41 - t69 * t42) * qJD(5) - t94 * qJD(4) + m(6) * (qJD(4) * t11 + t139 + t91) + m(5) * (t5 - t147) + t138) * t70 + (-m(4) * t45 - mrSges(4,3) * qJD(1) + t148 * t70 + t152 * t73) * qJD(1) + (m(6) + t149) * (-g(1) * t92 + g(2) * t62); t76 + (-t29 * t95 - t136) * g(2) + (-t28 * t95 - t135) * g(1) - t137 * pkin(4) + t150 * t15 + t132 * t14 + t131 * pkin(7); t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t33 + Ifges(6,6) * t32 + Ifges(6,3) * qJDD(5) - g(3) * t47 - t83 * t41 - t7 * t42 + (-t143 + t24 * t130 - t38 / 0.2e1 + (-t69 * t85 / 0.2e1 + t123 * t130) * t66 + t90 * mrSges(6,3) - (t25 + t54) * t72 / 0.2e1) * t66 + (-t28 * g(1) - t29 * g(2)) * t86;];
tau = t8;
