% Calculate vector of inverse dynamics joint torques for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:51
% EndTime: 2019-12-31 16:28:55
% DurationCPUTime: 2.38s
% Computational Cost: add. (559->204), mult. (1326->264), div. (0->0), fcn. (678->6), ass. (0->97)
t143 = Ifges(5,2) + Ifges(4,2);
t145 = Ifges(5,4) + Ifges(4,4);
t154 = Ifges(5,1) + Ifges(4,1);
t144 = Ifges(5,5) + Ifges(4,5);
t142 = Ifges(5,6) + Ifges(4,6);
t55 = sin(qJ(3));
t115 = Ifges(5,4) * t55;
t117 = Ifges(4,4) * t55;
t57 = cos(qJ(3));
t153 = t143 * t57 + t115 + t117;
t46 = pkin(3) * t57 + pkin(2);
t72 = -mrSges(5,1) * t57 + mrSges(5,2) * t55;
t74 = -mrSges(4,1) * t57 + mrSges(4,2) * t55;
t152 = m(4) * pkin(2) + m(5) * t46 + mrSges(3,1) - t72 - t74;
t54 = -qJ(4) - pkin(5);
t151 = -m(4) * pkin(5) + m(5) * t54 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t100 = qJD(2) * qJD(3);
t32 = qJDD(2) * t57 - t100 * t55;
t150 = t32 / 0.2e1;
t104 = qJD(2) * t57;
t149 = t145 * t104;
t52 = sin(pkin(6));
t53 = cos(pkin(6));
t148 = g(1) * t53 + g(2) * t52;
t147 = -mrSges(4,2) - mrSges(5,2);
t141 = t153 * qJD(2) + t142 * qJD(3);
t105 = qJD(2) * t55;
t140 = t144 * qJD(3) + t154 * t105 + t149;
t30 = t72 * qJD(2);
t139 = -t74 * qJD(2) - t30;
t118 = mrSges(5,2) * t57;
t58 = cos(qJ(2));
t107 = qJD(1) * t58;
t22 = -qJD(2) * t46 + qJD(4) - t107;
t43 = -qJD(2) * pkin(2) - t107;
t73 = mrSges(4,1) * t55 + mrSges(4,2) * t57;
t138 = t43 * t73 + t22 * (mrSges(5,1) * t55 + t118);
t137 = -t142 * t55 + t144 * t57;
t136 = m(5) * pkin(3) + mrSges(5,1);
t135 = t145 * t57;
t103 = qJD(3) * t55;
t101 = qJD(1) * qJD(2);
t45 = t58 * t101;
t56 = sin(qJ(2));
t35 = t56 * qJDD(1) + t45;
t28 = qJDD(2) * pkin(5) + t35;
t108 = qJD(1) * t56;
t42 = qJD(2) * pkin(5) + t108;
t3 = -t103 * t42 + t57 * t28;
t102 = qJD(3) * t57;
t90 = t42 * t102;
t4 = -t28 * t55 - t90;
t77 = t3 * t57 - t4 * t55;
t133 = -mrSges(4,1) - t136;
t132 = t42 * (t55 ^ 2 + t57 ^ 2);
t33 = qJDD(2) * t55 + t100 * t57;
t99 = qJD(2) * qJD(4);
t1 = -t90 + qJDD(3) * pkin(3) - qJ(4) * t33 + (-t28 - t99) * t55;
t79 = qJ(4) * qJD(2) + t42;
t10 = t79 * t57;
t2 = qJ(4) * t32 + t57 * t99 + t3;
t9 = t79 * t55;
t8 = qJD(3) * pkin(3) - t9;
t131 = -t1 * t55 - t10 * t103 - t8 * t102 + t2 * t57;
t59 = qJD(2) ^ 2;
t113 = t55 * t58;
t112 = t57 * t58;
t94 = mrSges(5,3) * t105;
t36 = qJD(3) * mrSges(5,1) - t94;
t37 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t105;
t111 = -t36 - t37;
t89 = mrSges(5,3) * t104;
t38 = -qJD(3) * mrSges(5,2) + t89;
t39 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t104;
t110 = t38 + t39;
t96 = pkin(3) * t103;
t93 = t55 * t107;
t92 = t57 * t107;
t44 = t56 * t101;
t6 = -t32 * mrSges(5,1) + t33 * mrSges(5,2);
t82 = qJD(3) * t54;
t34 = qJDD(1) * t58 - t44;
t76 = t10 * t57 - t55 * t8;
t65 = t55 * (Ifges(4,1) * t57 - t117);
t64 = t55 * (Ifges(5,1) * t57 - t115);
t27 = -qJDD(2) * pkin(2) - t34;
t41 = t54 * t57;
t40 = t54 * t55;
t17 = -qJD(4) * t55 + t57 * t82;
t16 = qJD(4) * t57 + t55 * t82;
t15 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t33;
t14 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t33;
t13 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t32;
t12 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t32;
t7 = -mrSges(4,1) * t32 + mrSges(4,2) * t33;
t5 = -pkin(3) * t32 + qJDD(4) + t27;
t11 = [m(2) * qJDD(1) + (-m(2) - m(3) - m(4) - m(5)) * g(3) + (qJDD(2) * mrSges(3,1) - t59 * mrSges(3,2) - t6 - t7 + (t110 * t57 + t111 * t55) * qJD(2) + m(3) * t34 + m(5) * (t10 * t104 - t105 * t8 - t5) + m(4) * (qJD(2) * t132 - t27)) * t58 + (-t59 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + (t12 + t13) * t57 + (-t14 - t15) * t55 - t139 * qJD(2) + (-t110 * t55 + t111 * t57) * qJD(3) + m(3) * t35 + m(5) * (qJD(2) * t22 + t131) + m(4) * (qJD(2) * t43 + t77)) * t56; (t1 * t40 + t10 * t16 + t17 * t8 - t2 * t41 + t22 * t96 - t46 * t5 - (t22 * t56 + t58 * t76) * qJD(1)) * m(5) + t139 * t108 + t140 * t102 / 0.2e1 - t141 * t103 / 0.2e1 + (t142 * t57 + t144 * t55) * qJDD(3) + t77 * mrSges(4,3) + (-pkin(2) * t27 - (t132 * t58 + t43 * t56) * qJD(1)) * m(4) + t131 * mrSges(5,3) + t135 * t33 / 0.2e1 + (t45 - t35) * mrSges(3,2) - t39 * t92 + (t44 + t34) * mrSges(3,1) + (-t92 + t16) * t38 + (t93 + t17) * t36 + (t138 + t137 * qJD(3) / 0.2e1) * qJD(3) + (t143 * t32 + t145 * t33) * t57 / 0.2e1 + t30 * t96 + t37 * t93 + t5 * t72 + t27 * t74 + t40 * t14 - t41 * t12 - t46 * t6 + ((-t143 * t55 + t135) * t57 + t65 + t64) * t100 / 0.2e1 + t153 * t150 + (-pkin(5) * t15 + t145 * t150 + t154 * t33) * t55 + (m(4) * t77 - t102 * t37 - t103 * t39 + t57 * t13) * pkin(5) + t148 * (t151 * t58 + t152 * t56) + (t151 * t56 - t152 * t58) * g(3) + Ifges(3,3) * qJDD(2) - pkin(2) * t7; t4 * mrSges(4,1) - t3 * mrSges(4,2) - t2 * mrSges(5,2) + t9 * t38 + t8 * t89 + t141 * t105 / 0.2e1 - t137 * t100 / 0.2e1 + (-t65 / 0.2e1 - t64 / 0.2e1) * t59 + (t37 * t57 + t39 * t55) * t42 + t144 * t33 + t142 * t32 + (t136 * t55 + t118 + t73) * g(3) * t56 + (-m(5) * (-t8 - t9) + t94 + t36) * t10 + t136 * t1 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) - t138 * qJD(2) + (t147 * (-t112 * t52 + t53 * t55) + t133 * (-t113 * t52 - t53 * t57)) * g(2) + (t147 * (-t112 * t53 - t52 * t55) + t133 * (-t113 * t53 + t52 * t57)) * g(1) - (-t105 * t143 + t140 + t149) * t104 / 0.2e1 + (t14 + (-m(5) * t22 - t30) * t105) * pkin(3); (t55 * t36 - t57 * t38) * qJD(2) + (g(3) * t58 - t76 * qJD(2) - t148 * t56 + t5) * m(5) + t6;];
tau = t11;
