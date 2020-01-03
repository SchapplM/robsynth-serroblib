% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR10_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR10_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:01
% EndTime: 2019-12-31 18:04:03
% DurationCPUTime: 1.04s
% Computational Cost: add. (3398->122), mult. (6807->175), div. (0->0), fcn. (7448->6), ass. (0->71)
t106 = sin(qJ(5));
t107 = cos(qJ(5));
t104 = sin(pkin(8));
t105 = cos(pkin(8));
t108 = cos(qJ(4));
t151 = sin(qJ(4));
t83 = -t104 * t151 - t105 * t108;
t84 = t104 * t108 - t105 * t151;
t113 = t106 * t83 + t107 * t84;
t115 = -t106 * t84 + t107 * t83;
t178 = t113 * mrSges(6,1) + t115 * mrSges(6,2);
t136 = -pkin(6) + qJ(2);
t90 = t136 * t104;
t91 = t136 * t105;
t68 = t108 * t90 - t151 * t91;
t111 = -t84 * pkin(7) + t68;
t69 = t108 * t91 + t151 * t90;
t45 = t83 * pkin(7) + t69;
t163 = -t106 * t45 + t107 * t111;
t28 = t106 * t111 + t107 * t45;
t3 = -t28 * mrSges(6,1) - t163 * mrSges(6,2) + Ifges(6,5) * t115 - Ifges(6,6) * t113;
t177 = t3 * qJD(5);
t167 = -t115 / 0.2e1;
t153 = t84 * pkin(4);
t176 = m(6) * t153;
t86 = -t106 * t151 + t107 * t108;
t87 = -t106 * t108 - t107 * t151;
t116 = t87 * mrSges(6,1) - t86 * mrSges(6,2);
t174 = t116 * qJD(5);
t150 = Ifges(6,4) * t113;
t156 = t113 / 0.2e1;
t157 = -t113 / 0.2e1;
t158 = t115 / 0.2e1;
t121 = -t105 * pkin(2) - t104 * qJ(3) - pkin(1);
t79 = t105 * pkin(3) - t121;
t67 = -t83 * pkin(4) + t79;
t2 = (Ifges(6,2) * t115 + t150) * t157 + t67 * t178 + (Ifges(6,1) * t115 - t150) * t156 + (0.2e1 * Ifges(6,4) * t115 + (Ifges(6,1) - Ifges(6,2)) * t113) * t158;
t159 = m(6) * pkin(4);
t162 = t84 ^ 2;
t161 = m(5) / 0.2e1;
t160 = m(6) / 0.2e1;
t117 = t158 + t167;
t118 = t156 + t157;
t6 = (t117 * t86 + t118 * t87) * mrSges(6,3);
t155 = t6 * qJD(4);
t137 = t84 * mrSges(5,1);
t133 = t106 * t115;
t130 = t107 * t113;
t129 = t2 * qJD(1);
t128 = t6 * qJD(1);
t7 = (-t113 ^ 2 - t115 ^ 2) * mrSges(6,3) + (-t83 ^ 2 - t162) * mrSges(5,3) + m(6) * (t113 * t163 - t115 * t28) + m(5) * (t68 * t84 - t69 * t83) + (mrSges(3,3) + mrSges(4,2) + 0.4e1 * (m(4) / 0.4e1 + m(3) / 0.4e1) * qJ(2)) * (t104 ^ 2 + t105 ^ 2);
t127 = t7 * qJD(1);
t12 = t83 * mrSges(5,2) + t137 + (-t133 / 0.2e1 + t130 / 0.2e1 + t84 / 0.2e1) * t159 + t178;
t126 = t12 * qJD(1);
t30 = -mrSges(6,1) * t115 + mrSges(6,2) * t113;
t15 = (-m(4) * t121 + m(5) * t79 + m(6) * t67 + t105 * mrSges(4,1) - t83 * mrSges(5,1) + t84 * mrSges(5,2) + t104 * mrSges(4,3) + t30) * t104;
t125 = t15 * qJD(1);
t109 = (t113 * t86 + t115 * t87) * t160 + (t108 * t84 - t151 * t83) * t161;
t16 = (m(4) + t160 + t161) * t104 + t109;
t124 = t16 * qJD(1);
t18 = 0.2e1 * t157 * mrSges(6,1) + 0.2e1 * t167 * mrSges(6,2);
t123 = t18 * qJD(1);
t1 = t30 * t153 + t67 * t176 + t79 * t137 - t162 * Ifges(5,4) + (t79 * mrSges(5,2) + Ifges(5,4) * t83 + (Ifges(5,1) - Ifges(5,2)) * t84) * t83 + t2;
t114 = t1 * qJD(1) + t6 * qJD(3);
t4 = t117 * Ifges(6,5) + t118 * Ifges(6,6);
t88 = (mrSges(6,1) * t106 + mrSges(6,2) * t107) * pkin(4);
t110 = -t4 * qJD(1) + t88 * qJD(4);
t85 = t88 * qJD(5);
t17 = t109 - (m(5) + m(6)) * t104 / 0.2e1;
t13 = -t176 / 0.2e1 + (t130 - t133) * t159 / 0.2e1;
t5 = [t7 * qJD(2) + t15 * qJD(3) + t1 * qJD(4) + t2 * qJD(5), t17 * qJD(3) + t13 * qJD(4) + t127, t17 * qJD(2) + t125 + t155, t13 * qJD(2) + t177 + t114 + (-t69 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,5) * t83 - Ifges(5,6) * t84 + (m(6) * (t106 * t163 - t107 * t28) + (-t106 * t113 - t107 * t115) * mrSges(6,3)) * pkin(4) + t3) * qJD(4), t3 * qJD(4) + t129 + t177; -t16 * qJD(3) - t12 * qJD(4) + t18 * qJD(5) - t127, 0, -t124, -t126, t123; t16 * qJD(2) - t125 + t155, t124, 0, t128 + (-t151 * mrSges(5,1) - t108 * mrSges(5,2) + t116 + (t106 * t86 + t107 * t87) * t159) * qJD(4) + t174, qJD(4) * t116 + t174; t12 * qJD(2) + t4 * qJD(5) - t114, t126, -t128, -t85, -t85 - t110; -t18 * qJD(2) - t4 * qJD(4) - t129, -t123, 0, t110, 0;];
Cq = t5;
