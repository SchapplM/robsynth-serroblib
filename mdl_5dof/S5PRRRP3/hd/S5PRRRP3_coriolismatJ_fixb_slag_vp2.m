% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:35
% EndTime: 2019-12-05 16:43:38
% DurationCPUTime: 1.00s
% Computational Cost: add. (2521->118), mult. (5210->156), div. (0->0), fcn. (5006->4), ass. (0->73)
t142 = cos(qJ(4));
t143 = cos(qJ(3));
t93 = sin(qJ(4));
t94 = sin(qJ(3));
t81 = t142 * t143 - t93 * t94;
t91 = -t143 * pkin(3) - pkin(2);
t62 = -t81 * pkin(4) + t91;
t82 = -t142 * t94 - t93 * t143;
t170 = m(6) * t62 - t81 * mrSges(6,1) - t82 * mrSges(6,2);
t104 = (Ifges(5,6) + Ifges(6,6)) * t82 + (Ifges(5,5) + Ifges(6,5)) * t81;
t86 = (-pkin(7) - pkin(6)) * t94;
t114 = t143 * pkin(6);
t87 = t143 * pkin(7) + t114;
t153 = t142 * t86 - t93 * t87;
t157 = t153 * mrSges(5,2);
t156 = t82 * qJ(5) + t153;
t161 = t156 * mrSges(6,2);
t58 = t142 * t87 + t93 * t86;
t162 = t58 * mrSges(5,1);
t45 = t81 * qJ(5) + t58;
t166 = t45 * mrSges(6,1);
t169 = t104 - t157 - t161 - t162 - t166;
t168 = -t157 / 0.2e1 - t161 / 0.2e1 - t162 / 0.2e1 - t166 / 0.2e1;
t113 = t142 * pkin(3);
t90 = t113 + pkin(4);
t101 = t113 - t90;
t163 = m(6) * t101;
t158 = Ifges(6,4) + Ifges(5,4);
t159 = mrSges(6,3) * t156 + t158 * t81;
t152 = (mrSges(5,2) + mrSges(6,2)) * t142;
t141 = mrSges(6,3) * t81;
t149 = m(6) / 0.2e1;
t151 = (-t45 * t149 - t141 / 0.2e1) * pkin(4) + t168;
t148 = m(6) * pkin(4);
t146 = m(6) * t82;
t145 = pkin(3) * t93;
t144 = t94 * pkin(3);
t132 = t82 * mrSges(6,1);
t13 = m(6) * (t156 * t82 + t45 * t81) + (t81 ^ 2 + t82 ^ 2) * mrSges(6,3);
t122 = t13 * qJD(2);
t73 = t81 * mrSges(6,2);
t106 = -t73 + t132;
t69 = t81 * t145;
t47 = t90 * t82 + t69;
t68 = -t82 * pkin(4) + t144;
t19 = 0.2e1 * (t47 / 0.4e1 - t68 / 0.4e1) * m(6) + t106;
t120 = t19 * qJD(2);
t111 = mrSges(6,1) + t148;
t24 = t111 * t82 - t73;
t119 = t24 * qJD(2);
t118 = t82 * t145;
t117 = t90 * t141;
t115 = t47 * t149;
t108 = t143 * mrSges(4,2);
t107 = t158 * t82;
t49 = -t82 * mrSges(5,1) + t81 * mrSges(5,2);
t105 = -Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2);
t103 = t81 * t113;
t102 = -t49 - t73;
t98 = -t106 * t62 - t141 * t156 + t91 * t49;
t1 = -pkin(2) * t108 + t143 ^ 2 * Ifges(4,4) + t159 * t81 + (-pkin(2) * mrSges(4,1) - Ifges(4,4) * t94 + (m(5) * t91 - t81 * mrSges(5,1)) * pkin(3) + (Ifges(4,1) - Ifges(4,2)) * t143) * t94 + (-mrSges(5,2) * t144 + t105 * t81 - t107) * t82 + t98 + t170 * t68;
t100 = t1 * qJD(2);
t2 = (-t170 * pkin(4) - t107) * t82 + (t105 * t82 + t159) * t81 + t98;
t99 = t2 * qJD(2);
t18 = (pkin(4) / 0.2e1 + t113 / 0.2e1 - t90 / 0.2e1) * t146;
t34 = t152 * pkin(3) + (mrSges(5,1) + mrSges(6,1) - t163) * t145;
t95 = t117 / 0.2e1 - mrSges(6,3) * t103 / 0.2e1 - t45 * t163 / 0.2e1 - t168;
t5 = t95 + t151;
t97 = t18 * qJD(1) + t5 * qJD(2) + t34 * qJD(3);
t21 = t68 * t149 + t115;
t14 = -t101 * t146 / 0.2e1 + (mrSges(6,1) + t148 / 0.2e1) * t82 + t102;
t4 = -t95 + t104 + t151;
t3 = [0, 0, t14 * qJD(4) + (-t94 * mrSges(4,1) + t102 - t108 + t132 + m(5) * (t82 * t113 + t69) + 0.2e1 * t115) * qJD(3), t14 * qJD(3) + (t24 - t49) * qJD(4), 0; 0, t1 * qJD(3) + t2 * qJD(4) + t13 * qJD(5), (mrSges(6,3) * t118 - t117 - mrSges(4,1) * t114 + Ifges(4,5) * t143 + m(6) * (t145 * t156 - t90 * t45) + m(5) * (-t142 * t58 + t153 * t93) * pkin(3) + (pkin(6) * mrSges(4,2) - Ifges(4,6)) * t94 + (t118 - t103) * mrSges(5,3) + t169) * qJD(3) + t4 * qJD(4) + t21 * qJD(5) + t100, t4 * qJD(3) + ((-m(6) * t45 - t141) * pkin(4) + t169) * qJD(4) + t99, t21 * qJD(3) + t122; -t18 * qJD(4), -t5 * qJD(4) + t19 * qJD(5) - t100, -t34 * qJD(4), ((-mrSges(5,1) - t111) * t93 - t152) * qJD(4) * pkin(3) - t97, t120; t18 * qJD(3), t5 * qJD(3) + t24 * qJD(5) - t99, t97, 0, t119; 0, -t19 * qJD(3) - t24 * qJD(4) - t122, -t120, -t119, 0;];
Cq = t3;
