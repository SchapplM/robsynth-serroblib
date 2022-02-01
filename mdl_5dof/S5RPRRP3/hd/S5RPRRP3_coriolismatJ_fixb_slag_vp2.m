% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:37
% EndTime: 2022-01-23 09:29:39
% DurationCPUTime: 0.99s
% Computational Cost: add. (3137->120), mult. (5826->158), div. (0->0), fcn. (5622->6), ass. (0->75)
t145 = cos(qJ(4));
t146 = cos(qJ(3));
t95 = sin(qJ(4));
t96 = sin(qJ(3));
t86 = t145 * t146 - t95 * t96;
t114 = -cos(pkin(8)) * pkin(1) - pkin(2);
t88 = -t146 * pkin(3) + t114;
t59 = -t86 * pkin(4) + t88;
t87 = -t145 * t96 - t95 * t146;
t173 = m(6) * t59 - mrSges(6,1) * t86 - mrSges(6,2) * t87;
t106 = (Ifges(5,6) + Ifges(6,6)) * t87 + (Ifges(5,5) + Ifges(6,5)) * t86;
t89 = sin(pkin(8)) * pkin(1) + pkin(6);
t76 = (-pkin(7) - t89) * t96;
t110 = t146 * t89;
t77 = t146 * pkin(7) + t110;
t156 = t145 * t76 - t95 * t77;
t160 = t156 * mrSges(5,2);
t159 = t87 * qJ(5) + t156;
t164 = t159 * mrSges(6,2);
t55 = t145 * t77 + t95 * t76;
t165 = t55 * mrSges(5,1);
t33 = t86 * qJ(5) + t55;
t169 = t33 * mrSges(6,1);
t172 = t106 - t160 - t164 - t165 - t169;
t171 = -t160 / 0.2e1 - t164 / 0.2e1 - t165 / 0.2e1 - t169 / 0.2e1;
t117 = t145 * pkin(3);
t92 = t117 + pkin(4);
t103 = t117 - t92;
t166 = m(6) * t103;
t161 = Ifges(6,4) + Ifges(5,4);
t162 = mrSges(6,3) * t159 + t161 * t86;
t155 = (mrSges(5,2) + mrSges(6,2)) * t145;
t144 = mrSges(6,3) * t86;
t152 = m(6) / 0.2e1;
t154 = (-t144 / 0.2e1 - t33 * t152) * pkin(4) + t171;
t151 = m(6) * pkin(4);
t149 = m(6) * t87;
t148 = pkin(3) * t95;
t147 = t96 * pkin(3);
t135 = t87 * mrSges(6,1);
t13 = m(6) * (t159 * t87 + t33 * t86) + (t86 ^ 2 + t87 ^ 2) * mrSges(6,3);
t125 = t13 * qJD(1);
t78 = t86 * mrSges(6,2);
t108 = -t78 + t135;
t72 = t86 * t148;
t51 = t92 * t87 + t72;
t68 = -t87 * pkin(4) + t147;
t19 = 0.2e1 * (t51 / 0.4e1 - t68 / 0.4e1) * m(6) + t108;
t123 = t19 * qJD(1);
t115 = mrSges(6,1) + t151;
t34 = t115 * t87 - t78;
t122 = t34 * qJD(1);
t121 = t87 * t148;
t120 = t92 * t144;
t118 = t51 * t152;
t111 = t146 * mrSges(4,2);
t109 = t161 * t87;
t57 = -t87 * mrSges(5,1) + t86 * mrSges(5,2);
t107 = -Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2);
t105 = t86 * t117;
t104 = -t57 - t78;
t100 = -t108 * t59 - t144 * t159 + t88 * t57;
t1 = t146 ^ 2 * Ifges(4,4) + t114 * t111 + t162 * t86 + (t114 * mrSges(4,1) - Ifges(4,4) * t96 + (m(5) * t88 - mrSges(5,1) * t86) * pkin(3) + (Ifges(4,1) - Ifges(4,2)) * t146) * t96 + (-mrSges(5,2) * t147 + t107 * t86 - t109) * t87 + t100 + t173 * t68;
t102 = t1 * qJD(1);
t2 = (-t173 * pkin(4) - t109) * t87 + (t107 * t87 + t162) * t86 + t100;
t101 = t2 * qJD(1);
t18 = (pkin(4) / 0.2e1 + t117 / 0.2e1 - t92 / 0.2e1) * t149;
t43 = t155 * pkin(3) + (mrSges(5,1) + mrSges(6,1) - t166) * t148;
t97 = t120 / 0.2e1 - mrSges(6,3) * t105 / 0.2e1 - t33 * t166 / 0.2e1 - t171;
t5 = t97 + t154;
t99 = t5 * qJD(1) + t18 * qJD(2) + t43 * qJD(3);
t22 = t68 * t152 + t118;
t14 = -t103 * t149 / 0.2e1 + (mrSges(6,1) + t151 / 0.2e1) * t87 + t104;
t4 = -t97 + t106 + t154;
t3 = [qJD(3) * t1 + qJD(4) * t2 + qJD(5) * t13, 0, (mrSges(6,3) * t121 - t120 - mrSges(4,1) * t110 + m(6) * (t148 * t159 - t92 * t33) + m(5) * (-t145 * t55 + t156 * t95) * pkin(3) + Ifges(4,5) * t146 + (t89 * mrSges(4,2) - Ifges(4,6)) * t96 + (t121 - t105) * mrSges(5,3) + t172) * qJD(3) + t4 * qJD(4) + t22 * qJD(5) + t102, t4 * qJD(3) + ((-m(6) * t33 - t144) * pkin(4) + t172) * qJD(4) + t101, t22 * qJD(3) + t125; 0, 0, t14 * qJD(4) + (-t96 * mrSges(4,1) + t104 - t111 + t135 + m(5) * (t87 * t117 + t72) + 0.2e1 * t118) * qJD(3), t14 * qJD(3) + (t34 - t57) * qJD(4), 0; -qJD(4) * t5 + qJD(5) * t19 - t102, -qJD(4) * t18, -t43 * qJD(4), ((-mrSges(5,1) - t115) * t95 - t155) * qJD(4) * pkin(3) - t99, t123; qJD(3) * t5 + qJD(5) * t34 - t101, qJD(3) * t18, t99, 0, t122; -t19 * qJD(3) - t34 * qJD(4) - t125, 0, -t123, -t122, 0;];
Cq = t3;
