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
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:03:13
% EndTime: 2019-12-05 18:03:16
% DurationCPUTime: 0.94s
% Computational Cost: add. (3137->122), mult. (5826->159), div. (0->0), fcn. (5622->6), ass. (0->72)
t145 = sin(qJ(3));
t146 = cos(qJ(3));
t95 = sin(qJ(4));
t96 = cos(qJ(4));
t86 = -t95 * t145 + t96 * t146;
t114 = -cos(pkin(8)) * pkin(1) - pkin(2);
t88 = -t146 * pkin(3) + t114;
t59 = -t86 * pkin(4) + t88;
t87 = -t96 * t145 - t95 * t146;
t169 = m(6) * t59 - t86 * mrSges(6,1) - t87 * mrSges(6,2);
t89 = sin(pkin(8)) * pkin(1) + pkin(6);
t112 = t145 * t89;
t76 = -t145 * pkin(7) - t112;
t113 = t146 * t89;
t77 = t146 * pkin(7) + t113;
t154 = t96 * t76 - t95 * t77;
t158 = t87 * qJ(5) + t154;
t55 = t95 * t76 + t96 * t77;
t29 = -t86 * qJ(5) - t55;
t168 = (Ifges(5,6) + Ifges(6,6)) * t87 + (Ifges(5,5) + Ifges(6,5)) * t86 - t154 * mrSges(5,2) + t29 * mrSges(6,1) - t158 * mrSges(6,2) - t55 * mrSges(5,1);
t144 = mrSges(6,3) * t86;
t162 = m(6) * t29 - t144;
t160 = Ifges(6,4) + Ifges(5,4);
t161 = mrSges(6,3) * t158 + t160 * t86;
t147 = t96 * pkin(3);
t92 = pkin(4) + t147;
t157 = t92 - t147;
t156 = pkin(3) * t95;
t152 = m(6) / 0.2e1;
t151 = m(6) * pkin(4);
t150 = -t29 / 0.2e1;
t148 = m(6) * t87;
t135 = t87 * mrSges(6,1);
t134 = t29 * t92;
t132 = t95 * t87;
t131 = mrSges(5,2) + mrSges(6,2);
t13 = m(6) * (t158 * t87 - t29 * t86) + (t86 ^ 2 + t87 ^ 2) * mrSges(6,3);
t123 = t13 * qJD(1);
t78 = t86 * mrSges(6,2);
t110 = -t78 + t135;
t72 = t86 * t156;
t51 = t92 * t87 + t72;
t116 = t145 * pkin(3);
t68 = -t87 * pkin(4) + t116;
t19 = 0.2e1 * (t51 / 0.4e1 - t68 / 0.4e1) * m(6) + t110;
t121 = t19 * qJD(1);
t115 = mrSges(6,1) + t151;
t34 = t115 * t87 - t78;
t120 = t34 * qJD(1);
t119 = t151 / 0.2e1;
t118 = t51 * t152;
t111 = t160 * t87;
t57 = -t87 * mrSges(5,1) + t86 * mrSges(5,2);
t109 = -Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2);
t107 = -t57 - t78;
t106 = t147 / 0.2e1 - t92 / 0.2e1;
t100 = -t110 * t59 - t144 * t158 + t88 * t57;
t99 = t145 * mrSges(4,1) + t146 * mrSges(4,2);
t1 = t114 * t99 + m(5) * t88 * t116 + (-mrSges(5,1) * t116 + t161) * t86 + (-mrSges(5,2) * t116 + t109 * t86 - t111) * t87 + t100 + (-Ifges(4,2) + Ifges(4,1)) * t146 * t145 + (-t145 ^ 2 + t146 ^ 2) * Ifges(4,4) + t169 * t68;
t103 = t1 * qJD(1);
t2 = (-t169 * pkin(4) - t111) * t87 + (t109 * t87 + t161) * t86 + t100;
t102 = t2 * qJD(1);
t101 = pkin(4) / 0.2e1 + t106;
t98 = t29 * t147;
t18 = t101 * t148;
t43 = t131 * t147 + (t157 * m(6) + mrSges(5,1) + mrSges(6,1)) * t156;
t5 = (t150 + t29 / 0.2e1) * mrSges(6,1) + (t134 / 0.2e1 - t98 / 0.2e1 + pkin(4) * t150) * m(6) + t101 * t144;
t97 = t5 * qJD(1) - t18 * qJD(2) - t43 * qJD(3);
t22 = t68 * t152 + t118;
t14 = t157 * t148 / 0.2e1 + (mrSges(6,1) + t119) * t87 + t107;
t4 = (-t98 + t134) * t152 + t29 * t119 + (-pkin(4) / 0.2e1 + t106) * t144 + t168;
t3 = [t1 * qJD(3) + t2 * qJD(4) + t13 * qJD(5), 0, t4 * qJD(4) + t22 * qJD(5) + t103 + (-mrSges(4,1) * t113 + mrSges(4,2) * t112 + Ifges(4,5) * t146 - Ifges(4,6) * t145 + (mrSges(6,3) * t132 + (-t96 * t86 + t132) * mrSges(5,3) + m(6) * t158 * t95 + m(5) * (t154 * t95 - t55 * t96)) * pkin(3) + t162 * t92 + t168) * qJD(3), t4 * qJD(3) + (pkin(4) * t162 + t168) * qJD(4) + t102, t22 * qJD(3) + t123; 0, 0, t14 * qJD(4) + (t107 - t99 + t135 + (t87 * t147 + t72) * m(5) + 0.2e1 * t118) * qJD(3), t14 * qJD(3) + (t34 - t57) * qJD(4), 0; t5 * qJD(4) + t19 * qJD(5) - t103, -t18 * qJD(4), -t43 * qJD(4), (-t131 * t96 + (-mrSges(5,1) - t115) * t95) * qJD(4) * pkin(3) + t97, t121; -t5 * qJD(3) + t34 * qJD(5) - t102, t18 * qJD(3), -t97, 0, t120; -t19 * qJD(3) - t34 * qJD(4) - t123, 0, -t121, -t120, 0;];
Cq = t3;
