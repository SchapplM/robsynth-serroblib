% Calculate kinetic energy for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:26
% EndTime: 2019-03-08 21:12:26
% DurationCPUTime: 0.41s
% Computational Cost: add. (352->95), mult. (783->139), div. (0->0), fcn. (515->10), ass. (0->39)
t121 = sin(qJ(2));
t117 = sin(pkin(6));
t132 = qJD(1) * t117;
t108 = qJD(2) * pkin(8) + t121 * t132;
t120 = sin(qJ(3));
t123 = cos(qJ(3));
t118 = cos(pkin(6));
t131 = qJD(1) * t118;
t102 = t123 * t108 + t120 * t131;
t100 = qJD(3) * qJ(4) + t102;
t124 = cos(qJ(2));
t128 = t124 * t132;
t103 = -t128 + (-pkin(3) * t123 - qJ(4) * t120 - pkin(2)) * qJD(2);
t116 = sin(pkin(11));
t133 = cos(pkin(11));
t93 = t133 * t100 + t116 * t103;
t101 = -t120 * t108 + t123 * t131;
t130 = qJD(2) * t120;
t129 = qJD(2) * t123;
t92 = -t116 * t100 + t133 * t103;
t90 = -qJ(5) * t129 + t93;
t127 = qJD(3) * pkin(3) - qJD(4) + t101;
t89 = pkin(4) * t129 + qJD(5) - t92;
t107 = t116 * qJD(3) + t133 * t130;
t126 = qJ(5) * t107 + t127;
t122 = cos(qJ(6));
t119 = sin(qJ(6));
t112 = qJD(6) + t129;
t109 = -qJD(2) * pkin(2) - t128;
t106 = -t133 * qJD(3) + t116 * t130;
t95 = t106 * t119 + t107 * t122;
t94 = t106 * t122 - t107 * t119;
t91 = pkin(4) * t106 - t126;
t88 = (-pkin(4) - pkin(5)) * t106 + t126;
t87 = pkin(9) * t106 + t90;
t86 = pkin(5) * t129 - t107 * pkin(9) + t89;
t85 = t119 * t86 + t122 * t87;
t84 = -t119 * t87 + t122 * t86;
t1 = m(6) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t127 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t84 ^ 2 + t85 ^ 2 + t88 ^ 2) / 0.2e1 + (t88 * mrSges(7,2) - t84 * mrSges(7,3) + Ifges(7,1) * t95 / 0.2e1) * t95 + (t101 * mrSges(4,1) - t102 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t88 * mrSges(7,1) + t85 * mrSges(7,3) + Ifges(7,4) * t95 + Ifges(7,2) * t94 / 0.2e1) * t94 + (m(3) * (t118 ^ 2 + (t121 ^ 2 + t124 ^ 2) * t117 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t84 * mrSges(7,1) - t85 * mrSges(7,2) + Ifges(7,5) * t95 + Ifges(7,6) * t94 + Ifges(7,3) * t112 / 0.2e1) * t112 + (-t127 * mrSges(5,2) + t89 * mrSges(6,2) - t92 * mrSges(5,3) - t91 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t107) * t107 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t124 - mrSges(3,2) * t121) * t132 + (t109 * mrSges(4,2) - t101 * mrSges(4,3) + Ifges(4,5) * qJD(3) + Ifges(4,1) * t130 / 0.2e1) * t120 + (-t109 * mrSges(4,1) - t92 * mrSges(5,1) + t89 * mrSges(6,1) + t93 * mrSges(5,2) + t102 * mrSges(4,3) - t90 * mrSges(6,3) + Ifges(4,6) * qJD(3) + (-Ifges(6,4) - Ifges(5,5)) * t107 + (Ifges(4,4) * t120 + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t123) * qJD(2)) * t123) * qJD(2) + (-t127 * mrSges(5,1) + t91 * mrSges(6,1) - t90 * mrSges(6,2) - t93 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t106 + (-Ifges(5,4) + Ifges(6,5)) * t107 + (Ifges(5,6) - Ifges(6,6)) * t129) * t106;
T  = t1;
