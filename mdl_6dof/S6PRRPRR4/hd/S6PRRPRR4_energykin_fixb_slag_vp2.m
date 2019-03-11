% Calculate kinetic energy for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:11
% EndTime: 2019-03-08 22:10:11
% DurationCPUTime: 0.39s
% Computational Cost: add. (340->94), mult. (703->142), div. (0->0), fcn. (439->10), ass. (0->39)
t122 = sin(qJ(5));
t123 = sin(qJ(3));
t126 = cos(qJ(5));
t127 = cos(qJ(3));
t104 = (t122 * t123 + t126 * t127) * qJD(2);
t124 = sin(qJ(2));
t118 = sin(pkin(6));
t135 = qJD(1) * t118;
t108 = qJD(2) * pkin(8) + t124 * t135;
t119 = cos(pkin(6));
t134 = qJD(1) * t119;
t100 = -t123 * t108 + t127 * t134;
t130 = qJD(4) - t100;
t132 = t123 * qJD(2);
t92 = -pkin(9) * t132 + (-pkin(3) - pkin(4)) * qJD(3) + t130;
t133 = qJD(2) * t127;
t101 = t127 * t108 + t123 * t134;
t97 = qJD(3) * qJ(4) + t101;
t94 = -pkin(9) * t133 + t97;
t89 = t122 * t92 + t126 * t94;
t128 = cos(qJ(2));
t109 = -qJD(2) * pkin(2) - t128 * t135;
t102 = -pkin(3) * t133 - qJ(4) * t132 + t109;
t88 = -t122 * t94 + t126 * t92;
t95 = pkin(4) * t133 - t102;
t125 = cos(qJ(6));
t121 = sin(qJ(6));
t116 = -qJD(3) + qJD(5);
t105 = (-t122 * t127 + t123 * t126) * qJD(2);
t103 = qJD(6) + t104;
t99 = t105 * t125 + t116 * t121;
t98 = -t105 * t121 + t116 * t125;
t96 = -qJD(3) * pkin(3) + t130;
t90 = pkin(5) * t104 - pkin(10) * t105 + t95;
t87 = pkin(10) * t116 + t89;
t86 = -pkin(5) * t116 - t88;
t85 = t121 * t90 + t125 * t87;
t84 = -t121 * t87 + t125 * t90;
t1 = m(4) * (t100 ^ 2 + t101 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t95 ^ 2) / 0.2e1 + m(7) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + (t86 * mrSges(7,2) - t84 * mrSges(7,3) + Ifges(7,1) * t99 / 0.2e1) * t99 + (t88 * mrSges(6,1) - t89 * mrSges(6,2) + Ifges(6,3) * t116 / 0.2e1) * t116 + (-t86 * mrSges(7,1) + t85 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,2) * t98 / 0.2e1) * t98 + (t95 * mrSges(6,2) - t88 * mrSges(6,3) + Ifges(6,5) * t116 + Ifges(6,1) * t105 / 0.2e1) * t105 + (m(2) / 0.2e1 + m(3) * (t119 ^ 2 + (t124 ^ 2 + t128 ^ 2) * t118 ^ 2) / 0.2e1) * qJD(1) ^ 2 - (-t95 * mrSges(6,1) + t89 * mrSges(6,3) + Ifges(6,4) * t105 + Ifges(6,6) * t116 - Ifges(6,2) * t104 / 0.2e1) * t104 + (t84 * mrSges(7,1) - t85 * mrSges(7,2) + Ifges(7,5) * t99 + Ifges(7,6) * t98 + Ifges(7,3) * t103 / 0.2e1) * t103 + (t100 * mrSges(4,1) - t96 * mrSges(5,1) - t101 * mrSges(4,2) + t97 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(3)) * qJD(3) + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t128 - mrSges(3,2) * t124) * t135 + (-t109 * mrSges(4,1) - t102 * mrSges(5,1) + t97 * mrSges(5,2) + t101 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t133 + (Ifges(4,6) - Ifges(5,6)) * qJD(3)) * t127 + (t109 * mrSges(4,2) + t96 * mrSges(5,2) - t100 * mrSges(4,3) - t102 * mrSges(5,3) + (Ifges(5,4) + Ifges(4,5)) * qJD(3) + ((Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t123 + (Ifges(4,4) - Ifges(5,5)) * t127) * qJD(2)) * t123) * qJD(2);
T  = t1;
