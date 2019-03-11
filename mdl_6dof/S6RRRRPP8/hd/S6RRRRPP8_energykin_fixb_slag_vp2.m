% Calculate kinetic energy for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:30:13
% EndTime: 2019-03-09 21:30:13
% DurationCPUTime: 0.53s
% Computational Cost: add. (738->107), mult. (1624->149), div. (0->0), fcn. (1224->8), ass. (0->42)
t134 = -pkin(4) - pkin(5);
t133 = cos(qJ(3));
t132 = cos(qJ(4));
t119 = sin(qJ(4));
t122 = cos(qJ(2));
t121 = sin(qJ(2));
t117 = sin(pkin(6));
t130 = t117 * qJD(1);
t128 = t121 * t130;
t131 = cos(pkin(6)) * qJD(1);
t129 = pkin(1) * t131;
t108 = -pkin(8) * t128 + t122 * t129;
t116 = qJD(2) + t131;
t101 = -pkin(2) * t116 - t108;
t120 = sin(qJ(3));
t106 = t133 * t116 - t120 * t128;
t107 = t120 * t116 + t133 * t128;
t88 = -pkin(3) * t106 - pkin(10) * t107 + t101;
t127 = t122 * t130;
t112 = qJD(3) - t127;
t109 = pkin(8) * t127 + t121 * t129;
t102 = pkin(9) * t116 + t109;
t104 = (-pkin(2) * t122 - pkin(9) * t121 - pkin(1)) * t130;
t94 = t133 * t102 + t120 * t104;
t92 = pkin(10) * t112 + t94;
t86 = t119 * t88 + t132 * t92;
t93 = -t120 * t102 + t133 * t104;
t105 = qJD(4) - t106;
t83 = t105 * qJ(5) + t86;
t126 = pkin(3) * t112 + t93;
t85 = -t119 * t92 + t132 * t88;
t125 = qJD(5) - t85;
t96 = t132 * t107 + t119 * t112;
t124 = qJ(5) * t96 + t126;
t123 = qJD(1) ^ 2;
t95 = t107 * t119 - t132 * t112;
t84 = pkin(4) * t95 - t124;
t82 = -t105 * pkin(4) + t125;
t81 = t134 * t95 + qJD(6) + t124;
t80 = qJ(6) * t95 + t83;
t79 = -t96 * qJ(6) + t134 * t105 + t125;
t1 = t123 * Ifges(2,3) / 0.2e1 + m(3) * (t117 ^ 2 * t123 * pkin(1) ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t126 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + (t93 * mrSges(4,1) - t94 * mrSges(4,2) + Ifges(4,3) * t112 / 0.2e1) * t112 + (t101 * mrSges(4,2) - t93 * mrSges(4,3) + Ifges(4,5) * t112 + Ifges(4,1) * t107 / 0.2e1) * t107 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t122 / 0.2e1) * t122 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t122 + Ifges(3,1) * t121 / 0.2e1) * t121) * t130 + (-t108 * t121 + t109 * t122) * mrSges(3,3)) * t130 + (t108 * mrSges(3,1) - t109 * mrSges(3,2) + Ifges(3,3) * t116 / 0.2e1 + (Ifges(3,5) * t121 + Ifges(3,6) * t122) * t130) * t116 + (-t101 * mrSges(4,1) + t94 * mrSges(4,3) + Ifges(4,4) * t107 + Ifges(4,6) * t112 + Ifges(4,2) * t106 / 0.2e1) * t106 + (-t126 * mrSges(5,2) + t82 * mrSges(6,2) + t81 * mrSges(7,2) - t85 * mrSges(5,3) - t84 * mrSges(6,3) - t79 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t96) * t96 + (-t126 * mrSges(5,1) + t84 * mrSges(6,1) - t81 * mrSges(7,1) - t83 * mrSges(6,2) - t86 * mrSges(5,3) + t80 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t95 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t96) * t95 + (t85 * mrSges(5,1) - t82 * mrSges(6,1) - t79 * mrSges(7,1) - t86 * mrSges(5,2) + t80 * mrSges(7,2) + t83 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t105 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t96 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t95) * t105;
T  = t1;
