% Calculate kinetic energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:30
% EndTime: 2019-03-09 16:03:30
% DurationCPUTime: 0.42s
% Computational Cost: add. (584->109), mult. (1292->149), div. (0->0), fcn. (928->8), ass. (0->43)
t132 = -pkin(4) - pkin(10);
t118 = sin(qJ(3));
t121 = cos(qJ(3));
t119 = sin(qJ(2));
t122 = cos(qJ(2));
t115 = sin(pkin(6));
t130 = t115 * qJD(1);
t127 = t122 * t130;
t131 = cos(pkin(6)) * qJD(1);
t129 = pkin(1) * t131;
t105 = pkin(8) * t127 + t119 * t129;
t114 = qJD(2) + t131;
t97 = pkin(9) * t114 + t105;
t98 = (-pkin(2) * t122 - pkin(9) * t119 - pkin(1)) * t130;
t90 = t118 * t98 + t121 * t97;
t128 = t119 * t130;
t104 = -pkin(8) * t128 + t122 * t129;
t108 = -qJD(3) + t127;
t88 = -t108 * qJ(4) + t90;
t96 = -t114 * pkin(2) - t104;
t89 = -t118 * t97 + t121 * t98;
t126 = qJD(4) - t89;
t102 = -t121 * t114 + t118 * t128;
t103 = t114 * t118 + t121 * t128;
t86 = t102 * pkin(3) - t103 * qJ(4) + t96;
t85 = -qJ(5) * t102 - t88;
t125 = qJD(5) - t86;
t124 = -qJ(5) * t103 + t126;
t123 = qJD(1) ^ 2;
t120 = cos(qJ(6));
t117 = sin(qJ(6));
t101 = qJD(6) + t103;
t92 = t102 * t120 + t108 * t117;
t91 = -t102 * t117 + t108 * t120;
t87 = pkin(3) * t108 + t126;
t84 = -pkin(4) * t102 + t125;
t83 = -pkin(5) * t108 - t85;
t82 = (pkin(3) + pkin(4)) * t108 + t124;
t81 = (pkin(3) - t132) * t108 + t124;
t80 = pkin(5) * t103 + t102 * t132 + t125;
t79 = t117 * t80 + t120 * t81;
t78 = -t117 * t81 + t120 * t80;
t1 = m(3) * (pkin(1) ^ 2 * t115 ^ 2 * t123 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + t123 * Ifges(2,3) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t78 ^ 2 + t79 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t82 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + (t83 * mrSges(7,2) - t78 * mrSges(7,3) + Ifges(7,1) * t92 / 0.2e1) * t92 + (-t83 * mrSges(7,1) + t79 * mrSges(7,3) + Ifges(7,4) * t92 + Ifges(7,2) * t91 / 0.2e1) * t91 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t122 / 0.2e1) * t122 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t122 + Ifges(3,1) * t119 / 0.2e1) * t119) * t130 + (-t104 * t119 + t105 * t122) * mrSges(3,3)) * t130 + (t104 * mrSges(3,1) - t105 * mrSges(3,2) + Ifges(3,3) * t114 / 0.2e1 + (Ifges(3,5) * t119 + Ifges(3,6) * t122) * t130) * t114 + (t78 * mrSges(7,1) - t79 * mrSges(7,2) + Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(7,3) * t101 / 0.2e1) * t101 + (-t89 * mrSges(4,1) + t87 * mrSges(5,1) + t85 * mrSges(6,1) + t90 * mrSges(4,2) - t82 * mrSges(6,2) - t88 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t108) * t108 + (t84 * mrSges(6,1) + t96 * mrSges(4,2) + t87 * mrSges(5,2) - t89 * mrSges(4,3) - t86 * mrSges(5,3) - t82 * mrSges(6,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t103 + (-Ifges(5,4) - Ifges(4,5) - Ifges(6,6)) * t108) * t103 + (t96 * mrSges(4,1) + t86 * mrSges(5,1) - t88 * mrSges(5,2) + t84 * mrSges(6,2) - t90 * mrSges(4,3) - t85 * mrSges(6,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t102 + (Ifges(6,5) + Ifges(4,6) - Ifges(5,6)) * t108 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t103) * t102;
T  = t1;
