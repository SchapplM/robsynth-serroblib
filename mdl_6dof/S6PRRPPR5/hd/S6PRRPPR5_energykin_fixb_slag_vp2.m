% Calculate kinetic energy for
% S6PRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:05
% EndTime: 2019-03-08 21:18:06
% DurationCPUTime: 0.35s
% Computational Cost: add. (342->96), mult. (729->140), div. (0->0), fcn. (451->10), ass. (0->41)
t130 = -pkin(3) - qJ(5);
t112 = sin(pkin(11));
t114 = cos(pkin(11));
t118 = sin(qJ(2));
t113 = sin(pkin(6));
t129 = qJD(1) * t113;
t105 = qJD(2) * pkin(8) + t118 * t129;
t117 = sin(qJ(3));
t120 = cos(qJ(3));
t115 = cos(pkin(6));
t128 = qJD(1) * t115;
t98 = -t117 * t105 + t120 * t128;
t123 = qJD(4) - t98;
t127 = qJD(2) * t117;
t90 = pkin(4) * t127 + t130 * qJD(3) + t123;
t124 = -qJ(4) * t117 - pkin(2);
t121 = cos(qJ(2));
t125 = t121 * t129;
t95 = -t125 + (t130 * t120 + t124) * qJD(2);
t86 = t112 * t90 + t114 * t95;
t99 = t120 * t105 + t117 * t128;
t126 = qJD(2) * t120;
t97 = -qJD(3) * qJ(4) - t99;
t85 = -t112 * t95 + t114 * t90;
t93 = pkin(4) * t126 + qJD(5) - t97;
t119 = cos(qJ(6));
t116 = sin(qJ(6));
t108 = qJD(6) + t127;
t106 = -qJD(2) * pkin(2) - t125;
t104 = qJD(3) * t114 - t112 * t126;
t103 = -qJD(3) * t112 - t114 * t126;
t100 = -t125 + (-pkin(3) * t120 + t124) * qJD(2);
t96 = -qJD(3) * pkin(3) + t123;
t92 = t103 * t116 + t104 * t119;
t91 = t103 * t119 - t104 * t116;
t87 = -pkin(5) * t103 + t93;
t84 = pkin(9) * t103 + t86;
t83 = pkin(5) * t127 - pkin(9) * t104 + t85;
t82 = t116 * t83 + t119 * t84;
t81 = -t116 * t84 + t119 * t83;
t1 = m(4) * (t106 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + (t87 * mrSges(7,2) - t81 * mrSges(7,3) + Ifges(7,1) * t92 / 0.2e1) * t92 + (t93 * mrSges(6,2) - t85 * mrSges(6,3) + Ifges(6,1) * t104 / 0.2e1) * t104 + (-t87 * mrSges(7,1) + t82 * mrSges(7,3) + Ifges(7,4) * t92 + Ifges(7,2) * t91 / 0.2e1) * t91 + (-t93 * mrSges(6,1) + t86 * mrSges(6,3) + Ifges(6,4) * t104 + Ifges(6,2) * t103 / 0.2e1) * t103 + (m(3) * (t115 ^ 2 + (t118 ^ 2 + t121 ^ 2) * t113 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t81 * mrSges(7,1) - t82 * mrSges(7,2) + Ifges(7,5) * t92 + Ifges(7,6) * t91 + Ifges(7,3) * t108 / 0.2e1) * t108 + (t98 * mrSges(4,1) - t99 * mrSges(4,2) + t96 * mrSges(5,2) - t97 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3)) * qJD(3) + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t121 - mrSges(3,2) * t118) * t129 + (-t106 * mrSges(4,1) - t97 * mrSges(5,1) + t100 * mrSges(5,2) + t99 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t126 + (-Ifges(5,5) + Ifges(4,6)) * qJD(3)) * t120 + (t96 * mrSges(5,1) + t85 * mrSges(6,1) + t106 * mrSges(4,2) - t86 * mrSges(6,2) - t98 * mrSges(4,3) - t100 * mrSges(5,3) + Ifges(6,5) * t104 + Ifges(6,6) * t103 + (-Ifges(5,4) + Ifges(4,5)) * qJD(3) + ((Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,3) / 0.2e1) * t117 + (Ifges(4,4) + Ifges(5,6)) * t120) * qJD(2)) * t117) * qJD(2);
T  = t1;
