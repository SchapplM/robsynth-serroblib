% Calculate kinetic energy for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:07
% EndTime: 2019-03-09 07:08:08
% DurationCPUTime: 0.74s
% Computational Cost: add. (921->104), mult. (2211->164), div. (0->0), fcn. (1724->10), ass. (0->46)
t122 = cos(pkin(11));
t137 = t122 ^ 2;
t136 = qJD(1) * (pkin(7) + qJ(2));
t135 = m(3) / 0.2e1;
t124 = sin(qJ(5));
t128 = cos(qJ(5));
t120 = qJD(3) + qJD(4);
t125 = sin(qJ(4));
t129 = cos(qJ(4));
t121 = sin(pkin(11));
t114 = t121 * t136;
t115 = t122 * t136;
t126 = sin(qJ(3));
t130 = cos(qJ(3));
t105 = -t130 * t114 - t115 * t126;
t113 = (t121 * t130 + t122 * t126) * qJD(1);
t97 = qJD(3) * pkin(3) - pkin(8) * t113 + t105;
t106 = -t126 * t114 + t130 * t115;
t112 = (-t121 * t126 + t122 * t130) * qJD(1);
t98 = pkin(8) * t112 + t106;
t90 = t125 * t97 + t129 * t98;
t86 = pkin(9) * t120 + t90;
t103 = t112 * t129 - t125 * t113;
t104 = t112 * t125 + t113 * t129;
t116 = qJD(2) + (-pkin(2) * t122 - pkin(1)) * qJD(1);
t107 = -pkin(3) * t112 + t116;
t91 = -pkin(4) * t103 - pkin(9) * t104 + t107;
t82 = t124 * t91 + t128 * t86;
t81 = -t124 * t86 + t128 * t91;
t89 = -t125 * t98 + t129 * t97;
t102 = qJD(5) - t103;
t85 = -pkin(4) * t120 - t89;
t127 = cos(qJ(6));
t123 = sin(qJ(6));
t117 = -qJD(1) * pkin(1) + qJD(2);
t101 = qJD(6) + t102;
t100 = t104 * t128 + t120 * t124;
t99 = -t104 * t124 + t120 * t128;
t93 = t100 * t127 + t123 * t99;
t92 = -t100 * t123 + t127 * t99;
t83 = -pkin(5) * t99 + t85;
t80 = pkin(10) * t99 + t82;
t79 = pkin(5) * t102 - pkin(10) * t100 + t81;
t78 = t123 * t79 + t127 * t80;
t77 = -t123 * t80 + t127 * t79;
t1 = t117 ^ 2 * t135 + m(4) * (t105 ^ 2 + t106 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t107 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t78 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + (-t85 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,2) * t99 / 0.2e1) * t99 + (t83 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,1) * t93 / 0.2e1) * t93 + (t89 * mrSges(5,1) - t90 * mrSges(5,2) + Ifges(5,3) * t120 / 0.2e1) * t120 + (t116 * mrSges(4,2) - t105 * mrSges(4,3) + Ifges(4,1) * t113 / 0.2e1) * t113 + (-t83 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t93 + Ifges(7,2) * t92 / 0.2e1) * t92 + (-t116 * mrSges(4,1) + t106 * mrSges(4,3) + Ifges(4,4) * t113 + Ifges(4,2) * t112 / 0.2e1) * t112 + (t107 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,5) * t120 + Ifges(5,1) * t104 / 0.2e1) * t104 + (t81 * mrSges(6,1) - t82 * mrSges(6,2) + Ifges(6,6) * t99 + Ifges(6,3) * t102 / 0.2e1) * t102 + (-t107 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t104 + Ifges(5,6) * t120 + Ifges(5,2) * t103 / 0.2e1) * t103 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,5) * t93 + Ifges(7,6) * t92 + Ifges(7,3) * t101 / 0.2e1) * t101 + (t85 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,5) * t102 + Ifges(6,1) * t100 / 0.2e1) * t100 + (t105 * mrSges(4,1) - t106 * mrSges(4,2) + Ifges(4,5) * t113 + Ifges(4,6) * t112 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t117 * (-mrSges(3,1) * t122 + mrSges(3,2) * t121) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t135 + mrSges(3,3)) * (t121 ^ 2 + t137) * qJ(2) + Ifges(3,2) * t137 / 0.2e1 + (Ifges(3,4) * t122 + Ifges(3,1) * t121 / 0.2e1) * t121) * qJD(1)) * qJD(1);
T  = t1;
