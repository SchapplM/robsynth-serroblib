% Calculate kinetic energy for
% S6RRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:56:23
% EndTime: 2019-03-10 03:56:24
% DurationCPUTime: 0.72s
% Computational Cost: add. (1506->114), mult. (3378->183), div. (0->0), fcn. (2746->12), ass. (0->53)
t136 = sin(qJ(5));
t141 = cos(qJ(5));
t149 = cos(pkin(6)) * qJD(1);
t132 = qJD(2) + t149;
t138 = sin(qJ(3));
t143 = cos(qJ(3));
t139 = sin(qJ(2));
t133 = sin(pkin(6));
t150 = qJD(1) * t133;
t147 = t139 * t150;
t122 = t132 * t143 - t138 * t147;
t123 = t132 * t138 + t143 * t147;
t137 = sin(qJ(4));
t142 = cos(qJ(4));
t114 = t122 * t137 + t123 * t142;
t144 = cos(qJ(2));
t146 = t144 * t150;
t128 = qJD(3) - t146;
t127 = qJD(4) + t128;
t148 = pkin(1) * t149;
t125 = pkin(8) * t146 + t139 * t148;
t120 = pkin(9) * t132 + t125;
t121 = (-pkin(2) * t144 - pkin(9) * t139 - pkin(1)) * t150;
t111 = -t120 * t138 + t143 * t121;
t108 = pkin(3) * t128 - pkin(10) * t123 + t111;
t112 = t143 * t120 + t138 * t121;
t110 = pkin(10) * t122 + t112;
t98 = t142 * t108 - t110 * t137;
t95 = pkin(4) * t127 - pkin(11) * t114 + t98;
t113 = t122 * t142 - t123 * t137;
t99 = t137 * t108 + t142 * t110;
t97 = pkin(11) * t113 + t99;
t92 = t136 * t95 + t141 * t97;
t91 = -t136 * t97 + t141 * t95;
t103 = t113 * t141 - t114 * t136;
t124 = -pkin(8) * t147 + t144 * t148;
t119 = -pkin(2) * t132 - t124;
t115 = -pkin(3) * t122 + t119;
t105 = -pkin(4) * t113 + t115;
t145 = qJD(1) ^ 2;
t140 = cos(qJ(6));
t135 = sin(qJ(6));
t126 = qJD(5) + t127;
t104 = t113 * t136 + t114 * t141;
t102 = qJD(6) - t103;
t101 = t104 * t140 + t126 * t135;
t100 = -t104 * t135 + t126 * t140;
t93 = -pkin(5) * t103 - pkin(12) * t104 + t105;
t90 = pkin(12) * t126 + t92;
t89 = -pkin(5) * t126 - t91;
t88 = t135 * t93 + t140 * t90;
t87 = -t135 * t90 + t140 * t93;
t1 = m(6) * (t105 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t115 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t119 ^ 2) / 0.2e1 + t145 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t133 ^ 2 * t145 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(7) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + (t111 * mrSges(4,1) - t112 * mrSges(4,2) + Ifges(4,3) * t128 / 0.2e1) * t128 + (t98 * mrSges(5,1) - t99 * mrSges(5,2) + Ifges(5,3) * t127 / 0.2e1) * t127 + (t91 * mrSges(6,1) - t92 * mrSges(6,2) + Ifges(6,3) * t126 / 0.2e1) * t126 + (t87 * mrSges(7,1) - t88 * mrSges(7,2) + Ifges(7,3) * t102 / 0.2e1) * t102 + (t119 * mrSges(4,2) - t111 * mrSges(4,3) + Ifges(4,5) * t128 + Ifges(4,1) * t123 / 0.2e1) * t123 + (t115 * mrSges(5,2) - t98 * mrSges(5,3) + Ifges(5,5) * t127 + Ifges(5,1) * t114 / 0.2e1) * t114 + (t105 * mrSges(6,2) - t91 * mrSges(6,3) + Ifges(6,5) * t126 + Ifges(6,1) * t104 / 0.2e1) * t104 + (t89 * mrSges(7,2) - t87 * mrSges(7,3) + Ifges(7,5) * t102 + Ifges(7,1) * t101 / 0.2e1) * t101 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t144 / 0.2e1) * t144 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t144 + Ifges(3,1) * t139 / 0.2e1) * t139) * t150 + (-t124 * t139 + t125 * t144) * mrSges(3,3)) * t150 + (t124 * mrSges(3,1) - t125 * mrSges(3,2) + Ifges(3,3) * t132 / 0.2e1 + (Ifges(3,5) * t139 + Ifges(3,6) * t144) * t150) * t132 + (-t119 * mrSges(4,1) + t112 * mrSges(4,3) + Ifges(4,4) * t123 + Ifges(4,6) * t128 + Ifges(4,2) * t122 / 0.2e1) * t122 + (-t115 * mrSges(5,1) + t99 * mrSges(5,3) + Ifges(5,4) * t114 + Ifges(5,6) * t127 + Ifges(5,2) * t113 / 0.2e1) * t113 + (-t105 * mrSges(6,1) + t92 * mrSges(6,3) + Ifges(6,4) * t104 + Ifges(6,6) * t126 + Ifges(6,2) * t103 / 0.2e1) * t103 + (-t89 * mrSges(7,1) + t88 * mrSges(7,3) + Ifges(7,4) * t101 + Ifges(7,6) * t102 + Ifges(7,2) * t100 / 0.2e1) * t100;
T  = t1;
