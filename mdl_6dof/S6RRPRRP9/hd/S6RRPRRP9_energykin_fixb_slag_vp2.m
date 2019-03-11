% Calculate kinetic energy for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:25
% EndTime: 2019-03-09 12:27:26
% DurationCPUTime: 0.56s
% Computational Cost: add. (1044->110), mult. (2428->163), div. (0->0), fcn. (1922->10), ass. (0->45)
t128 = sin(qJ(5));
t131 = cos(qJ(5));
t133 = cos(qJ(2));
t125 = sin(pkin(6));
t138 = t125 * qJD(1);
t135 = t133 * t138;
t119 = qJD(4) - t135;
t130 = sin(qJ(2));
t139 = cos(pkin(6)) * qJD(1);
t137 = pkin(1) * t139;
t118 = pkin(8) * t135 + t130 * t137;
t123 = qJD(2) + t139;
t113 = qJ(3) * t123 + t118;
t114 = (-pkin(2) * t133 - qJ(3) * t130 - pkin(1)) * t138;
t124 = sin(pkin(11));
t126 = cos(pkin(11));
t104 = t126 * t113 + t124 * t114;
t136 = t130 * t138;
t115 = t123 * t126 - t124 * t136;
t100 = pkin(9) * t115 + t104;
t129 = sin(qJ(4));
t132 = cos(qJ(4));
t103 = -t113 * t124 + t126 * t114;
t116 = t123 * t124 + t126 * t136;
t97 = -pkin(3) * t135 - pkin(9) * t116 + t103;
t92 = t132 * t100 + t129 * t97;
t90 = pkin(10) * t119 + t92;
t117 = -pkin(8) * t136 + t133 * t137;
t110 = -pkin(2) * t123 + qJD(3) - t117;
t106 = -pkin(3) * t115 + t110;
t107 = t115 * t132 - t116 * t129;
t108 = t115 * t129 + t116 * t132;
t95 = -pkin(4) * t107 - pkin(10) * t108 + t106;
t86 = t128 * t95 + t131 * t90;
t85 = -t128 * t90 + t131 * t95;
t91 = -t129 * t100 + t132 * t97;
t89 = -pkin(4) * t119 - t91;
t134 = qJD(1) ^ 2;
t105 = qJD(5) - t107;
t102 = t108 * t131 + t119 * t128;
t101 = -t108 * t128 + t119 * t131;
t87 = -pkin(5) * t101 + qJD(6) + t89;
t84 = qJ(6) * t101 + t86;
t83 = pkin(5) * t105 - qJ(6) * t102 + t85;
t1 = t134 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t125 ^ 2 * t134 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t83 ^ 2 + t84 ^ 2 + t87 ^ 2) / 0.2e1 + m(5) * (t106 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t110 ^ 2) / 0.2e1 + (t117 * mrSges(3,1) - t118 * mrSges(3,2) + Ifges(3,3) * t123 / 0.2e1) * t123 + (t91 * mrSges(5,1) - t92 * mrSges(5,2) + Ifges(5,3) * t119 / 0.2e1) * t119 + (t110 * mrSges(4,2) - t103 * mrSges(4,3) + Ifges(4,1) * t116 / 0.2e1) * t116 + (t106 * mrSges(5,2) - t91 * mrSges(5,3) + Ifges(5,5) * t119 + Ifges(5,1) * t108 / 0.2e1) * t108 + ((-t117 * mrSges(3,3) + Ifges(3,5) * t123 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) * t130 / 0.2e1) * t138) * t130 + (-t103 * mrSges(4,1) + t104 * mrSges(4,2) + t118 * mrSges(3,3) - Ifges(4,5) * t116 + Ifges(3,6) * t123 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t130 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t133) * t138) * t133) * t138 + (-Ifges(4,6) * t135 - t110 * mrSges(4,1) + t104 * mrSges(4,3) + Ifges(4,4) * t116 + Ifges(4,2) * t115 / 0.2e1) * t115 + (-t106 * mrSges(5,1) + t92 * mrSges(5,3) + Ifges(5,4) * t108 + Ifges(5,6) * t119 + Ifges(5,2) * t107 / 0.2e1) * t107 + (t85 * mrSges(6,1) + t83 * mrSges(7,1) - t86 * mrSges(6,2) - t84 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t105) * t105 + (t89 * mrSges(6,2) + t87 * mrSges(7,2) - t85 * mrSges(6,3) - t83 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t102 + (Ifges(6,5) + Ifges(7,5)) * t105) * t102 + (-t89 * mrSges(6,1) - t87 * mrSges(7,1) + t86 * mrSges(6,3) + t84 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t101 + (Ifges(6,6) + Ifges(7,6)) * t105 + (Ifges(6,4) + Ifges(7,4)) * t102) * t101;
T  = t1;
