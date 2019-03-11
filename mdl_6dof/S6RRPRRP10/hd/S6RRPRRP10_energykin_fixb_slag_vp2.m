% Calculate kinetic energy for
% S6RRPRRP10
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:27
% EndTime: 2019-03-09 12:36:27
% DurationCPUTime: 0.56s
% Computational Cost: add. (1042->110), mult. (2420->163), div. (0->0), fcn. (1914->10), ass. (0->45)
t139 = cos(qJ(5));
t128 = sin(qJ(5));
t132 = cos(qJ(2));
t125 = sin(pkin(6));
t137 = t125 * qJD(1);
t134 = t132 * t137;
t119 = qJD(4) - t134;
t129 = sin(qJ(4));
t131 = cos(qJ(4));
t130 = sin(qJ(2));
t138 = cos(pkin(6)) * qJD(1);
t136 = pkin(1) * t138;
t117 = pkin(8) * t134 + t130 * t136;
t123 = qJD(2) + t138;
t112 = qJ(3) * t123 + t117;
t113 = (-pkin(2) * t132 - qJ(3) * t130 - pkin(1)) * t137;
t124 = sin(pkin(11));
t126 = cos(pkin(11));
t102 = -t112 * t124 + t126 * t113;
t135 = t130 * t137;
t115 = t123 * t124 + t126 * t135;
t96 = -pkin(3) * t134 - pkin(9) * t115 + t102;
t103 = t126 * t112 + t124 * t113;
t114 = t123 * t126 - t124 * t135;
t99 = pkin(9) * t114 + t103;
t92 = t129 * t96 + t131 * t99;
t90 = pkin(10) * t119 + t92;
t116 = -pkin(8) * t135 + t132 * t136;
t109 = -pkin(2) * t123 + qJD(3) - t116;
t105 = -pkin(3) * t114 + t109;
t106 = t114 * t131 - t115 * t129;
t107 = t114 * t129 + t115 * t131;
t94 = -pkin(4) * t106 - pkin(10) * t107 + t105;
t86 = t128 * t94 + t139 * t90;
t91 = -t129 * t99 + t131 * t96;
t89 = -pkin(4) * t119 - t91;
t85 = -t128 * t90 + t139 * t94;
t133 = qJD(1) ^ 2;
t104 = qJD(5) - t106;
t101 = t139 * t107 + t128 * t119;
t100 = t107 * t128 - t139 * t119;
t87 = pkin(5) * t100 - qJ(6) * t101 + t89;
t84 = qJ(6) * t104 + t86;
t83 = -t104 * pkin(5) + qJD(6) - t85;
t1 = m(7) * (t83 ^ 2 + t84 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + t133 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t125 ^ 2 * t133 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + (t116 * mrSges(3,1) - t117 * mrSges(3,2) + Ifges(3,3) * t123 / 0.2e1) * t123 + (t91 * mrSges(5,1) - t92 * mrSges(5,2) + Ifges(5,3) * t119 / 0.2e1) * t119 + (t109 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,1) * t115 / 0.2e1) * t115 + (t105 * mrSges(5,2) - t91 * mrSges(5,3) + Ifges(5,5) * t119 + Ifges(5,1) * t107 / 0.2e1) * t107 + ((-t116 * mrSges(3,3) + Ifges(3,5) * t123 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) * t130 / 0.2e1) * t137) * t130 + (-t102 * mrSges(4,1) + t103 * mrSges(4,2) + t117 * mrSges(3,3) - Ifges(4,5) * t115 + Ifges(3,6) * t123 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t130 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t132) * t137) * t132) * t137 + (-Ifges(4,6) * t134 - t109 * mrSges(4,1) + t103 * mrSges(4,3) + Ifges(4,4) * t115 + Ifges(4,2) * t114 / 0.2e1) * t114 + (-t105 * mrSges(5,1) + t92 * mrSges(5,3) + Ifges(5,4) * t107 + Ifges(5,6) * t119 + Ifges(5,2) * t106 / 0.2e1) * t106 + (t85 * mrSges(6,1) - t83 * mrSges(7,1) - t86 * mrSges(6,2) + t84 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t104) * t104 + (t89 * mrSges(6,2) + t83 * mrSges(7,2) - t85 * mrSges(6,3) - t87 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t101 + (Ifges(7,4) + Ifges(6,5)) * t104) * t101 + (t89 * mrSges(6,1) + t87 * mrSges(7,1) - t84 * mrSges(7,2) - t86 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t100 + (-Ifges(6,6) + Ifges(7,6)) * t104 + (-Ifges(6,4) + Ifges(7,5)) * t101) * t100;
T  = t1;
