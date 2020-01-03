% Calculate kinetic energy for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:46:30
% EndTime: 2019-12-31 22:46:31
% DurationCPUTime: 0.57s
% Computational Cost: add. (1012->99), mult. (2653->166), div. (0->0), fcn. (2168->12), ass. (0->51)
t130 = sin(qJ(2));
t134 = cos(qJ(2));
t124 = sin(pkin(5));
t142 = t124 * qJD(1);
t138 = t134 * t142;
t143 = cos(pkin(5)) * qJD(1);
t141 = pkin(1) * t143;
t118 = pkin(8) * t138 + t130 * t141;
t125 = cos(pkin(6));
t122 = qJD(2) + t143;
t123 = sin(pkin(6));
t146 = t122 * t123;
t108 = (t125 * t138 + t146) * pkin(9) + t118;
t115 = (-pkin(9) * t123 * t130 - pkin(2) * t134 - pkin(1)) * t142;
t129 = sin(qJ(3));
t133 = cos(qJ(3));
t121 = t134 * t141;
t139 = t130 * t142;
t110 = pkin(2) * t122 + t121 + (-pkin(9) * t125 - pkin(8)) * t139;
t147 = t110 * t125;
t97 = -t129 * t108 + (t115 * t123 + t147) * t133;
t144 = t125 * t134;
t111 = (-t129 * t130 + t133 * t144) * t142 + t133 * t146;
t128 = sin(qJ(4));
t132 = cos(qJ(4));
t101 = -t110 * t123 + t125 * t115;
t145 = t123 * t129;
t112 = t122 * t145 + (t129 * t144 + t130 * t133) * t142;
t93 = -pkin(3) * t111 - pkin(10) * t112 + t101;
t116 = t122 * t125 - t123 * t138 + qJD(3);
t98 = t133 * t108 + t115 * t145 + t129 * t147;
t96 = pkin(10) * t116 + t98;
t90 = t128 * t93 + t132 * t96;
t89 = -t128 * t96 + t132 * t93;
t103 = -t112 * t128 + t116 * t132;
t95 = -pkin(3) * t116 - t97;
t135 = qJD(1) ^ 2;
t131 = cos(qJ(5));
t127 = sin(qJ(5));
t117 = -pkin(8) * t139 + t121;
t109 = qJD(4) - t111;
t104 = t112 * t132 + t116 * t128;
t102 = qJD(5) - t103;
t100 = t104 * t131 + t109 * t127;
t99 = -t104 * t127 + t109 * t131;
t91 = -pkin(4) * t103 - pkin(11) * t104 + t95;
t88 = pkin(11) * t109 + t90;
t87 = -pkin(4) * t109 - t89;
t86 = t127 * t91 + t131 * t88;
t85 = -t127 * t88 + t131 * t91;
t1 = t135 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t124 ^ 2 * t135 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + (-t87 * mrSges(6,1) + t86 * mrSges(6,3) + Ifges(6,2) * t99 / 0.2e1) * t99 + (t97 * mrSges(4,1) - t98 * mrSges(4,2) + Ifges(4,3) * t116 / 0.2e1) * t116 + (t89 * mrSges(5,1) - t90 * mrSges(5,2) + Ifges(5,3) * t109 / 0.2e1) * t109 + (t101 * mrSges(4,2) - t97 * mrSges(4,3) + Ifges(4,5) * t116 + Ifges(4,1) * t112 / 0.2e1) * t112 + (t95 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,5) * t109 + Ifges(5,1) * t104 / 0.2e1) * t104 + (t85 * mrSges(6,1) - t86 * mrSges(6,2) + Ifges(6,6) * t99 + Ifges(6,3) * t102 / 0.2e1) * t102 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t134 / 0.2e1) * t134 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t134 + Ifges(3,1) * t130 / 0.2e1) * t130) * t142 + (-t117 * t130 + t118 * t134) * mrSges(3,3)) * t142 + (t117 * mrSges(3,1) - t118 * mrSges(3,2) + Ifges(3,3) * t122 / 0.2e1 + (Ifges(3,5) * t130 + Ifges(3,6) * t134) * t142) * t122 + (-t101 * mrSges(4,1) + t98 * mrSges(4,3) + Ifges(4,4) * t112 + Ifges(4,6) * t116 + Ifges(4,2) * t111 / 0.2e1) * t111 + (-t95 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t104 + Ifges(5,6) * t109 + Ifges(5,2) * t103 / 0.2e1) * t103 + (t87 * mrSges(6,2) - t85 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,5) * t102 + Ifges(6,1) * t100 / 0.2e1) * t100;
T = t1;
