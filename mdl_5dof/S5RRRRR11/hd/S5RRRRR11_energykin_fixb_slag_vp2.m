% Calculate kinetic energy for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:44
% EndTime: 2019-12-31 22:38:44
% DurationCPUTime: 0.45s
% Computational Cost: add. (658->90), mult. (1479->148), div. (0->0), fcn. (1132->10), ass. (0->43)
t116 = sin(qJ(4));
t120 = cos(qJ(4));
t127 = cos(pkin(5)) * qJD(1);
t112 = qJD(2) + t127;
t117 = sin(qJ(3));
t121 = cos(qJ(3));
t118 = sin(qJ(2));
t113 = sin(pkin(5));
t128 = qJD(1) * t113;
t125 = t118 * t128;
t103 = t112 * t121 - t117 * t125;
t104 = t112 * t117 + t121 * t125;
t122 = cos(qJ(2));
t126 = pkin(1) * t127;
t105 = -pkin(7) * t125 + t122 * t126;
t98 = -pkin(2) * t112 - t105;
t87 = -pkin(3) * t103 - pkin(9) * t104 + t98;
t124 = t122 * t128;
t108 = qJD(3) - t124;
t101 = (-pkin(2) * t122 - pkin(8) * t118 - pkin(1)) * t128;
t106 = pkin(7) * t124 + t118 * t126;
t99 = pkin(8) * t112 + t106;
t92 = t117 * t101 + t121 * t99;
t90 = pkin(9) * t108 + t92;
t81 = t116 * t87 + t120 * t90;
t80 = -t116 * t90 + t120 * t87;
t91 = t101 * t121 - t117 * t99;
t102 = qJD(4) - t103;
t89 = -pkin(3) * t108 - t91;
t123 = qJD(1) ^ 2;
t119 = cos(qJ(5));
t115 = sin(qJ(5));
t100 = qJD(5) + t102;
t94 = t104 * t120 + t108 * t116;
t93 = -t104 * t116 + t108 * t120;
t84 = t115 * t93 + t119 * t94;
t83 = -t115 * t94 + t119 * t93;
t82 = -pkin(4) * t93 + t89;
t79 = pkin(10) * t93 + t81;
t78 = pkin(4) * t102 - pkin(10) * t94 + t80;
t77 = t115 * t78 + t119 * t79;
t76 = -t115 * t79 + t119 * t78;
t1 = m(6) * (t76 ^ 2 + t77 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t98 ^ 2) / 0.2e1 + t123 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t113 ^ 2 * t123 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + (t89 * mrSges(5,2) - t80 * mrSges(5,3) + Ifges(5,1) * t94 / 0.2e1) * t94 + (t82 * mrSges(6,2) - t76 * mrSges(6,3) + Ifges(6,1) * t84 / 0.2e1) * t84 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + Ifges(4,3) * t108 / 0.2e1) * t108 + (-t89 * mrSges(5,1) + t81 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,2) * t93 / 0.2e1) * t93 + (-t82 * mrSges(6,1) + t77 * mrSges(6,3) + Ifges(6,4) * t84 + Ifges(6,2) * t83 / 0.2e1) * t83 + (t98 * mrSges(4,2) - t91 * mrSges(4,3) + Ifges(4,5) * t108 + Ifges(4,1) * t104 / 0.2e1) * t104 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t122 / 0.2e1) * t122 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t122 + Ifges(3,1) * t118 / 0.2e1) * t118) * t128 + (-t105 * t118 + t106 * t122) * mrSges(3,3)) * t128 + (t105 * mrSges(3,1) - t106 * mrSges(3,2) + Ifges(3,3) * t112 / 0.2e1 + (Ifges(3,5) * t118 + Ifges(3,6) * t122) * t128) * t112 + (-t98 * mrSges(4,1) + t92 * mrSges(4,3) + Ifges(4,4) * t104 + Ifges(4,6) * t108 + Ifges(4,2) * t103 / 0.2e1) * t103 + (t80 * mrSges(5,1) - t81 * mrSges(5,2) + Ifges(5,5) * t94 + Ifges(5,6) * t93 + Ifges(5,3) * t102 / 0.2e1) * t102 + (t76 * mrSges(6,1) - t77 * mrSges(6,2) + Ifges(6,5) * t84 + Ifges(6,6) * t83 + Ifges(6,3) * t100 / 0.2e1) * t100;
T = t1;
