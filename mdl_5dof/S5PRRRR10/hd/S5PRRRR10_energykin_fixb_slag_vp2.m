% Calculate kinetic energy for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:23:10
% EndTime: 2019-12-05 17:23:11
% DurationCPUTime: 0.37s
% Computational Cost: add. (368->77), mult. (874->134), div. (0->0), fcn. (665->12), ass. (0->41)
t125 = cos(qJ(2));
t115 = sin(pkin(5));
t132 = qJD(1) * t115;
t108 = qJD(2) * pkin(2) + t125 * t132;
t114 = sin(pkin(6));
t116 = cos(pkin(6));
t117 = cos(pkin(5));
t131 = qJD(1) * t117;
t134 = t108 * t116 + t114 * t131;
t121 = sin(qJ(2));
t130 = qJD(2) * t114;
t107 = pkin(8) * t130 + t121 * t132;
t120 = sin(qJ(3));
t124 = cos(qJ(3));
t94 = -t120 * t107 + t134 * t124;
t119 = sin(qJ(4));
t123 = cos(qJ(4));
t112 = qJD(2) * t116 + qJD(3);
t95 = t124 * t107 + t134 * t120;
t93 = pkin(9) * t112 + t95;
t111 = t116 * t131;
t97 = t111 + (-t108 + (-pkin(3) * t124 - pkin(9) * t120) * qJD(2)) * t114;
t89 = t119 * t97 + t123 * t93;
t128 = t120 * t130;
t88 = -t119 * t93 + t123 * t97;
t102 = t112 * t123 - t119 * t128;
t92 = -pkin(3) * t112 - t94;
t122 = cos(qJ(5));
t118 = sin(qJ(5));
t110 = -t124 * t130 + qJD(4);
t103 = t112 * t119 + t123 * t128;
t101 = qJD(5) - t102;
t100 = -t108 * t114 + t111;
t99 = t103 * t122 + t110 * t118;
t98 = -t103 * t118 + t110 * t122;
t90 = -pkin(4) * t102 - pkin(10) * t103 + t92;
t87 = pkin(10) * t110 + t89;
t86 = -pkin(4) * t110 - t88;
t85 = t118 * t90 + t122 * t87;
t84 = -t118 * t87 + t122 * t90;
t1 = m(6) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t88 ^ 2 + t89 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + (t86 * mrSges(6,2) - t84 * mrSges(6,3) + Ifges(6,1) * t99 / 0.2e1) * t99 + (t94 * mrSges(4,1) - t95 * mrSges(4,2) + Ifges(4,3) * t112 / 0.2e1) * t112 + (t88 * mrSges(5,1) - t89 * mrSges(5,2) + Ifges(5,3) * t110 / 0.2e1) * t110 + (-t86 * mrSges(6,1) + t85 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,2) * t98 / 0.2e1) * t98 + (t92 * mrSges(5,2) - t88 * mrSges(5,3) + Ifges(5,5) * t110 + Ifges(5,1) * t103 / 0.2e1) * t103 + (m(3) * (t117 ^ 2 + (t121 ^ 2 + t125 ^ 2) * t115 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (-t92 * mrSges(5,1) + t89 * mrSges(5,3) + Ifges(5,4) * t103 + Ifges(5,6) * t110 + Ifges(5,2) * t102 / 0.2e1) * t102 + (t84 * mrSges(6,1) - t85 * mrSges(6,2) + Ifges(6,5) * t99 + Ifges(6,6) * t98 + Ifges(6,3) * t101 / 0.2e1) * t101 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t125 - mrSges(3,2) * t121) * t132 + (t100 * (-mrSges(4,1) * t124 + mrSges(4,2) * t120) + (Ifges(4,2) * t124 ^ 2 / 0.2e1 + (Ifges(4,4) * t124 + Ifges(4,1) * t120 / 0.2e1) * t120) * t130 + (-t94 * t120 + t95 * t124) * mrSges(4,3) + t112 * (Ifges(4,5) * t120 + Ifges(4,6) * t124)) * t114) * qJD(2);
T = t1;
