% Calculate kinetic energy for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:07
% EndTime: 2019-03-09 14:29:08
% DurationCPUTime: 0.55s
% Computational Cost: add. (652->108), mult. (1291->158), div. (0->0), fcn. (822->8), ass. (0->42)
t123 = -pkin(2) - pkin(8);
t122 = pkin(7) * mrSges(3,3);
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t112 = sin(qJ(2));
t107 = t112 * qJD(1);
t103 = t107 + qJD(4);
t111 = sin(qJ(4));
t115 = cos(qJ(4));
t116 = cos(qJ(2));
t119 = -qJ(3) * t112 - pkin(1);
t93 = (t116 * t123 + t119) * qJD(1);
t120 = pkin(7) * t107 + qJD(3);
t94 = pkin(3) * t107 + qJD(2) * t123 + t120;
t85 = -t111 * t93 + t115 * t94;
t121 = qJD(1) * t116;
t98 = qJD(2) * t115 - t111 * t121;
t81 = pkin(4) * t103 - pkin(9) * t98 + t85;
t86 = t111 * t94 + t115 * t93;
t97 = -qJD(2) * t111 - t115 * t121;
t84 = pkin(9) * t97 + t86;
t76 = t110 * t81 + t114 * t84;
t101 = -pkin(7) * t121 - qJD(2) * qJ(3);
t95 = pkin(3) * t121 - t101;
t102 = qJD(5) + t103;
t75 = -t110 * t84 + t114 * t81;
t89 = -pkin(4) * t97 + t95;
t113 = cos(qJ(6));
t109 = sin(qJ(6));
t100 = qJD(6) + t102;
t99 = -qJD(2) * pkin(2) + t120;
t96 = (-pkin(2) * t116 + t119) * qJD(1);
t88 = t110 * t97 + t114 * t98;
t87 = -t110 * t98 + t114 * t97;
t83 = -pkin(5) * t87 + t89;
t80 = t109 * t87 + t113 * t88;
t79 = -t109 * t88 + t113 * t87;
t74 = pkin(10) * t87 + t76;
t73 = pkin(5) * t102 - pkin(10) * t88 + t75;
t72 = t109 * t73 + t113 * t74;
t71 = -t109 * t74 + t113 * t73;
t1 = m(7) * (t71 ^ 2 + t72 ^ 2 + t83 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t96 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t89 ^ 2) / 0.2e1 + (t95 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t98 / 0.2e1) * t98 + (t89 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,1) * t88 / 0.2e1) * t88 + (t83 * mrSges(7,2) - t71 * mrSges(7,3) + Ifges(7,1) * t80 / 0.2e1) * t80 + (-t95 * mrSges(5,1) + t86 * mrSges(5,3) + Ifges(5,4) * t98 + Ifges(5,2) * t97 / 0.2e1) * t97 + (-t89 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t88 + Ifges(6,2) * t87 / 0.2e1) * t87 + (-t83 * mrSges(7,1) + t72 * mrSges(7,3) + Ifges(7,4) * t80 + Ifges(7,2) * t79 / 0.2e1) * t79 + (t85 * mrSges(5,1) - t86 * mrSges(5,2) + Ifges(5,5) * t98 + Ifges(5,6) * t97 + Ifges(5,3) * t103 / 0.2e1) * t103 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,5) * t88 + Ifges(6,6) * t87 + Ifges(6,3) * t102 / 0.2e1) * t102 + (t71 * mrSges(7,1) - t72 * mrSges(7,2) + Ifges(7,5) * t80 + Ifges(7,6) * t79 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t99 * mrSges(4,2) - t101 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(2)) * qJD(2) + ((-t101 * mrSges(4,1) + t96 * mrSges(4,2) + (-pkin(7) * mrSges(3,2) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t116 + (t99 * mrSges(4,1) - t96 * mrSges(4,3) + (-pkin(7) * mrSges(3,1) - Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t112 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t112 ^ 2 + t116 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t122 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t116) * t116 + (-pkin(1) * mrSges(3,2) + (t122 + Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t112 + (Ifges(3,4) + Ifges(4,6)) * t116) * t112) * qJD(1)) * qJD(1);
T  = t1;
