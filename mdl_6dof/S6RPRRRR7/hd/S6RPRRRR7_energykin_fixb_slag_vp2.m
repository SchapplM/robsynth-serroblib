% Calculate kinetic energy for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:53
% EndTime: 2019-03-09 07:16:53
% DurationCPUTime: 0.48s
% Computational Cost: add. (585->94), mult. (1134->149), div. (0->0), fcn. (752->8), ass. (0->41)
t98 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t116 = t98 * mrSges(4,3);
t115 = qJD(1) / 0.2e1;
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t102 = qJD(3) + qJD(4);
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t112 = cos(qJ(3));
t114 = -pkin(8) * qJD(1) + t98;
t92 = qJD(3) * pkin(3) + t114 * t112;
t108 = sin(qJ(3));
t93 = t114 * t108;
t83 = -t107 * t93 + t111 * t92;
t95 = (-t107 * t108 + t111 * t112) * qJD(1);
t78 = pkin(4) * t102 - pkin(9) * t95 + t83;
t84 = t107 * t92 + t111 * t93;
t94 = (-t107 * t112 - t108 * t111) * qJD(1);
t80 = pkin(9) * t94 + t84;
t75 = t106 * t78 + t110 * t80;
t104 = qJD(1) * qJ(2);
t96 = t108 * qJD(1) * pkin(3) + t104;
t88 = -pkin(4) * t94 + t96;
t74 = -t106 * t80 + t110 * t78;
t86 = -t106 * t95 + t110 * t94;
t113 = qJD(1) ^ 2;
t109 = cos(qJ(6));
t105 = sin(qJ(6));
t103 = t113 * qJ(2) ^ 2;
t101 = qJD(5) + t102;
t100 = -qJD(1) * pkin(1) + qJD(2);
t87 = t106 * t94 + t110 * t95;
t85 = qJD(6) - t86;
t82 = t101 * t105 + t109 * t87;
t81 = t101 * t109 - t105 * t87;
t76 = -pkin(5) * t86 - pkin(10) * t87 + t88;
t73 = pkin(10) * t101 + t75;
t72 = -pkin(5) * t101 - t74;
t71 = t105 * t76 + t109 * t73;
t70 = -t105 * t73 + t109 * t76;
t1 = m(4) * (t103 + (t108 ^ 2 + t112 ^ 2) * t98 ^ 2) / 0.2e1 + m(3) * (t100 ^ 2 + t103) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + (t96 * mrSges(5,2) - t83 * mrSges(5,3) + Ifges(5,1) * t95 / 0.2e1) * t95 + (t88 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,1) * t87 / 0.2e1) * t87 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,3) * t85 / 0.2e1) * t85 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t113 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t112 * mrSges(4,1) - t108 * mrSges(4,2)) * t98) * qJD(3) + (t100 * mrSges(3,2) + (mrSges(4,2) * t104 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t115 - t116) * t112) * t112 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t112) * qJD(1) + (Ifges(4,2) * t115 - t116) * t108) * t108) * qJD(1) + (-t96 * mrSges(5,1) + t84 * mrSges(5,3) + Ifges(5,4) * t95 + Ifges(5,2) * t94 / 0.2e1) * t94 + (-t88 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t87 + Ifges(6,2) * t86 / 0.2e1) * t86 + (t72 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,5) * t85 + Ifges(7,1) * t82 / 0.2e1) * t82 + (-t72 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t82 + Ifges(7,6) * t85 + Ifges(7,2) * t81 / 0.2e1) * t81 + (t83 * mrSges(5,1) - t84 * mrSges(5,2) + Ifges(5,5) * t95 + Ifges(5,6) * t94 + Ifges(5,3) * t102 / 0.2e1) * t102 + (t74 * mrSges(6,1) - t75 * mrSges(6,2) + Ifges(6,5) * t87 + Ifges(6,6) * t86 + Ifges(6,3) * t101 / 0.2e1) * t101;
T  = t1;
