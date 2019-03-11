% Calculate kinetic energy for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:34
% EndTime: 2019-03-08 19:46:34
% DurationCPUTime: 0.32s
% Computational Cost: add. (304->83), mult. (624->129), div. (0->0), fcn. (392->10), ass. (0->40)
t115 = sin(qJ(2));
t110 = sin(pkin(6));
t126 = qJD(1) * t110;
t122 = t115 * t126;
t102 = qJD(2) * qJ(3) + t122;
t127 = t102 ^ 2;
t109 = sin(pkin(11));
t111 = cos(pkin(11));
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t112 = cos(pkin(6));
t125 = qJD(1) * t112;
t118 = cos(qJ(2));
t121 = -t118 * t126 + qJD(3);
t98 = (-pkin(2) - pkin(8)) * qJD(2) + t121;
t93 = t114 * t98 + t117 * t125;
t91 = qJD(4) * qJ(5) + t93;
t96 = t122 + (pkin(4) * t114 - qJ(5) * t117 + qJ(3)) * qJD(2);
t85 = t109 * t96 + t111 * t91;
t124 = qJD(2) * t114;
t123 = qJD(2) * t117;
t84 = -t109 * t91 + t111 * t96;
t92 = -t114 * t125 + t117 * t98;
t90 = -qJD(4) * pkin(4) + qJD(5) - t92;
t120 = qJD(1) ^ 2;
t116 = cos(qJ(6));
t113 = sin(qJ(6));
t107 = t112 ^ 2 * t120;
t106 = qJD(6) + t124;
t101 = -qJD(2) * pkin(2) + t121;
t100 = qJD(4) * t109 + t111 * t123;
t99 = qJD(4) * t111 - t109 * t123;
t89 = t100 * t116 + t113 * t99;
t88 = -t100 * t113 + t116 * t99;
t86 = -pkin(5) * t99 + t90;
t83 = pkin(9) * t99 + t85;
t82 = pkin(5) * t124 - pkin(9) * t100 + t84;
t81 = t113 * t82 + t116 * t83;
t80 = -t113 * t83 + t116 * t82;
t1 = m(2) * t120 / 0.2e1 + m(3) * (t107 + (t115 ^ 2 + t118 ^ 2) * t120 * t110 ^ 2) / 0.2e1 + m(5) * (t92 ^ 2 + t93 ^ 2 + t127) / 0.2e1 + m(4) * (t101 ^ 2 + t107 + t127) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t90 ^ 2) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t86 ^ 2) / 0.2e1 + (-t90 * mrSges(6,1) + t85 * mrSges(6,3) + Ifges(6,2) * t99 / 0.2e1) * t99 + (t86 * mrSges(7,2) - t80 * mrSges(7,3) + Ifges(7,1) * t89 / 0.2e1) * t89 + (t92 * mrSges(5,1) - t93 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (-t86 * mrSges(7,1) + t81 * mrSges(7,3) + Ifges(7,4) * t89 + Ifges(7,2) * t88 / 0.2e1) * t88 + (t90 * mrSges(6,2) - t84 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,1) * t100 / 0.2e1) * t100 + (t80 * mrSges(7,1) - t81 * mrSges(7,2) + Ifges(7,5) * t89 + Ifges(7,6) * t88 + Ifges(7,3) * t106 / 0.2e1) * t106 + (t101 * mrSges(4,2) + t102 * mrSges(4,3) + (mrSges(3,1) * t118 - mrSges(3,2) * t115) * t126 + (t102 * mrSges(5,2) - t92 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t123 / 0.2e1) * t117 + (t102 * mrSges(5,1) + t84 * mrSges(6,1) - t85 * mrSges(6,2) - t93 * mrSges(5,3) + Ifges(6,5) * t100 - Ifges(5,6) * qJD(4) + Ifges(6,6) * t99) * t114 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t117 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t114) * t114) * qJD(2)) * qJD(2);
T  = t1;
