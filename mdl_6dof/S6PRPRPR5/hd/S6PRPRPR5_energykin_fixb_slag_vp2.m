% Calculate kinetic energy for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:42:52
% EndTime: 2019-03-08 19:42:53
% DurationCPUTime: 0.34s
% Computational Cost: add. (333->90), mult. (792->136), div. (0->0), fcn. (551->10), ass. (0->41)
t130 = pkin(4) + pkin(9);
t129 = cos(qJ(4));
t118 = sin(qJ(4));
t119 = sin(qJ(2));
t114 = sin(pkin(6));
t128 = qJD(1) * t114;
t108 = qJD(2) * qJ(3) + t119 * t128;
t115 = cos(pkin(11));
t116 = cos(pkin(6));
t127 = qJD(1) * t116;
t110 = t115 * t127;
t113 = sin(pkin(11));
t96 = t110 + (-pkin(8) * qJD(2) - t108) * t113;
t101 = t115 * t108 + t113 * t127;
t126 = qJD(2) * t115;
t97 = pkin(8) * t126 + t101;
t91 = t118 * t96 + t129 * t97;
t90 = -t118 * t97 + t129 * t96;
t89 = -qJD(4) * qJ(5) - t91;
t121 = cos(qJ(2));
t125 = -t121 * t128 + qJD(3);
t124 = qJD(5) - t90;
t105 = (t113 * t129 + t115 * t118) * qJD(2);
t102 = (-pkin(3) * t115 - pkin(2)) * qJD(2) + t125;
t123 = -qJ(5) * t105 + t102;
t120 = cos(qJ(6));
t117 = sin(qJ(6));
t107 = -qJD(2) * pkin(2) + t125;
t104 = qJD(2) * t113 * t118 - t126 * t129;
t103 = qJD(6) + t105;
t100 = -t108 * t113 + t110;
t99 = qJD(4) * t120 + t104 * t117;
t98 = -qJD(4) * t117 + t104 * t120;
t92 = pkin(4) * t104 + t123;
t88 = -qJD(4) * pkin(4) + t124;
t87 = t104 * t130 + t123;
t86 = -pkin(5) * t104 - t89;
t85 = t105 * pkin(5) - qJD(4) * t130 + t124;
t84 = t117 * t85 + t120 * t87;
t83 = -t117 * t87 + t120 * t85;
t1 = m(4) * (t100 ^ 2 + t101 ^ 2 + t107 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(7) * (t83 ^ 2 + t84 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t88 ^ 2 + t89 ^ 2 + t92 ^ 2) / 0.2e1 + (t86 * mrSges(7,2) - t83 * mrSges(7,3) + Ifges(7,1) * t99 / 0.2e1) * t99 + (-t86 * mrSges(7,1) + t84 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,2) * t98 / 0.2e1) * t98 + (m(2) / 0.2e1 + m(3) * (t116 ^ 2 + (t119 ^ 2 + t121 ^ 2) * t114 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t83 * mrSges(7,1) - t84 * mrSges(7,2) + Ifges(7,5) * t99 + Ifges(7,6) * t98 + Ifges(7,3) * t103 / 0.2e1) * t103 + (t88 * mrSges(6,1) + t102 * mrSges(5,2) - t90 * mrSges(5,3) - t92 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t105) * t105 + (t102 * mrSges(5,1) + t89 * mrSges(6,1) - t92 * mrSges(6,2) - t91 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t104 + (-Ifges(5,4) - Ifges(6,6)) * t105) * t104 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + t88 * mrSges(6,2) - t89 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(4) + (-Ifges(6,4) + Ifges(5,5)) * t105 + (Ifges(6,5) - Ifges(5,6)) * t104) * qJD(4) + (t107 * (-mrSges(4,1) * t115 + mrSges(4,2) * t113) + (Ifges(4,2) * t115 ^ 2 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t115 + Ifges(4,1) * t113 / 0.2e1) * t113) * qJD(2) + (mrSges(3,1) * t121 - mrSges(3,2) * t119) * t128 + (-t100 * t113 + t101 * t115) * mrSges(4,3)) * qJD(2);
T  = t1;
