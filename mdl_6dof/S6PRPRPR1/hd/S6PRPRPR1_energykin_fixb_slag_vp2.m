% Calculate kinetic energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:34
% EndTime: 2019-03-08 19:24:34
% DurationCPUTime: 0.39s
% Computational Cost: add. (345->85), mult. (777->137), div. (0->0), fcn. (540->12), ass. (0->40)
t115 = sin(pkin(12));
t118 = cos(pkin(12));
t122 = sin(qJ(4));
t125 = cos(qJ(4));
t105 = (t115 * t122 - t118 * t125) * qJD(2);
t120 = cos(pkin(6));
t113 = t120 * qJD(1) + qJD(3);
t112 = t125 * t113;
t132 = qJ(5) * qJD(2);
t126 = cos(qJ(2));
t117 = sin(pkin(6));
t131 = qJD(1) * t117;
t108 = qJD(2) * pkin(2) + t126 * t131;
t116 = sin(pkin(11));
t119 = cos(pkin(11));
t123 = sin(qJ(2));
t130 = t123 * t131;
t101 = t116 * t108 + t119 * t130;
t99 = qJD(2) * pkin(8) + t101;
t92 = qJD(4) * pkin(4) + t112 + (-t99 - t132) * t122;
t95 = t122 * t113 + t125 * t99;
t93 = t125 * t132 + t95;
t88 = t115 * t92 + t118 * t93;
t100 = t119 * t108 - t116 * t130;
t87 = -t115 * t93 + t118 * t92;
t96 = qJD(5) + (-pkin(4) * t125 - pkin(3)) * qJD(2) - t100;
t124 = cos(qJ(6));
t121 = sin(qJ(6));
t106 = (t115 * t125 + t118 * t122) * qJD(2);
t104 = qJD(6) + t105;
t103 = t121 * qJD(4) + t124 * t106;
t102 = t124 * qJD(4) - t121 * t106;
t98 = -qJD(2) * pkin(3) - t100;
t94 = -t122 * t99 + t112;
t89 = t105 * pkin(5) - t106 * pkin(9) + t96;
t86 = qJD(4) * pkin(9) + t88;
t85 = -qJD(4) * pkin(5) - t87;
t84 = t121 * t89 + t124 * t86;
t83 = -t121 * t86 + t124 * t89;
t1 = m(7) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t87 ^ 2 + t88 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t94 ^ 2 + t95 ^ 2 + t98 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t113 ^ 2) / 0.2e1 + (t96 * mrSges(6,2) - t87 * mrSges(6,3) + Ifges(6,1) * t106 / 0.2e1) * t106 + (t83 * mrSges(7,1) - t84 * mrSges(7,2) + Ifges(7,3) * t104 / 0.2e1) * t104 - (-t96 * mrSges(6,1) + t88 * mrSges(6,3) + Ifges(6,4) * t106 - Ifges(6,2) * t105 / 0.2e1) * t105 + (t85 * mrSges(7,2) - t83 * mrSges(7,3) + Ifges(7,5) * t104 + Ifges(7,1) * t103 / 0.2e1) * t103 + (m(2) / 0.2e1 + m(3) * (t120 ^ 2 + (t123 ^ 2 + t126 ^ 2) * t117 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (-t85 * mrSges(7,1) + t84 * mrSges(7,3) + Ifges(7,4) * t103 + Ifges(7,6) * t104 + Ifges(7,2) * t102 / 0.2e1) * t102 + (t94 * mrSges(5,1) + t87 * mrSges(6,1) - t95 * mrSges(5,2) - t88 * mrSges(6,2) + Ifges(6,5) * t106 - Ifges(6,6) * t105 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4)) * qJD(4) + (t100 * mrSges(4,1) - t101 * mrSges(4,2) + t98 * (-mrSges(5,1) * t125 + mrSges(5,2) * t122) + (mrSges(3,1) * t126 - mrSges(3,2) * t123) * t131 + (-t94 * t122 + t95 * t125) * mrSges(5,3) + qJD(4) * (Ifges(5,5) * t122 + Ifges(5,6) * t125) + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,2) * t125 ^ 2 / 0.2e1 + (Ifges(5,4) * t125 + Ifges(5,1) * t122 / 0.2e1) * t122) * qJD(2)) * qJD(2);
T  = t1;
