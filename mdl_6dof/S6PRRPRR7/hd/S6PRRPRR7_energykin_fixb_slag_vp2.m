% Calculate kinetic energy for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:09
% EndTime: 2019-03-08 22:30:10
% DurationCPUTime: 0.36s
% Computational Cost: add. (354->96), mult. (729->142), div. (0->0), fcn. (451->10), ass. (0->42)
t132 = -pkin(3) - pkin(9);
t118 = sin(qJ(5));
t122 = cos(qJ(5));
t119 = sin(qJ(3));
t112 = t119 * qJD(2);
t120 = sin(qJ(2));
t115 = sin(pkin(6));
t131 = qJD(1) * t115;
t106 = qJD(2) * pkin(8) + t120 * t131;
t123 = cos(qJ(3));
t116 = cos(pkin(6));
t130 = qJD(1) * t116;
t99 = -t106 * t119 + t123 * t130;
t126 = qJD(4) - t99;
t91 = pkin(4) * t112 + qJD(3) * t132 + t126;
t127 = -qJ(4) * t119 - pkin(2);
t124 = cos(qJ(2));
t128 = t124 * t131;
t96 = -t128 + (t123 * t132 + t127) * qJD(2);
t87 = t118 * t91 + t122 * t96;
t100 = t106 * t123 + t119 * t130;
t129 = qJD(2) * t123;
t110 = t112 + qJD(5);
t98 = -qJD(3) * qJ(4) - t100;
t86 = -t118 * t96 + t122 * t91;
t94 = pkin(4) * t129 - t98;
t121 = cos(qJ(6));
t117 = sin(qJ(6));
t108 = qJD(6) + t110;
t107 = -qJD(2) * pkin(2) - t128;
t105 = qJD(3) * t122 - t118 * t129;
t104 = -qJD(3) * t118 - t122 * t129;
t101 = -t128 + (-pkin(3) * t123 + t127) * qJD(2);
t97 = -qJD(3) * pkin(3) + t126;
t93 = t104 * t117 + t105 * t121;
t92 = t104 * t121 - t105 * t117;
t88 = -pkin(5) * t104 + t94;
t85 = pkin(10) * t104 + t87;
t84 = pkin(5) * t110 - pkin(10) * t105 + t86;
t83 = t117 * t84 + t121 * t85;
t82 = -t117 * t85 + t121 * t84;
t1 = m(4) * (t100 ^ 2 + t107 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t86 ^ 2 + t87 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t101 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t82 ^ 2 + t83 ^ 2 + t88 ^ 2) / 0.2e1 + (t88 * mrSges(7,2) - t82 * mrSges(7,3) + Ifges(7,1) * t93 / 0.2e1) * t93 + (t86 * mrSges(6,1) - t87 * mrSges(6,2) + Ifges(6,3) * t110 / 0.2e1) * t110 + (-t88 * mrSges(7,1) + t83 * mrSges(7,3) + Ifges(7,4) * t93 + Ifges(7,2) * t92 / 0.2e1) * t92 + (t94 * mrSges(6,2) - t86 * mrSges(6,3) + Ifges(6,5) * t110 + Ifges(6,1) * t105 / 0.2e1) * t105 + (m(3) * (t116 ^ 2 + (t120 ^ 2 + t124 ^ 2) * t115 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t82 * mrSges(7,1) - t83 * mrSges(7,2) + Ifges(7,5) * t93 + Ifges(7,6) * t92 + Ifges(7,3) * t108 / 0.2e1) * t108 + (-t94 * mrSges(6,1) + t87 * mrSges(6,3) + Ifges(6,4) * t105 + Ifges(6,6) * t110 + Ifges(6,2) * t104 / 0.2e1) * t104 + (t99 * mrSges(4,1) - t100 * mrSges(4,2) + t97 * mrSges(5,2) - t98 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3)) * qJD(3) + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t124 - mrSges(3,2) * t120) * t131 + (-t107 * mrSges(4,1) - t98 * mrSges(5,1) + t101 * mrSges(5,2) + t100 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t129 + (-Ifges(5,5) + Ifges(4,6)) * qJD(3)) * t123 + (t97 * mrSges(5,1) + t107 * mrSges(4,2) - t99 * mrSges(4,3) - t101 * mrSges(5,3) + (-Ifges(5,4) + Ifges(4,5)) * qJD(3) + ((Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t119 + (Ifges(4,4) + Ifges(5,6)) * t123) * qJD(2)) * t119) * qJD(2);
T  = t1;
