% Calculate kinetic energy for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:04
% EndTime: 2019-03-09 09:17:05
% DurationCPUTime: 0.44s
% Computational Cost: add. (492->110), mult. (1088->154), div. (0->0), fcn. (720->8), ass. (0->43)
t127 = -pkin(2) - pkin(3);
t115 = sin(qJ(5));
t118 = cos(qJ(5));
t113 = cos(pkin(6));
t125 = qJD(1) * t113;
t110 = qJD(2) + t125;
t116 = sin(qJ(2));
t112 = sin(pkin(6));
t126 = qJD(1) * t112;
t123 = t116 * t126;
t105 = pkin(8) * t123;
t119 = cos(qJ(2));
t121 = qJD(3) + t105 + (-pkin(1) * t113 * t119 - qJ(4) * t112 * t116) * qJD(1);
t86 = (-pkin(9) + t127) * t110 + t121;
t122 = t119 * t126;
t96 = -pkin(1) * t126 - pkin(2) * t122 - qJ(3) * t123;
t93 = pkin(3) * t122 + qJD(4) - t96;
t87 = (pkin(4) * t116 + pkin(9) * t119) * t126 + t93;
t82 = t115 * t87 + t118 * t86;
t124 = pkin(1) * t125;
t101 = pkin(8) * t122 + t116 * t124;
t95 = t110 * qJ(3) + t101;
t81 = -t115 * t86 + t118 * t87;
t100 = t119 * t124 - t105;
t98 = -t110 * t118 + t115 * t122;
t92 = qJ(4) * t122 - t95;
t89 = pkin(4) * t110 - t92;
t120 = qJD(1) ^ 2;
t117 = cos(qJ(6));
t114 = sin(qJ(6));
t102 = qJD(5) + t123;
t99 = -t110 * t115 - t118 * t122;
t97 = qJD(6) - t98;
t94 = -pkin(2) * t110 + qJD(3) - t100;
t91 = t102 * t114 + t117 * t99;
t90 = t102 * t117 - t114 * t99;
t88 = t127 * t110 + t121;
t83 = -pkin(5) * t98 - pkin(10) * t99 + t89;
t80 = pkin(10) * t102 + t82;
t79 = -pkin(5) * t102 - t81;
t78 = t114 * t83 + t117 * t80;
t77 = -t114 * t80 + t117 * t83;
t1 = m(5) * (t88 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + Ifges(2,3) * t120 / 0.2e1 + m(3) * (pkin(1) ^ 2 * t112 ^ 2 * t120 + t100 ^ 2 + t101 ^ 2) / 0.2e1 + (t89 * mrSges(6,2) - t81 * mrSges(6,3) + Ifges(6,1) * t99 / 0.2e1) * t99 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,3) * t97 / 0.2e1) * t97 + (-t89 * mrSges(6,1) + t82 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,2) * t98 / 0.2e1) * t98 + (t79 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,5) * t97 + Ifges(7,1) * t91 / 0.2e1) * t91 + (-t79 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t91 + Ifges(7,6) * t97 + Ifges(7,2) * t90 / 0.2e1) * t90 + (t81 * mrSges(6,1) - t82 * mrSges(6,2) + Ifges(6,5) * t99 + Ifges(6,6) * t98 + Ifges(6,3) * t102 / 0.2e1) * t102 + ((-t96 * mrSges(4,1) + t95 * mrSges(4,2) - t93 * mrSges(5,2) + t101 * mrSges(3,3) + t92 * mrSges(5,3) + (pkin(1) * mrSges(3,1) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t119) * t126) * t119 + (t93 * mrSges(5,1) + t94 * mrSges(4,2) - t100 * mrSges(3,3) - t96 * mrSges(4,3) - t88 * mrSges(5,3) + (-pkin(1) * mrSges(3,2) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t116 + (Ifges(3,4) + Ifges(5,4) - Ifges(4,5)) * t119) * t126) * t116) * t126 + (t100 * mrSges(3,1) - t94 * mrSges(4,1) - t92 * mrSges(5,1) - t101 * mrSges(3,2) + t88 * mrSges(5,2) + t95 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t110 + ((Ifges(5,5) + Ifges(3,6) - Ifges(4,6)) * t119 + (Ifges(4,4) + Ifges(3,5) + Ifges(5,6)) * t116) * t126) * t110;
T  = t1;
