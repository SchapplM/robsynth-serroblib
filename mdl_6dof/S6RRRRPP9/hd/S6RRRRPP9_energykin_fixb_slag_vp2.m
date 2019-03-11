% Calculate kinetic energy for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:41:52
% EndTime: 2019-03-09 21:41:52
% DurationCPUTime: 0.42s
% Computational Cost: add. (738->107), mult. (1624->149), div. (0->0), fcn. (1224->8), ass. (0->42)
t126 = cos(qJ(4));
t125 = pkin(4) + qJ(6);
t112 = sin(qJ(4));
t124 = cos(pkin(6)) * qJD(1);
t109 = qJD(2) + t124;
t113 = sin(qJ(3));
t115 = cos(qJ(3));
t114 = sin(qJ(2));
t110 = sin(pkin(6));
t123 = t110 * qJD(1);
t121 = t114 * t123;
t100 = t109 * t115 - t113 * t121;
t101 = t109 * t113 + t115 * t121;
t116 = cos(qJ(2));
t122 = pkin(1) * t124;
t102 = -pkin(8) * t121 + t116 * t122;
t96 = -pkin(2) * t109 - t102;
t84 = -pkin(3) * t100 - pkin(10) * t101 + t96;
t120 = t116 * t123;
t105 = qJD(3) - t120;
t103 = pkin(8) * t120 + t114 * t122;
t97 = pkin(9) * t109 + t103;
t98 = (-pkin(2) * t116 - pkin(9) * t114 - pkin(1)) * t123;
t90 = t113 * t98 + t115 * t97;
t88 = pkin(10) * t105 + t90;
t82 = t112 * t84 + t126 * t88;
t89 = -t113 * t97 + t115 * t98;
t99 = qJD(4) - t100;
t79 = -qJ(5) * t99 - t82;
t81 = -t112 * t88 + t126 * t84;
t119 = qJD(5) - t81;
t87 = -pkin(3) * t105 - t89;
t92 = t101 * t126 + t112 * t105;
t118 = -qJ(5) * t92 + t87;
t117 = qJD(1) ^ 2;
t91 = t101 * t112 - t105 * t126;
t80 = pkin(4) * t91 + t118;
t78 = -t99 * pkin(4) + t119;
t77 = t125 * t91 + t118;
t76 = -pkin(5) * t91 + qJD(6) - t79;
t75 = t92 * pkin(5) - t125 * t99 + t119;
t1 = m(3) * (pkin(1) ^ 2 * t110 ^ 2 * t117 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + t117 * Ifges(2,3) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t87 ^ 2) / 0.2e1 + m(7) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + (t89 * mrSges(4,1) - t90 * mrSges(4,2) + Ifges(4,3) * t105 / 0.2e1) * t105 + (t96 * mrSges(4,2) - t89 * mrSges(4,3) + Ifges(4,5) * t105 + Ifges(4,1) * t101 / 0.2e1) * t101 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t116 / 0.2e1) * t116 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t116 + Ifges(3,1) * t114 / 0.2e1) * t114) * t123 + (-t102 * t114 + t103 * t116) * mrSges(3,3)) * t123 + (t102 * mrSges(3,1) - t103 * mrSges(3,2) + Ifges(3,3) * t109 / 0.2e1 + (Ifges(3,5) * t114 + Ifges(3,6) * t116) * t123) * t109 + (-t96 * mrSges(4,1) + t90 * mrSges(4,3) + Ifges(4,4) * t101 + Ifges(4,6) * t105 + Ifges(4,2) * t100 / 0.2e1) * t100 + (t81 * mrSges(5,1) - t82 * mrSges(5,2) + t78 * mrSges(6,2) + t76 * mrSges(7,2) - t79 * mrSges(6,3) - t75 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t99) * t99 + (t78 * mrSges(6,1) + t75 * mrSges(7,1) + t87 * mrSges(5,2) - t77 * mrSges(7,2) - t81 * mrSges(5,3) - t80 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t92 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t99) * t92 + (t87 * mrSges(5,1) + t79 * mrSges(6,1) - t76 * mrSges(7,1) - t80 * mrSges(6,2) - t82 * mrSges(5,3) + t77 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t91 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t99 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t92) * t91;
T  = t1;
