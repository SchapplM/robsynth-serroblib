% Calculate kinetic energy for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:09
% EndTime: 2019-03-09 05:00:10
% DurationCPUTime: 0.54s
% Computational Cost: add. (606->97), mult. (1269->154), div. (0->0), fcn. (850->10), ass. (0->41)
t125 = m(3) / 0.2e1;
t111 = sin(pkin(11));
t113 = cos(pkin(11));
t116 = sin(qJ(4));
t119 = cos(qJ(4));
t117 = sin(qJ(3));
t124 = qJD(1) * t117;
t105 = qJD(3) * t116 + t119 * t124;
t120 = cos(qJ(3));
t109 = -t120 * qJD(1) + qJD(4);
t112 = sin(pkin(10));
t106 = (pkin(1) * t112 + pkin(7)) * qJD(1);
t101 = t117 * qJD(2) + t120 * t106;
t98 = qJD(3) * pkin(8) + t101;
t114 = cos(pkin(10));
t123 = -pkin(1) * t114 - pkin(2);
t99 = (-pkin(3) * t120 - pkin(8) * t117 + t123) * qJD(1);
t89 = -t116 * t98 + t119 * t99;
t85 = pkin(4) * t109 - qJ(5) * t105 + t89;
t104 = qJD(3) * t119 - t116 * t124;
t90 = t116 * t99 + t119 * t98;
t88 = qJ(5) * t104 + t90;
t80 = t111 * t85 + t113 * t88;
t79 = -t111 * t88 + t113 * t85;
t100 = qJD(2) * t120 - t117 * t106;
t97 = -qJD(3) * pkin(3) - t100;
t91 = -pkin(4) * t104 + qJD(5) + t97;
t118 = cos(qJ(6));
t115 = sin(qJ(6));
t108 = qJD(6) + t109;
t107 = t123 * qJD(1);
t93 = t104 * t111 + t105 * t113;
t92 = t104 * t113 - t105 * t111;
t86 = -pkin(5) * t92 + t91;
t84 = t115 * t92 + t118 * t93;
t83 = -t115 * t93 + t118 * t92;
t78 = pkin(9) * t92 + t80;
t77 = pkin(5) * t109 - pkin(9) * t93 + t79;
t76 = t115 * t77 + t118 * t78;
t75 = -t115 * t78 + t118 * t77;
t1 = qJD(2) ^ 2 * t125 + m(4) * (t100 ^ 2 + t101 ^ 2 + t107 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t97 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t91 ^ 2) / 0.2e1 + m(7) * (t75 ^ 2 + t76 ^ 2 + t86 ^ 2) / 0.2e1 + (t91 * mrSges(6,2) - t79 * mrSges(6,3) + Ifges(6,1) * t93 / 0.2e1) * t93 + (t86 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,1) * t84 / 0.2e1) * t84 + (t97 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,1) * t105 / 0.2e1) * t105 + (-t91 * mrSges(6,1) + t80 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,2) * t92 / 0.2e1) * t92 + (-t86 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t84 + Ifges(7,2) * t83 / 0.2e1) * t83 + (-t97 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,2) * t104 / 0.2e1) * t104 + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + Ifges(7,5) * t84 + Ifges(7,6) * t83 + Ifges(7,3) * t108 / 0.2e1) * t108 + (t100 * mrSges(4,1) - t101 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t89 * mrSges(5,1) + t79 * mrSges(6,1) - t90 * mrSges(5,2) - t80 * mrSges(6,2) + Ifges(5,5) * t105 + Ifges(6,5) * t93 + Ifges(5,6) * t104 + Ifges(6,6) * t92 + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t109) * t109 + (t107 * (-mrSges(4,1) * t120 + mrSges(4,2) * t117) + (-t100 * t117 + t101 * t120) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t117 + Ifges(4,6) * t120) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t114 * mrSges(3,1) - t112 * mrSges(3,2) + (t112 ^ 2 + t114 ^ 2) * t125 * pkin(1)) * pkin(1) + Ifges(4,2) * t120 ^ 2 / 0.2e1 + (Ifges(4,4) * t120 + Ifges(4,1) * t117 / 0.2e1) * t117) * qJD(1)) * qJD(1);
T  = t1;
