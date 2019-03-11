% Calculate kinetic energy for
% S6RPRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:42:33
% EndTime: 2019-03-09 04:42:33
% DurationCPUTime: 0.48s
% Computational Cost: add. (455->97), mult. (1071->130), div. (0->0), fcn. (738->6), ass. (0->37)
t102 = cos(pkin(9));
t118 = t102 ^ 2;
t117 = m(3) / 0.2e1;
t116 = -pkin(4) - pkin(5);
t115 = cos(qJ(3));
t114 = cos(qJ(4));
t113 = pkin(7) + qJ(2);
t103 = sin(qJ(4));
t104 = sin(qJ(3));
t110 = t102 * qJD(1);
t101 = sin(pkin(9));
t111 = t101 * qJD(1);
t91 = -t104 * t111 + t115 * t110;
t92 = (t115 * t101 + t102 * t104) * qJD(1);
t95 = qJD(2) + (-pkin(2) * t102 - pkin(1)) * qJD(1);
t76 = -pkin(3) * t91 - pkin(8) * t92 + t95;
t93 = t113 * t111;
t94 = t113 * t110;
t82 = -t104 * t93 + t115 * t94;
t80 = qJD(3) * pkin(8) + t82;
t73 = t103 * t76 + t114 * t80;
t81 = -t104 * t94 - t115 * t93;
t86 = qJD(4) - t91;
t71 = t86 * qJ(5) + t73;
t109 = qJD(3) * pkin(3) + t81;
t72 = -t103 * t80 + t114 * t76;
t108 = qJD(5) - t72;
t84 = t103 * qJD(3) + t114 * t92;
t107 = qJ(5) * t84 + t109;
t97 = -qJD(1) * pkin(1) + qJD(2);
t83 = -t114 * qJD(3) + t103 * t92;
t74 = pkin(4) * t83 - t107;
t70 = -t86 * pkin(4) + t108;
t69 = t116 * t83 + qJD(6) + t107;
t68 = qJ(6) * t83 + t71;
t67 = -t84 * qJ(6) + t116 * t86 + t108;
t1 = t97 ^ 2 * t117 + m(4) * (t81 ^ 2 + t82 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t74 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + (t95 * mrSges(4,2) - t81 * mrSges(4,3) + Ifges(4,1) * t92 / 0.2e1) * t92 + (-t95 * mrSges(4,1) + t82 * mrSges(4,3) + Ifges(4,4) * t92 + Ifges(4,2) * t91 / 0.2e1) * t91 + (t81 * mrSges(4,1) - t82 * mrSges(4,2) + Ifges(4,5) * t92 + Ifges(4,6) * t91 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t72 * mrSges(5,1) - t70 * mrSges(6,1) - t67 * mrSges(7,1) - t73 * mrSges(5,2) + t68 * mrSges(7,2) + t71 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t86) * t86 + (-t109 * mrSges(5,2) + t70 * mrSges(6,2) + t69 * mrSges(7,2) - t72 * mrSges(5,3) - t74 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t84 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t86) * t84 + (-t109 * mrSges(5,1) + t74 * mrSges(6,1) - t69 * mrSges(7,1) - t71 * mrSges(6,2) - t73 * mrSges(5,3) + t68 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t83 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t86 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t84) * t83 + (t97 * (-mrSges(3,1) * t102 + mrSges(3,2) * t101) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t117 + mrSges(3,3)) * (t101 ^ 2 + t118) * qJ(2) + Ifges(3,2) * t118 / 0.2e1 + (Ifges(3,4) * t102 + Ifges(3,1) * t101 / 0.2e1) * t101) * qJD(1)) * qJD(1);
T  = t1;
