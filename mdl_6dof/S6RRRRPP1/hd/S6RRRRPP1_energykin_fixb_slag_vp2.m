% Calculate kinetic energy for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:43:48
% EndTime: 2019-03-09 20:43:48
% DurationCPUTime: 0.56s
% Computational Cost: add. (804->106), mult. (1641->153), div. (0->0), fcn. (1186->8), ass. (0->37)
t117 = qJD(1) * (-pkin(8) - pkin(7));
t106 = sin(qJ(3));
t107 = sin(qJ(2));
t109 = cos(qJ(3));
t110 = cos(qJ(2));
t96 = (t106 * t107 - t109 * t110) * qJD(1);
t115 = pkin(7) * mrSges(3,3);
t104 = sin(pkin(10));
t114 = cos(pkin(10));
t105 = sin(qJ(4));
t108 = cos(qJ(4));
t102 = (-pkin(2) * t110 - pkin(1)) * qJD(1);
t97 = (t106 * t110 + t107 * t109) * qJD(1);
t85 = pkin(3) * t96 - pkin(9) * t97 + t102;
t103 = qJD(2) + qJD(3);
t100 = qJD(2) * pkin(2) + t107 * t117;
t101 = t110 * t117;
t90 = t106 * t100 - t109 * t101;
t88 = pkin(9) * t103 + t90;
t78 = -t105 * t88 + t108 * t85;
t93 = t103 * t105 + t108 * t97;
t95 = qJD(4) + t96;
t75 = pkin(4) * t95 - qJ(5) * t93 + t78;
t79 = t105 * t85 + t108 * t88;
t92 = t103 * t108 - t105 * t97;
t77 = qJ(5) * t92 + t79;
t72 = t104 * t75 + t114 * t77;
t89 = t100 * t109 + t106 * t101;
t87 = -pkin(3) * t103 - t89;
t71 = -t104 * t77 + t114 * t75;
t80 = -pkin(4) * t92 + qJD(5) + t87;
t82 = t104 * t92 + t114 * t93;
t81 = t104 * t93 - t114 * t92;
t73 = pkin(5) * t81 - qJ(6) * t82 + t80;
t70 = qJ(6) * t95 + t72;
t69 = -t95 * pkin(5) + qJD(6) - t71;
t1 = m(6) * (t71 ^ 2 + t72 ^ 2 + t80 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t87 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t102 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + (t102 * mrSges(4,2) - t89 * mrSges(4,3) + Ifges(4,1) * t97 / 0.2e1) * t97 + (t87 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t93 / 0.2e1) * t93 - (-t102 * mrSges(4,1) + t90 * mrSges(4,3) + Ifges(4,4) * t97 - Ifges(4,2) * t96 / 0.2e1) * t96 + (-t87 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t92 / 0.2e1) * t92 + (t89 * mrSges(4,1) - t90 * mrSges(4,2) + Ifges(4,5) * t97 - Ifges(4,6) * t96 + Ifges(4,3) * t103 / 0.2e1) * t103 + (t80 * mrSges(6,2) + t69 * mrSges(7,2) - t71 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t82) * t82 + (t80 * mrSges(6,1) + t73 * mrSges(7,1) - t70 * mrSges(7,2) - t72 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t81 + (-Ifges(6,4) + Ifges(7,5)) * t82) * t81 + (t78 * mrSges(5,1) + t71 * mrSges(6,1) - t69 * mrSges(7,1) - t79 * mrSges(5,2) - t72 * mrSges(6,2) + t70 * mrSges(7,3) + Ifges(5,5) * t93 + Ifges(5,6) * t92 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t95 + (Ifges(7,4) + Ifges(6,5)) * t82 + (-Ifges(6,6) + Ifges(7,6)) * t81) * t95 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t107 ^ 2 + t110 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t115) * t110) * t110 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t110 + (Ifges(3,1) / 0.2e1 + t115) * t107) * t107) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t110 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t107) * qJD(2)) * qJD(1);
T  = t1;
