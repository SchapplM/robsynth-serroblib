% Calculate kinetic energy for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:07
% EndTime: 2019-03-09 02:53:07
% DurationCPUTime: 0.36s
% Computational Cost: add. (511->94), mult. (1060->146), div. (0->0), fcn. (690->8), ass. (0->39)
t99 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t115 = t99 * mrSges(4,3);
t114 = qJD(1) / 0.2e1;
t104 = sin(pkin(10));
t106 = cos(pkin(10));
t105 = sin(pkin(9));
t107 = cos(pkin(9));
t111 = cos(qJ(3));
t113 = -qJ(4) * qJD(1) + t99;
t92 = qJD(3) * pkin(3) + t113 * t111;
t109 = sin(qJ(3));
t93 = t113 * t109;
t85 = t105 * t92 + t107 * t93;
t81 = qJD(3) * qJ(5) + t85;
t95 = (t105 * t111 + t107 * t109) * qJD(1);
t96 = (-t105 * t109 + t107 * t111) * qJD(1);
t103 = qJD(1) * qJ(2);
t97 = t109 * qJD(1) * pkin(3) + qJD(4) + t103;
t86 = pkin(4) * t95 - qJ(5) * t96 + t97;
t75 = t104 * t86 + t106 * t81;
t74 = -t104 * t81 + t106 * t86;
t84 = -t105 * t93 + t107 * t92;
t80 = -qJD(3) * pkin(4) + qJD(5) - t84;
t112 = qJD(1) ^ 2;
t110 = cos(qJ(6));
t108 = sin(qJ(6));
t102 = t112 * qJ(2) ^ 2;
t101 = -qJD(1) * pkin(1) + qJD(2);
t94 = qJD(6) + t95;
t88 = qJD(3) * t104 + t106 * t96;
t87 = qJD(3) * t106 - t104 * t96;
t78 = t108 * t87 + t110 * t88;
t77 = -t108 * t88 + t110 * t87;
t76 = -pkin(5) * t87 + t80;
t73 = pkin(8) * t87 + t75;
t72 = pkin(5) * t95 - pkin(8) * t88 + t74;
t71 = t108 * t72 + t110 * t73;
t70 = -t108 * t73 + t110 * t72;
t1 = m(4) * (t102 + (t109 ^ 2 + t111 ^ 2) * t99 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t101 ^ 2 + t102) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t80 ^ 2) / 0.2e1 + (t97 * mrSges(5,2) - t84 * mrSges(5,3) + Ifges(5,1) * t96 / 0.2e1) * t96 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,3) * t94 / 0.2e1) * t94 + (t80 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,1) * t88 / 0.2e1) * t88 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t112 + (-t80 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t88 + Ifges(6,2) * t87 / 0.2e1) * t87 + (t76 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,5) * t94 + Ifges(7,1) * t78 / 0.2e1) * t78 + (-t76 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t78 + Ifges(7,6) * t94 + Ifges(7,2) * t77 / 0.2e1) * t77 + (t101 * mrSges(3,2) + (mrSges(4,2) * t103 + (Ifges(4,1) * t114 - t115) * t111) * t111 + ((qJ(2) * mrSges(4,1) - Ifges(4,4) * t111) * qJD(1) + (Ifges(4,2) * t114 - t115) * t109) * t109) * qJD(1) + (t84 * mrSges(5,1) - t85 * mrSges(5,2) + Ifges(5,5) * t96 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (t111 * mrSges(4,1) - t109 * mrSges(4,2)) * t99 + (Ifges(4,5) * t111 - Ifges(4,6) * t109) * qJD(1)) * qJD(3) + (t97 * mrSges(5,1) + t74 * mrSges(6,1) - t75 * mrSges(6,2) - t85 * mrSges(5,3) - Ifges(5,4) * t96 + Ifges(6,5) * t88 - Ifges(5,6) * qJD(3) + Ifges(6,6) * t87 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t95) * t95;
T  = t1;
