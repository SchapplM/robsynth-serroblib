% Calculate kinetic energy for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:32:54
% EndTime: 2019-03-09 02:32:54
% DurationCPUTime: 0.44s
% Computational Cost: add. (524->87), mult. (1098->139), div. (0->0), fcn. (748->8), ass. (0->39)
t106 = cos(pkin(10));
t117 = t106 ^ 2;
t116 = m(3) / 0.2e1;
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t109 = sin(qJ(4));
t112 = cos(qJ(4));
t105 = sin(pkin(10));
t97 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t114 = -pkin(7) * qJD(1) + t97;
t91 = t114 * t105;
t92 = t114 * t106;
t82 = -t109 * t91 + t112 * t92;
t94 = (-t105 * t109 + t106 * t112) * qJD(1);
t78 = qJD(4) * pkin(4) - pkin(8) * t94 + t82;
t83 = t109 * t92 + t112 * t91;
t93 = (-t105 * t112 - t106 * t109) * qJD(1);
t79 = pkin(8) * t93 + t83;
t74 = t108 * t78 + t111 * t79;
t115 = t105 ^ 2 + t117;
t99 = qJD(1) * qJ(2) + qJD(3);
t95 = t105 * qJD(1) * pkin(3) + t99;
t73 = -t108 * t79 + t111 * t78;
t85 = -t108 * t94 + t111 * t93;
t87 = -pkin(4) * t93 + t95;
t110 = cos(qJ(6));
t107 = sin(qJ(6));
t103 = qJD(4) + qJD(5);
t100 = -qJD(1) * pkin(1) + qJD(2);
t86 = t108 * t93 + t111 * t94;
t84 = qJD(6) - t85;
t81 = t103 * t107 + t110 * t86;
t80 = t103 * t110 - t107 * t86;
t75 = -pkin(5) * t85 - pkin(9) * t86 + t87;
t72 = pkin(9) * t103 + t74;
t71 = -pkin(5) * t103 - t73;
t70 = t107 * t75 + t110 * t72;
t69 = -t107 * t72 + t110 * t75;
t1 = t100 ^ 2 * t116 + m(5) * (t82 ^ 2 + t83 ^ 2 + t95 ^ 2) / 0.2e1 + m(4) * (t115 * t97 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t87 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + (t95 * mrSges(5,2) - t82 * mrSges(5,3) + Ifges(5,1) * t94 / 0.2e1) * t94 + (t87 * mrSges(6,2) - t73 * mrSges(6,3) + Ifges(6,1) * t86 / 0.2e1) * t86 + (t69 * mrSges(7,1) - t70 * mrSges(7,2) + Ifges(7,3) * t84 / 0.2e1) * t84 + (-t95 * mrSges(5,1) + t83 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,2) * t93 / 0.2e1) * t93 + (-t87 * mrSges(6,1) + t74 * mrSges(6,3) + Ifges(6,4) * t86 + Ifges(6,2) * t85 / 0.2e1) * t85 + (t71 * mrSges(7,2) - t69 * mrSges(7,3) + Ifges(7,5) * t84 + Ifges(7,1) * t81 / 0.2e1) * t81 + (-t71 * mrSges(7,1) + t70 * mrSges(7,3) + Ifges(7,4) * t81 + Ifges(7,6) * t84 + Ifges(7,2) * t80 / 0.2e1) * t80 + (t73 * mrSges(6,1) - t74 * mrSges(6,2) + Ifges(6,5) * t86 + Ifges(6,6) * t85 + Ifges(6,3) * t103 / 0.2e1) * t103 + (t82 * mrSges(5,1) - t83 * mrSges(5,2) + Ifges(5,5) * t94 + Ifges(5,6) * t93 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t99 * (mrSges(4,1) * t105 + mrSges(4,2) * t106) + t100 * mrSges(3,2) - t115 * t97 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t116 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t117 / 0.2e1 + (-Ifges(4,4) * t106 + Ifges(4,2) * t105 / 0.2e1) * t105) * qJD(1)) * qJD(1);
T  = t1;
