% Calculate kinetic energy for
% S6RPPRRR4
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:28
% EndTime: 2019-03-09 02:25:29
% DurationCPUTime: 0.32s
% Computational Cost: add. (385->83), mult. (682->128), div. (0->0), fcn. (362->8), ass. (0->37)
t117 = m(3) / 0.2e1;
t109 = sin(qJ(5));
t112 = cos(qJ(5));
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t105 = sin(pkin(10));
t106 = cos(pkin(10));
t116 = qJ(2) * qJD(1);
t98 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t93 = t105 * t98 + t106 * t116;
t91 = -qJD(1) * pkin(7) + t93;
t87 = t110 * qJD(3) + t113 * t91;
t82 = qJD(4) * pkin(8) + t87;
t92 = -t105 * t116 + t106 * t98;
t90 = qJD(1) * pkin(3) - t92;
t83 = (pkin(4) * t113 + pkin(8) * t110) * qJD(1) + t90;
t76 = t109 * t83 + t112 * t82;
t115 = qJD(1) * t110;
t104 = t113 * qJD(1);
t99 = t104 + qJD(5);
t75 = -t109 * t82 + t112 * t83;
t86 = t113 * qJD(3) - t110 * t91;
t81 = -qJD(4) * pkin(4) - t86;
t111 = cos(qJ(6));
t108 = sin(qJ(6));
t102 = -qJD(1) * pkin(1) + qJD(2);
t97 = qJD(6) + t99;
t95 = t109 * qJD(4) - t112 * t115;
t94 = t112 * qJD(4) + t109 * t115;
t85 = t108 * t94 + t111 * t95;
t84 = -t108 * t95 + t111 * t94;
t77 = -t94 * pkin(5) + t81;
t74 = t94 * pkin(9) + t76;
t73 = t99 * pkin(5) - t95 * pkin(9) + t75;
t72 = t108 * t73 + t111 * t74;
t71 = -t108 * t74 + t111 * t73;
t1 = t102 ^ 2 * t117 + m(5) * (t86 ^ 2 + t87 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t71 ^ 2 + t72 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t99 / 0.2e1) * t99 + (t71 * mrSges(7,1) - t72 * mrSges(7,2) + Ifges(7,3) * t97 / 0.2e1) * t97 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t81 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t99 + Ifges(6,1) * t95 / 0.2e1) * t95 + (t77 * mrSges(7,2) - t71 * mrSges(7,3) + Ifges(7,5) * t97 + Ifges(7,1) * t85 / 0.2e1) * t85 + (-t81 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t95 + Ifges(6,6) * t99 + Ifges(6,2) * t94 / 0.2e1) * t94 + (-t77 * mrSges(7,1) + t72 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,6) * t97 + Ifges(7,2) * t84 / 0.2e1) * t84 + (-t102 * mrSges(3,1) - t92 * mrSges(4,1) + t93 * mrSges(4,2) + (t90 * mrSges(5,1) - t87 * mrSges(5,3) - Ifges(5,6) * qJD(4) + Ifges(5,2) * t104 / 0.2e1) * t113 + (-t90 * mrSges(5,2) + t86 * mrSges(5,3) - Ifges(5,5) * qJD(4)) * t110 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (qJ(2) * t117 + mrSges(3,3)) * qJ(2) + (Ifges(5,4) * t113 + Ifges(5,1) * t110 / 0.2e1) * t110) * qJD(1)) * qJD(1);
T  = t1;
