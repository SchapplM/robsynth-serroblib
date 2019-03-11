% Calculate kinetic energy for
% S6RPPRRR8
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
% Datum: 2019-03-09 02:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:35:20
% EndTime: 2019-03-09 02:35:20
% DurationCPUTime: 0.45s
% Computational Cost: add. (498->87), mult. (1028->139), div. (0->0), fcn. (686->8), ass. (0->40)
t111 = cos(pkin(10));
t123 = t111 ^ 2;
t122 = m(3) / 0.2e1;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t110 = sin(pkin(10));
t101 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t119 = -pkin(7) * qJD(1) + t101;
t94 = t119 * t110;
t95 = t119 * t111;
t87 = t114 * t95 + t117 * t94;
t84 = qJD(4) * pkin(8) + t87;
t120 = qJD(1) * t110;
t97 = -t114 * t111 * qJD(1) - t117 * t120;
t98 = (-t110 * t114 + t111 * t117) * qJD(1);
t105 = qJD(1) * qJ(2) + qJD(3);
t99 = pkin(3) * t120 + t105;
t85 = -pkin(4) * t97 - pkin(8) * t98 + t99;
t76 = t113 * t85 + t116 * t84;
t121 = t110 ^ 2 + t123;
t96 = qJD(5) - t97;
t75 = -t113 * t84 + t116 * t85;
t86 = -t114 * t94 + t117 * t95;
t83 = -qJD(4) * pkin(4) - t86;
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t106 = -qJD(1) * pkin(1) + qJD(2);
t93 = qJD(6) + t96;
t89 = qJD(4) * t113 + t116 * t98;
t88 = qJD(4) * t116 - t113 * t98;
t79 = t112 * t88 + t115 * t89;
t78 = -t112 * t89 + t115 * t88;
t77 = -pkin(5) * t88 + t83;
t74 = pkin(9) * t88 + t76;
t73 = pkin(5) * t96 - pkin(9) * t89 + t75;
t72 = t112 * t73 + t115 * t74;
t71 = -t112 * t74 + t115 * t73;
t1 = m(5) * (t86 ^ 2 + t87 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t121 * t101 ^ 2 + t105 ^ 2) / 0.2e1 + t106 ^ 2 * t122 + m(7) * (t71 ^ 2 + t72 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t83 ^ 2) / 0.2e1 + (t99 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,1) * t98 / 0.2e1) * t98 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t96 / 0.2e1) * t96 + (t71 * mrSges(7,1) - t72 * mrSges(7,2) + Ifges(7,3) * t93 / 0.2e1) * t93 + (-t99 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,4) * t98 + Ifges(5,2) * t97 / 0.2e1) * t97 + (t83 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t96 + Ifges(6,1) * t89 / 0.2e1) * t89 + (t77 * mrSges(7,2) - t71 * mrSges(7,3) + Ifges(7,5) * t93 + Ifges(7,1) * t79 / 0.2e1) * t79 + (-t83 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t89 + Ifges(6,6) * t96 + Ifges(6,2) * t88 / 0.2e1) * t88 + (-t77 * mrSges(7,1) + t72 * mrSges(7,3) + Ifges(7,4) * t79 + Ifges(7,6) * t93 + Ifges(7,2) * t78 / 0.2e1) * t78 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,5) * t98 + Ifges(5,6) * t97 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t105 * (mrSges(4,1) * t110 + mrSges(4,2) * t111) + t106 * mrSges(3,2) - t121 * t101 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + (qJ(2) * t122 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t123 / 0.2e1 + (-Ifges(4,4) * t111 + Ifges(4,2) * t110 / 0.2e1) * t110) * qJD(1)) * qJD(1);
T  = t1;
