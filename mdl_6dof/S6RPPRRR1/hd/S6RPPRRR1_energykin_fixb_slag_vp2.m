% Calculate kinetic energy for
% S6RPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:43
% EndTime: 2019-03-09 02:17:44
% DurationCPUTime: 0.47s
% Computational Cost: add. (533->93), mult. (1228->148), div. (0->0), fcn. (868->10), ass. (0->42)
t123 = m(3) / 0.2e1;
t114 = sin(qJ(5));
t117 = cos(qJ(5));
t109 = sin(pkin(11));
t111 = cos(pkin(11));
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t102 = (t109 * t118 + t111 * t115) * qJD(1);
t110 = sin(pkin(10));
t105 = (pkin(1) * t110 + qJ(3)) * qJD(1);
t107 = t111 * qJD(2);
t122 = pkin(7) * qJD(1);
t96 = t107 + (-t105 - t122) * t109;
t99 = t109 * qJD(2) + t111 * t105;
t97 = t111 * t122 + t99;
t85 = -t115 * t97 + t118 * t96;
t83 = qJD(4) * pkin(4) - pkin(8) * t102 + t85;
t101 = (-t109 * t115 + t111 * t118) * qJD(1);
t86 = t115 * t96 + t118 * t97;
t84 = pkin(8) * t101 + t86;
t79 = t114 * t83 + t117 * t84;
t112 = cos(pkin(10));
t121 = -pkin(1) * t112 - pkin(2);
t78 = -t114 * t84 + t117 * t83;
t90 = t101 * t117 - t102 * t114;
t100 = qJD(3) + (-pkin(3) * t111 + t121) * qJD(1);
t92 = -pkin(4) * t101 + t100;
t116 = cos(qJ(6));
t113 = sin(qJ(6));
t108 = qJD(4) + qJD(5);
t104 = t121 * qJD(1) + qJD(3);
t98 = -t105 * t109 + t107;
t91 = t101 * t114 + t102 * t117;
t89 = qJD(6) - t90;
t88 = t108 * t113 + t116 * t91;
t87 = t108 * t116 - t113 * t91;
t80 = -pkin(5) * t90 - pkin(9) * t91 + t92;
t77 = pkin(9) * t108 + t79;
t76 = -pkin(5) * t108 - t78;
t75 = t113 * t80 + t116 * t77;
t74 = -t113 * t77 + t116 * t80;
t1 = qJD(2) ^ 2 * t123 + m(4) * (t104 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + (t92 * mrSges(6,2) - t78 * mrSges(6,3) + Ifges(6,1) * t91 / 0.2e1) * t91 + (t74 * mrSges(7,1) - t75 * mrSges(7,2) + Ifges(7,3) * t89 / 0.2e1) * t89 + (t100 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t102 / 0.2e1) * t102 + (-t92 * mrSges(6,1) + t79 * mrSges(6,3) + Ifges(6,4) * t91 + Ifges(6,2) * t90 / 0.2e1) * t90 + (t76 * mrSges(7,2) - t74 * mrSges(7,3) + Ifges(7,5) * t89 + Ifges(7,1) * t88 / 0.2e1) * t88 + (-t100 * mrSges(5,1) + t86 * mrSges(5,3) + Ifges(5,4) * t102 + Ifges(5,2) * t101 / 0.2e1) * t101 + (-t76 * mrSges(7,1) + t75 * mrSges(7,3) + Ifges(7,4) * t88 + Ifges(7,6) * t89 + Ifges(7,2) * t87 / 0.2e1) * t87 + (t78 * mrSges(6,1) - t79 * mrSges(6,2) + Ifges(6,5) * t91 + Ifges(6,6) * t90 + Ifges(6,3) * t108 / 0.2e1) * t108 + (t85 * mrSges(5,1) - t86 * mrSges(5,2) + Ifges(5,5) * t102 + Ifges(5,6) * t101 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t104 * (-mrSges(4,1) * t111 + mrSges(4,2) * t109) + (-t98 * t109 + t99 * t111) * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t112 * mrSges(3,1) - t110 * mrSges(3,2) + (t110 ^ 2 + t112 ^ 2) * t123 * pkin(1)) * pkin(1) + Ifges(4,2) * t111 ^ 2 / 0.2e1 + (Ifges(4,4) * t111 + Ifges(4,1) * t109 / 0.2e1) * t109) * qJD(1)) * qJD(1);
T  = t1;
