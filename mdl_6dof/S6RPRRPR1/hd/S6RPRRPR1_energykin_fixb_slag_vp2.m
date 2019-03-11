% Calculate kinetic energy for
% S6RPRRPR1
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:18
% EndTime: 2019-03-09 04:57:18
% DurationCPUTime: 0.54s
% Computational Cost: add. (576->98), mult. (1263->155), div. (0->0), fcn. (870->10), ass. (0->42)
t123 = m(3) / 0.2e1;
t109 = sin(pkin(11));
t111 = cos(pkin(11));
t114 = sin(qJ(4));
t115 = sin(qJ(3));
t117 = cos(qJ(4));
t118 = cos(qJ(3));
t101 = (t114 * t118 + t115 * t117) * qJD(1);
t108 = qJD(3) + qJD(4);
t110 = sin(pkin(10));
t104 = (pkin(1) * t110 + pkin(7)) * qJD(1);
t107 = t118 * qJD(2);
t122 = pkin(8) * qJD(1);
t96 = qJD(3) * pkin(3) + t107 + (-t104 - t122) * t115;
t99 = t115 * qJD(2) + t118 * t104;
t97 = t118 * t122 + t99;
t85 = -t114 * t97 + t117 * t96;
t82 = pkin(4) * t108 - qJ(5) * t101 + t85;
t100 = (-t114 * t115 + t117 * t118) * qJD(1);
t86 = t114 * t96 + t117 * t97;
t84 = qJ(5) * t100 + t86;
t79 = t109 * t82 + t111 * t84;
t112 = cos(pkin(10));
t121 = -pkin(1) * t112 - pkin(2);
t78 = -t109 * t84 + t111 * t82;
t90 = t100 * t111 - t101 * t109;
t102 = (-pkin(3) * t118 + t121) * qJD(1);
t92 = -pkin(4) * t100 + qJD(5) + t102;
t116 = cos(qJ(6));
t113 = sin(qJ(6));
t105 = t121 * qJD(1);
t98 = -t104 * t115 + t107;
t91 = t100 * t109 + t101 * t111;
t89 = qJD(6) - t90;
t88 = t108 * t113 + t116 * t91;
t87 = t108 * t116 - t113 * t91;
t80 = -pkin(5) * t90 - pkin(9) * t91 + t92;
t77 = pkin(9) * t108 + t79;
t76 = -pkin(5) * t108 - t78;
t75 = t113 * t80 + t116 * t77;
t74 = -t113 * t77 + t116 * t80;
t1 = qJD(2) ^ 2 * t123 + m(5) * (t102 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t105 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + (t92 * mrSges(6,2) - t78 * mrSges(6,3) + Ifges(6,1) * t91 / 0.2e1) * t91 + (t74 * mrSges(7,1) - t75 * mrSges(7,2) + Ifges(7,3) * t89 / 0.2e1) * t89 + (t102 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t101 / 0.2e1) * t101 + (-t92 * mrSges(6,1) + t79 * mrSges(6,3) + Ifges(6,4) * t91 + Ifges(6,2) * t90 / 0.2e1) * t90 + (t76 * mrSges(7,2) - t74 * mrSges(7,3) + Ifges(7,5) * t89 + Ifges(7,1) * t88 / 0.2e1) * t88 + (-t102 * mrSges(5,1) + t86 * mrSges(5,3) + Ifges(5,4) * t101 + Ifges(5,2) * t100 / 0.2e1) * t100 + (-t76 * mrSges(7,1) + t75 * mrSges(7,3) + Ifges(7,4) * t88 + Ifges(7,6) * t89 + Ifges(7,2) * t87 / 0.2e1) * t87 + (t98 * mrSges(4,1) - t99 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t85 * mrSges(5,1) + t78 * mrSges(6,1) - t86 * mrSges(5,2) - t79 * mrSges(6,2) + Ifges(5,5) * t101 + Ifges(6,5) * t91 + Ifges(5,6) * t100 + Ifges(6,6) * t90 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t108) * t108 + (t105 * (-mrSges(4,1) * t118 + mrSges(4,2) * t115) + (-t98 * t115 + t99 * t118) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t115 + Ifges(4,6) * t118) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t112 * mrSges(3,1) - t110 * mrSges(3,2) + (t110 ^ 2 + t112 ^ 2) * t123 * pkin(1)) * pkin(1) + Ifges(4,2) * t118 ^ 2 / 0.2e1 + (Ifges(4,4) * t118 + Ifges(4,1) * t115 / 0.2e1) * t115) * qJD(1)) * qJD(1);
T  = t1;
