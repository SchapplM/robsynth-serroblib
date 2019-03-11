% Calculate kinetic energy for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:40:38
% EndTime: 2019-03-08 20:40:38
% DurationCPUTime: 0.35s
% Computational Cost: add. (320->82), mult. (626->132), div. (0->0), fcn. (396->10), ass. (0->40)
t112 = sin(qJ(5));
t113 = sin(qJ(4));
t116 = cos(qJ(5));
t117 = cos(qJ(4));
t99 = (t112 * t117 + t113 * t116) * qJD(2);
t114 = sin(qJ(2));
t109 = sin(pkin(6));
t126 = qJD(1) * t109;
t102 = qJD(2) * qJ(3) + t114 * t126;
t127 = t102 ^ 2;
t123 = qJD(2) * t117;
t110 = cos(pkin(6));
t125 = qJD(1) * t110;
t118 = cos(qJ(2));
t122 = -t118 * t126 + qJD(3);
t98 = (-pkin(2) - pkin(8)) * qJD(2) + t122;
t90 = -t113 * t125 + t117 * t98;
t88 = qJD(4) * pkin(4) - pkin(9) * t123 + t90;
t124 = qJD(2) * t113;
t91 = t113 * t98 + t117 * t125;
t89 = -pkin(9) * t124 + t91;
t84 = t112 * t88 + t116 * t89;
t83 = -t112 * t89 + t116 * t88;
t96 = pkin(4) * t124 + t102;
t120 = qJD(1) ^ 2;
t115 = cos(qJ(6));
t111 = sin(qJ(6));
t108 = qJD(4) + qJD(5);
t106 = t110 ^ 2 * t120;
t101 = -qJD(2) * pkin(2) + t122;
t100 = (-t112 * t113 + t116 * t117) * qJD(2);
t97 = qJD(6) + t99;
t93 = t115 * t100 + t111 * t108;
t92 = -t111 * t100 + t115 * t108;
t85 = t99 * pkin(5) - t100 * pkin(10) + t96;
t82 = t108 * pkin(10) + t84;
t81 = -t108 * pkin(5) - t83;
t80 = t111 * t85 + t115 * t82;
t79 = -t111 * t82 + t115 * t85;
t1 = m(4) * (t101 ^ 2 + t106 + t127) / 0.2e1 + m(2) * t120 / 0.2e1 + m(3) * (t106 + (t114 ^ 2 + t118 ^ 2) * t120 * t109 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t90 ^ 2 + t91 ^ 2 + t127) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 - (-t96 * mrSges(6,1) + t84 * mrSges(6,3) - Ifges(6,2) * t99 / 0.2e1) * t99 + (t79 * mrSges(7,1) - t80 * mrSges(7,2) + Ifges(7,3) * t97 / 0.2e1) * t97 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t81 * mrSges(7,2) - t79 * mrSges(7,3) + Ifges(7,5) * t97 + Ifges(7,1) * t93 / 0.2e1) * t93 + (t83 * mrSges(6,1) - t84 * mrSges(6,2) - Ifges(6,6) * t99 + Ifges(6,3) * t108 / 0.2e1) * t108 + (-t81 * mrSges(7,1) + t80 * mrSges(7,3) + Ifges(7,4) * t93 + Ifges(7,6) * t97 + Ifges(7,2) * t92 / 0.2e1) * t92 + (t96 * mrSges(6,2) - t83 * mrSges(6,3) - Ifges(6,4) * t99 + Ifges(6,5) * t108 + Ifges(6,1) * t100 / 0.2e1) * t100 + (t101 * mrSges(4,2) + t102 * mrSges(4,3) + (mrSges(3,1) * t118 - mrSges(3,2) * t114) * t126 + (t102 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t123 / 0.2e1) * t117 + (t102 * mrSges(5,1) - t91 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t113 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t117 + Ifges(5,2) * t113 / 0.2e1) * t113) * qJD(2)) * qJD(2);
T  = t1;
