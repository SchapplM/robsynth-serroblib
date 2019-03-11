% Calculate kinetic energy for
% S6PRPRRR6
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:05
% EndTime: 2019-03-08 20:45:06
% DurationCPUTime: 0.33s
% Computational Cost: add. (316->83), mult. (624->131), div. (0->0), fcn. (392->10), ass. (0->40)
t117 = sin(qJ(2));
t112 = sin(pkin(6));
t128 = qJD(1) * t112;
t125 = t117 * t128;
t103 = qJD(2) * qJ(3) + t125;
t129 = t103 ^ 2;
t115 = sin(qJ(5));
t119 = cos(qJ(5));
t116 = sin(qJ(4));
t120 = cos(qJ(4));
t113 = cos(pkin(6));
t127 = qJD(1) * t113;
t121 = cos(qJ(2));
t124 = -t121 * t128 + qJD(3);
t99 = (-pkin(2) - pkin(8)) * qJD(2) + t124;
t94 = t116 * t99 + t120 * t127;
t92 = qJD(4) * pkin(9) + t94;
t97 = t125 + (pkin(4) * t116 - pkin(9) * t120 + qJ(3)) * qJD(2);
t86 = t115 * t97 + t119 * t92;
t126 = qJD(2) * t120;
t108 = t116 * qJD(2) + qJD(5);
t85 = -t115 * t92 + t119 * t97;
t93 = -t116 * t127 + t120 * t99;
t91 = -qJD(4) * pkin(4) - t93;
t123 = qJD(1) ^ 2;
t118 = cos(qJ(6));
t114 = sin(qJ(6));
t109 = t113 ^ 2 * t123;
t105 = qJD(6) + t108;
t102 = qJD(4) * t115 + t119 * t126;
t101 = qJD(4) * t119 - t115 * t126;
t100 = -qJD(2) * pkin(2) + t124;
t90 = t101 * t114 + t102 * t118;
t89 = t101 * t118 - t102 * t114;
t87 = -t101 * pkin(5) + t91;
t84 = pkin(10) * t101 + t86;
t83 = t108 * pkin(5) - t102 * pkin(10) + t85;
t82 = t114 * t83 + t118 * t84;
t81 = -t114 * t84 + t118 * t83;
t1 = m(5) * (t93 ^ 2 + t94 ^ 2 + t129) / 0.2e1 + m(7) * (t81 ^ 2 + t82 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t85 ^ 2 + t86 ^ 2 + t91 ^ 2) / 0.2e1 + m(3) * (t109 + (t117 ^ 2 + t121 ^ 2) * t123 * t112 ^ 2) / 0.2e1 + m(2) * t123 / 0.2e1 + m(4) * (t100 ^ 2 + t109 + t129) / 0.2e1 + (t87 * mrSges(7,2) - t81 * mrSges(7,3) + Ifges(7,1) * t90 / 0.2e1) * t90 + (t85 * mrSges(6,1) - t86 * mrSges(6,2) + Ifges(6,3) * t108 / 0.2e1) * t108 + (t93 * mrSges(5,1) - t94 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (-t87 * mrSges(7,1) + t82 * mrSges(7,3) + Ifges(7,4) * t90 + Ifges(7,2) * t89 / 0.2e1) * t89 + (t91 * mrSges(6,2) - t85 * mrSges(6,3) + Ifges(6,5) * t108 + Ifges(6,1) * t102 / 0.2e1) * t102 + (t81 * mrSges(7,1) - t82 * mrSges(7,2) + Ifges(7,5) * t90 + Ifges(7,6) * t89 + Ifges(7,3) * t105 / 0.2e1) * t105 + (-t91 * mrSges(6,1) + t86 * mrSges(6,3) + Ifges(6,4) * t102 + Ifges(6,6) * t108 + Ifges(6,2) * t101 / 0.2e1) * t101 + (t100 * mrSges(4,2) + t103 * mrSges(4,3) + (mrSges(3,1) * t121 - mrSges(3,2) * t117) * t128 + (t103 * mrSges(5,2) - t93 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t126 / 0.2e1) * t120 + (t103 * mrSges(5,1) - t94 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t116 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t120 + Ifges(5,2) * t116 / 0.2e1) * t116) * qJD(2)) * qJD(2);
T  = t1;
