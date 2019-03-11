% Calculate kinetic energy for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:31:04
% EndTime: 2019-03-08 23:31:04
% DurationCPUTime: 0.43s
% Computational Cost: add. (380->95), mult. (783->141), div. (0->0), fcn. (515->10), ass. (0->42)
t137 = -pkin(4) - pkin(5);
t136 = cos(qJ(4));
t123 = sin(qJ(2));
t118 = sin(pkin(6));
t135 = qJD(1) * t118;
t109 = qJD(2) * pkin(8) + t123 * t135;
t122 = sin(qJ(3));
t125 = cos(qJ(3));
t119 = cos(pkin(6));
t134 = qJD(1) * t119;
t103 = t125 * t109 + t122 * t134;
t100 = qJD(3) * pkin(9) + t103;
t126 = cos(qJ(2));
t131 = t126 * t135;
t104 = -t131 + (-pkin(3) * t125 - pkin(9) * t122 - pkin(2)) * qJD(2);
t121 = sin(qJ(4));
t93 = t136 * t100 + t121 * t104;
t102 = -t122 * t109 + t125 * t134;
t133 = qJD(2) * t122;
t132 = qJD(2) * t125;
t115 = qJD(4) - t132;
t91 = t115 * qJ(5) + t93;
t130 = qJD(3) * pkin(3) + t102;
t92 = -t121 * t100 + t136 * t104;
t129 = qJD(5) - t92;
t108 = t121 * qJD(3) + t136 * t133;
t128 = qJ(5) * t108 + t130;
t124 = cos(qJ(6));
t120 = sin(qJ(6));
t112 = qJD(6) - t115;
t110 = -qJD(2) * pkin(2) - t131;
t107 = -t136 * qJD(3) + t121 * t133;
t96 = t107 * t120 + t108 * t124;
t95 = t107 * t124 - t108 * t120;
t94 = pkin(4) * t107 - t128;
t90 = -t115 * pkin(4) + t129;
t89 = t137 * t107 + t128;
t88 = pkin(10) * t107 + t91;
t87 = -t108 * pkin(10) + t137 * t115 + t129;
t86 = t120 * t87 + t124 * t88;
t85 = -t120 * t88 + t124 * t87;
t1 = m(6) * (t90 ^ 2 + t91 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t130 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t103 ^ 2 + t110 ^ 2) / 0.2e1 + (t89 * mrSges(7,2) - t85 * mrSges(7,3) + Ifges(7,1) * t96 / 0.2e1) * t96 + (t102 * mrSges(4,1) - t103 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t89 * mrSges(7,1) + t86 * mrSges(7,3) + Ifges(7,4) * t96 + Ifges(7,2) * t95 / 0.2e1) * t95 + (m(3) * (t119 ^ 2 + (t123 ^ 2 + t126 ^ 2) * t118 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t85 * mrSges(7,1) - t86 * mrSges(7,2) + Ifges(7,5) * t96 + Ifges(7,6) * t95 + Ifges(7,3) * t112 / 0.2e1) * t112 + (t92 * mrSges(5,1) - t90 * mrSges(6,1) - t93 * mrSges(5,2) + t91 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t115) * t115 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t126 - mrSges(3,2) * t123) * t135 + (-t110 * mrSges(4,1) + t103 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t132 / 0.2e1) * t125 + (t110 * mrSges(4,2) - t102 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t125 + Ifges(4,1) * t122 / 0.2e1) * qJD(2)) * t122) * qJD(2) + (-t130 * mrSges(5,2) + t90 * mrSges(6,2) - t92 * mrSges(5,3) - t94 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t108 + (Ifges(6,4) + Ifges(5,5)) * t115) * t108 + (-t130 * mrSges(5,1) + t94 * mrSges(6,1) - t91 * mrSges(6,2) - t93 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t107 + (-Ifges(5,6) + Ifges(6,6)) * t115 + (-Ifges(5,4) + Ifges(6,5)) * t108) * t107;
T  = t1;
