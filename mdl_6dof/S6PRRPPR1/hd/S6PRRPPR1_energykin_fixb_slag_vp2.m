% Calculate kinetic energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:31
% EndTime: 2019-03-08 20:58:31
% DurationCPUTime: 0.42s
% Computational Cost: add. (526->99), mult. (1231->158), div. (0->0), fcn. (907->12), ass. (0->44)
t132 = cos(qJ(3));
t133 = cos(qJ(2));
t125 = sin(pkin(6));
t138 = qJD(1) * t125;
t135 = t133 * t138;
t112 = -t135 + qJD(4) + (-pkin(3) * t132 - pkin(2)) * qJD(2);
t124 = sin(pkin(11));
t129 = sin(qJ(3));
t136 = qJD(2) * t132;
t139 = cos(pkin(11));
t114 = qJD(2) * t124 * t129 - t139 * t136;
t115 = (t124 * t132 + t139 * t129) * qJD(2);
t102 = pkin(4) * t114 - qJ(5) * t115 + t112;
t123 = sin(pkin(12));
t126 = cos(pkin(12));
t130 = sin(qJ(2));
t117 = qJD(2) * pkin(8) + t130 * t138;
t127 = cos(pkin(6));
t137 = qJD(1) * t127;
t121 = t132 * t137;
t106 = qJD(3) * pkin(3) + t121 + (-qJ(4) * qJD(2) - t117) * t129;
t111 = t132 * t117 + t129 * t137;
t107 = qJ(4) * t136 + t111;
t97 = t124 * t106 + t139 * t107;
t95 = qJD(3) * qJ(5) + t97;
t91 = t123 * t102 + t126 * t95;
t90 = t126 * t102 - t123 * t95;
t96 = t139 * t106 - t124 * t107;
t94 = -qJD(3) * pkin(4) + qJD(5) - t96;
t131 = cos(qJ(6));
t128 = sin(qJ(6));
t118 = -qJD(2) * pkin(2) - t135;
t113 = qJD(6) + t114;
t110 = -t117 * t129 + t121;
t109 = qJD(3) * t123 + t115 * t126;
t108 = qJD(3) * t126 - t115 * t123;
t99 = t108 * t128 + t109 * t131;
t98 = t108 * t131 - t109 * t128;
t92 = -t108 * pkin(5) + t94;
t89 = pkin(9) * t108 + t91;
t88 = pkin(5) * t114 - pkin(9) * t109 + t90;
t87 = t128 * t88 + t131 * t89;
t86 = -t128 * t89 + t131 * t88;
t1 = m(4) * (t110 ^ 2 + t111 ^ 2 + t118 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t86 ^ 2 + t87 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t90 ^ 2 + t91 ^ 2 + t94 ^ 2) / 0.2e1 + (t92 * mrSges(7,2) - t86 * mrSges(7,3) + Ifges(7,1) * t99 / 0.2e1) * t99 + (t112 * mrSges(5,2) - t96 * mrSges(5,3) + Ifges(5,1) * t115 / 0.2e1) * t115 + (t94 * mrSges(6,2) - t90 * mrSges(6,3) + Ifges(6,1) * t109 / 0.2e1) * t109 + (-t92 * mrSges(7,1) + t87 * mrSges(7,3) + Ifges(7,4) * t99 + Ifges(7,2) * t98 / 0.2e1) * t98 + (-t94 * mrSges(6,1) + t91 * mrSges(6,3) + Ifges(6,4) * t109 + Ifges(6,2) * t108 / 0.2e1) * t108 + (m(2) / 0.2e1 + m(3) * (t127 ^ 2 + (t130 ^ 2 + t133 ^ 2) * t125 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t86 * mrSges(7,1) - t87 * mrSges(7,2) + Ifges(7,5) * t99 + Ifges(7,6) * t98 + Ifges(7,3) * t113 / 0.2e1) * t113 + (t112 * mrSges(5,1) + t90 * mrSges(6,1) - t91 * mrSges(6,2) - t97 * mrSges(5,3) - Ifges(5,4) * t115 + Ifges(6,5) * t109 + Ifges(6,6) * t108 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t114) * t114 + (t118 * (-mrSges(4,1) * t132 + mrSges(4,2) * t129) + (Ifges(4,2) * t132 ^ 2 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t132 + Ifges(4,1) * t129 / 0.2e1) * t129) * qJD(2) + (mrSges(3,1) * t133 - mrSges(3,2) * t130) * t138 + (-t110 * t129 + t111 * t132) * mrSges(4,3)) * qJD(2) + (t110 * mrSges(4,1) + t96 * mrSges(5,1) - t111 * mrSges(4,2) - t97 * mrSges(5,2) + Ifges(5,5) * t115 - Ifges(5,6) * t114 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (Ifges(4,5) * t129 + Ifges(4,6) * t132) * qJD(2)) * qJD(3);
T  = t1;
