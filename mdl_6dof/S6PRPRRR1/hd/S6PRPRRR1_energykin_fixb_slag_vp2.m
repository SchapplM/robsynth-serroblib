% Calculate kinetic energy for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:20
% EndTime: 2019-03-08 20:22:20
% DurationCPUTime: 0.36s
% Computational Cost: add. (363->85), mult. (777->138), div. (0->0), fcn. (540->12), ass. (0->41)
t123 = sin(qJ(5));
t124 = sin(qJ(4));
t127 = cos(qJ(5));
t128 = cos(qJ(4));
t107 = (t123 * t124 - t127 * t128) * qJD(2);
t129 = cos(qJ(2));
t119 = sin(pkin(6));
t135 = qJD(1) * t119;
t110 = qJD(2) * pkin(2) + t129 * t135;
t118 = sin(pkin(12));
t120 = cos(pkin(12));
t125 = sin(qJ(2));
t133 = t125 * t135;
t103 = t110 * t118 + t120 * t133;
t101 = qJD(2) * pkin(8) + t103;
t121 = cos(pkin(6));
t115 = qJD(1) * t121 + qJD(3);
t114 = t128 * t115;
t94 = qJD(4) * pkin(4) + t114 + (-pkin(9) * qJD(2) - t101) * t124;
t134 = qJD(2) * t128;
t97 = t101 * t128 + t115 * t124;
t95 = pkin(9) * t134 + t97;
t90 = t123 * t94 + t127 * t95;
t102 = t110 * t120 - t118 * t133;
t89 = -t123 * t95 + t127 * t94;
t98 = (-pkin(4) * t128 - pkin(3)) * qJD(2) - t102;
t126 = cos(qJ(6));
t122 = sin(qJ(6));
t117 = qJD(4) + qJD(5);
t108 = (t123 * t128 + t124 * t127) * qJD(2);
t106 = qJD(6) + t107;
t105 = t108 * t126 + t117 * t122;
t104 = -t108 * t122 + t117 * t126;
t100 = -qJD(2) * pkin(3) - t102;
t96 = -t101 * t124 + t114;
t91 = pkin(5) * t107 - pkin(10) * t108 + t98;
t88 = pkin(10) * t117 + t90;
t87 = -pkin(5) * t117 - t89;
t86 = t122 * t91 + t126 * t88;
t85 = -t122 * t88 + t126 * t91;
t1 = m(4) * (t102 ^ 2 + t103 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t98 ^ 2) / 0.2e1 + (t89 * mrSges(6,1) - t90 * mrSges(6,2) + Ifges(6,3) * t117 / 0.2e1) * t117 + (t85 * mrSges(7,1) - t86 * mrSges(7,2) + Ifges(7,3) * t106 / 0.2e1) * t106 + (t96 * mrSges(5,1) - t97 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t98 * mrSges(6,2) - t89 * mrSges(6,3) + Ifges(6,5) * t117 + Ifges(6,1) * t108 / 0.2e1) * t108 + (t87 * mrSges(7,2) - t85 * mrSges(7,3) + Ifges(7,5) * t106 + Ifges(7,1) * t105 / 0.2e1) * t105 + (m(3) * (t121 ^ 2 + (t125 ^ 2 + t129 ^ 2) * t119 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 - (-t98 * mrSges(6,1) + t90 * mrSges(6,3) + Ifges(6,4) * t108 + Ifges(6,6) * t117 - Ifges(6,2) * t107 / 0.2e1) * t107 + (-t87 * mrSges(7,1) + t86 * mrSges(7,3) + Ifges(7,4) * t105 + Ifges(7,6) * t106 + Ifges(7,2) * t104 / 0.2e1) * t104 + (t102 * mrSges(4,1) - t103 * mrSges(4,2) + (mrSges(3,1) * t129 - mrSges(3,2) * t125) * t135 + (-t100 * mrSges(5,1) + t97 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t134 / 0.2e1) * t128 + (t100 * mrSges(5,2) - t96 * mrSges(5,3) + Ifges(5,5) * qJD(4)) * t124 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(5,4) * t128 + Ifges(5,1) * t124 / 0.2e1) * t124) * qJD(2)) * qJD(2);
T  = t1;
