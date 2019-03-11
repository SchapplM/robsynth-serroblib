% Calculate kinetic energy for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:27
% EndTime: 2019-03-08 19:29:27
% DurationCPUTime: 0.30s
% Computational Cost: add. (355->84), mult. (789->134), div. (0->0), fcn. (544->12), ass. (0->40)
t118 = sin(pkin(12));
t121 = cos(pkin(12));
t129 = cos(qJ(2));
t120 = sin(pkin(6));
t135 = qJD(1) * t120;
t111 = qJD(2) * pkin(2) + t129 * t135;
t119 = sin(pkin(11));
t122 = cos(pkin(11));
t126 = sin(qJ(2));
t132 = t126 * t135;
t107 = t119 * t111 + t122 * t132;
t105 = qJD(2) * pkin(8) + t107;
t123 = cos(pkin(6));
t115 = t123 * qJD(1) + qJD(3);
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t98 = t128 * t105 + t125 * t115;
t94 = qJD(4) * qJ(5) + t98;
t106 = t122 * t111 - t119 * t132;
t99 = (-pkin(4) * t128 - qJ(5) * t125 - pkin(3)) * qJD(2) - t106;
t90 = t118 * t99 + t121 * t94;
t134 = qJD(2) * t125;
t133 = t128 * qJD(2);
t89 = -t118 * t94 + t121 * t99;
t97 = -t125 * t105 + t128 * t115;
t93 = -qJD(4) * pkin(4) + qJD(5) - t97;
t127 = cos(qJ(6));
t124 = sin(qJ(6));
t116 = qJD(6) - t133;
t110 = t118 * qJD(4) + t121 * t134;
t109 = t121 * qJD(4) - t118 * t134;
t104 = -qJD(2) * pkin(3) - t106;
t101 = t124 * t109 + t127 * t110;
t100 = t127 * t109 - t124 * t110;
t91 = -t109 * pkin(5) + t93;
t88 = t109 * pkin(9) + t90;
t87 = -pkin(5) * t133 - t110 * pkin(9) + t89;
t86 = t124 * t87 + t127 * t88;
t85 = -t124 * t88 + t127 * t87;
t1 = m(4) * (t106 ^ 2 + t107 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(7) * (t85 ^ 2 + t86 ^ 2 + t91 ^ 2) / 0.2e1 + m(6) * (t89 ^ 2 + t90 ^ 2 + t93 ^ 2) / 0.2e1 + (t85 * mrSges(7,1) - t86 * mrSges(7,2) + Ifges(7,3) * t116 / 0.2e1) * t116 + (t93 * mrSges(6,2) - t89 * mrSges(6,3) + Ifges(6,1) * t110 / 0.2e1) * t110 + (t97 * mrSges(5,1) - t98 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (-t93 * mrSges(6,1) + t90 * mrSges(6,3) + Ifges(6,4) * t110 + Ifges(6,2) * t109 / 0.2e1) * t109 + (t91 * mrSges(7,2) - t85 * mrSges(7,3) + Ifges(7,5) * t116 + Ifges(7,1) * t101 / 0.2e1) * t101 + (m(2) / 0.2e1 + m(3) * (t123 ^ 2 + (t126 ^ 2 + t129 ^ 2) * t120 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (-t91 * mrSges(7,1) + t86 * mrSges(7,3) + Ifges(7,4) * t101 + Ifges(7,6) * t116 + Ifges(7,2) * t100 / 0.2e1) * t100 + (t106 * mrSges(4,1) - t107 * mrSges(4,2) + (mrSges(3,1) * t129 - mrSges(3,2) * t126) * t135 + (t104 * mrSges(5,2) - t97 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t134 / 0.2e1) * t125 + (-t104 * mrSges(5,1) - t89 * mrSges(6,1) + t90 * mrSges(6,2) + t98 * mrSges(5,3) - Ifges(6,5) * t110 + Ifges(5,6) * qJD(4) - Ifges(6,6) * t109) * t128 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(5,4) * t125 + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t128) * t128) * qJD(2)) * qJD(2);
T  = t1;
