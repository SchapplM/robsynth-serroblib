% Calculate kinetic energy for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR14_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:48
% EndTime: 2019-12-31 19:16:49
% DurationCPUTime: 0.54s
% Computational Cost: add. (834->100), mult. (2651->172), div. (0->0), fcn. (2168->12), ass. (0->49)
t123 = sin(pkin(11));
t125 = sin(pkin(5));
t131 = sin(qJ(3));
t134 = cos(qJ(3));
t126 = cos(pkin(11));
t127 = cos(pkin(6));
t142 = t126 * t127;
t124 = sin(pkin(6));
t128 = cos(pkin(5));
t144 = t124 * t128;
t111 = ((-t123 * t131 + t134 * t142) * t125 + t134 * t144) * qJD(1);
t140 = qJD(1) * t125;
t137 = qJ(2) * t140;
t139 = pkin(1) * qJD(1) * t128;
t118 = t123 * t139 + t126 * t137;
t108 = (t125 * t142 + t144) * qJD(1) * pkin(8) + t118;
t121 = t126 * t139;
t110 = t121 + (pkin(2) * t128 + (-pkin(8) * t127 - qJ(2)) * t125 * t123) * qJD(1);
t115 = qJD(2) + (-pkin(8) * t123 * t124 - pkin(2) * t126 - pkin(1)) * t140;
t97 = -t131 * t108 + (t110 * t127 + t115 * t124) * t134;
t130 = sin(qJ(4));
t133 = cos(qJ(4));
t101 = -t110 * t124 + t127 * t115;
t141 = t127 * t131;
t143 = t124 * t131;
t112 = (t128 * t143 + (t123 * t134 + t126 * t141) * t125) * qJD(1);
t93 = -pkin(3) * t111 - pkin(9) * t112 + t101;
t116 = qJD(3) + (-t124 * t125 * t126 + t127 * t128) * qJD(1);
t98 = t134 * t108 + t110 * t141 + t115 * t143;
t96 = pkin(9) * t116 + t98;
t90 = t130 * t93 + t133 * t96;
t89 = -t130 * t96 + t133 * t93;
t103 = -t112 * t130 + t116 * t133;
t95 = -t116 * pkin(3) - t97;
t132 = cos(qJ(5));
t129 = sin(qJ(5));
t122 = -pkin(1) * t140 + qJD(2);
t117 = -t123 * t137 + t121;
t109 = qJD(4) - t111;
t104 = t112 * t133 + t116 * t130;
t102 = qJD(5) - t103;
t100 = t104 * t132 + t109 * t129;
t99 = -t104 * t129 + t109 * t132;
t91 = -t103 * pkin(4) - t104 * pkin(10) + t95;
t88 = pkin(10) * t109 + t90;
t87 = -pkin(4) * t109 - t89;
t86 = t129 * t91 + t132 * t88;
t85 = -t129 * t88 + t132 * t91;
t1 = m(6) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t95 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t97 ^ 2 + t98 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t118 ^ 2 + t122 ^ 2) / 0.2e1 + (-t87 * mrSges(6,1) + t86 * mrSges(6,3) + Ifges(6,2) * t99 / 0.2e1) * t99 + (t97 * mrSges(4,1) - t98 * mrSges(4,2) + Ifges(4,3) * t116 / 0.2e1) * t116 + (t89 * mrSges(5,1) - t90 * mrSges(5,2) + Ifges(5,3) * t109 / 0.2e1) * t109 + (t101 * mrSges(4,2) - t97 * mrSges(4,3) + Ifges(4,5) * t116 + Ifges(4,1) * t112 / 0.2e1) * t112 + (t95 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,5) * t109 + Ifges(5,1) * t104 / 0.2e1) * t104 + (t85 * mrSges(6,1) - t86 * mrSges(6,2) + Ifges(6,6) * t99 + Ifges(6,3) * t102 / 0.2e1) * t102 + (-t101 * mrSges(4,1) + t98 * mrSges(4,3) + Ifges(4,4) * t112 + Ifges(4,6) * t116 + Ifges(4,2) * t111 / 0.2e1) * t111 + (-t95 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t104 + Ifges(5,6) * t109 + Ifges(5,2) * t103 / 0.2e1) * t103 + (t87 * mrSges(6,2) - t85 * mrSges(6,3) + Ifges(6,4) * t99 + Ifges(6,5) * t102 + Ifges(6,1) * t100 / 0.2e1) * t100 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t122 * (-mrSges(3,1) * t126 + mrSges(3,2) * t123) + (Ifges(3,2) * t126 ^ 2 / 0.2e1 + (Ifges(3,4) * t126 + Ifges(3,1) * t123 / 0.2e1) * t123) * t140 + (-t117 * t123 + t118 * t126) * mrSges(3,3)) * t125 + (t117 * mrSges(3,1) - t118 * mrSges(3,2) + (Ifges(3,3) * t128 / 0.2e1 + (Ifges(3,5) * t123 + Ifges(3,6) * t126) * t125) * qJD(1)) * t128) * qJD(1);
T = t1;
