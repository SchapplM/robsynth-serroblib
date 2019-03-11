% Calculate kinetic energy for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:14
% EndTime: 2019-03-09 14:08:15
% DurationCPUTime: 0.70s
% Computational Cost: add. (1446->114), mult. (3378->180), div. (0->0), fcn. (2746->12), ass. (0->52)
t139 = sin(qJ(5));
t143 = cos(qJ(5));
t141 = sin(qJ(2));
t145 = cos(qJ(2));
t135 = sin(pkin(6));
t150 = t135 * qJD(1);
t147 = t145 * t150;
t151 = cos(pkin(6)) * qJD(1);
t149 = pkin(1) * t151;
t127 = pkin(8) * t147 + t141 * t149;
t133 = qJD(2) + t151;
t122 = qJ(3) * t133 + t127;
t123 = (-pkin(2) * t145 - qJ(3) * t141 - pkin(1)) * t150;
t134 = sin(pkin(12));
t136 = cos(pkin(12));
t113 = -t122 * t134 + t136 * t123;
t148 = t141 * t150;
t125 = t133 * t134 + t136 * t148;
t110 = -pkin(3) * t147 - pkin(9) * t125 + t113;
t114 = t136 * t122 + t134 * t123;
t124 = t133 * t136 - t134 * t148;
t112 = pkin(9) * t124 + t114;
t140 = sin(qJ(4));
t144 = cos(qJ(4));
t100 = t144 * t110 - t112 * t140;
t117 = t124 * t140 + t125 * t144;
t129 = qJD(4) - t147;
t97 = pkin(4) * t129 - pkin(10) * t117 + t100;
t101 = t140 * t110 + t144 * t112;
t116 = t124 * t144 - t125 * t140;
t99 = pkin(10) * t116 + t101;
t94 = t139 * t97 + t143 * t99;
t93 = -t139 * t99 + t143 * t97;
t105 = t116 * t143 - t117 * t139;
t126 = -pkin(8) * t148 + t145 * t149;
t119 = -pkin(2) * t133 + qJD(3) - t126;
t115 = -pkin(3) * t124 + t119;
t107 = -pkin(4) * t116 + t115;
t146 = qJD(1) ^ 2;
t142 = cos(qJ(6));
t138 = sin(qJ(6));
t128 = qJD(5) + t129;
t106 = t116 * t139 + t117 * t143;
t104 = qJD(6) - t105;
t103 = t106 * t142 + t128 * t138;
t102 = -t106 * t138 + t128 * t142;
t95 = -pkin(5) * t105 - pkin(11) * t106 + t107;
t92 = pkin(11) * t128 + t94;
t91 = -pkin(5) * t128 - t93;
t90 = t138 * t95 + t142 * t92;
t89 = -t138 * t92 + t142 * t95;
t1 = m(4) * (t113 ^ 2 + t114 ^ 2 + t119 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t115 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t146 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t135 ^ 2 * t146 + t126 ^ 2 + t127 ^ 2) / 0.2e1 + (t126 * mrSges(3,1) - t127 * mrSges(3,2) + Ifges(3,3) * t133 / 0.2e1) * t133 + (t100 * mrSges(5,1) - t101 * mrSges(5,2) + Ifges(5,3) * t129 / 0.2e1) * t129 + (t93 * mrSges(6,1) - t94 * mrSges(6,2) + Ifges(6,3) * t128 / 0.2e1) * t128 + (t119 * mrSges(4,2) - t113 * mrSges(4,3) + Ifges(4,1) * t125 / 0.2e1) * t125 + (t89 * mrSges(7,1) - t90 * mrSges(7,2) + Ifges(7,3) * t104 / 0.2e1) * t104 + (t115 * mrSges(5,2) - t100 * mrSges(5,3) + Ifges(5,5) * t129 + Ifges(5,1) * t117 / 0.2e1) * t117 + (t107 * mrSges(6,2) - t93 * mrSges(6,3) + Ifges(6,5) * t128 + Ifges(6,1) * t106 / 0.2e1) * t106 + (t91 * mrSges(7,2) - t89 * mrSges(7,3) + Ifges(7,5) * t104 + Ifges(7,1) * t103 / 0.2e1) * t103 + ((-t126 * mrSges(3,3) + Ifges(3,5) * t133 + (-pkin(1) * mrSges(3,2) + Ifges(3,1) * t141 / 0.2e1) * t150) * t141 + (-t113 * mrSges(4,1) + t114 * mrSges(4,2) + t127 * mrSges(3,3) - Ifges(4,5) * t125 + Ifges(3,6) * t133 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t141 + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t145) * t150) * t145) * t150 + (-Ifges(4,6) * t147 - t119 * mrSges(4,1) + t114 * mrSges(4,3) + Ifges(4,4) * t125 + Ifges(4,2) * t124 / 0.2e1) * t124 + (-t115 * mrSges(5,1) + t101 * mrSges(5,3) + Ifges(5,4) * t117 + Ifges(5,6) * t129 + Ifges(5,2) * t116 / 0.2e1) * t116 + (-t107 * mrSges(6,1) + t94 * mrSges(6,3) + Ifges(6,4) * t106 + Ifges(6,6) * t128 + Ifges(6,2) * t105 / 0.2e1) * t105 + (-t91 * mrSges(7,1) + t90 * mrSges(7,3) + Ifges(7,4) * t103 + Ifges(7,6) * t104 + Ifges(7,2) * t102 / 0.2e1) * t102;
T  = t1;
