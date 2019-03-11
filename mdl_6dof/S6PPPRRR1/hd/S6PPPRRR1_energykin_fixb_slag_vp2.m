% Calculate kinetic energy for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPPRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:04
% EndTime: 2019-03-08 18:39:05
% DurationCPUTime: 0.24s
% Computational Cost: add. (417->67), mult. (1024->117), div. (0->0), fcn. (907->16), ass. (0->42)
t134 = sin(pkin(14));
t135 = sin(pkin(13));
t139 = cos(pkin(14));
t138 = sin(pkin(6));
t154 = qJD(1) * t138;
t140 = cos(pkin(13));
t142 = cos(pkin(7));
t155 = t140 * t142;
t131 = cos(pkin(6)) * qJD(1) + qJD(2);
t137 = sin(pkin(7));
t156 = t131 * t137;
t125 = t139 * t156 + (-t134 * t135 + t139 * t155) * t154;
t128 = -t137 * t140 * t154 + t131 * t142 + qJD(3);
t136 = sin(pkin(8));
t141 = cos(pkin(8));
t159 = t125 * t141 + t128 * t136;
t126 = t134 * t156 + (t134 * t155 + t135 * t139) * t154;
t145 = sin(qJ(4));
t148 = cos(qJ(4));
t118 = -t126 * t145 + t148 * t159;
t119 = t148 * t126 + t145 * t159;
t117 = qJD(4) * pkin(10) + t119;
t121 = -t125 * t136 + t128 * t141;
t144 = sin(qJ(5));
t147 = cos(qJ(5));
t113 = t117 * t147 + t121 * t144;
t153 = qJD(4) * t144;
t152 = t147 * qJD(4);
t112 = -t117 * t144 + t121 * t147;
t149 = qJD(1) ^ 2;
t146 = cos(qJ(6));
t143 = sin(qJ(6));
t132 = qJD(6) - t152;
t130 = qJD(5) * t143 + t146 * t153;
t129 = qJD(5) * t146 - t143 * t153;
t116 = -qJD(4) * pkin(4) - t118;
t114 = (-pkin(5) * t147 - pkin(11) * t144 - pkin(4)) * qJD(4) - t118;
t111 = qJD(5) * pkin(11) + t113;
t110 = -qJD(5) * pkin(5) - t112;
t109 = t111 * t146 + t114 * t143;
t108 = -t111 * t143 + t114 * t146;
t1 = m(6) * (t112 ^ 2 + t113 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t118 ^ 2 + t119 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t131 ^ 2 + (t135 ^ 2 + t140 ^ 2) * t149 * t138 ^ 2) / 0.2e1 + m(4) * (t125 ^ 2 + t126 ^ 2 + t128 ^ 2) / 0.2e1 + m(2) * t149 / 0.2e1 + (t108 * mrSges(7,1) - t109 * mrSges(7,2) + Ifges(7,3) * t132 / 0.2e1) * t132 + (t112 * mrSges(6,1) - t113 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t110 * mrSges(7,2) - t108 * mrSges(7,3) + Ifges(7,5) * t132 + Ifges(7,1) * t130 / 0.2e1) * t130 + (-t110 * mrSges(7,1) + t109 * mrSges(7,3) + Ifges(7,4) * t130 + Ifges(7,6) * t132 + Ifges(7,2) * t129 / 0.2e1) * t129 + (-t119 * mrSges(5,2) + t118 * mrSges(5,1) + Ifges(5,3) * qJD(4) / 0.2e1 + (-t116 * mrSges(6,1) + t113 * mrSges(6,3) + Ifges(6,6) * qJD(5) + Ifges(6,2) * t152 / 0.2e1) * t147 + (t116 * mrSges(6,2) - t112 * mrSges(6,3) + Ifges(6,5) * qJD(5) + (Ifges(6,4) * t147 + Ifges(6,1) * t144 / 0.2e1) * qJD(4)) * t144) * qJD(4);
T  = t1;
