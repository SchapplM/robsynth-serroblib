% Calculate kinetic energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:13
% EndTime: 2019-03-08 19:07:13
% DurationCPUTime: 0.41s
% Computational Cost: add. (606->85), mult. (1463->147), div. (0->0), fcn. (1250->16), ass. (0->51)
t145 = cos(pkin(6)) * qJD(1) + qJD(2);
t149 = sin(pkin(7));
t153 = cos(pkin(7));
t151 = cos(pkin(14));
t150 = sin(pkin(6));
t169 = qJD(1) * t150;
t166 = t151 * t169;
t174 = t145 * t149 + t153 * t166;
t157 = sin(qJ(3));
t161 = cos(qJ(3));
t147 = sin(pkin(14));
t167 = t147 * t169;
t130 = -t157 * t167 + t174 * t161;
t127 = qJD(3) * pkin(3) + t130;
t135 = t153 * t145 - t149 * t166;
t148 = sin(pkin(8));
t152 = cos(pkin(8));
t173 = t127 * t152 + t135 * t148;
t131 = t174 * t157 + t161 * t167;
t168 = t148 * qJD(3);
t126 = pkin(10) * t168 + t131;
t156 = sin(qJ(4));
t160 = cos(qJ(4));
t118 = -t156 * t126 + t173 * t160;
t119 = t160 * t126 + t173 * t156;
t144 = t152 * qJD(3) + qJD(4);
t117 = t144 * pkin(11) + t119;
t133 = t152 * t135;
t121 = t133 + (-t127 + (-pkin(4) * t160 - pkin(11) * t156) * qJD(3)) * t148;
t155 = sin(qJ(5));
t159 = cos(qJ(5));
t113 = t159 * t117 + t155 * t121;
t165 = t156 * t168;
t112 = -t155 * t117 + t159 * t121;
t136 = t159 * t144 - t155 * t165;
t116 = -t144 * pkin(4) - t118;
t162 = qJD(1) ^ 2;
t158 = cos(qJ(6));
t154 = sin(qJ(6));
t143 = -t160 * t168 + qJD(5);
t137 = t155 * t144 + t159 * t165;
t134 = qJD(6) - t136;
t129 = t158 * t137 + t154 * t143;
t128 = -t154 * t137 + t158 * t143;
t122 = -t148 * t127 + t133;
t114 = -t136 * pkin(5) - t137 * pkin(12) + t116;
t111 = t143 * pkin(12) + t113;
t110 = -t143 * pkin(5) - t112;
t109 = t158 * t111 + t154 * t114;
t108 = -t154 * t111 + t158 * t114;
t1 = m(3) * (t145 ^ 2 + (t147 ^ 2 + t151 ^ 2) * t162 * t150 ^ 2) / 0.2e1 + m(2) * t162 / 0.2e1 + m(4) * (t130 ^ 2 + t131 ^ 2 + t135 ^ 2) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t116 ^ 2) / 0.2e1 + m(5) * (t118 ^ 2 + t119 ^ 2 + t122 ^ 2) / 0.2e1 + m(7) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + (t118 * mrSges(5,1) - t119 * mrSges(5,2) + Ifges(5,3) * t144 / 0.2e1) * t144 + (t112 * mrSges(6,1) - t113 * mrSges(6,2) + Ifges(6,3) * t143 / 0.2e1) * t143 + (t108 * mrSges(7,1) - t109 * mrSges(7,2) + Ifges(7,3) * t134 / 0.2e1) * t134 + (t116 * mrSges(6,2) - t112 * mrSges(6,3) + Ifges(6,5) * t143 + Ifges(6,1) * t137 / 0.2e1) * t137 + (t110 * mrSges(7,2) - t108 * mrSges(7,3) + Ifges(7,5) * t134 + Ifges(7,1) * t129 / 0.2e1) * t129 + (-t116 * mrSges(6,1) + t113 * mrSges(6,3) + Ifges(6,4) * t137 + Ifges(6,6) * t143 + Ifges(6,2) * t136 / 0.2e1) * t136 + (-t110 * mrSges(7,1) + t109 * mrSges(7,3) + Ifges(7,4) * t129 + Ifges(7,6) * t134 + Ifges(7,2) * t128 / 0.2e1) * t128 + (t130 * mrSges(4,1) - t131 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (t122 * (-mrSges(5,1) * t160 + mrSges(5,2) * t156) + (Ifges(5,2) * t160 ^ 2 / 0.2e1 + (Ifges(5,4) * t160 + Ifges(5,1) * t156 / 0.2e1) * t156) * t168 + (-t118 * t156 + t119 * t160) * mrSges(5,3) + t144 * (Ifges(5,5) * t156 + Ifges(5,6) * t160)) * t148) * qJD(3);
T  = t1;
