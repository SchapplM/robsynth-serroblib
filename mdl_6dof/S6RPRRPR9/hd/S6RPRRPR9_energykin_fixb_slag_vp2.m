% Calculate kinetic energy for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:26:55
% EndTime: 2019-03-09 05:26:56
% DurationCPUTime: 0.76s
% Computational Cost: add. (1790->124), mult. (5642->205), div. (0->0), fcn. (4776->14), ass. (0->58)
t146 = sin(pkin(12));
t148 = sin(pkin(6));
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t150 = cos(pkin(12));
t151 = cos(pkin(7));
t166 = t150 * t151;
t147 = sin(pkin(7));
t152 = cos(pkin(6));
t168 = t147 * t152;
t133 = (t148 * (-t146 * t155 + t158 * t166) + t158 * t168) * qJD(1);
t164 = qJD(1) * t148;
t161 = qJ(2) * t164;
t163 = pkin(1) * qJD(1) * t152;
t140 = t146 * t163 + t150 * t161;
t130 = (t148 * t166 + t168) * qJD(1) * pkin(9) + t140;
t143 = t150 * t163;
t132 = t143 + (pkin(2) * t152 + (-pkin(9) * t151 - qJ(2)) * t148 * t146) * qJD(1);
t137 = qJD(2) + (-pkin(9) * t146 * t147 - pkin(2) * t150 - pkin(1)) * t164;
t122 = -t155 * t130 + (t132 * t151 + t137 * t147) * t158;
t167 = t147 * t155;
t165 = t151 * t155;
t124 = -t132 * t147 + t151 * t137;
t134 = (t152 * t167 + (t146 * t158 + t150 * t165) * t148) * qJD(1);
t115 = -pkin(3) * t133 - pkin(10) * t134 + t124;
t123 = t158 * t130 + t132 * t165 + t137 * t167;
t138 = qJD(3) + (-t147 * t148 * t150 + t151 * t152) * qJD(1);
t118 = pkin(10) * t138 + t123;
t154 = sin(qJ(4));
t157 = cos(qJ(4));
t108 = t157 * t115 - t118 * t154;
t126 = t134 * t157 + t138 * t154;
t131 = qJD(4) - t133;
t105 = pkin(4) * t131 - qJ(5) * t126 + t108;
t109 = t154 * t115 + t157 * t118;
t125 = -t134 * t154 + t138 * t157;
t107 = qJ(5) * t125 + t109;
t145 = sin(pkin(13));
t149 = cos(pkin(13));
t102 = t145 * t105 + t149 * t107;
t101 = t105 * t149 - t107 * t145;
t120 = t125 * t149 - t126 * t145;
t117 = -t138 * pkin(3) - t122;
t110 = -t125 * pkin(4) + qJD(5) + t117;
t156 = cos(qJ(6));
t153 = sin(qJ(6));
t144 = -pkin(1) * t164 + qJD(2);
t139 = -t146 * t161 + t143;
t121 = t125 * t145 + t126 * t149;
t119 = qJD(6) - t120;
t112 = t121 * t156 + t131 * t153;
t111 = -t121 * t153 + t131 * t156;
t103 = -t120 * pkin(5) - t121 * pkin(11) + t110;
t100 = pkin(11) * t131 + t102;
t99 = -pkin(5) * t131 - t101;
t98 = t100 * t156 + t103 * t153;
t97 = -t100 * t153 + t103 * t156;
t1 = m(3) * (t139 ^ 2 + t140 ^ 2 + t144 ^ 2) / 0.2e1 + m(4) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + (t122 * mrSges(4,1) - t123 * mrSges(4,2) + Ifges(4,3) * t138 / 0.2e1) * t138 + (t117 * mrSges(5,2) - t108 * mrSges(5,3) + Ifges(5,1) * t126 / 0.2e1) * t126 + (t110 * mrSges(6,2) - t101 * mrSges(6,3) + Ifges(6,1) * t121 / 0.2e1) * t121 + (t97 * mrSges(7,1) - t98 * mrSges(7,2) + Ifges(7,3) * t119 / 0.2e1) * t119 + (t124 * mrSges(4,2) - t122 * mrSges(4,3) + Ifges(4,5) * t138 + Ifges(4,1) * t134 / 0.2e1) * t134 + (-t117 * mrSges(5,1) + t109 * mrSges(5,3) + Ifges(5,4) * t126 + Ifges(5,2) * t125 / 0.2e1) * t125 + (-t110 * mrSges(6,1) + t102 * mrSges(6,3) + Ifges(6,4) * t121 + Ifges(6,2) * t120 / 0.2e1) * t120 + (t99 * mrSges(7,2) - t97 * mrSges(7,3) + Ifges(7,5) * t119 + Ifges(7,1) * t112 / 0.2e1) * t112 + (-t124 * mrSges(4,1) + t123 * mrSges(4,3) + Ifges(4,4) * t134 + Ifges(4,6) * t138 + Ifges(4,2) * t133 / 0.2e1) * t133 + (-t99 * mrSges(7,1) + t98 * mrSges(7,3) + Ifges(7,4) * t112 + Ifges(7,6) * t119 + Ifges(7,2) * t111 / 0.2e1) * t111 + (t108 * mrSges(5,1) + t101 * mrSges(6,1) - t109 * mrSges(5,2) - t102 * mrSges(6,2) + Ifges(5,5) * t126 + Ifges(6,5) * t121 + Ifges(5,6) * t125 + Ifges(6,6) * t120 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t131) * t131 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t144 * (-mrSges(3,1) * t150 + mrSges(3,2) * t146) + (Ifges(3,2) * t150 ^ 2 / 0.2e1 + (Ifges(3,4) * t150 + Ifges(3,1) * t146 / 0.2e1) * t146) * t164 + (-t139 * t146 + t140 * t150) * mrSges(3,3)) * t148 + (t139 * mrSges(3,1) - t140 * mrSges(3,2) + (Ifges(3,3) * t152 / 0.2e1 + (Ifges(3,5) * t146 + Ifges(3,6) * t150) * t148) * qJD(1)) * t152) * qJD(1);
T  = t1;
