% Calculate kinetic energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:23
% EndTime: 2019-03-08 18:45:23
% DurationCPUTime: 0.30s
% Computational Cost: add. (412->82), mult. (975->134), div. (0->0), fcn. (776->14), ass. (0->44)
t129 = cos(pkin(6)) * qJD(1) + qJD(2);
t134 = sin(pkin(7));
t138 = cos(pkin(7));
t137 = cos(pkin(12));
t135 = sin(pkin(6));
t153 = qJD(1) * t135;
t149 = t137 * t153;
t155 = t129 * t134 + t138 * t149;
t141 = sin(qJ(3));
t144 = cos(qJ(3));
t133 = sin(pkin(12));
t150 = t133 * t153;
t117 = -t141 * t150 + t155 * t144;
t118 = t155 * t141 + t144 * t150;
t116 = qJD(3) * pkin(9) + t118;
t122 = t129 * t138 - t134 * t149;
t140 = sin(qJ(4));
t143 = cos(qJ(4));
t109 = t143 * t116 + t140 * t122;
t107 = qJD(4) * qJ(5) + t109;
t112 = (-pkin(4) * t143 - qJ(5) * t140 - pkin(3)) * qJD(3) - t117;
t132 = sin(pkin(13));
t136 = cos(pkin(13));
t103 = t136 * t107 + t132 * t112;
t152 = qJD(3) * t140;
t151 = qJD(3) * t143;
t102 = -t107 * t132 + t136 * t112;
t108 = -t140 * t116 + t122 * t143;
t106 = -qJD(4) * pkin(4) + qJD(5) - t108;
t145 = qJD(1) ^ 2;
t142 = cos(qJ(6));
t139 = sin(qJ(6));
t130 = qJD(6) - t151;
t125 = qJD(4) * t132 + t136 * t152;
t124 = qJD(4) * t136 - t132 * t152;
t120 = t124 * t139 + t125 * t142;
t119 = t124 * t142 - t125 * t139;
t115 = -qJD(3) * pkin(3) - t117;
t104 = -pkin(5) * t124 + t106;
t101 = pkin(10) * t124 + t103;
t100 = -pkin(5) * t151 - pkin(10) * t125 + t102;
t99 = t100 * t139 + t101 * t142;
t98 = t100 * t142 - t101 * t139;
t1 = m(3) * (t129 ^ 2 + (t133 ^ 2 + t137 ^ 2) * t145 * t135 ^ 2) / 0.2e1 + m(2) * t145 / 0.2e1 + m(4) * (t117 ^ 2 + t118 ^ 2 + t122 ^ 2) / 0.2e1 + m(5) * (t108 ^ 2 + t109 ^ 2 + t115 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t106 ^ 2) / 0.2e1 + (t98 * mrSges(7,1) - t99 * mrSges(7,2) + Ifges(7,3) * t130 / 0.2e1) * t130 + (t106 * mrSges(6,2) - t102 * mrSges(6,3) + Ifges(6,1) * t125 / 0.2e1) * t125 + (t108 * mrSges(5,1) - t109 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (-t106 * mrSges(6,1) + t103 * mrSges(6,3) + Ifges(6,4) * t125 + Ifges(6,2) * t124 / 0.2e1) * t124 + (t104 * mrSges(7,2) - t98 * mrSges(7,3) + Ifges(7,5) * t130 + Ifges(7,1) * t120 / 0.2e1) * t120 + (-t104 * mrSges(7,1) + t99 * mrSges(7,3) + Ifges(7,4) * t120 + Ifges(7,6) * t130 + Ifges(7,2) * t119 / 0.2e1) * t119 + (t117 * mrSges(4,1) - t118 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (t115 * mrSges(5,2) - t108 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t152 / 0.2e1) * t140 + (-t115 * mrSges(5,1) - t102 * mrSges(6,1) + t103 * mrSges(6,2) + t109 * mrSges(5,3) - Ifges(6,5) * t125 + Ifges(5,6) * qJD(4) - Ifges(6,6) * t124 + (Ifges(5,4) * t140 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t143) * qJD(3)) * t143) * qJD(3);
T  = t1;
