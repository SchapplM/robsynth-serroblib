% Calculate kinetic energy for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:11
% EndTime: 2019-03-08 19:03:11
% DurationCPUTime: 0.31s
% Computational Cost: add. (424->82), mult. (975->136), div. (0->0), fcn. (776->14), ass. (0->45)
t131 = cos(pkin(6)) * qJD(1) + qJD(2);
t135 = sin(pkin(7));
t138 = cos(pkin(7));
t137 = cos(pkin(13));
t136 = sin(pkin(6));
t155 = qJD(1) * t136;
t151 = t137 * t155;
t157 = t131 * t135 + t138 * t151;
t142 = sin(qJ(3));
t146 = cos(qJ(3));
t134 = sin(pkin(13));
t152 = t134 * t155;
t118 = -t142 * t152 + t157 * t146;
t119 = t157 * t142 + t146 * t152;
t117 = qJD(3) * pkin(9) + t119;
t123 = t138 * t131 - t135 * t151;
t141 = sin(qJ(4));
t145 = cos(qJ(4));
t110 = t145 * t117 + t141 * t123;
t108 = qJD(4) * pkin(10) + t110;
t113 = (-pkin(4) * t145 - pkin(10) * t141 - pkin(3)) * qJD(3) - t118;
t140 = sin(qJ(5));
t144 = cos(qJ(5));
t104 = t144 * t108 + t140 * t113;
t154 = qJD(3) * t141;
t153 = t145 * qJD(3);
t103 = -t140 * t108 + t144 * t113;
t109 = -t141 * t117 + t145 * t123;
t132 = qJD(5) - t153;
t107 = -qJD(4) * pkin(4) - t109;
t147 = qJD(1) ^ 2;
t143 = cos(qJ(6));
t139 = sin(qJ(6));
t130 = qJD(6) + t132;
t127 = t140 * qJD(4) + t144 * t154;
t126 = t144 * qJD(4) - t140 * t154;
t121 = t139 * t126 + t143 * t127;
t120 = t143 * t126 - t139 * t127;
t116 = -qJD(3) * pkin(3) - t118;
t105 = -t126 * pkin(5) + t107;
t102 = t126 * pkin(11) + t104;
t101 = t132 * pkin(5) - t127 * pkin(11) + t103;
t100 = t139 * t101 + t143 * t102;
t99 = t143 * t101 - t139 * t102;
t1 = m(3) * (t131 ^ 2 + (t134 ^ 2 + t137 ^ 2) * t147 * t136 ^ 2) / 0.2e1 + m(2) * t147 / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t116 ^ 2) / 0.2e1 + m(4) * (t118 ^ 2 + t119 ^ 2 + t123 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t105 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t103 ^ 2 + t104 ^ 2 + t107 ^ 2) / 0.2e1 + (t103 * mrSges(6,1) - t104 * mrSges(6,2) + Ifges(6,3) * t132 / 0.2e1) * t132 + (t99 * mrSges(7,1) - t100 * mrSges(7,2) + Ifges(7,3) * t130 / 0.2e1) * t130 + (t109 * mrSges(5,1) - t110 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t107 * mrSges(6,2) - t103 * mrSges(6,3) + Ifges(6,5) * t132 + Ifges(6,1) * t127 / 0.2e1) * t127 + (t105 * mrSges(7,2) - t99 * mrSges(7,3) + Ifges(7,5) * t130 + Ifges(7,1) * t121 / 0.2e1) * t121 + (-t107 * mrSges(6,1) + t104 * mrSges(6,3) + Ifges(6,4) * t127 + Ifges(6,6) * t132 + Ifges(6,2) * t126 / 0.2e1) * t126 + (-t105 * mrSges(7,1) + t100 * mrSges(7,3) + Ifges(7,4) * t121 + Ifges(7,6) * t130 + Ifges(7,2) * t120 / 0.2e1) * t120 + (t118 * mrSges(4,1) - t119 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t116 * mrSges(5,1) + t110 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t153 / 0.2e1) * t145 + (t116 * mrSges(5,2) - t109 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t145 + Ifges(5,1) * t141 / 0.2e1) * qJD(3)) * t141) * qJD(3);
T  = t1;
