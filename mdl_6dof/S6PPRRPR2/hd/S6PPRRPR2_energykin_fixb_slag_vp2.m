% Calculate kinetic energy for
% S6PPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-03-08 18:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:49:07
% EndTime: 2019-03-08 18:49:07
% DurationCPUTime: 0.25s
% Computational Cost: add. (272->80), mult. (633->120), div. (0->0), fcn. (464->12), ass. (0->41)
t116 = cos(pkin(6)) * qJD(1) + qJD(2);
t120 = sin(pkin(7));
t123 = cos(pkin(7));
t122 = cos(pkin(12));
t121 = sin(pkin(6));
t140 = qJD(1) * t121;
t136 = t122 * t140;
t143 = t116 * t120 + t123 * t136;
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t119 = sin(pkin(12));
t137 = t119 * t140;
t106 = -t126 * t137 + t143 * t129;
t142 = -pkin(4) - pkin(10);
t107 = t143 * t126 + t129 * t137;
t105 = qJD(3) * pkin(9) + t107;
t109 = t123 * t116 - t120 * t136;
t125 = sin(qJ(4));
t128 = cos(qJ(4));
t99 = t128 * t105 + t125 * t109;
t139 = qJD(3) * t128;
t138 = t125 * qJD(3);
t135 = -qJ(5) * t125 - pkin(3);
t98 = -t125 * t105 + t128 * t109;
t97 = -qJD(4) * qJ(5) - t99;
t133 = qJD(5) - t98;
t130 = qJD(1) ^ 2;
t127 = cos(qJ(6));
t124 = sin(qJ(6));
t117 = qJD(6) + t138;
t113 = t127 * qJD(4) - t124 * t139;
t112 = -t124 * qJD(4) - t127 * t139;
t104 = -qJD(3) * pkin(3) - t106;
t101 = (-pkin(4) * t128 + t135) * qJD(3) - t106;
t100 = (t142 * t128 + t135) * qJD(3) - t106;
t96 = -qJD(4) * pkin(4) + t133;
t95 = pkin(5) * t139 - t97;
t94 = pkin(5) * t138 + t142 * qJD(4) + t133;
t93 = t127 * t100 + t124 * t94;
t92 = -t124 * t100 + t127 * t94;
t1 = m(3) * (t116 ^ 2 + (t119 ^ 2 + t122 ^ 2) * t130 * t121 ^ 2) / 0.2e1 + m(2) * t130 / 0.2e1 + m(6) * (t101 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t106 ^ 2 + t107 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t92 ^ 2 + t93 ^ 2 + t95 ^ 2) / 0.2e1 + (t92 * mrSges(7,1) - t93 * mrSges(7,2) + Ifges(7,3) * t117 / 0.2e1) * t117 + (t95 * mrSges(7,2) - t92 * mrSges(7,3) + Ifges(7,5) * t117 + Ifges(7,1) * t113 / 0.2e1) * t113 + (-t95 * mrSges(7,1) + t93 * mrSges(7,3) + Ifges(7,4) * t113 + Ifges(7,6) * t117 + Ifges(7,2) * t112 / 0.2e1) * t112 + (t98 * mrSges(5,1) - t99 * mrSges(5,2) + t96 * mrSges(6,2) - t97 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(4)) * qJD(4) + (t106 * mrSges(4,1) - t107 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t104 * mrSges(5,1) - t97 * mrSges(6,1) + t101 * mrSges(6,2) + t99 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t139 + (-Ifges(6,5) + Ifges(5,6)) * qJD(4)) * t128 + (t96 * mrSges(6,1) + t104 * mrSges(5,2) - t98 * mrSges(5,3) - t101 * mrSges(6,3) + (-Ifges(6,4) + Ifges(5,5)) * qJD(4) + ((Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t125 + (Ifges(5,4) + Ifges(6,6)) * t128) * qJD(3)) * t125) * qJD(3);
T  = t1;
