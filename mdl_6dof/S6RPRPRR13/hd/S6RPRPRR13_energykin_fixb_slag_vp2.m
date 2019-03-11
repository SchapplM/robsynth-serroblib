% Calculate kinetic energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:32
% EndTime: 2019-03-09 04:21:32
% DurationCPUTime: 0.66s
% Computational Cost: add. (1090->120), mult. (3446->187), div. (0->0), fcn. (2832->12), ass. (0->55)
t136 = sin(pkin(7));
t139 = cos(pkin(7));
t140 = cos(pkin(6));
t137 = sin(pkin(6));
t138 = cos(pkin(12));
t154 = t137 * t138;
t160 = qJD(1) * (t136 * t140 + t139 * t154);
t158 = pkin(3) + pkin(10);
t157 = cos(qJ(3));
t143 = sin(qJ(3));
t135 = sin(pkin(12));
t152 = qJD(1) * t137;
t149 = t135 * t152;
t120 = t143 * t149 - t157 * t160;
t151 = pkin(1) * qJD(1) * t140;
t133 = t138 * t151;
t119 = t133 + (pkin(2) * t140 + (-pkin(9) * t139 - qJ(2)) * t137 * t135) * qJD(1);
t124 = qJD(2) + (-pkin(9) * t135 * t136 - pkin(2) * t138 - pkin(1)) * t152;
t110 = -t119 * t136 + t139 * t124;
t153 = t139 * t143;
t155 = t136 * t143;
t121 = (t140 * t155 + (t157 * t135 + t138 * t153) * t137) * qJD(1);
t147 = -qJ(4) * t121 + t110;
t101 = t158 * t120 + t147;
t142 = sin(qJ(5));
t145 = cos(qJ(5));
t126 = qJD(3) + (-t136 * t154 + t139 * t140) * qJD(1);
t128 = t138 * qJ(2) * t152 + t135 * t151;
t117 = pkin(9) * t160 + t128;
t106 = -t143 * t117 + (t119 * t139 + t124 * t136) * t157;
t146 = qJD(4) - t106;
t99 = t121 * pkin(4) - t158 * t126 + t146;
t96 = t145 * t101 + t142 * t99;
t107 = t157 * t117 + t119 * t153 + t124 * t155;
t105 = -t126 * qJ(4) - t107;
t95 = -t142 * t101 + t145 * t99;
t112 = t120 * t145 - t142 * t126;
t102 = -pkin(4) * t120 - t105;
t144 = cos(qJ(6));
t141 = sin(qJ(6));
t134 = -pkin(1) * t152 + qJD(2);
t127 = -qJ(2) * t149 + t133;
t118 = qJD(5) + t121;
t113 = t142 * t120 + t126 * t145;
t111 = qJD(6) - t112;
t109 = t113 * t144 + t118 * t141;
t108 = -t113 * t141 + t118 * t144;
t104 = -t126 * pkin(3) + t146;
t103 = pkin(3) * t120 + t147;
t97 = -pkin(5) * t112 - pkin(11) * t113 + t102;
t94 = pkin(11) * t118 + t96;
t93 = -t118 * pkin(5) - t95;
t92 = t141 * t97 + t144 * t94;
t91 = -t141 * t94 + t144 * t97;
t1 = m(4) * (t106 ^ 2 + t107 ^ 2 + t110 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(3) * (t127 ^ 2 + t128 ^ 2 + t134 ^ 2) / 0.2e1 + m(7) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + (t95 * mrSges(6,1) - t96 * mrSges(6,2) + Ifges(6,3) * t118 / 0.2e1) * t118 + (t91 * mrSges(7,1) - t92 * mrSges(7,2) + Ifges(7,3) * t111 / 0.2e1) * t111 + (t102 * mrSges(6,2) - t95 * mrSges(6,3) + Ifges(6,5) * t118 + Ifges(6,1) * t113 / 0.2e1) * t113 + (t93 * mrSges(7,2) - t91 * mrSges(7,3) + Ifges(7,5) * t111 + Ifges(7,1) * t109 / 0.2e1) * t109 + (-t102 * mrSges(6,1) + t96 * mrSges(6,3) + Ifges(6,4) * t113 + Ifges(6,6) * t118 + Ifges(6,2) * t112 / 0.2e1) * t112 + (-t93 * mrSges(7,1) + t92 * mrSges(7,3) + Ifges(7,4) * t109 + Ifges(7,6) * t111 + Ifges(7,2) * t108 / 0.2e1) * t108 + (t106 * mrSges(4,1) - t107 * mrSges(4,2) + t104 * mrSges(5,2) - t105 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t126) * t126 + (t104 * mrSges(5,1) + t110 * mrSges(4,2) - t106 * mrSges(4,3) - t103 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t121 + (-Ifges(5,4) + Ifges(4,5)) * t126) * t121 + (t110 * mrSges(4,1) + t105 * mrSges(5,1) - t103 * mrSges(5,2) - t107 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t120 + (Ifges(5,5) - Ifges(4,6)) * t126 + (-Ifges(4,4) - Ifges(5,6)) * t121) * t120 + (Ifges(2,3) * qJD(1) / 0.2e1 + (t134 * (-mrSges(3,1) * t138 + mrSges(3,2) * t135) + (Ifges(3,2) * t138 ^ 2 / 0.2e1 + (Ifges(3,4) * t138 + Ifges(3,1) * t135 / 0.2e1) * t135) * t152 + (-t127 * t135 + t128 * t138) * mrSges(3,3)) * t137 + (t127 * mrSges(3,1) - t128 * mrSges(3,2) + (Ifges(3,3) * t140 / 0.2e1 + (Ifges(3,5) * t135 + Ifges(3,6) * t138) * t137) * qJD(1)) * t140) * qJD(1);
T  = t1;
