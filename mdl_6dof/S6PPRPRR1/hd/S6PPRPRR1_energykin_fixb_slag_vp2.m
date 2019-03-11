% Calculate kinetic energy for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRPRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:41:56
% EndTime: 2019-03-08 18:41:56
% DurationCPUTime: 0.23s
% Computational Cost: add. (321->69), mult. (779->117), div. (0->0), fcn. (635->14), ass. (0->40)
t121 = cos(pkin(6)) * qJD(1) + qJD(2);
t126 = sin(pkin(7));
t143 = t121 * t126;
t125 = sin(pkin(12));
t130 = cos(pkin(7));
t133 = sin(qJ(3));
t136 = cos(qJ(3));
t129 = cos(pkin(12));
t127 = sin(pkin(6));
t142 = qJD(1) * t127;
t139 = t129 * t142;
t111 = -t133 * t125 * t142 + (t130 * t139 + t143) * t136;
t110 = qJD(3) * pkin(3) + t111;
t112 = t133 * t143 + (t129 * t130 * t133 + t125 * t136) * t142;
t124 = sin(pkin(13));
t128 = cos(pkin(13));
t106 = t124 * t110 + t128 * t112;
t104 = qJD(3) * pkin(9) + t106;
t115 = t130 * t121 - t126 * t139;
t114 = qJD(4) + t115;
t132 = sin(qJ(5));
t135 = cos(qJ(5));
t100 = t135 * t104 + t132 * t114;
t141 = qJD(3) * t132;
t140 = t135 * qJD(3);
t105 = t128 * t110 - t124 * t112;
t99 = -t132 * t104 + t135 * t114;
t138 = qJD(1) ^ 2;
t134 = cos(qJ(6));
t131 = sin(qJ(6));
t122 = qJD(6) - t140;
t119 = t131 * qJD(5) + t134 * t141;
t118 = t134 * qJD(5) - t131 * t141;
t103 = -qJD(3) * pkin(4) - t105;
t101 = (-pkin(5) * t135 - pkin(10) * t132 - pkin(4)) * qJD(3) - t105;
t98 = qJD(5) * pkin(10) + t100;
t97 = -qJD(5) * pkin(5) - t99;
t96 = t131 * t101 + t134 * t98;
t95 = t134 * t101 - t131 * t98;
t1 = m(3) * (t121 ^ 2 + (t125 ^ 2 + t129 ^ 2) * t138 * t127 ^ 2) / 0.2e1 + m(7) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t106 ^ 2 + t114 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t103 ^ 2 + t99 ^ 2) / 0.2e1 + m(4) * (t111 ^ 2 + t112 ^ 2 + t115 ^ 2) / 0.2e1 + m(2) * t138 / 0.2e1 + (t95 * mrSges(7,1) - t96 * mrSges(7,2) + Ifges(7,3) * t122 / 0.2e1) * t122 + (t99 * mrSges(6,1) - t100 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t97 * mrSges(7,2) - t95 * mrSges(7,3) + Ifges(7,5) * t122 + Ifges(7,1) * t119 / 0.2e1) * t119 + (-t97 * mrSges(7,1) + t96 * mrSges(7,3) + Ifges(7,4) * t119 + Ifges(7,6) * t122 + Ifges(7,2) * t118 / 0.2e1) * t118 + (t111 * mrSges(4,1) + t105 * mrSges(5,1) - t112 * mrSges(4,2) - t106 * mrSges(5,2) + (-t103 * mrSges(6,1) + t100 * mrSges(6,3) + Ifges(6,6) * qJD(5) + Ifges(6,2) * t140 / 0.2e1) * t135 + (t103 * mrSges(6,2) - t99 * mrSges(6,3) + Ifges(6,5) * qJD(5)) * t132 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(6,4) * t135 + Ifges(6,1) * t132 / 0.2e1) * t132) * qJD(3)) * qJD(3);
T  = t1;
