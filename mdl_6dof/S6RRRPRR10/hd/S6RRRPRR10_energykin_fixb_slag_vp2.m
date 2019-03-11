% Calculate kinetic energy for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:32
% EndTime: 2019-03-09 19:15:32
% DurationCPUTime: 0.57s
% Computational Cost: add. (654->108), mult. (1297->156), div. (0->0), fcn. (868->8), ass. (0->40)
t126 = pkin(7) * mrSges(3,3);
t125 = cos(qJ(3));
t115 = sin(qJ(5));
t119 = cos(qJ(5));
t116 = sin(qJ(3));
t117 = sin(qJ(2));
t124 = t117 * qJD(1);
t101 = t116 * qJD(2) + t124 * t125;
t120 = cos(qJ(2));
t112 = t120 * qJD(1);
t109 = -t112 + qJD(3);
t106 = pkin(7) * t112 + qJD(2) * pkin(8);
t97 = (-pkin(2) * t120 - pkin(8) * t117 - pkin(1)) * qJD(1);
t94 = -t116 * t106 + t125 * t97;
t123 = qJD(4) - t94;
t85 = -t101 * pkin(9) + (-pkin(3) - pkin(4)) * t109 + t123;
t100 = -qJD(2) * t125 + t116 * t124;
t95 = t125 * t106 + t116 * t97;
t90 = t109 * qJ(4) + t95;
t87 = pkin(9) * t100 + t90;
t79 = t115 * t85 + t119 * t87;
t105 = -qJD(2) * pkin(2) + pkin(7) * t124;
t108 = qJD(5) - t109;
t78 = -t115 * t87 + t119 * t85;
t91 = t100 * pkin(3) - t101 * qJ(4) + t105;
t88 = -pkin(4) * t100 - t91;
t118 = cos(qJ(6));
t114 = sin(qJ(6));
t104 = qJD(6) + t108;
t93 = t100 * t115 + t101 * t119;
t92 = t100 * t119 - t101 * t115;
t89 = -t109 * pkin(3) + t123;
t82 = t114 * t92 + t118 * t93;
t81 = -t114 * t93 + t118 * t92;
t80 = -pkin(5) * t92 + t88;
t77 = pkin(10) * t92 + t79;
t76 = pkin(5) * t108 - pkin(10) * t93 + t78;
t75 = t114 * t76 + t118 * t77;
t74 = -t114 * t77 + t118 * t76;
t1 = m(7) * (t74 ^ 2 + t75 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t88 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t105 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + (t88 * mrSges(6,2) - t78 * mrSges(6,3) + Ifges(6,1) * t93 / 0.2e1) * t93 + (t80 * mrSges(7,2) - t74 * mrSges(7,3) + Ifges(7,1) * t82 / 0.2e1) * t82 + (-t88 * mrSges(6,1) + t79 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,2) * t92 / 0.2e1) * t92 + (-t80 * mrSges(7,1) + t75 * mrSges(7,3) + Ifges(7,4) * t82 + Ifges(7,2) * t81 / 0.2e1) * t81 + (t78 * mrSges(6,1) - t79 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(6,6) * t92 + Ifges(6,3) * t108 / 0.2e1) * t108 + (t74 * mrSges(7,1) - t75 * mrSges(7,2) + Ifges(7,5) * t82 + Ifges(7,6) * t81 + Ifges(7,3) * t104 / 0.2e1) * t104 + (t94 * mrSges(4,1) - t89 * mrSges(5,1) - t95 * mrSges(4,2) + t90 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t109) * t109 + (t105 * mrSges(4,2) + t89 * mrSges(5,2) - t94 * mrSges(4,3) - t91 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t101 + (Ifges(5,4) + Ifges(4,5)) * t109) * t101 + (t105 * mrSges(4,1) + t91 * mrSges(5,1) - t90 * mrSges(5,2) - t95 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t100 + (-Ifges(4,6) + Ifges(5,6)) * t109 + (-Ifges(4,4) + Ifges(5,5)) * t101) * t100 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t117 ^ 2 + t120 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t126 + Ifges(3,2) / 0.2e1) * t120) * t120 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t120 + (t126 + Ifges(3,1) / 0.2e1) * t117) * t117) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t120 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t117) * qJD(2)) * qJD(1);
T  = t1;
