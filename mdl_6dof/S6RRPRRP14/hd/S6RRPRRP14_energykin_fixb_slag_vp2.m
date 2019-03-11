% Calculate kinetic energy for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP14_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:03:40
% EndTime: 2019-03-09 13:03:40
% DurationCPUTime: 0.45s
% Computational Cost: add. (620->109), mult. (1372->154), div. (0->0), fcn. (974->8), ass. (0->43)
t125 = -pkin(2) - pkin(9);
t124 = cos(qJ(5));
t112 = sin(qJ(5));
t114 = sin(qJ(2));
t110 = sin(pkin(6));
t123 = qJD(1) * t110;
t119 = t114 * t123;
t103 = qJD(4) + t119;
t113 = sin(qJ(4));
t115 = cos(qJ(4));
t105 = pkin(8) * t119;
t111 = cos(pkin(6));
t122 = qJD(1) * t111;
t109 = qJD(2) + t122;
t116 = cos(qJ(2));
t87 = qJD(3) + t105 + t125 * t109 + (-pkin(1) * t111 * t116 + pkin(3) * t110 * t114) * qJD(1);
t118 = -qJ(3) * t114 - pkin(1);
t93 = (t116 * t125 + t118) * t123;
t83 = t113 * t87 + t115 * t93;
t81 = pkin(10) * t103 + t83;
t120 = t116 * t123;
t121 = pkin(1) * t122;
t101 = pkin(8) * t120 + t114 * t121;
t95 = -t109 * qJ(3) - t101;
t92 = pkin(3) * t120 - t95;
t98 = -t109 * t113 - t115 * t120;
t99 = t109 * t115 - t113 * t120;
t85 = -pkin(4) * t98 - pkin(10) * t99 + t92;
t77 = t112 * t85 + t124 * t81;
t82 = -t113 * t93 + t115 * t87;
t100 = t116 * t121 - t105;
t80 = -pkin(4) * t103 - t82;
t76 = -t112 * t81 + t124 * t85;
t117 = qJD(1) ^ 2;
t97 = qJD(5) - t98;
t96 = (-pkin(2) * t116 + t118) * t123;
t94 = -pkin(2) * t109 + qJD(3) - t100;
t89 = t112 * t103 + t124 * t99;
t88 = -t103 * t124 + t112 * t99;
t78 = pkin(5) * t88 - qJ(6) * t89 + t80;
t75 = qJ(6) * t97 + t77;
t74 = -t97 * pkin(5) + qJD(6) - t76;
t1 = m(7) * (t74 ^ 2 + t75 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) / 0.2e1 + m(5) * (t82 ^ 2 + t83 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t117 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t110 ^ 2 * t117 + t100 ^ 2 + t101 ^ 2) / 0.2e1 + (t92 * mrSges(5,2) - t82 * mrSges(5,3) + Ifges(5,1) * t99 / 0.2e1) * t99 + (-t92 * mrSges(5,1) + t83 * mrSges(5,3) + Ifges(5,4) * t99 + Ifges(5,2) * t98 / 0.2e1) * t98 + (t82 * mrSges(5,1) - t83 * mrSges(5,2) + Ifges(5,5) * t99 + Ifges(5,6) * t98 + Ifges(5,3) * t103 / 0.2e1) * t103 + (t76 * mrSges(6,1) - t74 * mrSges(7,1) - t77 * mrSges(6,2) + t75 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t97) * t97 + (t80 * mrSges(6,2) + t74 * mrSges(7,2) - t76 * mrSges(6,3) - t78 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t89 + (Ifges(7,4) + Ifges(6,5)) * t97) * t89 + (t80 * mrSges(6,1) + t78 * mrSges(7,1) - t75 * mrSges(7,2) - t77 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t88 + (-Ifges(6,6) + Ifges(7,6)) * t97 + (-Ifges(6,4) + Ifges(7,5)) * t89) * t88 + ((-t95 * mrSges(4,1) + t96 * mrSges(4,2) + t101 * mrSges(3,3) + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t116) * t123) * t116 + (t94 * mrSges(4,1) - t100 * mrSges(3,3) - t96 * mrSges(4,3) + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t114 + (Ifges(3,4) + Ifges(4,6)) * t116) * t123) * t114) * t123 + (t100 * mrSges(3,1) - t101 * mrSges(3,2) + t94 * mrSges(4,2) - t95 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t109 + ((-Ifges(4,5) + Ifges(3,6)) * t116 + (-Ifges(4,4) + Ifges(3,5)) * t114) * t123) * t109;
T  = t1;
