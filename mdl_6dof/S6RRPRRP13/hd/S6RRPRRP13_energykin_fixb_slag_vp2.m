% Calculate kinetic energy for
% S6RRPRRP13
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
% Datum: 2019-03-09 13:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:55:29
% EndTime: 2019-03-09 12:55:29
% DurationCPUTime: 0.44s
% Computational Cost: add. (622->109), mult. (1380->154), div. (0->0), fcn. (982->8), ass. (0->43)
t125 = -pkin(2) - pkin(9);
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t114 = sin(qJ(2));
t110 = sin(pkin(6));
t124 = t110 * qJD(1);
t120 = t114 * t124;
t103 = qJD(4) + t120;
t113 = sin(qJ(4));
t116 = cos(qJ(4));
t105 = pkin(8) * t120;
t111 = cos(pkin(6));
t123 = t111 * qJD(1);
t109 = qJD(2) + t123;
t117 = cos(qJ(2));
t88 = qJD(3) + t105 + t125 * t109 + (-pkin(1) * t111 * t117 + pkin(3) * t110 * t114) * qJD(1);
t119 = -qJ(3) * t114 - pkin(1);
t94 = (t125 * t117 + t119) * t124;
t83 = t113 * t88 + t116 * t94;
t81 = t103 * pkin(10) + t83;
t121 = t117 * t124;
t100 = t116 * t109 - t113 * t121;
t122 = pkin(1) * t123;
t102 = pkin(8) * t121 + t114 * t122;
t96 = -t109 * qJ(3) - t102;
t93 = pkin(3) * t121 - t96;
t99 = -t113 * t109 - t116 * t121;
t86 = -t99 * pkin(4) - t100 * pkin(10) + t93;
t77 = t112 * t86 + t115 * t81;
t76 = -t112 * t81 + t115 * t86;
t82 = -t113 * t94 + t116 * t88;
t101 = t117 * t122 - t105;
t80 = -t103 * pkin(4) - t82;
t118 = qJD(1) ^ 2;
t98 = qJD(5) - t99;
t97 = (-pkin(2) * t117 + t119) * t124;
t95 = -t109 * pkin(2) + qJD(3) - t101;
t90 = t115 * t100 + t112 * t103;
t89 = -t112 * t100 + t115 * t103;
t78 = -t89 * pkin(5) + qJD(6) + t80;
t75 = t89 * qJ(6) + t77;
t74 = t98 * pkin(5) - t90 * qJ(6) + t76;
t1 = m(5) * (t82 ^ 2 + t83 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t96 ^ 2 + t97 ^ 2) / 0.2e1 + m(7) * (t74 ^ 2 + t75 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t76 ^ 2 + t77 ^ 2 + t80 ^ 2) / 0.2e1 + m(3) * (t110 ^ 2 * t118 * pkin(1) ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + t118 * Ifges(2,3) / 0.2e1 + (-t93 * mrSges(5,1) + t83 * mrSges(5,3) + Ifges(5,2) * t99 / 0.2e1) * t99 + (t82 * mrSges(5,1) - t83 * mrSges(5,2) + Ifges(5,6) * t99 + Ifges(5,3) * t103 / 0.2e1) * t103 + (t93 * mrSges(5,2) - t82 * mrSges(5,3) + Ifges(5,4) * t99 + Ifges(5,5) * t103 + Ifges(5,1) * t100 / 0.2e1) * t100 + (t76 * mrSges(6,1) + t74 * mrSges(7,1) - t77 * mrSges(6,2) - t75 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t98) * t98 + (t80 * mrSges(6,2) + t78 * mrSges(7,2) - t76 * mrSges(6,3) - t74 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t90 + (Ifges(6,5) + Ifges(7,5)) * t98) * t90 + (-t80 * mrSges(6,1) - t78 * mrSges(7,1) + t77 * mrSges(6,3) + t75 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t89 + (Ifges(6,6) + Ifges(7,6)) * t98 + (Ifges(6,4) + Ifges(7,4)) * t90) * t89 + ((-t96 * mrSges(4,1) + t97 * mrSges(4,2) + t102 * mrSges(3,3) + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t117) * t124) * t117 + (t95 * mrSges(4,1) - t101 * mrSges(3,3) - t97 * mrSges(4,3) + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t114 + (Ifges(3,4) + Ifges(4,6)) * t117) * t124) * t114) * t124 + (t101 * mrSges(3,1) - t102 * mrSges(3,2) + t95 * mrSges(4,2) - t96 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t109 + ((-Ifges(4,5) + Ifges(3,6)) * t117 + (-Ifges(4,4) + Ifges(3,5)) * t114) * t124) * t109;
T  = t1;
