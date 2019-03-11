% Calculate kinetic energy for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:50:23
% EndTime: 2019-03-09 13:50:24
% DurationCPUTime: 0.55s
% Computational Cost: add. (590->108), mult. (1211->159), div. (0->0), fcn. (788->8), ass. (0->40)
t125 = pkin(7) * mrSges(3,3);
t113 = sin(qJ(5));
t117 = cos(qJ(5));
t109 = -qJD(2) + qJD(4);
t114 = sin(qJ(4));
t118 = cos(qJ(4));
t115 = sin(qJ(2));
t124 = qJD(1) * t115;
t122 = pkin(7) * t124 + qJD(3);
t95 = -pkin(8) * t124 + (-pkin(2) - pkin(3)) * qJD(2) + t122;
t119 = cos(qJ(2));
t123 = qJD(1) * t119;
t102 = pkin(7) * t123 + qJD(2) * qJ(3);
t99 = -pkin(8) * t123 + t102;
t86 = -t114 * t99 + t118 * t95;
t98 = (-t114 * t119 + t115 * t118) * qJD(1);
t81 = pkin(4) * t109 - pkin(9) * t98 + t86;
t87 = t114 * t95 + t118 * t99;
t97 = (-t114 * t115 - t118 * t119) * qJD(1);
t83 = pkin(9) * t97 + t87;
t78 = t113 * t81 + t117 * t83;
t100 = -qJD(1) * pkin(1) - pkin(2) * t123 - qJ(3) * t124;
t94 = pkin(3) * t123 - t100;
t77 = -t113 * t83 + t117 * t81;
t89 = -t113 * t98 + t117 * t97;
t91 = -pkin(4) * t97 + t94;
t116 = cos(qJ(6));
t112 = sin(qJ(6));
t108 = qJD(5) + t109;
t101 = -qJD(2) * pkin(2) + t122;
t90 = t113 * t97 + t117 * t98;
t88 = qJD(6) - t89;
t85 = t108 * t112 + t116 * t90;
t84 = t108 * t116 - t112 * t90;
t79 = -pkin(5) * t89 - pkin(10) * t90 + t91;
t76 = pkin(10) * t108 + t78;
t75 = -pkin(5) * t108 - t77;
t74 = t112 * t79 + t116 * t76;
t73 = -t112 * t76 + t116 * t79;
t1 = m(4) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + (t94 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,1) * t98 / 0.2e1) * t98 + (t91 * mrSges(6,2) - t77 * mrSges(6,3) + Ifges(6,1) * t90 / 0.2e1) * t90 + (t73 * mrSges(7,1) - t74 * mrSges(7,2) + Ifges(7,3) * t88 / 0.2e1) * t88 + (-t94 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,4) * t98 + Ifges(5,2) * t97 / 0.2e1) * t97 + (-t91 * mrSges(6,1) + t78 * mrSges(6,3) + Ifges(6,4) * t90 + Ifges(6,2) * t89 / 0.2e1) * t89 + (t75 * mrSges(7,2) - t73 * mrSges(7,3) + Ifges(7,5) * t88 + Ifges(7,1) * t85 / 0.2e1) * t85 + (-t75 * mrSges(7,1) + t74 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,6) * t88 + Ifges(7,2) * t84 / 0.2e1) * t84 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,5) * t98 + Ifges(5,6) * t97 + Ifges(5,3) * t109 / 0.2e1) * t109 + (t77 * mrSges(6,1) - t78 * mrSges(6,2) + Ifges(6,5) * t90 + Ifges(6,6) * t89 + Ifges(6,3) * t108 / 0.2e1) * t108 + (-t101 * mrSges(4,1) + t102 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t100 * mrSges(4,1) + t102 * mrSges(4,2) + (-pkin(7) * mrSges(3,2) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t119 + (t101 * mrSges(4,2) - t100 * mrSges(4,3) + (-pkin(7) * mrSges(3,1) + Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t115 + (m(3) * (pkin(1) ^ 2 + (t115 ^ 2 + t119 ^ 2) * pkin(7) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + t125) * t119) * t119 + (-pkin(1) * mrSges(3,2) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + t125) * t115 + (Ifges(3,4) - Ifges(4,5)) * t119) * t115) * qJD(1)) * qJD(1);
T  = t1;
