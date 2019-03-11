% Calculate kinetic energy for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:44
% EndTime: 2019-03-09 03:36:44
% DurationCPUTime: 0.53s
% Computational Cost: add. (538->98), mult. (1185->156), div. (0->0), fcn. (798->10), ass. (0->42)
t126 = m(3) / 0.2e1;
t117 = sin(qJ(5));
t120 = cos(qJ(5));
t112 = sin(pkin(11));
t114 = cos(pkin(11));
t113 = sin(pkin(10));
t107 = (pkin(1) * t113 + pkin(7)) * qJD(1);
t121 = cos(qJ(3));
t111 = t121 * qJD(2);
t118 = sin(qJ(3));
t95 = qJD(3) * pkin(3) + t111 + (-qJ(4) * qJD(1) - t107) * t118;
t101 = t118 * qJD(2) + t121 * t107;
t125 = qJD(1) * t121;
t98 = qJ(4) * t125 + t101;
t86 = t112 * t95 + t114 * t98;
t84 = qJD(3) * pkin(8) + t86;
t115 = cos(pkin(10));
t124 = -pkin(1) * t115 - pkin(2);
t103 = qJD(4) + (-pkin(3) * t121 + t124) * qJD(1);
t104 = -t112 * t118 * qJD(1) + t114 * t125;
t105 = (t112 * t121 + t114 * t118) * qJD(1);
t91 = -pkin(4) * t104 - pkin(8) * t105 + t103;
t80 = t117 * t91 + t120 * t84;
t85 = -t112 * t98 + t114 * t95;
t79 = -t117 * t84 + t120 * t91;
t83 = -qJD(3) * pkin(4) - t85;
t102 = qJD(5) - t104;
t119 = cos(qJ(6));
t116 = sin(qJ(6));
t108 = t124 * qJD(1);
t100 = -t107 * t118 + t111;
t99 = qJD(6) + t102;
t97 = qJD(3) * t117 + t105 * t120;
t96 = qJD(3) * t120 - t105 * t117;
t88 = t116 * t96 + t119 * t97;
t87 = -t116 * t97 + t119 * t96;
t81 = -pkin(5) * t96 + t83;
t78 = pkin(9) * t96 + t80;
t77 = pkin(5) * t102 - pkin(9) * t97 + t79;
t76 = t116 * t77 + t119 * t78;
t75 = -t116 * t78 + t119 * t77;
t1 = m(5) * (t103 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t126 + m(4) * (t100 ^ 2 + t101 ^ 2 + t108 ^ 2) / 0.2e1 + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + Ifges(7,3) * t99 / 0.2e1) * t99 + (t83 * mrSges(6,2) - t79 * mrSges(6,3) + Ifges(6,1) * t97 / 0.2e1) * t97 + (t103 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t105 / 0.2e1) * t105 + (-t83 * mrSges(6,1) + t80 * mrSges(6,3) + Ifges(6,4) * t97 + Ifges(6,2) * t96 / 0.2e1) * t96 + (t81 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,5) * t99 + Ifges(7,1) * t88 / 0.2e1) * t88 + (-t103 * mrSges(5,1) + t86 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,2) * t104 / 0.2e1) * t104 + (-t81 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t88 + Ifges(7,6) * t99 + Ifges(7,2) * t87 / 0.2e1) * t87 + (t79 * mrSges(6,1) - t80 * mrSges(6,2) + Ifges(6,5) * t97 + Ifges(6,6) * t96 + Ifges(6,3) * t102 / 0.2e1) * t102 + (t100 * mrSges(4,1) + t85 * mrSges(5,1) - t101 * mrSges(4,2) - t86 * mrSges(5,2) + Ifges(5,5) * t105 + Ifges(5,6) * t104 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3)) * qJD(3) + (t108 * (-mrSges(4,1) * t121 + mrSges(4,2) * t118) + (-t100 * t118 + t101 * t121) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t118 + Ifges(4,6) * t121) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t115 * mrSges(3,1) - t113 * mrSges(3,2) + (t113 ^ 2 + t115 ^ 2) * t126 * pkin(1)) * pkin(1) + Ifges(4,2) * t121 ^ 2 / 0.2e1 + (Ifges(4,4) * t121 + Ifges(4,1) * t118 / 0.2e1) * t118) * qJD(1)) * qJD(1);
T  = t1;
