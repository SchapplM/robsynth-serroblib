% Calculate kinetic energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:46
% EndTime: 2019-03-09 08:13:46
% DurationCPUTime: 0.42s
% Computational Cost: add. (362->105), mult. (721->139), div. (0->0), fcn. (364->6), ass. (0->34)
t112 = -pkin(2) - pkin(3);
t111 = pkin(7) * mrSges(3,3);
t102 = sin(qJ(2));
t104 = cos(qJ(2));
t108 = qJD(1) * t104;
t109 = qJD(1) * t102;
t85 = -qJD(1) * pkin(1) - pkin(2) * t108 - qJ(3) * t109;
t81 = pkin(3) * t108 + qJD(4) - t85;
t78 = (pkin(4) * t102 + qJ(5) * t104) * qJD(1) + t81;
t110 = pkin(7) * t109 + qJD(3);
t107 = -qJ(4) * t109 + t110;
t80 = (-qJ(5) + t112) * qJD(2) + t107;
t97 = sin(pkin(9));
t98 = cos(pkin(9));
t72 = t78 * t97 + t80 * t98;
t89 = pkin(7) * t108 + qJD(2) * qJ(3);
t71 = t78 * t98 - t80 * t97;
t84 = qJ(4) * t108 - t89;
t82 = qJD(2) * pkin(4) + qJD(5) - t84;
t103 = cos(qJ(6));
t101 = sin(qJ(6));
t90 = qJD(6) + t109;
t88 = -qJD(2) * pkin(2) + t110;
t87 = -qJD(2) * t97 - t108 * t98;
t86 = -qJD(2) * t98 + t108 * t97;
t83 = qJD(2) * t112 + t107;
t77 = t101 * t86 + t103 * t87;
t76 = -t101 * t87 + t103 * t86;
t75 = -pkin(5) * t86 + t82;
t70 = pkin(8) * t86 + t72;
t69 = pkin(5) * t109 - pkin(8) * t87 + t71;
t68 = t101 * t69 + t103 * t70;
t67 = -t101 * t70 + t103 * t69;
t1 = m(4) * (t85 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t75 ^ 2) / 0.2e1 + (t67 * mrSges(7,1) - t68 * mrSges(7,2) + Ifges(7,3) * t90 / 0.2e1) * t90 + (t82 * mrSges(6,2) - t71 * mrSges(6,3) + Ifges(6,1) * t87 / 0.2e1) * t87 + (-t82 * mrSges(6,1) + t72 * mrSges(6,3) + Ifges(6,4) * t87 + Ifges(6,2) * t86 / 0.2e1) * t86 + (t75 * mrSges(7,2) - t67 * mrSges(7,3) + Ifges(7,5) * t90 + Ifges(7,1) * t77 / 0.2e1) * t77 + (-t75 * mrSges(7,1) + t68 * mrSges(7,3) + Ifges(7,4) * t77 + Ifges(7,6) * t90 + Ifges(7,2) * t76 / 0.2e1) * t76 + (-t88 * mrSges(4,1) - t84 * mrSges(5,1) + t83 * mrSges(5,2) + t89 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t85 * mrSges(4,1) + t89 * mrSges(4,2) - t81 * mrSges(5,2) + t84 * mrSges(5,3) + (-mrSges(3,2) * pkin(7) + Ifges(5,5) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t104 + (t81 * mrSges(5,1) + t71 * mrSges(6,1) + t88 * mrSges(4,2) - t72 * mrSges(6,2) - t85 * mrSges(4,3) - t83 * mrSges(5,3) + Ifges(6,5) * t87 + Ifges(6,6) * t86 + (-mrSges(3,1) * pkin(7) + Ifges(4,4) + Ifges(3,5) + Ifges(5,6)) * qJD(2)) * t102 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t102 ^ 2 + t104 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t111 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t104) * t104 + (-pkin(1) * mrSges(3,2) + (t111 + Ifges(6,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t102 + (Ifges(3,4) + Ifges(5,4) - Ifges(4,5)) * t104) * t102) * qJD(1)) * qJD(1);
T  = t1;
