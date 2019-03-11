% Calculate kinetic energy for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:11
% EndTime: 2019-03-09 08:17:11
% DurationCPUTime: 0.45s
% Computational Cost: add. (374->105), mult. (777->139), div. (0->0), fcn. (424->6), ass. (0->37)
t111 = -pkin(4) - pkin(5);
t110 = pkin(7) * mrSges(3,3);
t109 = -pkin(2) - qJ(4);
t100 = cos(qJ(2));
t98 = sin(qJ(2));
t105 = -qJ(3) * t98 - pkin(1);
t81 = (t109 * t100 + t105) * qJD(1);
t108 = qJD(1) * t98;
t107 = pkin(7) * t108 + qJD(3);
t82 = pkin(3) * t108 + t109 * qJD(2) + t107;
t95 = sin(pkin(9));
t96 = cos(pkin(9));
t74 = t96 * t81 + t95 * t82;
t106 = qJD(1) * t100;
t88 = -pkin(7) * t106 - qJD(2) * qJ(3);
t72 = qJ(5) * t108 + t74;
t73 = -t95 * t81 + t82 * t96;
t83 = pkin(3) * t106 + qJD(4) - t88;
t104 = qJD(5) - t73;
t86 = qJD(2) * t96 - t95 * t106;
t103 = qJ(5) * t86 - t83;
t99 = cos(qJ(6));
t97 = sin(qJ(6));
t89 = qJD(6) - t108;
t87 = -qJD(2) * pkin(2) + t107;
t85 = qJD(2) * t95 + t96 * t106;
t84 = (-pkin(2) * t100 + t105) * qJD(1);
t77 = t85 * t97 + t86 * t99;
t76 = t85 * t99 - t86 * t97;
t75 = pkin(4) * t85 - t103;
t71 = -pkin(4) * t108 + t104;
t70 = t111 * t85 + t103;
t69 = pkin(8) * t85 + t72;
t68 = -pkin(8) * t86 + t111 * t108 + t104;
t67 = t68 * t97 + t69 * t99;
t66 = t68 * t99 - t69 * t97;
t1 = m(4) * (t84 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t75 ^ 2) / 0.2e1 + m(7) * (t66 ^ 2 + t67 ^ 2 + t70 ^ 2) / 0.2e1 + (t66 * mrSges(7,1) - t67 * mrSges(7,2) + Ifges(7,3) * t89 / 0.2e1) * t89 + (t70 * mrSges(7,2) - t66 * mrSges(7,3) + Ifges(7,5) * t89 + Ifges(7,1) * t77 / 0.2e1) * t77 + (-t70 * mrSges(7,1) + t67 * mrSges(7,3) + Ifges(7,4) * t77 + Ifges(7,6) * t89 + Ifges(7,2) * t76 / 0.2e1) * t76 + (t87 * mrSges(4,2) - t88 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(2)) * qJD(2) + (t83 * mrSges(5,2) + t71 * mrSges(6,2) - t73 * mrSges(5,3) - t75 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t86) * t86 + (t83 * mrSges(5,1) + t75 * mrSges(6,1) - t72 * mrSges(6,2) - t74 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t85 + (-Ifges(5,4) + Ifges(6,5)) * t86 + (-Ifges(5,6) + Ifges(6,6)) * t108) * t85 + ((-t88 * mrSges(4,1) + t84 * mrSges(4,2) + (-pkin(7) * mrSges(3,2) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t100 + (t87 * mrSges(4,1) + t73 * mrSges(5,1) - t71 * mrSges(6,1) - t74 * mrSges(5,2) - t84 * mrSges(4,3) + t72 * mrSges(6,3) + (Ifges(6,4) + Ifges(5,5)) * t86 + (-pkin(7) * mrSges(3,1) - Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t98 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t100 ^ 2 + t98 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t110 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t100) * t100 + (-pkin(1) * mrSges(3,2) + (t110 + Ifges(3,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t98 + (Ifges(3,4) + Ifges(4,6)) * t100) * t98) * qJD(1)) * qJD(1);
T  = t1;
