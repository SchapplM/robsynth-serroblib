% Calculate kinetic energy for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR11_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:11:36
% EndTime: 2019-03-09 11:11:36
% DurationCPUTime: 0.55s
% Computational Cost: add. (640->108), mult. (1291->156), div. (0->0), fcn. (822->8), ass. (0->41)
t121 = -pkin(2) - pkin(8);
t120 = pkin(7) * mrSges(3,3);
t107 = sin(pkin(10));
t108 = cos(pkin(10));
t111 = sin(qJ(2));
t105 = t111 * qJD(1);
t101 = t105 + qJD(4);
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t114 = cos(qJ(2));
t117 = -qJ(3) * t111 - pkin(1);
t92 = (t114 * t121 + t117) * qJD(1);
t118 = pkin(7) * t105 + qJD(3);
t93 = pkin(3) * t105 + qJD(2) * t121 + t118;
t84 = -t110 * t92 + t113 * t93;
t119 = qJD(1) * t114;
t97 = qJD(2) * t113 - t110 * t119;
t80 = pkin(4) * t101 - qJ(5) * t97 + t84;
t85 = t110 * t93 + t113 * t92;
t96 = -qJD(2) * t110 - t113 * t119;
t83 = qJ(5) * t96 + t85;
t75 = t107 * t80 + t108 * t83;
t99 = -pkin(7) * t119 - qJD(2) * qJ(3);
t94 = pkin(3) * t119 - t99;
t74 = -t107 * t83 + t108 * t80;
t88 = -pkin(4) * t96 + qJD(5) + t94;
t112 = cos(qJ(6));
t109 = sin(qJ(6));
t100 = qJD(6) + t101;
t98 = -qJD(2) * pkin(2) + t118;
t95 = (-pkin(2) * t114 + t117) * qJD(1);
t87 = t107 * t96 + t108 * t97;
t86 = -t107 * t97 + t108 * t96;
t81 = -pkin(5) * t86 + t88;
t77 = t109 * t86 + t112 * t87;
t76 = -t109 * t87 + t112 * t86;
t73 = pkin(9) * t86 + t75;
t72 = pkin(5) * t101 - pkin(9) * t87 + t74;
t71 = t109 * t72 + t112 * t73;
t70 = -t109 * t73 + t112 * t72;
t1 = m(5) * (t84 ^ 2 + t85 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t95 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t81 ^ 2) / 0.2e1 + (t94 * mrSges(5,2) - t84 * mrSges(5,3) + Ifges(5,1) * t97 / 0.2e1) * t97 + (t88 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,1) * t87 / 0.2e1) * t87 + (t81 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,1) * t77 / 0.2e1) * t77 + (-t94 * mrSges(5,1) + t85 * mrSges(5,3) + Ifges(5,4) * t97 + Ifges(5,2) * t96 / 0.2e1) * t96 + (-t88 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t87 + Ifges(6,2) * t86 / 0.2e1) * t86 + (-t81 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t77 + Ifges(7,2) * t76 / 0.2e1) * t76 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,5) * t77 + Ifges(7,6) * t76 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t98 * mrSges(4,2) - t99 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + (t84 * mrSges(5,1) + t74 * mrSges(6,1) - t85 * mrSges(5,2) - t75 * mrSges(6,2) + Ifges(5,5) * t97 + Ifges(6,5) * t87 + Ifges(5,6) * t96 + Ifges(6,6) * t86 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t101) * t101 + ((-t99 * mrSges(4,1) + t95 * mrSges(4,2) + (-pkin(7) * mrSges(3,2) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t114 + (t98 * mrSges(4,1) - t95 * mrSges(4,3) + (-pkin(7) * mrSges(3,1) - Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t111 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t111 ^ 2 + t114 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t120) * t114) * t114 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + t120) * t111 + (Ifges(3,4) + Ifges(4,6)) * t114) * t111) * qJD(1)) * qJD(1);
T  = t1;
