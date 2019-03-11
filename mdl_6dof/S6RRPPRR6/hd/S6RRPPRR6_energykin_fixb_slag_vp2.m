% Calculate kinetic energy for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:12:57
% EndTime: 2019-03-09 09:12:58
% DurationCPUTime: 0.54s
% Computational Cost: add. (558->108), mult. (1211->158), div. (0->0), fcn. (788->8), ass. (0->39)
t123 = pkin(7) * mrSges(3,3);
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t109 = sin(pkin(10));
t110 = cos(pkin(10));
t114 = sin(qJ(2));
t122 = qJD(1) * t114;
t120 = pkin(7) * t122 + qJD(3);
t94 = -qJ(4) * t122 + (-pkin(2) - pkin(3)) * qJD(2) + t120;
t117 = cos(qJ(2));
t121 = qJD(1) * t117;
t101 = pkin(7) * t121 + qJD(2) * qJ(3);
t98 = -qJ(4) * t121 + t101;
t85 = -t109 * t98 + t110 * t94;
t97 = (-t109 * t117 + t110 * t114) * qJD(1);
t81 = -qJD(2) * pkin(4) - pkin(8) * t97 + t85;
t86 = t109 * t94 + t110 * t98;
t96 = (-t109 * t114 - t110 * t117) * qJD(1);
t82 = pkin(8) * t96 + t86;
t77 = t113 * t81 + t116 * t82;
t99 = -qJD(1) * pkin(1) - pkin(2) * t121 - qJ(3) * t122;
t76 = -t113 * t82 + t116 * t81;
t89 = -t113 * t97 + t116 * t96;
t93 = pkin(3) * t121 + qJD(4) - t99;
t88 = -pkin(4) * t96 + t93;
t115 = cos(qJ(6));
t112 = sin(qJ(6));
t107 = -qJD(2) + qJD(5);
t100 = -qJD(2) * pkin(2) + t120;
t90 = t113 * t96 + t116 * t97;
t87 = qJD(6) - t89;
t84 = t107 * t112 + t115 * t90;
t83 = t107 * t115 - t112 * t90;
t78 = -pkin(5) * t89 - pkin(9) * t90 + t88;
t75 = pkin(9) * t107 + t77;
t74 = -pkin(5) * t107 - t76;
t73 = t112 * t78 + t115 * t75;
t72 = -t112 * t75 + t115 * t78;
t1 = m(4) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t76 ^ 2 + t77 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + (t93 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t97 / 0.2e1) * t97 + (t88 * mrSges(6,2) - t76 * mrSges(6,3) + Ifges(6,1) * t90 / 0.2e1) * t90 + (t72 * mrSges(7,1) - t73 * mrSges(7,2) + Ifges(7,3) * t87 / 0.2e1) * t87 + (-t93 * mrSges(5,1) + t86 * mrSges(5,3) + Ifges(5,4) * t97 + Ifges(5,2) * t96 / 0.2e1) * t96 + (-t88 * mrSges(6,1) + t77 * mrSges(6,3) + Ifges(6,4) * t90 + Ifges(6,2) * t89 / 0.2e1) * t89 + (t74 * mrSges(7,2) - t72 * mrSges(7,3) + Ifges(7,5) * t87 + Ifges(7,1) * t84 / 0.2e1) * t84 + (-t74 * mrSges(7,1) + t73 * mrSges(7,3) + Ifges(7,4) * t84 + Ifges(7,6) * t87 + Ifges(7,2) * t83 / 0.2e1) * t83 + (t76 * mrSges(6,1) - t77 * mrSges(6,2) + Ifges(6,5) * t90 + Ifges(6,6) * t89 + Ifges(6,3) * t107 / 0.2e1) * t107 + (-t100 * mrSges(4,1) - t85 * mrSges(5,1) + t86 * mrSges(5,2) + t101 * mrSges(4,3) - Ifges(5,5) * t97 - Ifges(5,6) * t96 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t99 * mrSges(4,1) + t101 * mrSges(4,2)) * t117 + (t100 * mrSges(4,2) - t99 * mrSges(4,3)) * t114 + ((-mrSges(3,2) * pkin(7) + Ifges(3,6) - Ifges(4,6)) * t117 + (-mrSges(3,1) * pkin(7) + Ifges(4,4) + Ifges(3,5)) * t114) * qJD(2) + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t114 ^ 2 + t117 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + t123) * t117) * t117 + (-pkin(1) * mrSges(3,2) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + t123) * t114 + (Ifges(3,4) - Ifges(4,5)) * t117) * t114) * qJD(1)) * qJD(1);
T  = t1;
