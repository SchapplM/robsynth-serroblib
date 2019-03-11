% Calculate kinetic energy for
% S6RRRPRR5
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:21
% EndTime: 2019-03-09 18:20:21
% DurationCPUTime: 0.59s
% Computational Cost: add. (646->107), mult. (1305->155), div. (0->0), fcn. (892->8), ass. (0->43)
t123 = pkin(3) + pkin(9);
t122 = -pkin(8) - pkin(7);
t121 = pkin(7) * mrSges(3,3);
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t106 = qJD(2) + qJD(3);
t110 = sin(qJ(2));
t119 = t110 * qJD(1);
t101 = qJD(2) * pkin(2) + t119 * t122;
t114 = cos(qJ(2));
t120 = qJD(1) * t114;
t102 = t122 * t120;
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t90 = t101 * t113 + t109 * t102;
t118 = qJD(4) - t90;
t98 = (t109 * t114 + t110 * t113) * qJD(1);
t82 = pkin(4) * t98 - t106 * t123 + t118;
t103 = (-pkin(2) * t114 - pkin(1)) * qJD(1);
t117 = -qJ(4) * t98 + t103;
t97 = t109 * t119 - t113 * t120;
t83 = t123 * t97 + t117;
t77 = t108 * t82 + t112 * t83;
t91 = t109 * t101 - t113 * t102;
t89 = -t106 * qJ(4) - t91;
t76 = -t108 * t83 + t112 * t82;
t86 = -pkin(4) * t97 - t89;
t96 = qJD(5) + t98;
t111 = cos(qJ(6));
t107 = sin(qJ(6));
t94 = qJD(6) + t96;
t93 = t106 * t112 + t108 * t97;
t92 = -t106 * t108 + t112 * t97;
t88 = -pkin(3) * t106 + t118;
t87 = pkin(3) * t97 + t117;
t85 = t107 * t92 + t111 * t93;
t84 = -t107 * t93 + t111 * t92;
t78 = -pkin(5) * t92 + t86;
t75 = pkin(10) * t92 + t77;
t74 = pkin(5) * t96 - pkin(10) * t93 + t76;
t73 = t107 * t74 + t111 * t75;
t72 = -t107 * t75 + t111 * t74;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t103 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(6) * (t76 ^ 2 + t77 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t72 ^ 2 + t73 ^ 2 + t78 ^ 2) / 0.2e1 + (t76 * mrSges(6,1) - t77 * mrSges(6,2) + Ifges(6,3) * t96 / 0.2e1) * t96 + (t72 * mrSges(7,1) - t73 * mrSges(7,2) + Ifges(7,3) * t94 / 0.2e1) * t94 + (t86 * mrSges(6,2) - t76 * mrSges(6,3) + Ifges(6,5) * t96 + Ifges(6,1) * t93 / 0.2e1) * t93 + (t78 * mrSges(7,2) - t72 * mrSges(7,3) + Ifges(7,5) * t94 + Ifges(7,1) * t85 / 0.2e1) * t85 + (-t86 * mrSges(6,1) + t77 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,6) * t96 + Ifges(6,2) * t92 / 0.2e1) * t92 + (-t78 * mrSges(7,1) + t73 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,6) * t94 + Ifges(7,2) * t84 / 0.2e1) * t84 + (t88 * mrSges(5,1) + t103 * mrSges(4,2) - t90 * mrSges(4,3) - t87 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t98) * t98 + (t103 * mrSges(4,1) + t89 * mrSges(5,1) - t87 * mrSges(5,2) - t91 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t97 + (-Ifges(4,4) - Ifges(5,6)) * t98) * t97 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + t88 * mrSges(5,2) - t89 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t106 + (-Ifges(5,4) + Ifges(4,5)) * t98 + (Ifges(5,5) - Ifges(4,6)) * t97) * t106 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t110 ^ 2 + t114 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t121 + Ifges(3,2) / 0.2e1) * t114) * t114 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t114 + (t121 + Ifges(3,1) / 0.2e1) * t110) * t110) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t114 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t110) * qJD(2)) * qJD(1);
T  = t1;
