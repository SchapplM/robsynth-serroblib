% Calculate kinetic energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:51
% EndTime: 2019-03-09 15:24:52
% DurationCPUTime: 0.56s
% Computational Cost: add. (662->107), mult. (1535->153), div. (0->0), fcn. (1110->8), ass. (0->40)
t118 = qJD(1) * (-pkin(8) - pkin(7));
t117 = pkin(4) + pkin(9);
t115 = pkin(7) * mrSges(3,3);
t103 = sin(pkin(10));
t114 = cos(pkin(10));
t102 = qJD(2) + qJD(3);
t109 = cos(qJ(2));
t100 = t109 * t118;
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t106 = sin(qJ(2));
t99 = qJD(2) * pkin(2) + t106 * t118;
t90 = t100 * t105 + t108 * t99;
t97 = (t105 * t109 + t106 * t108) * qJD(1);
t81 = pkin(3) * t102 - qJ(4) * t97 + t90;
t91 = -t108 * t100 + t105 * t99;
t96 = (-t105 * t106 + t108 * t109) * qJD(1);
t84 = qJ(4) * t96 + t91;
t78 = t103 * t81 + t114 * t84;
t76 = -qJ(5) * t102 - t78;
t101 = (-pkin(2) * t109 - pkin(1)) * qJD(1);
t77 = -t103 * t84 + t114 * t81;
t113 = qJD(5) - t77;
t89 = t103 * t96 + t114 * t97;
t92 = -pkin(3) * t96 + qJD(4) + t101;
t112 = -qJ(5) * t89 + t92;
t107 = cos(qJ(6));
t104 = sin(qJ(6));
t88 = t103 * t97 - t114 * t96;
t87 = qJD(6) + t89;
t86 = t102 * t107 + t104 * t88;
t85 = -t102 * t104 + t107 * t88;
t79 = pkin(4) * t88 + t112;
t75 = -t102 * pkin(4) + t113;
t74 = t117 * t88 + t112;
t73 = -pkin(5) * t88 - t76;
t72 = t89 * pkin(5) - t117 * t102 + t113;
t71 = t104 * t72 + t107 * t74;
t70 = -t104 * t74 + t107 * t72;
t1 = m(7) * (t70 ^ 2 + t71 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t92 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + (t101 * mrSges(4,2) - t90 * mrSges(4,3) + Ifges(4,1) * t97 / 0.2e1) * t97 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,3) * t87 / 0.2e1) * t87 + (-t101 * mrSges(4,1) + t91 * mrSges(4,3) + Ifges(4,4) * t97 + Ifges(4,2) * t96 / 0.2e1) * t96 + (t73 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,5) * t87 + Ifges(7,1) * t86 / 0.2e1) * t86 + (-t73 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,6) * t87 + Ifges(7,2) * t85 / 0.2e1) * t85 + (t75 * mrSges(6,1) + t92 * mrSges(5,2) - t77 * mrSges(5,3) - t79 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t89) * t89 + (t92 * mrSges(5,1) + t76 * mrSges(6,1) - t79 * mrSges(6,2) - t78 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t88 + (-Ifges(5,4) - Ifges(6,6)) * t89) * t88 + (t90 * mrSges(4,1) + t77 * mrSges(5,1) - t91 * mrSges(4,2) - t78 * mrSges(5,2) + t75 * mrSges(6,2) - t76 * mrSges(6,3) + Ifges(4,5) * t97 + Ifges(4,6) * t96 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t102 + (-Ifges(6,4) + Ifges(5,5)) * t89 + (Ifges(6,5) - Ifges(5,6)) * t88) * t102 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t106 ^ 2 + t109 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t115) * t109) * t109 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t109 + (Ifges(3,1) / 0.2e1 + t115) * t106) * t106) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t109 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t106) * qJD(2)) * qJD(1);
T  = t1;
