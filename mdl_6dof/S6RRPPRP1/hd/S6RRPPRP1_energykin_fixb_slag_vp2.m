% Calculate kinetic energy for
% S6RRPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-03-09 08:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:25:20
% EndTime: 2019-03-09 08:25:21
% DurationCPUTime: 0.53s
% Computational Cost: add. (694->106), mult. (1641->152), div. (0->0), fcn. (1186->8), ass. (0->38)
t116 = pkin(7) * mrSges(3,3);
t115 = cos(qJ(5));
t114 = pkin(7) + qJ(3);
t107 = sin(qJ(5));
t103 = sin(pkin(10));
t105 = cos(pkin(10));
t109 = cos(qJ(2));
t101 = qJD(3) + (-pkin(2) * t109 - pkin(1)) * qJD(1);
t104 = sin(pkin(9));
t106 = cos(pkin(9));
t108 = sin(qJ(2));
t112 = t108 * qJD(1);
t113 = qJD(1) * t109;
t95 = t104 * t112 - t106 * t113;
t96 = (t104 * t109 + t106 * t108) * qJD(1);
t84 = pkin(3) * t95 - qJ(4) * t96 + t101;
t100 = t114 * t113;
t99 = qJD(2) * pkin(2) - t112 * t114;
t89 = t106 * t100 + t104 * t99;
t87 = qJD(2) * qJ(4) + t89;
t77 = -t103 * t87 + t105 * t84;
t92 = qJD(2) * t103 + t105 * t96;
t74 = pkin(4) * t95 - pkin(8) * t92 + t77;
t78 = t103 * t84 + t105 * t87;
t91 = qJD(2) * t105 - t103 * t96;
t76 = pkin(8) * t91 + t78;
t71 = t107 * t74 + t115 * t76;
t88 = -t104 * t100 + t106 * t99;
t70 = -t107 * t76 + t115 * t74;
t86 = -qJD(2) * pkin(3) + qJD(4) - t88;
t79 = -pkin(4) * t91 + t86;
t94 = qJD(5) + t95;
t81 = t107 * t91 + t115 * t92;
t80 = t107 * t92 - t115 * t91;
t72 = pkin(5) * t80 - qJ(6) * t81 + t79;
t69 = qJ(6) * t94 + t71;
t68 = -t94 * pkin(5) + qJD(6) - t70;
t1 = m(4) * (t101 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t79 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t72 ^ 2) / 0.2e1 + (t101 * mrSges(4,2) - t88 * mrSges(4,3) + Ifges(4,1) * t96 / 0.2e1) * t96 + (t86 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,1) * t92 / 0.2e1) * t92 + (-t86 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t92 + Ifges(5,2) * t91 / 0.2e1) * t91 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + Ifges(4,5) * t96 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t108 + Ifges(3,6) * t109 + (-mrSges(3,1) * t108 - mrSges(3,2) * t109) * pkin(7)) * qJD(1)) * qJD(2) + (t70 * mrSges(6,1) - t68 * mrSges(7,1) - t71 * mrSges(6,2) + t69 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t94) * t94 + (t79 * mrSges(6,2) + t68 * mrSges(7,2) - t70 * mrSges(6,3) - t72 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t81 + (Ifges(7,4) + Ifges(6,5)) * t94) * t81 + (t101 * mrSges(4,1) + t77 * mrSges(5,1) - t78 * mrSges(5,2) - t89 * mrSges(4,3) - Ifges(4,4) * t96 + Ifges(5,5) * t92 - Ifges(4,6) * qJD(2) + Ifges(5,6) * t91 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t95) * t95 + (t79 * mrSges(6,1) + t72 * mrSges(7,1) - t69 * mrSges(7,2) - t71 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t80 + (-Ifges(6,6) + Ifges(7,6)) * t94 + (-Ifges(6,4) + Ifges(7,5)) * t81) * t80 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t108 ^ 2 + t109 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t116 + Ifges(3,2) / 0.2e1) * t109) * t109 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t109 + (t116 + Ifges(3,1) / 0.2e1) * t108) * t108) * qJD(1) ^ 2;
T  = t1;
