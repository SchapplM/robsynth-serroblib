% Calculate kinetic energy for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:22
% EndTime: 2019-03-08 22:54:22
% DurationCPUTime: 0.30s
% Computational Cost: add. (298->91), mult. (623->124), div. (0->0), fcn. (385->8), ass. (0->35)
t116 = cos(qJ(4));
t115 = pkin(4) + qJ(6);
t102 = sin(qJ(4));
t103 = sin(qJ(3));
t105 = cos(qJ(3));
t101 = cos(pkin(6));
t113 = qJD(1) * t101;
t104 = sin(qJ(2));
t100 = sin(pkin(6));
t114 = qJD(1) * t100;
t94 = qJD(2) * pkin(8) + t104 * t114;
t88 = t103 * t113 + t105 * t94;
t85 = qJD(3) * pkin(9) + t88;
t106 = cos(qJ(2));
t110 = t106 * t114;
t89 = -t110 + (-pkin(3) * t105 - pkin(9) * t103 - pkin(2)) * qJD(2);
t80 = t102 * t89 + t116 * t85;
t112 = qJD(2) * t103;
t111 = t105 * qJD(2);
t97 = qJD(4) - t111;
t78 = -t97 * qJ(5) - t80;
t79 = -t102 * t85 + t116 * t89;
t87 = -t103 * t94 + t105 * t113;
t109 = qJD(5) - t79;
t84 = -qJD(3) * pkin(3) - t87;
t93 = t102 * qJD(3) + t116 * t112;
t108 = -t93 * qJ(5) + t84;
t95 = -qJD(2) * pkin(2) - t110;
t92 = -t116 * qJD(3) + t102 * t112;
t81 = t92 * pkin(4) + t108;
t77 = -t97 * pkin(4) + t109;
t76 = t115 * t92 + t108;
t75 = -t92 * pkin(5) + qJD(6) - t78;
t74 = t93 * pkin(5) - t115 * t97 + t109;
t1 = m(4) * (t87 ^ 2 + t88 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t84 ^ 2) / 0.2e1 + m(7) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t81 ^ 2) / 0.2e1 + (t87 * mrSges(4,1) - t88 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (m(3) * (t101 ^ 2 + (t104 ^ 2 + t106 ^ 2) * t100 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t106 - mrSges(3,2) * t104) * t114 + (-t95 * mrSges(4,1) + t88 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t111 / 0.2e1) * t105 + (t95 * mrSges(4,2) - t87 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t105 + Ifges(4,1) * t103 / 0.2e1) * qJD(2)) * t103) * qJD(2) + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + t77 * mrSges(6,2) + t75 * mrSges(7,2) - t78 * mrSges(6,3) - t74 * mrSges(7,3) + (Ifges(5,3) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t97) * t97 + (t77 * mrSges(6,1) + t74 * mrSges(7,1) + t84 * mrSges(5,2) - t76 * mrSges(7,2) - t79 * mrSges(5,3) - t81 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t93 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t97) * t93 + (t84 * mrSges(5,1) + t78 * mrSges(6,1) - t75 * mrSges(7,1) - t81 * mrSges(6,2) - t80 * mrSges(5,3) + t76 * mrSges(7,3) + (Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t92 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t97 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t93) * t92;
T  = t1;
