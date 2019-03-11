% Calculate kinetic energy for
% S6PRRRPP2
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:39
% EndTime: 2019-03-08 22:48:40
% DurationCPUTime: 0.37s
% Computational Cost: add. (298->91), mult. (623->124), div. (0->0), fcn. (385->8), ass. (0->35)
t123 = pkin(4) + pkin(5);
t122 = cos(qJ(4));
t108 = sin(qJ(4));
t109 = sin(qJ(3));
t111 = cos(qJ(3));
t107 = cos(pkin(6));
t120 = qJD(1) * t107;
t110 = sin(qJ(2));
t106 = sin(pkin(6));
t121 = qJD(1) * t106;
t98 = qJD(2) * pkin(8) + t110 * t121;
t92 = t109 * t120 + t111 * t98;
t89 = qJD(3) * pkin(9) + t92;
t112 = cos(qJ(2));
t117 = t112 * t121;
t93 = -t117 + (-pkin(3) * t111 - pkin(9) * t109 - pkin(2)) * qJD(2);
t84 = t108 * t93 + t122 * t89;
t91 = -t109 * t98 + t111 * t120;
t119 = qJD(2) * t109;
t118 = qJD(2) * t111;
t103 = -qJD(4) + t118;
t82 = -t103 * qJ(5) + t84;
t116 = qJD(3) * pkin(3) + t91;
t83 = -t108 * t89 + t122 * t93;
t115 = qJD(5) - t83;
t97 = t108 * qJD(3) + t122 * t119;
t114 = qJ(5) * t97 + t116;
t99 = -qJD(2) * pkin(2) - t117;
t96 = -t122 * qJD(3) + t108 * t119;
t85 = pkin(4) * t96 - t114;
t81 = t103 * pkin(4) + t115;
t80 = -t123 * t96 + qJD(6) + t114;
t79 = qJ(6) * t96 + t82;
t78 = -t97 * qJ(6) + t123 * t103 + t115;
t1 = m(7) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t85 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t99 ^ 2) / 0.2e1 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (m(3) * (t107 ^ 2 + (t110 ^ 2 + t112 ^ 2) * t106 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t112 - mrSges(3,2) * t110) * t121 + (-t99 * mrSges(4,1) + t92 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t118 / 0.2e1) * t111 + (t99 * mrSges(4,2) - t91 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t111 + Ifges(4,1) * t109 / 0.2e1) * qJD(2)) * t109) * qJD(2) + (-t116 * mrSges(5,2) + t81 * mrSges(6,2) + t80 * mrSges(7,2) - t83 * mrSges(5,3) - t85 * mrSges(6,3) - t78 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t97) * t97 + (-t116 * mrSges(5,1) + t85 * mrSges(6,1) - t80 * mrSges(7,1) - t82 * mrSges(6,2) - t84 * mrSges(5,3) + t79 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t96 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t97) * t96 + (-t83 * mrSges(5,1) + t81 * mrSges(6,1) + t78 * mrSges(7,1) + t84 * mrSges(5,2) - t79 * mrSges(7,2) - t82 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t103 + (-Ifges(6,4) - Ifges(5,5) + Ifges(7,5)) * t97 + (Ifges(5,6) - Ifges(6,6) + Ifges(7,6)) * t96) * t103;
T  = t1;
