% Calculate kinetic energy for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:46
% EndTime: 2019-03-08 22:42:46
% DurationCPUTime: 0.36s
% Computational Cost: add. (442->94), mult. (951->139), div. (0->0), fcn. (661->10), ass. (0->38)
t111 = sin(pkin(11));
t126 = cos(pkin(11));
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t115 = sin(qJ(3));
t123 = t115 * qJD(2);
t105 = qJD(3) * t114 + t117 * t123;
t118 = cos(qJ(3));
t122 = t118 * qJD(2);
t109 = qJD(4) - t122;
t119 = cos(qJ(2));
t112 = sin(pkin(6));
t125 = qJD(1) * t112;
t121 = t119 * t125;
t100 = -t121 + (-pkin(3) * t118 - pkin(9) * t115 - pkin(2)) * qJD(2);
t116 = sin(qJ(2));
t106 = qJD(2) * pkin(8) + t116 * t125;
t113 = cos(pkin(6));
t124 = qJD(1) * t113;
t99 = t118 * t106 + t115 * t124;
t95 = qJD(3) * pkin(9) + t99;
t88 = t117 * t100 - t114 * t95;
t85 = pkin(4) * t109 - qJ(5) * t105 + t88;
t104 = qJD(3) * t117 - t114 * t123;
t89 = t114 * t100 + t117 * t95;
t87 = qJ(5) * t104 + t89;
t82 = t111 * t85 + t126 * t87;
t98 = -t115 * t106 + t118 * t124;
t81 = -t111 * t87 + t126 * t85;
t94 = -qJD(3) * pkin(3) - t98;
t90 = -pkin(4) * t104 + qJD(5) + t94;
t107 = -qJD(2) * pkin(2) - t121;
t92 = t111 * t104 + t105 * t126;
t91 = -t104 * t126 + t105 * t111;
t83 = pkin(5) * t91 - qJ(6) * t92 + t90;
t80 = qJ(6) * t109 + t82;
t79 = -t109 * pkin(5) + qJD(6) - t81;
t1 = m(5) * (t88 ^ 2 + t89 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t82 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + (t94 * mrSges(5,2) - t88 * mrSges(5,3) + Ifges(5,1) * t105 / 0.2e1) * t105 + (t98 * mrSges(4,1) - t99 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t94 * mrSges(5,1) + t89 * mrSges(5,3) + Ifges(5,4) * t105 + Ifges(5,2) * t104 / 0.2e1) * t104 + (m(2) / 0.2e1 + m(3) * (t113 ^ 2 + (t116 ^ 2 + t119 ^ 2) * t112 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t90 * mrSges(6,2) + t79 * mrSges(7,2) - t81 * mrSges(6,3) - t83 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t92) * t92 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t119 - mrSges(3,2) * t116) * t125 + (-t107 * mrSges(4,1) + t99 * mrSges(4,3) + Ifges(4,6) * qJD(3) + Ifges(4,2) * t122 / 0.2e1) * t118 + (t107 * mrSges(4,2) - t98 * mrSges(4,3) + Ifges(4,5) * qJD(3) + (Ifges(4,4) * t118 + Ifges(4,1) * t115 / 0.2e1) * qJD(2)) * t115) * qJD(2) + (t90 * mrSges(6,1) + t83 * mrSges(7,1) - t80 * mrSges(7,2) - t82 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t91 + (-Ifges(6,4) + Ifges(7,5)) * t92) * t91 + (t88 * mrSges(5,1) + t81 * mrSges(6,1) - t79 * mrSges(7,1) - t89 * mrSges(5,2) - t82 * mrSges(6,2) + t80 * mrSges(7,3) + Ifges(5,5) * t105 + Ifges(5,6) * t104 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t109 + (Ifges(7,4) + Ifges(6,5)) * t92 + (-Ifges(6,6) + Ifges(7,6)) * t91) * t109;
T  = t1;
