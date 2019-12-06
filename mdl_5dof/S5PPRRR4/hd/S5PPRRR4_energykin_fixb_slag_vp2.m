% Calculate kinetic energy for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:18:32
% EndTime: 2019-12-05 15:18:32
% DurationCPUTime: 0.23s
% Computational Cost: add. (190->58), mult. (462->101), div. (0->0), fcn. (350->12), ass. (0->35)
t104 = cos(pkin(5)) * qJD(1) + qJD(2);
t108 = sin(pkin(6));
t111 = cos(pkin(6));
t110 = cos(pkin(11));
t109 = sin(pkin(5));
t126 = qJD(1) * t109;
t122 = t110 * t126;
t128 = t104 * t108 + t111 * t122;
t114 = sin(qJ(3));
t117 = cos(qJ(3));
t107 = sin(pkin(11));
t123 = t107 * t126;
t94 = -t114 * t123 + t128 * t117;
t113 = sin(qJ(4));
t116 = cos(qJ(4));
t95 = t128 * t114 + t117 * t123;
t93 = qJD(3) * pkin(8) + t95;
t97 = t111 * t104 - t108 * t122;
t89 = t113 * t97 + t116 * t93;
t125 = qJD(3) * t113;
t124 = t116 * qJD(3);
t88 = -t113 * t93 + t116 * t97;
t118 = qJD(1) ^ 2;
t115 = cos(qJ(5));
t112 = sin(qJ(5));
t105 = qJD(5) - t124;
t101 = t112 * qJD(4) + t115 * t125;
t100 = t115 * qJD(4) - t112 * t125;
t92 = -qJD(3) * pkin(3) - t94;
t90 = (-pkin(4) * t116 - pkin(9) * t113 - pkin(3)) * qJD(3) - t94;
t87 = qJD(4) * pkin(9) + t89;
t86 = -qJD(4) * pkin(4) - t88;
t85 = t112 * t90 + t115 * t87;
t84 = -t112 * t87 + t115 * t90;
t1 = m(5) * (t88 ^ 2 + t89 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(4) * (t94 ^ 2 + t95 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t104 ^ 2 + (t107 ^ 2 + t110 ^ 2) * t118 * t109 ^ 2) / 0.2e1 + m(2) * t118 / 0.2e1 + (t84 * mrSges(6,1) - t85 * mrSges(6,2) + Ifges(6,3) * t105 / 0.2e1) * t105 + (t88 * mrSges(5,1) - t89 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t86 * mrSges(6,2) - t84 * mrSges(6,3) + Ifges(6,5) * t105 + Ifges(6,1) * t101 / 0.2e1) * t101 + (-t86 * mrSges(6,1) + t85 * mrSges(6,3) + Ifges(6,4) * t101 + Ifges(6,6) * t105 + Ifges(6,2) * t100 / 0.2e1) * t100 + (-t95 * mrSges(4,2) + t94 * mrSges(4,1) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t92 * mrSges(5,1) + t89 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t124 / 0.2e1) * t116 + (t92 * mrSges(5,2) - t88 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t116 + Ifges(5,1) * t113 / 0.2e1) * qJD(3)) * t113) * qJD(3);
T = t1;
