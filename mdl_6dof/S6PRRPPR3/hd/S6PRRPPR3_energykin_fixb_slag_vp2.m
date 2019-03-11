% Calculate kinetic energy for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:04
% EndTime: 2019-03-08 21:08:04
% DurationCPUTime: 0.29s
% Computational Cost: add. (226->91), mult. (475->123), div. (0->0), fcn. (245->8), ass. (0->34)
t120 = -pkin(3) - pkin(4);
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t105 = cos(pkin(6));
t118 = qJD(1) * t105;
t109 = sin(qJ(2));
t104 = sin(pkin(6));
t119 = qJD(1) * t104;
t94 = qJD(2) * pkin(8) + t109 * t119;
t88 = t108 * t118 + t111 * t94;
t112 = cos(qJ(2));
t95 = -qJD(2) * pkin(2) - t112 * t119;
t117 = qJD(2) * t111;
t116 = t108 * qJD(2);
t86 = qJD(3) * qJ(4) + t88;
t89 = -pkin(3) * t117 - qJ(4) * t116 + t95;
t87 = -t108 * t94 + t111 * t118;
t115 = qJD(4) - t87;
t84 = pkin(4) * t117 + qJD(5) - t89;
t83 = qJ(5) * t117 - t86;
t114 = -qJ(5) * t116 + t115;
t110 = cos(qJ(6));
t107 = sin(qJ(6));
t98 = qJD(6) + t116;
t93 = -t107 * qJD(3) - t110 * t117;
t92 = -t110 * qJD(3) + t107 * t117;
t85 = -qJD(3) * pkin(3) + t115;
t82 = qJD(3) * pkin(5) - t83;
t81 = t120 * qJD(3) + t114;
t80 = (-pkin(9) + t120) * qJD(3) + t114;
t79 = (pkin(5) * t108 + pkin(9) * t111) * qJD(2) + t84;
t78 = t107 * t79 + t110 * t80;
t77 = -t107 * t80 + t110 * t79;
t1 = m(7) * (t77 ^ 2 + t78 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t81 ^ 2 + t83 ^ 2 + t84 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t95 ^ 2) / 0.2e1 + (t77 * mrSges(7,1) - t78 * mrSges(7,2) + Ifges(7,3) * t98 / 0.2e1) * t98 + (t82 * mrSges(7,2) - t77 * mrSges(7,3) + Ifges(7,5) * t98 + Ifges(7,1) * t93 / 0.2e1) * t93 + (-t82 * mrSges(7,1) + t78 * mrSges(7,3) + Ifges(7,4) * t93 + Ifges(7,6) * t98 + Ifges(7,2) * t92 / 0.2e1) * t92 + (m(3) * (t105 ^ 2 + (t109 ^ 2 + t112 ^ 2) * t104 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t87 * mrSges(4,1) - t85 * mrSges(5,1) - t83 * mrSges(6,1) - t88 * mrSges(4,2) + t81 * mrSges(6,2) + t86 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(3)) * qJD(3) + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t112 - mrSges(3,2) * t109) * t119 + (-t95 * mrSges(4,1) - t89 * mrSges(5,1) + t86 * mrSges(5,2) - t84 * mrSges(6,2) + t88 * mrSges(4,3) + t83 * mrSges(6,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t117 + (Ifges(6,5) + Ifges(4,6) - Ifges(5,6)) * qJD(3)) * t111 + (-t81 * mrSges(6,3) + t85 * mrSges(5,2) - t87 * mrSges(4,3) + t84 * mrSges(6,1) - t89 * mrSges(5,3) + t95 * mrSges(4,2) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t116 + (Ifges(4,4) + Ifges(6,4) - Ifges(5,5)) * t117 + (Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * qJD(3)) * t108) * qJD(2);
T  = t1;
