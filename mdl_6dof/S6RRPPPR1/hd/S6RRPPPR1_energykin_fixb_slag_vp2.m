% Calculate kinetic energy for
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:05:59
% EndTime: 2019-03-09 08:06:00
% DurationCPUTime: 0.60s
% Computational Cost: add. (578->107), mult. (1375->152), div. (0->0), fcn. (962->8), ass. (0->41)
t124 = -pkin(4) - pkin(5);
t123 = pkin(7) * mrSges(3,3);
t122 = pkin(7) + qJ(3);
t107 = sin(pkin(10));
t120 = cos(pkin(10));
t112 = cos(qJ(2));
t104 = qJD(3) + (-pkin(2) * t112 - pkin(1)) * qJD(1);
t108 = sin(pkin(9));
t118 = t112 * qJD(1);
t110 = sin(qJ(2));
t119 = t110 * qJD(1);
t121 = cos(pkin(9));
t98 = t108 * t119 - t121 * t118;
t99 = (t108 * t112 + t110 * t121) * qJD(1);
t85 = pkin(3) * t98 - qJ(4) * t99 + t104;
t102 = qJD(2) * pkin(2) - t119 * t122;
t103 = t122 * t118;
t91 = t108 * t102 + t121 * t103;
t89 = qJD(2) * qJ(4) + t91;
t81 = t107 * t85 + t120 * t89;
t90 = t121 * t102 - t108 * t103;
t78 = t98 * qJ(5) + t81;
t80 = -t107 * t89 + t120 * t85;
t117 = qJD(2) * pkin(3) - qJD(4) + t90;
t116 = qJD(5) - t80;
t93 = t107 * qJD(2) + t120 * t99;
t115 = qJ(5) * t93 + t117;
t111 = cos(qJ(6));
t109 = sin(qJ(6));
t97 = qJD(6) - t98;
t92 = -qJD(2) * t120 + t107 * t99;
t83 = t109 * t92 + t111 * t93;
t82 = -t109 * t93 + t111 * t92;
t79 = pkin(4) * t92 - t115;
t77 = -t98 * pkin(4) + t116;
t76 = t124 * t92 + t115;
t75 = pkin(8) * t92 + t78;
t74 = -t93 * pkin(8) + t124 * t98 + t116;
t73 = t109 * t74 + t111 * t75;
t72 = -t109 * t75 + t111 * t74;
t1 = m(5) * (t117 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(7) * (t72 ^ 2 + t73 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(4) * (t104 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + (t104 * mrSges(4,2) - t90 * mrSges(4,3) + Ifges(4,1) * t99 / 0.2e1) * t99 + (t72 * mrSges(7,1) - t73 * mrSges(7,2) + Ifges(7,3) * t97 / 0.2e1) * t97 + (t76 * mrSges(7,2) - t72 * mrSges(7,3) + Ifges(7,5) * t97 + Ifges(7,1) * t83 / 0.2e1) * t83 + (-t76 * mrSges(7,1) + t73 * mrSges(7,3) + Ifges(7,4) * t83 + Ifges(7,6) * t97 + Ifges(7,2) * t82 / 0.2e1) * t82 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + Ifges(4,5) * t99 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t110 + Ifges(3,6) * t112 + (-mrSges(3,1) * t110 - mrSges(3,2) * t112) * pkin(7)) * qJD(1)) * qJD(2) + (-t117 * mrSges(5,2) + t77 * mrSges(6,2) - t80 * mrSges(5,3) - t79 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t93) * t93 + (-t117 * mrSges(5,1) + t79 * mrSges(6,1) - t78 * mrSges(6,2) - t81 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t92 + (-Ifges(5,4) + Ifges(6,5)) * t93) * t92 + (t104 * mrSges(4,1) + t80 * mrSges(5,1) - t77 * mrSges(6,1) - t81 * mrSges(5,2) - t91 * mrSges(4,3) + t78 * mrSges(6,3) - Ifges(4,4) * t99 - Ifges(4,6) * qJD(2) + (Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t98 + (Ifges(6,4) + Ifges(5,5)) * t93 + (-Ifges(5,6) + Ifges(6,6)) * t92) * t98 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t110 ^ 2 + t112 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t123) * t112) * t112 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t112 + (Ifges(3,1) / 0.2e1 + t123) * t110) * t110) * qJD(1) ^ 2;
T  = t1;
