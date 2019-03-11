% Calculate kinetic energy for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:41
% EndTime: 2019-03-09 08:49:42
% DurationCPUTime: 0.64s
% Computational Cost: add. (942->110), mult. (2231->169), div. (0->0), fcn. (1678->10), ass. (0->45)
t130 = pkin(7) * mrSges(3,3);
t129 = pkin(7) + qJ(3);
t120 = sin(qJ(5));
t123 = cos(qJ(5));
t116 = sin(pkin(10));
t118 = cos(pkin(10));
t121 = sin(qJ(2));
t124 = cos(qJ(2));
t108 = (t116 * t124 + t118 * t121) * qJD(1);
t115 = sin(pkin(11));
t117 = cos(pkin(11));
t103 = qJD(2) * t115 + t108 * t117;
t127 = t121 * qJD(1);
t128 = qJD(1) * t124;
t107 = t116 * t127 - t118 * t128;
t113 = qJD(3) + (-pkin(2) * t124 - pkin(1)) * qJD(1);
t96 = pkin(3) * t107 - qJ(4) * t108 + t113;
t111 = qJD(2) * pkin(2) - t129 * t127;
t112 = t129 * t128;
t101 = t116 * t111 + t118 * t112;
t99 = qJD(2) * qJ(4) + t101;
t89 = -t115 * t99 + t117 * t96;
t83 = pkin(4) * t107 - pkin(8) * t103 + t89;
t102 = qJD(2) * t117 - t108 * t115;
t90 = t115 * t96 + t117 * t99;
t88 = pkin(8) * t102 + t90;
t80 = t120 * t83 + t123 * t88;
t79 = -t120 * t88 + t123 * t83;
t100 = t111 * t118 - t116 * t112;
t106 = qJD(5) + t107;
t98 = -qJD(2) * pkin(3) + qJD(4) - t100;
t91 = -pkin(4) * t102 + t98;
t122 = cos(qJ(6));
t119 = sin(qJ(6));
t104 = qJD(6) + t106;
t93 = t102 * t120 + t103 * t123;
t92 = t102 * t123 - t103 * t120;
t86 = t119 * t92 + t122 * t93;
t85 = -t119 * t93 + t122 * t92;
t84 = -pkin(5) * t92 + t91;
t78 = pkin(9) * t92 + t80;
t77 = pkin(5) * t106 - pkin(9) * t93 + t79;
t76 = t119 * t77 + t122 * t78;
t75 = -t119 * t78 + t122 * t77;
t1 = m(7) * (t75 ^ 2 + t76 ^ 2 + t84 ^ 2) / 0.2e1 + m(5) * (t89 ^ 2 + t90 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t79 ^ 2 + t80 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t101 ^ 2 + t113 ^ 2) / 0.2e1 + (t91 * mrSges(6,2) - t79 * mrSges(6,3) + Ifges(6,1) * t93 / 0.2e1) * t93 + (t84 * mrSges(7,2) - t75 * mrSges(7,3) + Ifges(7,1) * t86 / 0.2e1) * t86 + (t113 * mrSges(4,2) - t100 * mrSges(4,3) + Ifges(4,1) * t108 / 0.2e1) * t108 + (t98 * mrSges(5,2) - t89 * mrSges(5,3) + Ifges(5,1) * t103 / 0.2e1) * t103 + (-t91 * mrSges(6,1) + t80 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,2) * t92 / 0.2e1) * t92 + (-t84 * mrSges(7,1) + t76 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,2) * t85 / 0.2e1) * t85 + (-t98 * mrSges(5,1) + t90 * mrSges(5,3) + Ifges(5,4) * t103 + Ifges(5,2) * t102 / 0.2e1) * t102 + (t79 * mrSges(6,1) - t80 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(6,6) * t92 + Ifges(6,3) * t106 / 0.2e1) * t106 + (t75 * mrSges(7,1) - t76 * mrSges(7,2) + Ifges(7,5) * t86 + Ifges(7,6) * t85 + Ifges(7,3) * t104 / 0.2e1) * t104 + (t100 * mrSges(4,1) - t101 * mrSges(4,2) + Ifges(4,5) * t108 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t121 + Ifges(3,6) * t124 + (-mrSges(3,1) * t121 - mrSges(3,2) * t124) * pkin(7)) * qJD(1)) * qJD(2) + (t113 * mrSges(4,1) + t89 * mrSges(5,1) - t90 * mrSges(5,2) - t101 * mrSges(4,3) - Ifges(4,4) * t108 + Ifges(5,5) * t103 - Ifges(4,6) * qJD(2) + Ifges(5,6) * t102 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t107) * t107 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t121 ^ 2 + t124 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t130) * t124) * t124 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t124 + (Ifges(3,1) / 0.2e1 + t130) * t121) * t121) * qJD(1) ^ 2;
T  = t1;
