% Calculate kinetic energy for
% S5RRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:31:54
% EndTime: 2019-12-31 22:31:55
% DurationCPUTime: 0.45s
% Computational Cost: add. (662->90), mult. (1503->148), div. (0->0), fcn. (1154->10), ass. (0->43)
t112 = sin(qJ(4));
t116 = cos(qJ(4));
t124 = cos(pkin(5)) * qJD(1);
t108 = qJD(2) + t124;
t113 = sin(qJ(3));
t117 = cos(qJ(3));
t114 = sin(qJ(2));
t109 = sin(pkin(5));
t123 = t109 * qJD(1);
t121 = t114 * t123;
t100 = t108 * t113 + t117 * t121;
t118 = cos(qJ(2));
t120 = t118 * t123;
t104 = qJD(3) - t120;
t122 = pkin(1) * t124;
t102 = pkin(7) * t120 + t114 * t122;
t97 = pkin(8) * t108 + t102;
t98 = (-pkin(2) * t118 - pkin(8) * t114 - pkin(1)) * t123;
t87 = -t113 * t97 + t117 * t98;
t82 = pkin(3) * t104 - pkin(9) * t100 + t87;
t88 = t113 * t98 + t117 * t97;
t99 = t108 * t117 - t113 * t121;
t84 = pkin(9) * t99 + t88;
t79 = t112 * t82 + t116 * t84;
t78 = -t112 * t84 + t116 * t82;
t90 = -t100 * t112 + t116 * t99;
t101 = -pkin(7) * t121 + t118 * t122;
t96 = -pkin(2) * t108 - t101;
t92 = -pkin(3) * t99 + t96;
t119 = qJD(1) ^ 2;
t115 = cos(qJ(5));
t111 = sin(qJ(5));
t103 = qJD(4) + t104;
t91 = t100 * t116 + t112 * t99;
t89 = qJD(5) - t90;
t86 = t103 * t111 + t115 * t91;
t85 = t103 * t115 - t111 * t91;
t80 = -pkin(4) * t90 - pkin(10) * t91 + t92;
t77 = pkin(10) * t103 + t79;
t76 = -pkin(4) * t103 - t78;
t75 = t111 * t80 + t115 * t77;
t74 = -t111 * t77 + t115 * t80;
t1 = m(3) * (pkin(1) ^ 2 * t109 ^ 2 * t119 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + t119 * Ifges(2,3) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t92 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + (-t96 * mrSges(4,1) + t88 * mrSges(4,3) + Ifges(4,2) * t99 / 0.2e1) * t99 + (t92 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t91 / 0.2e1) * t91 + (t74 * mrSges(6,1) - t75 * mrSges(6,2) + Ifges(6,3) * t89 / 0.2e1) * t89 + (-t92 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t90 / 0.2e1) * t90 + (t76 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,5) * t89 + Ifges(6,1) * t86 / 0.2e1) * t86 + (t87 * mrSges(4,1) - t88 * mrSges(4,2) + Ifges(4,6) * t99 + Ifges(4,3) * t104 / 0.2e1) * t104 + (-t76 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t86 + Ifges(6,6) * t89 + Ifges(6,2) * t85 / 0.2e1) * t85 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t118 / 0.2e1) * t118 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t118 + Ifges(3,1) * t114 / 0.2e1) * t114) * t123 + (-t101 * t114 + t102 * t118) * mrSges(3,3)) * t123 + (t101 * mrSges(3,1) - t102 * mrSges(3,2) + Ifges(3,3) * t108 / 0.2e1 + (Ifges(3,5) * t114 + Ifges(3,6) * t118) * t123) * t108 + (t78 * mrSges(5,1) - t79 * mrSges(5,2) + Ifges(5,5) * t91 + Ifges(5,6) * t90 + Ifges(5,3) * t103 / 0.2e1) * t103 + (t96 * mrSges(4,2) - t87 * mrSges(4,3) + Ifges(4,4) * t99 + Ifges(4,5) * t104 + Ifges(4,1) * t100 / 0.2e1) * t100;
T = t1;
