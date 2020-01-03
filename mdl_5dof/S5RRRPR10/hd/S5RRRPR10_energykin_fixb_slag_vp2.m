% Calculate kinetic energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:33
% EndTime: 2019-12-31 21:26:33
% DurationCPUTime: 0.38s
% Computational Cost: add. (644->90), mult. (1503->146), div. (0->0), fcn. (1154->10), ass. (0->42)
t107 = sin(pkin(10));
t109 = cos(pkin(10));
t116 = cos(qJ(2));
t108 = sin(pkin(5));
t122 = qJD(1) * t108;
t118 = t116 * t122;
t102 = qJD(3) - t118;
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t113 = sin(qJ(2));
t121 = cos(pkin(5)) * qJD(1);
t120 = pkin(1) * t121;
t101 = pkin(7) * t118 + t113 * t120;
t106 = qJD(2) + t121;
t96 = pkin(8) * t106 + t101;
t97 = (-pkin(2) * t116 - pkin(8) * t113 - pkin(1)) * t122;
t86 = -t112 * t96 + t115 * t97;
t119 = t113 * t122;
t99 = t106 * t112 + t115 * t119;
t81 = pkin(3) * t102 - qJ(4) * t99 + t86;
t87 = t112 * t97 + t115 * t96;
t98 = t106 * t115 - t112 * t119;
t83 = qJ(4) * t98 + t87;
t78 = t107 * t81 + t109 * t83;
t77 = -t107 * t83 + t109 * t81;
t90 = -t107 * t99 + t109 * t98;
t100 = -pkin(7) * t119 + t116 * t120;
t95 = -pkin(2) * t106 - t100;
t89 = -pkin(3) * t98 + qJD(4) + t95;
t117 = qJD(1) ^ 2;
t114 = cos(qJ(5));
t111 = sin(qJ(5));
t91 = t107 * t98 + t109 * t99;
t88 = qJD(5) - t90;
t85 = t102 * t111 + t114 * t91;
t84 = t102 * t114 - t111 * t91;
t79 = -pkin(4) * t90 - pkin(9) * t91 + t89;
t76 = pkin(9) * t102 + t78;
t75 = -pkin(4) * t102 - t77;
t74 = t111 * t79 + t114 * t76;
t73 = -t111 * t76 + t114 * t79;
t1 = m(3) * (pkin(1) ^ 2 * t108 ^ 2 * t117 + t100 ^ 2 + t101 ^ 2) / 0.2e1 + t117 * Ifges(2,3) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + (t95 * mrSges(4,2) - t86 * mrSges(4,3) + Ifges(4,1) * t99 / 0.2e1) * t99 + (t89 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,1) * t91 / 0.2e1) * t91 + (t73 * mrSges(6,1) - t74 * mrSges(6,2) + Ifges(6,3) * t88 / 0.2e1) * t88 + (-t95 * mrSges(4,1) + t87 * mrSges(4,3) + Ifges(4,4) * t99 + Ifges(4,2) * t98 / 0.2e1) * t98 + (-t89 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t90 / 0.2e1) * t90 + (t75 * mrSges(6,2) - t73 * mrSges(6,3) + Ifges(6,5) * t88 + Ifges(6,1) * t85 / 0.2e1) * t85 + (-t75 * mrSges(6,1) + t74 * mrSges(6,3) + Ifges(6,4) * t85 + Ifges(6,6) * t88 + Ifges(6,2) * t84 / 0.2e1) * t84 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t116 / 0.2e1) * t116 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t116 + Ifges(3,1) * t113 / 0.2e1) * t113) * t122 + (-t100 * t113 + t101 * t116) * mrSges(3,3)) * t122 + (t100 * mrSges(3,1) - t101 * mrSges(3,2) + Ifges(3,3) * t106 / 0.2e1 + (Ifges(3,5) * t113 + Ifges(3,6) * t116) * t122) * t106 + (t86 * mrSges(4,1) + t77 * mrSges(5,1) - t87 * mrSges(4,2) - t78 * mrSges(5,2) + Ifges(4,5) * t99 + Ifges(5,5) * t91 + Ifges(4,6) * t98 + Ifges(5,6) * t90 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t102) * t102;
T = t1;
