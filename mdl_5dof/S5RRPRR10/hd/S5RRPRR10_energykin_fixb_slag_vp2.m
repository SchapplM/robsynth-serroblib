% Calculate kinetic energy for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:42
% EndTime: 2019-12-31 20:23:43
% DurationCPUTime: 0.44s
% Computational Cost: add. (544->91), mult. (1481->148), div. (0->0), fcn. (1136->10), ass. (0->43)
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t112 = sin(qJ(2));
t115 = cos(qJ(2));
t107 = sin(pkin(5));
t122 = qJD(1) * t107;
t97 = (t106 * t112 - t108 * t115) * t122;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t121 = cos(pkin(5)) * qJD(1);
t105 = qJD(2) + t121;
t120 = pkin(1) * t121;
t104 = t115 * t120;
t119 = t112 * t122;
t92 = pkin(2) * t105 + t104 + (-pkin(7) - qJ(3)) * t119;
t118 = t115 * t122;
t100 = pkin(7) * t118 + t112 * t120;
t95 = qJ(3) * t118 + t100;
t83 = t106 * t92 + t108 * t95;
t81 = pkin(8) * t105 + t83;
t101 = qJD(3) + (-pkin(2) * t115 - pkin(1)) * t122;
t98 = (t106 * t115 + t108 * t112) * t122;
t87 = pkin(3) * t97 - pkin(8) * t98 + t101;
t78 = t111 * t87 + t114 * t81;
t82 = -t106 * t95 + t108 * t92;
t77 = -t111 * t81 + t114 * t87;
t90 = t105 * t114 - t111 * t98;
t80 = -pkin(3) * t105 - t82;
t116 = qJD(1) ^ 2;
t113 = cos(qJ(5));
t110 = sin(qJ(5));
t99 = -pkin(7) * t119 + t104;
t96 = qJD(4) + t97;
t91 = t105 * t111 + t114 * t98;
t89 = qJD(5) - t90;
t85 = t110 * t96 + t113 * t91;
t84 = -t110 * t91 + t113 * t96;
t76 = -pkin(4) * t90 - pkin(9) * t91 + t80;
t75 = pkin(9) * t96 + t78;
t74 = -pkin(4) * t96 - t77;
t73 = t110 * t76 + t113 * t75;
t72 = -t110 * t75 + t113 * t76;
t1 = m(3) * (pkin(1) ^ 2 * t107 ^ 2 * t116 + t100 ^ 2 + t99 ^ 2) / 0.2e1 + t116 * Ifges(2,3) / 0.2e1 + m(4) * (t101 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + (t101 * mrSges(4,2) - t82 * mrSges(4,3) + Ifges(4,1) * t98 / 0.2e1) * t98 + (t77 * mrSges(5,1) - t78 * mrSges(5,2) + Ifges(5,3) * t96 / 0.2e1) * t96 + (t72 * mrSges(6,1) - t73 * mrSges(6,2) + Ifges(6,3) * t89 / 0.2e1) * t89 - (-t101 * mrSges(4,1) + t83 * mrSges(4,3) + Ifges(4,4) * t98 - Ifges(4,2) * t97 / 0.2e1) * t97 + (t80 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,5) * t96 + Ifges(5,1) * t91 / 0.2e1) * t91 + (t74 * mrSges(6,2) - t72 * mrSges(6,3) + Ifges(6,5) * t89 + Ifges(6,1) * t85 / 0.2e1) * t85 + (-t80 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,6) * t96 + Ifges(5,2) * t90 / 0.2e1) * t90 + (-t74 * mrSges(6,1) + t73 * mrSges(6,3) + Ifges(6,4) * t85 + Ifges(6,6) * t89 + Ifges(6,2) * t84 / 0.2e1) * t84 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t115 / 0.2e1) * t115 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t115 + Ifges(3,1) * t112 / 0.2e1) * t112) * t122 + (t100 * t115 - t112 * t99) * mrSges(3,3)) * t122 + (t99 * mrSges(3,1) + t82 * mrSges(4,1) - t100 * mrSges(3,2) - t83 * mrSges(4,2) + Ifges(4,5) * t98 - Ifges(4,6) * t97 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t105 + (Ifges(3,5) * t112 + Ifges(3,6) * t115) * t122) * t105;
T = t1;
