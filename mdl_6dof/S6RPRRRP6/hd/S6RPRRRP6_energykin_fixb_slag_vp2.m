% Calculate kinetic energy for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:09
% EndTime: 2019-03-09 06:14:10
% DurationCPUTime: 0.57s
% Computational Cost: add. (687->100), mult. (1619->147), div. (0->0), fcn. (1202->8), ass. (0->41)
t110 = cos(pkin(10));
t124 = t110 ^ 2;
t123 = m(3) / 0.2e1;
t122 = pkin(7) + qJ(2);
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t112 = sin(qJ(4));
t115 = cos(qJ(4));
t113 = sin(qJ(3));
t116 = cos(qJ(3));
t119 = qJD(1) * t110;
t109 = sin(pkin(10));
t120 = qJD(1) * t109;
t100 = -t113 * t120 + t116 * t119;
t101 = (t109 * t116 + t110 * t113) * qJD(1);
t104 = qJD(2) + (-pkin(2) * t110 - pkin(1)) * qJD(1);
t87 = -pkin(3) * t100 - pkin(8) * t101 + t104;
t102 = t122 * t120;
t103 = t122 * t119;
t92 = -t113 * t102 + t116 * t103;
t90 = qJD(3) * pkin(8) + t92;
t80 = -t112 * t90 + t115 * t87;
t94 = qJD(3) * t112 + t101 * t115;
t96 = qJD(4) - t100;
t76 = pkin(4) * t96 - pkin(9) * t94 + t80;
t81 = t112 * t87 + t115 * t90;
t93 = qJD(3) * t115 - t101 * t112;
t79 = pkin(9) * t93 + t81;
t73 = t111 * t76 + t114 * t79;
t72 = -t111 * t79 + t114 * t76;
t91 = -t102 * t116 - t113 * t103;
t89 = -qJD(3) * pkin(3) - t91;
t82 = -pkin(4) * t93 + t89;
t106 = -qJD(1) * pkin(1) + qJD(2);
t95 = qJD(5) + t96;
t84 = t111 * t93 + t114 * t94;
t83 = -t111 * t94 + t114 * t93;
t77 = -pkin(5) * t83 + qJD(6) + t82;
t71 = qJ(6) * t83 + t73;
t70 = pkin(5) * t95 - qJ(6) * t84 + t72;
t1 = m(5) * (t80 ^ 2 + t81 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t82 ^ 2) / 0.2e1 + t106 ^ 2 * t123 + m(4) * (t104 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t77 ^ 2) / 0.2e1 + (t80 * mrSges(5,1) - t81 * mrSges(5,2) + Ifges(5,3) * t96 / 0.2e1) * t96 + (t104 * mrSges(4,2) - t91 * mrSges(4,3) + Ifges(4,1) * t101 / 0.2e1) * t101 + (t89 * mrSges(5,2) - t80 * mrSges(5,3) + Ifges(5,5) * t96 + Ifges(5,1) * t94 / 0.2e1) * t94 + (-t104 * mrSges(4,1) + t92 * mrSges(4,3) + Ifges(4,4) * t101 + Ifges(4,2) * t100 / 0.2e1) * t100 + (-t89 * mrSges(5,1) + t81 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,6) * t96 + Ifges(5,2) * t93 / 0.2e1) * t93 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + Ifges(4,5) * t101 + Ifges(4,6) * t100 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t72 * mrSges(6,1) + t70 * mrSges(7,1) - t73 * mrSges(6,2) - t71 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t95) * t95 + (t82 * mrSges(6,2) + t77 * mrSges(7,2) - t72 * mrSges(6,3) - t70 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t84 + (Ifges(6,5) + Ifges(7,5)) * t95) * t84 + (-t82 * mrSges(6,1) - t77 * mrSges(7,1) + t73 * mrSges(6,3) + t71 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t83 + (Ifges(6,6) + Ifges(7,6)) * t95 + (Ifges(6,4) + Ifges(7,4)) * t84) * t83 + (t106 * (-mrSges(3,1) * t110 + mrSges(3,2) * t109) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t123 + mrSges(3,3)) * (t109 ^ 2 + t124) * qJ(2) + Ifges(3,2) * t124 / 0.2e1 + (Ifges(3,4) * t110 + Ifges(3,1) * t109 / 0.2e1) * t109) * qJD(1)) * qJD(1);
T  = t1;
