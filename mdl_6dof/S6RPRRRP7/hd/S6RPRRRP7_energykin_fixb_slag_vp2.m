% Calculate kinetic energy for
% S6RPRRRP7
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:18:33
% EndTime: 2019-03-09 06:18:34
% DurationCPUTime: 0.61s
% Computational Cost: add. (683->100), mult. (1599->147), div. (0->0), fcn. (1182->8), ass. (0->41)
t110 = cos(pkin(10));
t124 = t110 ^ 2;
t123 = m(3) / 0.2e1;
t122 = cos(qJ(5));
t121 = pkin(7) + qJ(2);
t111 = sin(qJ(5));
t112 = sin(qJ(4));
t114 = cos(qJ(4));
t113 = sin(qJ(3));
t115 = cos(qJ(3));
t118 = qJD(1) * t110;
t109 = sin(pkin(10));
t119 = qJD(1) * t109;
t100 = -t113 * t119 + t115 * t118;
t101 = (t109 * t115 + t110 * t113) * qJD(1);
t104 = qJD(2) + (-pkin(2) * t110 - pkin(1)) * qJD(1);
t86 = -pkin(3) * t100 - pkin(8) * t101 + t104;
t102 = t121 * t119;
t103 = t121 * t118;
t91 = -t113 * t102 + t115 * t103;
t89 = qJD(3) * pkin(8) + t91;
t79 = -t112 * t89 + t114 * t86;
t94 = qJD(3) * t112 + t101 * t114;
t96 = qJD(4) - t100;
t76 = pkin(4) * t96 - pkin(9) * t94 + t79;
t80 = t112 * t86 + t114 * t89;
t93 = qJD(3) * t114 - t101 * t112;
t78 = pkin(9) * t93 + t80;
t73 = t111 * t76 + t122 * t78;
t90 = -t102 * t115 - t113 * t103;
t88 = -qJD(3) * pkin(3) - t90;
t72 = -t111 * t78 + t122 * t76;
t81 = -pkin(4) * t93 + t88;
t106 = -qJD(1) * pkin(1) + qJD(2);
t95 = qJD(5) + t96;
t83 = t111 * t93 + t122 * t94;
t82 = t111 * t94 - t122 * t93;
t74 = pkin(5) * t82 - qJ(6) * t83 + t81;
t71 = qJ(6) * t95 + t73;
t70 = -t95 * pkin(5) + qJD(6) - t72;
t1 = t106 ^ 2 * t123 + m(4) * (t104 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t81 ^ 2) / 0.2e1 + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,3) * t96 / 0.2e1) * t96 + (t104 * mrSges(4,2) - t90 * mrSges(4,3) + Ifges(4,1) * t101 / 0.2e1) * t101 + (t88 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,5) * t96 + Ifges(5,1) * t94 / 0.2e1) * t94 + (-t104 * mrSges(4,1) + t91 * mrSges(4,3) + Ifges(4,4) * t101 + Ifges(4,2) * t100 / 0.2e1) * t100 + (-t88 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,6) * t96 + Ifges(5,2) * t93 / 0.2e1) * t93 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + Ifges(4,5) * t101 + Ifges(4,6) * t100 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t72 * mrSges(6,1) - t70 * mrSges(7,1) - t73 * mrSges(6,2) + t71 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t95) * t95 + (t81 * mrSges(6,2) + t70 * mrSges(7,2) - t72 * mrSges(6,3) - t74 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t83 + (Ifges(7,4) + Ifges(6,5)) * t95) * t83 + (t81 * mrSges(6,1) + t74 * mrSges(7,1) - t71 * mrSges(7,2) - t73 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t82 + (-Ifges(6,6) + Ifges(7,6)) * t95 + (-Ifges(6,4) + Ifges(7,5)) * t83) * t82 + (t106 * (-mrSges(3,1) * t110 + mrSges(3,2) * t109) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t123 + mrSges(3,3)) * (t109 ^ 2 + t124) * qJ(2) + Ifges(3,2) * t124 / 0.2e1 + (Ifges(3,4) * t110 + Ifges(3,1) * t109 / 0.2e1) * t109) * qJD(1)) * qJD(1);
T  = t1;
