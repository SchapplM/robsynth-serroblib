% Calculate kinetic energy for
% S6RPRRRP4
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
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:06:53
% EndTime: 2019-03-09 06:06:54
% DurationCPUTime: 0.61s
% Computational Cost: add. (699->100), mult. (1697->147), div. (0->0), fcn. (1284->8), ass. (0->39)
t108 = cos(pkin(10));
t121 = t108 ^ 2;
t120 = qJD(1) * (pkin(7) + qJ(2));
t119 = m(3) / 0.2e1;
t109 = sin(qJ(5));
t112 = cos(qJ(5));
t106 = qJD(3) + qJD(4);
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t107 = sin(pkin(10));
t100 = t107 * t120;
t101 = t108 * t120;
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t92 = -t114 * t100 - t101 * t111;
t99 = (t107 * t114 + t108 * t111) * qJD(1);
t85 = qJD(3) * pkin(3) - pkin(8) * t99 + t92;
t93 = -t111 * t100 + t114 * t101;
t98 = (-t107 * t111 + t108 * t114) * qJD(1);
t86 = pkin(8) * t98 + t93;
t80 = t110 * t85 + t113 * t86;
t76 = pkin(9) * t106 + t80;
t90 = -t110 * t99 + t113 * t98;
t91 = t110 * t98 + t113 * t99;
t102 = qJD(2) + (-pkin(2) * t108 - pkin(1)) * qJD(1);
t94 = -pkin(3) * t98 + t102;
t81 = -pkin(4) * t90 - pkin(9) * t91 + t94;
t72 = t109 * t81 + t112 * t76;
t71 = -t109 * t76 + t112 * t81;
t79 = -t110 * t86 + t113 * t85;
t75 = -pkin(4) * t106 - t79;
t103 = -qJD(1) * pkin(1) + qJD(2);
t89 = qJD(5) - t90;
t88 = t106 * t109 + t112 * t91;
t87 = t106 * t112 - t109 * t91;
t73 = -pkin(5) * t87 + qJD(6) + t75;
t70 = qJ(6) * t87 + t72;
t69 = pkin(5) * t89 - qJ(6) * t88 + t71;
t1 = t103 ^ 2 * t119 + m(4) * (t102 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t75 ^ 2) / 0.2e1 + (t102 * mrSges(4,2) - t92 * mrSges(4,3) + Ifges(4,1) * t99 / 0.2e1) * t99 + (t94 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,1) * t91 / 0.2e1) * t91 + (-t102 * mrSges(4,1) + t93 * mrSges(4,3) + Ifges(4,4) * t99 + Ifges(4,2) * t98 / 0.2e1) * t98 + (-t94 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t90 / 0.2e1) * t90 + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,5) * t91 + Ifges(5,6) * t90 + Ifges(5,3) * t106 / 0.2e1) * t106 + (t92 * mrSges(4,1) - t93 * mrSges(4,2) + Ifges(4,5) * t99 + Ifges(4,6) * t98 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t71 * mrSges(6,1) + t69 * mrSges(7,1) - t72 * mrSges(6,2) - t70 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t89) * t89 + (t75 * mrSges(6,2) + t73 * mrSges(7,2) - t71 * mrSges(6,3) - t69 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t88 + (Ifges(6,5) + Ifges(7,5)) * t89) * t88 + (-t75 * mrSges(6,1) - t73 * mrSges(7,1) + t72 * mrSges(6,3) + t70 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t87 + (Ifges(6,6) + Ifges(7,6)) * t89 + (Ifges(6,4) + Ifges(7,4)) * t88) * t87 + (t103 * (-mrSges(3,1) * t108 + mrSges(3,2) * t107) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t119 + mrSges(3,3)) * (t107 ^ 2 + t121) * qJ(2) + Ifges(3,2) * t121 / 0.2e1 + (Ifges(3,4) * t108 + Ifges(3,1) * t107 / 0.2e1) * t107) * qJD(1)) * qJD(1);
T  = t1;
