% Calculate kinetic energy for
% S6RPRRRP5
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
% Datum: 2019-03-09 06:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:10:31
% EndTime: 2019-03-09 06:10:32
% DurationCPUTime: 0.60s
% Computational Cost: add. (697->100), mult. (1693->147), div. (0->0), fcn. (1280->8), ass. (0->39)
t108 = cos(pkin(10));
t121 = t108 ^ 2;
t120 = qJD(1) * (pkin(7) + qJ(2));
t119 = m(3) / 0.2e1;
t118 = cos(qJ(5));
t109 = sin(qJ(5));
t106 = qJD(3) + qJD(4);
t110 = sin(qJ(4));
t112 = cos(qJ(4));
t100 = t108 * t120;
t111 = sin(qJ(3));
t113 = cos(qJ(3));
t107 = sin(pkin(10));
t99 = t107 * t120;
t91 = -t100 * t111 - t113 * t99;
t98 = (t107 * t113 + t108 * t111) * qJD(1);
t84 = qJD(3) * pkin(3) - pkin(8) * t98 + t91;
t92 = t113 * t100 - t111 * t99;
t97 = (-t107 * t111 + t108 * t113) * qJD(1);
t85 = pkin(8) * t97 + t92;
t79 = t110 * t84 + t112 * t85;
t76 = pkin(9) * t106 + t79;
t89 = -t110 * t98 + t112 * t97;
t90 = t110 * t97 + t112 * t98;
t101 = qJD(2) + (-pkin(2) * t108 - pkin(1)) * qJD(1);
t93 = -pkin(3) * t97 + t101;
t80 = -pkin(4) * t89 - pkin(9) * t90 + t93;
t72 = t109 * t80 + t118 * t76;
t78 = -t110 * t85 + t112 * t84;
t75 = -pkin(4) * t106 - t78;
t71 = -t109 * t76 + t118 * t80;
t103 = -qJD(1) * pkin(1) + qJD(2);
t88 = qJD(5) - t89;
t87 = t109 * t106 + t118 * t90;
t86 = -t118 * t106 + t109 * t90;
t73 = pkin(5) * t86 - qJ(6) * t87 + t75;
t70 = qJ(6) * t88 + t72;
t69 = -t88 * pkin(5) + qJD(6) - t71;
t1 = t103 ^ 2 * t119 + m(5) * (t78 ^ 2 + t79 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t75 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + (t101 * mrSges(4,2) - t91 * mrSges(4,3) + Ifges(4,1) * t98 / 0.2e1) * t98 + (t93 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t90 / 0.2e1) * t90 + (-t101 * mrSges(4,1) + t92 * mrSges(4,3) + Ifges(4,4) * t98 + Ifges(4,2) * t97 / 0.2e1) * t97 + (-t93 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t90 + Ifges(5,2) * t89 / 0.2e1) * t89 + (t78 * mrSges(5,1) - t79 * mrSges(5,2) + Ifges(5,5) * t90 + Ifges(5,6) * t89 + Ifges(5,3) * t106 / 0.2e1) * t106 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + Ifges(4,5) * t98 + Ifges(4,6) * t97 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t71 * mrSges(6,1) - t69 * mrSges(7,1) - t72 * mrSges(6,2) + t70 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t88) * t88 + (t75 * mrSges(6,2) + t69 * mrSges(7,2) - t71 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t87 + (Ifges(7,4) + Ifges(6,5)) * t88) * t87 + (t75 * mrSges(6,1) + t73 * mrSges(7,1) - t70 * mrSges(7,2) - t72 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t86 + (-Ifges(6,6) + Ifges(7,6)) * t88 + (-Ifges(6,4) + Ifges(7,5)) * t87) * t86 + (t103 * (-mrSges(3,1) * t108 + mrSges(3,2) * t107) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t119 + mrSges(3,3)) * (t107 ^ 2 + t121) * qJ(2) + Ifges(3,2) * t121 / 0.2e1 + (Ifges(3,4) * t108 + Ifges(3,1) * t107 / 0.2e1) * t107) * qJD(1)) * qJD(1);
T  = t1;
