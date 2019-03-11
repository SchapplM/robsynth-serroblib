% Calculate kinetic energy for
% S6RPRRRP2
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
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:58:42
% EndTime: 2019-03-09 05:58:43
% DurationCPUTime: 0.46s
% Computational Cost: add. (456->93), mult. (933->139), div. (0->0), fcn. (588->8), ass. (0->35)
t113 = m(3) / 0.2e1;
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t101 = sin(pkin(10));
t96 = (pkin(1) * t101 + pkin(7)) * qJD(1);
t91 = t105 * qJD(2) + t108 * t96;
t88 = qJD(3) * pkin(8) + t91;
t102 = cos(pkin(10));
t111 = -pkin(1) * t102 - pkin(2);
t89 = (-pkin(3) * t108 - pkin(8) * t105 + t111) * qJD(1);
t79 = -t104 * t88 + t107 * t89;
t112 = qJD(1) * t105;
t95 = t104 * qJD(3) + t107 * t112;
t99 = -t108 * qJD(1) + qJD(4);
t75 = t99 * pkin(4) - t95 * pkin(9) + t79;
t80 = t104 * t89 + t107 * t88;
t94 = t107 * qJD(3) - t104 * t112;
t78 = t94 * pkin(9) + t80;
t72 = t103 * t75 + t106 * t78;
t71 = -t103 * t78 + t106 * t75;
t90 = t108 * qJD(2) - t105 * t96;
t87 = -qJD(3) * pkin(3) - t90;
t81 = -t94 * pkin(4) + t87;
t98 = qJD(5) + t99;
t97 = t111 * qJD(1);
t83 = t103 * t94 + t106 * t95;
t82 = -t103 * t95 + t106 * t94;
t76 = -t82 * pkin(5) + qJD(6) + t81;
t70 = t82 * qJ(6) + t72;
t69 = t98 * pkin(5) - t83 * qJ(6) + t71;
t1 = qJD(2) ^ 2 * t113 + m(4) * (t90 ^ 2 + t91 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t81 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t76 ^ 2) / 0.2e1 + (t79 * mrSges(5,1) - t80 * mrSges(5,2) + Ifges(5,3) * t99 / 0.2e1) * t99 + (t87 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,5) * t99 + Ifges(5,1) * t95 / 0.2e1) * t95 + (-t87 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t95 + Ifges(5,6) * t99 + Ifges(5,2) * t94 / 0.2e1) * t94 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t71 * mrSges(6,1) + t69 * mrSges(7,1) - t72 * mrSges(6,2) - t70 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t98) * t98 + (t81 * mrSges(6,2) + t76 * mrSges(7,2) - t71 * mrSges(6,3) - t69 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t83 + (Ifges(6,5) + Ifges(7,5)) * t98) * t83 + (-t81 * mrSges(6,1) - t76 * mrSges(7,1) + t72 * mrSges(6,3) + t70 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t82 + (Ifges(6,6) + Ifges(7,6)) * t98 + (Ifges(6,4) + Ifges(7,4)) * t83) * t82 + (t97 * (-mrSges(4,1) * t108 + mrSges(4,2) * t105) + (-t90 * t105 + t91 * t108) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t105 + Ifges(4,6) * t108) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t102 * mrSges(3,1) - t101 * mrSges(3,2) + (t101 ^ 2 + t102 ^ 2) * t113 * pkin(1)) * pkin(1) + Ifges(4,2) * t108 ^ 2 / 0.2e1 + (Ifges(4,4) * t108 + Ifges(4,1) * t105 / 0.2e1) * t105) * qJD(1)) * qJD(1);
T  = t1;
