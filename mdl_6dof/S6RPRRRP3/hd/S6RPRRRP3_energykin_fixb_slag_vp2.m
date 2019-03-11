% Calculate kinetic energy for
% S6RPRRRP3
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:02:59
% EndTime: 2019-03-09 06:03:00
% DurationCPUTime: 0.45s
% Computational Cost: add. (454->93), mult. (919->139), div. (0->0), fcn. (574->8), ass. (0->35)
t113 = m(3) / 0.2e1;
t112 = cos(qJ(5));
t103 = sin(qJ(5));
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t101 = sin(pkin(10));
t96 = (pkin(1) * t101 + pkin(7)) * qJD(1);
t90 = t105 * qJD(2) + t107 * t96;
t87 = qJD(3) * pkin(8) + t90;
t102 = cos(pkin(10));
t110 = -pkin(1) * t102 - pkin(2);
t88 = (-pkin(3) * t107 - pkin(8) * t105 + t110) * qJD(1);
t78 = -t104 * t87 + t106 * t88;
t111 = qJD(1) * t105;
t95 = qJD(3) * t104 + t106 * t111;
t99 = -qJD(1) * t107 + qJD(4);
t75 = pkin(4) * t99 - pkin(9) * t95 + t78;
t79 = t104 * t88 + t106 * t87;
t94 = qJD(3) * t106 - t104 * t111;
t77 = pkin(9) * t94 + t79;
t72 = t103 * t75 + t112 * t77;
t89 = qJD(2) * t107 - t105 * t96;
t71 = -t103 * t77 + t112 * t75;
t86 = -qJD(3) * pkin(3) - t89;
t80 = -pkin(4) * t94 + t86;
t98 = qJD(5) + t99;
t97 = t110 * qJD(1);
t82 = t103 * t94 + t112 * t95;
t81 = t103 * t95 - t112 * t94;
t73 = pkin(5) * t81 - qJ(6) * t82 + t80;
t70 = qJ(6) * t98 + t72;
t69 = -t98 * pkin(5) + qJD(6) - t71;
t1 = qJD(2) ^ 2 * t113 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t80 ^ 2) / 0.2e1 + m(4) * (t89 ^ 2 + t90 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t86 ^ 2) / 0.2e1 + (t78 * mrSges(5,1) - t79 * mrSges(5,2) + Ifges(5,3) * t99 / 0.2e1) * t99 + (t86 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,5) * t99 + Ifges(5,1) * t95 / 0.2e1) * t95 + (-t86 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t95 + Ifges(5,6) * t99 + Ifges(5,2) * t94 / 0.2e1) * t94 + (t89 * mrSges(4,1) - t90 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t71 * mrSges(6,1) - t69 * mrSges(7,1) - t72 * mrSges(6,2) + t70 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t98) * t98 + (t80 * mrSges(6,2) + t69 * mrSges(7,2) - t71 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t82 + (Ifges(7,4) + Ifges(6,5)) * t98) * t82 + (t80 * mrSges(6,1) + t73 * mrSges(7,1) - t70 * mrSges(7,2) - t72 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t81 + (-Ifges(6,6) + Ifges(7,6)) * t98 + (-Ifges(6,4) + Ifges(7,5)) * t82) * t81 + (t97 * (-mrSges(4,1) * t107 + mrSges(4,2) * t105) + (-t89 * t105 + t90 * t107) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t105 + Ifges(4,6) * t107) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t102 * mrSges(3,1) - t101 * mrSges(3,2) + (t101 ^ 2 + t102 ^ 2) * t113 * pkin(1)) * pkin(1) + Ifges(4,2) * t107 ^ 2 / 0.2e1 + (Ifges(4,4) * t107 + Ifges(4,1) * t105 / 0.2e1) * t105) * qJD(1)) * qJD(1);
T  = t1;
