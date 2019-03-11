% Calculate kinetic energy for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:03
% EndTime: 2019-03-09 09:45:04
% DurationCPUTime: 0.55s
% Computational Cost: add. (706->106), mult. (1641->152), div. (0->0), fcn. (1186->8), ass. (0->36)
t115 = qJD(1) * (pkin(7) + qJ(3));
t102 = sin(pkin(9));
t103 = cos(pkin(9));
t105 = sin(qJ(2));
t107 = cos(qJ(2));
t94 = (t102 * t105 - t103 * t107) * qJD(1);
t113 = pkin(7) * mrSges(3,3);
t101 = sin(pkin(10));
t111 = cos(pkin(10));
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t100 = qJD(3) + (-pkin(2) * t107 - pkin(1)) * qJD(1);
t95 = (t102 * t107 + t103 * t105) * qJD(1);
t83 = pkin(3) * t94 - pkin(8) * t95 + t100;
t98 = qJD(2) * pkin(2) - t105 * t115;
t99 = t107 * t115;
t88 = t102 * t98 + t103 * t99;
t86 = qJD(2) * pkin(8) + t88;
t76 = -t104 * t86 + t106 * t83;
t91 = qJD(2) * t104 + t106 * t95;
t93 = qJD(4) + t94;
t73 = pkin(4) * t93 - qJ(5) * t91 + t76;
t77 = t104 * t83 + t106 * t86;
t90 = qJD(2) * t106 - t104 * t95;
t75 = qJ(5) * t90 + t77;
t70 = t101 * t73 + t111 * t75;
t87 = -t102 * t99 + t103 * t98;
t85 = -qJD(2) * pkin(3) - t87;
t69 = -t101 * t75 + t111 * t73;
t78 = -pkin(4) * t90 + qJD(5) + t85;
t80 = t101 * t90 + t111 * t91;
t79 = t101 * t91 - t111 * t90;
t71 = pkin(5) * t79 - qJ(6) * t80 + t78;
t68 = qJ(6) * t93 + t70;
t67 = -t93 * pkin(5) + qJD(6) - t69;
t1 = m(4) * (t100 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t76 ^ 2 + t77 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t78 ^ 2) / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + (t100 * mrSges(4,2) - t87 * mrSges(4,3) + Ifges(4,1) * t95 / 0.2e1) * t95 + (t85 * mrSges(5,2) - t76 * mrSges(5,3) + Ifges(5,1) * t91 / 0.2e1) * t91 - (-t100 * mrSges(4,1) + t88 * mrSges(4,3) + Ifges(4,4) * t95 - Ifges(4,2) * t94 / 0.2e1) * t94 + (-t85 * mrSges(5,1) + t77 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t90 / 0.2e1) * t90 + (t87 * mrSges(4,1) - t88 * mrSges(4,2) + Ifges(4,5) * t95 - Ifges(4,6) * t94 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t105 + Ifges(3,6) * t107 + (-mrSges(3,1) * t105 - mrSges(3,2) * t107) * pkin(7)) * qJD(1)) * qJD(2) + (t78 * mrSges(6,2) + t67 * mrSges(7,2) - t69 * mrSges(6,3) - t71 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t80) * t80 + (t78 * mrSges(6,1) + t71 * mrSges(7,1) - t68 * mrSges(7,2) - t70 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t79 + (-Ifges(6,4) + Ifges(7,5)) * t80) * t79 + (t76 * mrSges(5,1) + t69 * mrSges(6,1) - t67 * mrSges(7,1) - t77 * mrSges(5,2) - t70 * mrSges(6,2) + t68 * mrSges(7,3) + Ifges(5,5) * t91 + Ifges(5,6) * t90 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t93 + (Ifges(7,4) + Ifges(6,5)) * t80 + (-Ifges(6,6) + Ifges(7,6)) * t79) * t93 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t105 ^ 2 + t107 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t113 + Ifges(3,2) / 0.2e1) * t107) * t107 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t107 + (t113 + Ifges(3,1) / 0.2e1) * t105) * t105) * qJD(1) ^ 2;
T  = t1;
