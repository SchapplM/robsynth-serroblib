% Calculate kinetic energy for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:38
% EndTime: 2019-03-09 01:31:38
% DurationCPUTime: 0.31s
% Computational Cost: add. (287->76), mult. (577->121), div. (0->0), fcn. (328->8), ass. (0->34)
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t97 = sin(pkin(10));
t99 = cos(pkin(10));
t87 = (t102 * t99 + t104 * t97) * qJD(1);
t111 = m(3) / 0.2e1;
t100 = cos(pkin(9));
t109 = -pkin(1) * t100 - pkin(2);
t89 = qJD(3) + (-qJ(4) + t109) * qJD(1);
t81 = -t97 * qJD(2) + t99 * t89;
t77 = -t99 * qJD(1) * pkin(7) + t81;
t110 = qJD(1) * t97;
t82 = t99 * qJD(2) + t97 * t89;
t78 = -pkin(7) * t110 + t82;
t73 = t102 * t77 + t104 * t78;
t98 = sin(pkin(9));
t92 = (-pkin(1) * t98 - qJ(3)) * qJD(1);
t90 = qJD(4) - t92;
t86 = pkin(4) * t110 + t90;
t72 = -t102 * t78 + t104 * t77;
t105 = qJD(2) ^ 2;
t103 = cos(qJ(6));
t101 = sin(qJ(6));
t91 = qJD(1) * t109 + qJD(3);
t88 = (-t102 * t97 + t104 * t99) * qJD(1);
t85 = qJD(6) + t87;
t80 = t101 * qJD(5) + t103 * t88;
t79 = t103 * qJD(5) - t101 * t88;
t74 = t87 * pkin(5) - t88 * pkin(8) + t86;
t71 = qJD(5) * pkin(8) + t73;
t70 = -qJD(5) * pkin(5) - t72;
t69 = t101 * t74 + t103 * t71;
t68 = -t101 * t71 + t103 * t74;
t1 = m(4) * (t91 ^ 2 + t92 ^ 2 + t105) / 0.2e1 + t105 * t111 + m(6) * (t72 ^ 2 + t73 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t90 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + (t86 * mrSges(6,2) - t72 * mrSges(6,3) + Ifges(6,1) * t88 / 0.2e1) * t88 + (t68 * mrSges(7,1) - t69 * mrSges(7,2) + Ifges(7,3) * t85 / 0.2e1) * t85 - (-t86 * mrSges(6,1) + t73 * mrSges(6,3) + Ifges(6,4) * t88 - Ifges(6,2) * t87 / 0.2e1) * t87 + (t70 * mrSges(7,2) - t68 * mrSges(7,3) + Ifges(7,5) * t85 + Ifges(7,1) * t80 / 0.2e1) * t80 + (-t70 * mrSges(7,1) + t69 * mrSges(7,3) + Ifges(7,4) * t80 + Ifges(7,6) * t85 + Ifges(7,2) * t79 / 0.2e1) * t79 + (t72 * mrSges(6,1) - t73 * mrSges(6,2) + Ifges(6,5) * t88 - Ifges(6,6) * t87 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t90 * (mrSges(5,1) * t97 + mrSges(5,2) * t99) + t91 * mrSges(4,2) - t92 * mrSges(4,3) + (-t81 * t99 - t82 * t97) * mrSges(5,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + (t100 * mrSges(3,1) - t98 * mrSges(3,2) + (t100 ^ 2 + t98 ^ 2) * t111 * pkin(1)) * pkin(1) + Ifges(5,1) * t99 ^ 2 / 0.2e1 + (-Ifges(5,4) * t99 + Ifges(5,2) * t97 / 0.2e1) * t97) * qJD(1)) * qJD(1);
T  = t1;
