% Calculate kinetic energy for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:27
% EndTime: 2019-03-09 02:30:28
% DurationCPUTime: 0.29s
% Computational Cost: add. (290->79), mult. (513->121), div. (0->0), fcn. (266->6), ass. (0->33)
t101 = -pkin(1) - qJ(3);
t85 = -t101 * qJD(1) - qJD(2);
t105 = t85 ^ 2;
t104 = m(3) / 0.2e1;
t89 = qJD(1) * qJ(2) + qJD(3);
t84 = -qJD(1) * pkin(7) + t89;
t103 = t84 * mrSges(5,3);
t102 = qJD(1) / 0.2e1;
t95 = sin(qJ(4));
t98 = cos(qJ(4));
t77 = -qJD(2) + (pkin(4) * t95 - pkin(8) * t98 - t101) * qJD(1);
t79 = qJD(4) * pkin(8) + t84 * t95;
t94 = sin(qJ(5));
t97 = cos(qJ(5));
t71 = t94 * t77 + t97 * t79;
t100 = qJD(1) * t98;
t88 = t95 * qJD(1) + qJD(5);
t70 = t97 * t77 - t79 * t94;
t80 = -qJD(4) * pkin(4) - t84 * t98;
t96 = cos(qJ(6));
t93 = sin(qJ(6));
t90 = -qJD(1) * pkin(1) + qJD(2);
t87 = qJD(6) + t88;
t82 = qJD(4) * t94 + t97 * t100;
t81 = qJD(4) * t97 - t94 * t100;
t74 = -pkin(5) * t81 + t80;
t73 = t81 * t93 + t82 * t96;
t72 = t81 * t96 - t82 * t93;
t69 = pkin(9) * t81 + t71;
t68 = pkin(5) * t88 - pkin(9) * t82 + t70;
t67 = t68 * t93 + t69 * t96;
t66 = t68 * t96 - t69 * t93;
t1 = m(4) * (t89 ^ 2 + t105) / 0.2e1 + m(5) * (t105 + (t95 ^ 2 + t98 ^ 2) * t84 ^ 2) / 0.2e1 + t90 ^ 2 * t104 + m(7) * (t66 ^ 2 + t67 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t80 ^ 2) / 0.2e1 + (t70 * mrSges(6,1) - t71 * mrSges(6,2) + Ifges(6,3) * t88 / 0.2e1) * t88 + (t66 * mrSges(7,1) - t67 * mrSges(7,2) + Ifges(7,3) * t87 / 0.2e1) * t87 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t98 * mrSges(5,1) - t95 * mrSges(5,2)) * t84) * qJD(4) + (t80 * mrSges(6,2) - t70 * mrSges(6,3) + Ifges(6,5) * t88 + Ifges(6,1) * t82 / 0.2e1) * t82 + (t74 * mrSges(7,2) - t66 * mrSges(7,3) + Ifges(7,5) * t87 + Ifges(7,1) * t73 / 0.2e1) * t73 + (-t80 * mrSges(6,1) + t71 * mrSges(6,3) + Ifges(6,4) * t82 + Ifges(6,6) * t88 + Ifges(6,2) * t81 / 0.2e1) * t81 + (-t74 * mrSges(7,1) + t67 * mrSges(7,3) + Ifges(7,4) * t73 + Ifges(7,6) * t87 + Ifges(7,2) * t72 / 0.2e1) * t72 + (t90 * mrSges(3,2) + t89 * mrSges(4,2) + t85 * mrSges(4,3) + (t85 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t102 - t103) * t98) * t98 + (-Ifges(5,4) * t100 + t85 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t102 - t103) * t95) * t95 + (Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t104 + mrSges(3,3)) * qJ(2)) * qJD(1)) * qJD(1);
T  = t1;
