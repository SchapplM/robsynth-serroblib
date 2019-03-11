% Calculate kinetic energy for
% S6RPPRRR5
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
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:28:11
% EndTime: 2019-03-09 02:28:12
% DurationCPUTime: 0.32s
% Computational Cost: add. (280->77), mult. (505->121), div. (0->0), fcn. (266->6), ass. (0->32)
t91 = sin(qJ(5));
t92 = sin(qJ(4));
t94 = cos(qJ(5));
t95 = cos(qJ(4));
t78 = (t91 * t95 + t92 * t94) * qJD(1);
t83 = (pkin(1) + qJ(3)) * qJD(1) - qJD(2);
t103 = t83 ^ 2;
t102 = m(3) / 0.2e1;
t86 = qJD(1) * qJ(2) + qJD(3);
t82 = -qJD(1) * pkin(7) + t86;
t101 = t82 * mrSges(5,3);
t100 = qJD(1) / 0.2e1;
t98 = -pkin(8) * qJD(1) + t82;
t75 = qJD(4) * pkin(4) + t98 * t95;
t76 = t98 * t92;
t70 = t91 * t75 + t94 * t76;
t69 = t75 * t94 - t76 * t91;
t80 = t92 * qJD(1) * pkin(4) + t83;
t93 = cos(qJ(6));
t90 = sin(qJ(6));
t88 = qJD(4) + qJD(5);
t87 = -qJD(1) * pkin(1) + qJD(2);
t79 = (-t91 * t92 + t94 * t95) * qJD(1);
t77 = qJD(6) + t78;
t72 = t79 * t93 + t88 * t90;
t71 = -t79 * t90 + t88 * t93;
t68 = pkin(5) * t78 - pkin(9) * t79 + t80;
t67 = pkin(9) * t88 + t70;
t66 = -pkin(5) * t88 - t69;
t65 = t67 * t93 + t68 * t90;
t64 = -t67 * t90 + t68 * t93;
t1 = m(5) * (t103 + (t92 ^ 2 + t95 ^ 2) * t82 ^ 2) / 0.2e1 + t87 ^ 2 * t102 + m(6) * (t69 ^ 2 + t70 ^ 2 + t80 ^ 2) / 0.2e1 + m(4) * (t86 ^ 2 + t103) / 0.2e1 + m(7) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t69 * mrSges(6,1) - t70 * mrSges(6,2) + Ifges(6,3) * t88 / 0.2e1) * t88 + (t64 * mrSges(7,1) - t65 * mrSges(7,2) + Ifges(7,3) * t77 / 0.2e1) * t77 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t95 * mrSges(5,1) - t92 * mrSges(5,2)) * t82) * qJD(4) + (t80 * mrSges(6,2) - t69 * mrSges(6,3) + Ifges(6,5) * t88 + Ifges(6,1) * t79 / 0.2e1) * t79 + (t66 * mrSges(7,2) - t64 * mrSges(7,3) + Ifges(7,5) * t77 + Ifges(7,1) * t72 / 0.2e1) * t72 - (-t80 * mrSges(6,1) + t70 * mrSges(6,3) + Ifges(6,4) * t79 + Ifges(6,6) * t88 - Ifges(6,2) * t78 / 0.2e1) * t78 + (-t66 * mrSges(7,1) + t65 * mrSges(7,3) + Ifges(7,4) * t72 + Ifges(7,6) * t77 + Ifges(7,2) * t71 / 0.2e1) * t71 + (t87 * mrSges(3,2) + t86 * mrSges(4,2) + t83 * mrSges(4,3) + (t83 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t100 - t101) * t95) * t95 + (t83 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t100 - t101) * t92) * t92 + (Ifges(2,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t102 + mrSges(3,3)) * qJ(2) - Ifges(5,4) * t95 * t92) * qJD(1)) * qJD(1);
T  = t1;
