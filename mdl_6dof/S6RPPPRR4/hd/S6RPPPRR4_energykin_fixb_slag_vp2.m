% Calculate kinetic energy for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
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
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:14
% EndTime: 2019-03-09 01:35:15
% DurationCPUTime: 0.19s
% Computational Cost: add. (215->69), mult. (355->101), div. (0->0), fcn. (138->6), ass. (0->31)
t87 = cos(pkin(9));
t95 = qJ(2) * qJD(1);
t80 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t86 = sin(pkin(9));
t97 = t86 * t80;
t77 = t87 * t95 + t97;
t85 = qJD(1) * qJ(4);
t74 = t77 - t85;
t99 = t74 ^ 2;
t98 = m(3) / 0.2e1;
t76 = t80 * t87 - t86 * t95;
t73 = qJD(1) * pkin(3) + qJD(4) - t76;
t72 = qJD(1) * pkin(7) + t73;
t90 = sin(qJ(5));
t92 = cos(qJ(5));
t69 = qJD(3) * t92 + t72 * t90;
t96 = qJD(1) * t92;
t68 = -qJD(3) * t90 + t72 * t92;
t93 = qJD(3) ^ 2;
t91 = cos(qJ(6));
t89 = sin(qJ(6));
t83 = -pkin(1) * qJD(1) + qJD(2);
t81 = -qJD(1) * t90 + qJD(6);
t79 = qJD(5) * t89 - t91 * t96;
t78 = qJD(5) * t91 + t89 * t96;
t70 = t97 - t85 + (-pkin(5) * t90 + pkin(8) * t92 + qJ(2) * t87) * qJD(1);
t67 = qJD(5) * pkin(8) + t69;
t66 = -qJD(5) * pkin(5) - t68;
t65 = t67 * t91 + t70 * t89;
t64 = -t67 * t89 + t70 * t91;
t1 = t83 ^ 2 * t98 + m(5) * (t73 ^ 2 + t93 + t99) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t93) / 0.2e1 + m(7) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t69 ^ 2 + t99) / 0.2e1 + (t64 * mrSges(7,1) - t65 * mrSges(7,2) + Ifges(7,3) * t81 / 0.2e1) * t81 + (t68 * mrSges(6,1) - t69 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t66 * mrSges(7,2) - t64 * mrSges(7,3) + Ifges(7,5) * t81 + Ifges(7,1) * t79 / 0.2e1) * t79 + (-t66 * mrSges(7,1) + t65 * mrSges(7,3) + Ifges(7,4) * t79 + Ifges(7,6) * t81 + Ifges(7,2) * t78 / 0.2e1) * t78 + (-t83 * mrSges(3,1) - t76 * mrSges(4,1) + t77 * mrSges(4,2) - t73 * mrSges(5,2) - t74 * mrSges(5,3) + (-t74 * mrSges(6,2) + t68 * mrSges(6,3) - Ifges(6,5) * qJD(5) + Ifges(6,1) * t96 / 0.2e1) * t92 + (-t74 * mrSges(6,1) + t69 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t90 + (Ifges(2,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + (qJ(2) * t98 + mrSges(3,3)) * qJ(2) + (-Ifges(6,4) * t92 + Ifges(6,2) * t90 / 0.2e1) * t90) * qJD(1)) * qJD(1);
T  = t1;
