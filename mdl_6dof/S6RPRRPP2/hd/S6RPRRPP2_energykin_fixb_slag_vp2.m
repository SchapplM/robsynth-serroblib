% Calculate kinetic energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:29
% EndTime: 2019-03-09 04:31:29
% DurationCPUTime: 0.39s
% Computational Cost: add. (298->90), mult. (599->122), div. (0->0), fcn. (322->6), ass. (0->31)
t107 = m(3) / 0.2e1;
t106 = pkin(4) + pkin(5);
t105 = cos(qJ(4));
t93 = sin(pkin(9));
t86 = (pkin(1) * t93 + pkin(7)) * qJD(1);
t96 = sin(qJ(3));
t97 = cos(qJ(3));
t81 = t96 * qJD(2) + t97 * t86;
t78 = qJD(3) * pkin(8) + t81;
t94 = cos(pkin(9));
t103 = -pkin(1) * t94 - pkin(2);
t79 = (-pkin(3) * t97 - pkin(8) * t96 + t103) * qJD(1);
t95 = sin(qJ(4));
t73 = t105 * t78 + t95 * t79;
t80 = t97 * qJD(2) - t96 * t86;
t104 = qJD(1) * t96;
t89 = t97 * qJD(1) - qJD(4);
t70 = -t89 * qJ(5) + t73;
t102 = qJD(3) * pkin(3) + t80;
t72 = t105 * t79 - t95 * t78;
t101 = qJD(5) - t72;
t85 = t95 * qJD(3) + t105 * t104;
t100 = t85 * qJ(5) + t102;
t87 = t103 * qJD(1);
t84 = -t105 * qJD(3) + t95 * t104;
t71 = t84 * pkin(4) - t100;
t69 = t89 * pkin(4) + t101;
t68 = -t106 * t84 + qJD(6) + t100;
t67 = t84 * qJ(6) + t70;
t66 = -t85 * qJ(6) + t106 * t89 + t101;
t1 = qJD(2) ^ 2 * t107 + m(4) * (t80 ^ 2 + t81 ^ 2 + t87 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(7) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + (t80 * mrSges(4,1) - t81 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t72 * mrSges(5,1) + t69 * mrSges(6,1) + t66 * mrSges(7,1) + t73 * mrSges(5,2) - t67 * mrSges(7,2) - t70 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t89) * t89 + (-t102 * mrSges(5,2) + t69 * mrSges(6,2) + t68 * mrSges(7,2) - t72 * mrSges(5,3) - t71 * mrSges(6,3) - t66 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t85 + (-Ifges(6,4) - Ifges(5,5) + Ifges(7,5)) * t89) * t85 + (-t102 * mrSges(5,1) + t71 * mrSges(6,1) - t68 * mrSges(7,1) - t70 * mrSges(6,2) - t73 * mrSges(5,3) + t67 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t84 + (Ifges(5,6) - Ifges(6,6) + Ifges(7,6)) * t89 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t85) * t84 + (t87 * (-mrSges(4,1) * t97 + mrSges(4,2) * t96) + (-t80 * t96 + t81 * t97) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t96 + Ifges(4,6) * t97) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t94 * mrSges(3,1) - t93 * mrSges(3,2) + (t93 ^ 2 + t94 ^ 2) * t107 * pkin(1)) * pkin(1) + Ifges(4,2) * t97 ^ 2 / 0.2e1 + (Ifges(4,4) * t97 + Ifges(4,1) * t96 / 0.2e1) * t96) * qJD(1)) * qJD(1);
T  = t1;
