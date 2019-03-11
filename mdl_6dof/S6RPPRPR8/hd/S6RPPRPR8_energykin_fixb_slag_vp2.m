% Calculate kinetic energy for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:14
% EndTime: 2019-03-09 01:55:14
% DurationCPUTime: 0.38s
% Computational Cost: add. (324->84), mult. (668->122), div. (0->0), fcn. (392->6), ass. (0->36)
t95 = cos(pkin(9));
t108 = t95 ^ 2;
t107 = m(3) / 0.2e1;
t106 = pkin(4) + pkin(8);
t105 = sin(qJ(4));
t86 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t102 = -pkin(7) * qJD(1) + t86;
t94 = sin(pkin(9));
t79 = t102 * t94;
t80 = t102 * t95;
t98 = cos(qJ(4));
t73 = t105 * t80 + t98 * t79;
t104 = t94 ^ 2 + t108;
t103 = qJD(1) * t94;
t89 = qJD(1) * qJ(2) + qJD(3);
t84 = pkin(3) * t103 + t89;
t72 = -t105 * t79 + t80 * t98;
t83 = qJD(1) * t95 * t98 - t105 * t103;
t101 = qJD(5) - t72;
t70 = -qJD(4) * qJ(5) - t73;
t100 = -qJ(5) * t83 + t84;
t97 = cos(qJ(6));
t96 = sin(qJ(6));
t90 = -qJD(1) * pkin(1) + qJD(2);
t82 = (t105 * t95 + t94 * t98) * qJD(1);
t81 = qJD(6) + t83;
t75 = qJD(4) * t97 + t82 * t96;
t74 = -qJD(4) * t96 + t82 * t97;
t71 = pkin(4) * t82 + t100;
t69 = -qJD(4) * pkin(4) + t101;
t68 = t106 * t82 + t100;
t67 = -pkin(5) * t82 - t70;
t66 = pkin(5) * t83 - t106 * qJD(4) + t101;
t65 = t66 * t96 + t68 * t97;
t64 = t66 * t97 - t68 * t96;
t1 = t90 ^ 2 * t107 + m(5) * (t72 ^ 2 + t73 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t104 * t86 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(7) * (t64 ^ 2 + t65 ^ 2 + t67 ^ 2) / 0.2e1 + (t64 * mrSges(7,1) - t65 * mrSges(7,2) + Ifges(7,3) * t81 / 0.2e1) * t81 + (t67 * mrSges(7,2) - t64 * mrSges(7,3) + Ifges(7,5) * t81 + Ifges(7,1) * t75 / 0.2e1) * t75 + (-t67 * mrSges(7,1) + t65 * mrSges(7,3) + Ifges(7,4) * t75 + Ifges(7,6) * t81 + Ifges(7,2) * t74 / 0.2e1) * t74 + (t69 * mrSges(6,1) + t84 * mrSges(5,2) - t72 * mrSges(5,3) - t71 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t83) * t83 + (t84 * mrSges(5,1) + t70 * mrSges(6,1) - t71 * mrSges(6,2) - t73 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t82 + (-Ifges(5,4) - Ifges(6,6)) * t83) * t82 + (t72 * mrSges(5,1) - t73 * mrSges(5,2) + t69 * mrSges(6,2) - t70 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(4) + (-Ifges(6,4) + Ifges(5,5)) * t83 + (Ifges(6,5) - Ifges(5,6)) * t82) * qJD(4) + (t89 * (mrSges(4,1) * t94 + mrSges(4,2) * t95) + t90 * mrSges(3,2) - t104 * t86 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t107 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t108 / 0.2e1 + (-Ifges(4,4) * t95 + Ifges(4,2) * t94 / 0.2e1) * t94) * qJD(1)) * qJD(1);
T  = t1;
