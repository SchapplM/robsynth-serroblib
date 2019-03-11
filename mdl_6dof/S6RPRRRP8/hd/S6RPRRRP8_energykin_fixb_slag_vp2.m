% Calculate kinetic energy for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:22:56
% EndTime: 2019-03-09 06:22:56
% DurationCPUTime: 0.39s
% Computational Cost: add. (433->90), mult. (798->132), div. (0->0), fcn. (482->6), ass. (0->34)
t96 = sin(qJ(4));
t97 = sin(qJ(3));
t98 = cos(qJ(4));
t99 = cos(qJ(3));
t84 = (t96 * t99 + t97 * t98) * qJD(1);
t105 = cos(qJ(5));
t88 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t104 = t88 * mrSges(4,3);
t103 = qJD(1) / 0.2e1;
t102 = -pkin(8) * qJD(1) + t88;
t81 = qJD(3) * pkin(3) + t102 * t99;
t82 = t102 * t97;
t74 = t96 * t81 + t98 * t82;
t92 = qJD(3) + qJD(4);
t71 = pkin(9) * t92 + t74;
t85 = (-t96 * t97 + t98 * t99) * qJD(1);
t94 = qJD(1) * qJ(2);
t86 = t97 * qJD(1) * pkin(3) + t94;
t75 = pkin(4) * t84 - pkin(9) * t85 + t86;
t95 = sin(qJ(5));
t68 = t105 * t71 + t95 * t75;
t73 = t81 * t98 - t96 * t82;
t70 = -pkin(4) * t92 - t73;
t67 = t105 * t75 - t95 * t71;
t100 = qJD(1) ^ 2;
t93 = t100 * qJ(2) ^ 2;
t91 = -qJD(1) * pkin(1) + qJD(2);
t83 = qJD(5) + t84;
t77 = t105 * t85 + t95 * t92;
t76 = -t105 * t92 + t85 * t95;
t66 = pkin(5) * t76 - qJ(6) * t77 + t70;
t65 = qJ(6) * t83 + t68;
t64 = -t83 * pkin(5) + qJD(6) - t67;
t1 = m(3) * (t91 ^ 2 + t93) / 0.2e1 + m(4) * (t93 + (t97 ^ 2 + t99 ^ 2) * t88 ^ 2) / 0.2e1 + m(5) * (t73 ^ 2 + t74 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t70 ^ 2) / 0.2e1 + m(7) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (t73 * mrSges(5,1) - t74 * mrSges(5,2) + Ifges(5,3) * t92 / 0.2e1) * t92 + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + qJ(2) * mrSges(3,3)) * t100 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t99 * mrSges(4,1) - t97 * mrSges(4,2)) * t88) * qJD(3) + (t91 * mrSges(3,2) + (mrSges(4,2) * t94 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t103 - t104) * t99) * t99 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t99) * qJD(1) + (Ifges(4,2) * t103 - t104) * t97) * t97) * qJD(1) + (t86 * mrSges(5,2) - t73 * mrSges(5,3) + Ifges(5,5) * t92 + Ifges(5,1) * t85 / 0.2e1) * t85 - (-t86 * mrSges(5,1) + t74 * mrSges(5,3) + Ifges(5,4) * t85 + Ifges(5,6) * t92 - Ifges(5,2) * t84 / 0.2e1) * t84 + (t67 * mrSges(6,1) - t64 * mrSges(7,1) - t68 * mrSges(6,2) + t65 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t83) * t83 + (t70 * mrSges(6,2) + t64 * mrSges(7,2) - t67 * mrSges(6,3) - t66 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t77 + (Ifges(7,4) + Ifges(6,5)) * t83) * t77 + (t70 * mrSges(6,1) + t66 * mrSges(7,1) - t65 * mrSges(7,2) - t68 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t76 + (-Ifges(6,6) + Ifges(7,6)) * t83 + (-Ifges(6,4) + Ifges(7,5)) * t77) * t76;
T  = t1;
