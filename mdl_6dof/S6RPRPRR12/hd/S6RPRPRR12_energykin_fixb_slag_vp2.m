% Calculate kinetic energy for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR12_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:06
% EndTime: 2019-03-09 04:18:06
% DurationCPUTime: 0.36s
% Computational Cost: add. (339->94), mult. (618->136), div. (0->0), fcn. (314->6), ass. (0->35)
t89 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t106 = t89 * mrSges(4,3);
t101 = cos(qJ(5));
t102 = cos(qJ(3));
t77 = qJD(4) + (pkin(4) * qJD(1) - t89) * t102 + (-pkin(3) - pkin(8)) * qJD(3);
t99 = sin(qJ(3));
t104 = qJD(1) * t99;
t96 = qJD(1) * qJ(2);
t105 = pkin(3) * t104 + t96;
t79 = (pkin(8) * t99 - qJ(4) * t102) * qJD(1) + t105;
t98 = sin(qJ(5));
t71 = t101 * t79 + t98 * t77;
t83 = -qJD(3) * qJ(4) - t99 * t89;
t93 = t102 * qJD(1);
t90 = t93 + qJD(5);
t70 = t101 * t77 - t79 * t98;
t80 = -pkin(4) * t104 - t83;
t103 = qJD(1) ^ 2;
t100 = cos(qJ(6));
t97 = sin(qJ(6));
t95 = t103 * qJ(2) ^ 2;
t92 = -qJD(1) * pkin(1) + qJD(2);
t88 = qJD(6) + t90;
t85 = qJD(3) * t101 + t98 * t104;
t84 = -qJD(3) * t98 + t101 * t104;
t82 = -qJ(4) * t93 + t105;
t81 = -qJD(3) * pkin(3) - t102 * t89 + qJD(4);
t74 = t100 * t85 + t84 * t97;
t73 = t100 * t84 - t85 * t97;
t72 = -pkin(5) * t84 + t80;
t69 = pkin(9) * t84 + t71;
t68 = pkin(5) * t90 - pkin(9) * t85 + t70;
t67 = t100 * t69 + t68 * t97;
t66 = t100 * t68 - t69 * t97;
t1 = m(3) * (t92 ^ 2 + t95) / 0.2e1 + m(4) * (t95 + (t102 ^ 2 + t99 ^ 2) * t89 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t80 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(7) * (t66 ^ 2 + t67 ^ 2 + t72 ^ 2) / 0.2e1 + (t70 * mrSges(6,1) - t71 * mrSges(6,2) + Ifges(6,3) * t90 / 0.2e1) * t90 + (t66 * mrSges(7,1) - t67 * mrSges(7,2) + Ifges(7,3) * t88 / 0.2e1) * t88 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t103 + (t80 * mrSges(6,2) - t70 * mrSges(6,3) + Ifges(6,5) * t90 + Ifges(6,1) * t85 / 0.2e1) * t85 + (t72 * mrSges(7,2) - t66 * mrSges(7,3) + Ifges(7,5) * t88 + Ifges(7,1) * t74 / 0.2e1) * t74 + (-t80 * mrSges(6,1) + t71 * mrSges(6,3) + Ifges(6,4) * t85 + Ifges(6,6) * t90 + Ifges(6,2) * t84 / 0.2e1) * t84 + (-t72 * mrSges(7,1) + t67 * mrSges(7,3) + Ifges(7,4) * t74 + Ifges(7,6) * t88 + Ifges(7,2) * t73 / 0.2e1) * t73 + (t81 * mrSges(5,2) - t83 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3) + (mrSges(4,1) * t102 - mrSges(4,2) * t99) * t89) * qJD(3) + (t92 * mrSges(3,2) + (mrSges(4,1) * t96 + t83 * mrSges(5,1) - t82 * mrSges(5,2) + (-t106 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(1)) * t99 + (Ifges(5,5) - Ifges(4,6)) * qJD(3)) * t99 + (-t82 * mrSges(5,3) + t81 * mrSges(5,1) + (qJ(2) * mrSges(4,2) + (-Ifges(4,4) - Ifges(5,6)) * t99) * qJD(1) + (-t106 + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(1)) * t102 + (-Ifges(5,4) + Ifges(4,5)) * qJD(3)) * t102) * qJD(1);
T  = t1;
