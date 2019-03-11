% Calculate kinetic energy for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:24:25
% EndTime: 2019-03-09 03:24:25
% DurationCPUTime: 0.38s
% Computational Cost: add. (397->90), mult. (798->131), div. (0->0), fcn. (482->6), ass. (0->33)
t93 = sin(pkin(9));
t94 = cos(pkin(9));
t96 = sin(qJ(3));
t97 = cos(qJ(3));
t83 = (t93 * t97 + t94 * t96) * qJD(1);
t103 = cos(qJ(5));
t87 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t102 = t87 * mrSges(4,3);
t101 = qJD(1) / 0.2e1;
t100 = -qJ(4) * qJD(1) + t87;
t80 = qJD(3) * pkin(3) + t100 * t97;
t81 = t100 * t96;
t73 = t93 * t80 + t94 * t81;
t70 = qJD(3) * pkin(8) + t73;
t84 = (-t93 * t96 + t94 * t97) * qJD(1);
t92 = qJD(1) * qJ(2);
t85 = t96 * qJD(1) * pkin(3) + qJD(4) + t92;
t74 = t83 * pkin(4) - t84 * pkin(8) + t85;
t95 = sin(qJ(5));
t66 = t103 * t70 + t95 * t74;
t72 = t94 * t80 - t93 * t81;
t69 = -qJD(3) * pkin(4) - t72;
t65 = t103 * t74 - t95 * t70;
t98 = qJD(1) ^ 2;
t91 = t98 * qJ(2) ^ 2;
t89 = -qJD(1) * pkin(1) + qJD(2);
t82 = qJD(5) + t83;
t76 = t95 * qJD(3) + t103 * t84;
t75 = -t103 * qJD(3) + t95 * t84;
t67 = t75 * pkin(5) - t76 * qJ(6) + t69;
t64 = t82 * qJ(6) + t66;
t63 = -t82 * pkin(5) + qJD(6) - t65;
t1 = m(4) * (t91 + (t96 ^ 2 + t97 ^ 2) * t87 ^ 2) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t85 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t91) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t98 + (t85 * mrSges(5,2) - t72 * mrSges(5,3) + Ifges(5,1) * t84 / 0.2e1) * t84 - (-t85 * mrSges(5,1) + t73 * mrSges(5,3) + Ifges(5,4) * t84 - Ifges(5,2) * t83 / 0.2e1) * t83 + (t65 * mrSges(6,1) - t63 * mrSges(7,1) - t66 * mrSges(6,2) + t64 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t82) * t82 + (t89 * mrSges(3,2) + (mrSges(4,2) * t92 + (Ifges(4,1) * t101 - t102) * t97) * t97 + ((qJ(2) * mrSges(4,1) - Ifges(4,4) * t97) * qJD(1) + (Ifges(4,2) * t101 - t102) * t96) * t96) * qJD(1) + (t72 * mrSges(5,1) - t73 * mrSges(5,2) + Ifges(5,5) * t84 - Ifges(5,6) * t83 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (t97 * mrSges(4,1) - t96 * mrSges(4,2)) * t87 + (Ifges(4,5) * t97 - Ifges(4,6) * t96) * qJD(1)) * qJD(3) + (t69 * mrSges(6,2) + t63 * mrSges(7,2) - t65 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t76 + (Ifges(7,4) + Ifges(6,5)) * t82) * t76 + (t69 * mrSges(6,1) + t67 * mrSges(7,1) - t64 * mrSges(7,2) - t66 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t75 + (-Ifges(6,6) + Ifges(7,6)) * t82 + (-Ifges(6,4) + Ifges(7,5)) * t76) * t75;
T  = t1;
