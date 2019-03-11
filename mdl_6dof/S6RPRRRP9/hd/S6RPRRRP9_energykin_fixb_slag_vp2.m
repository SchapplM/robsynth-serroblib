% Calculate kinetic energy for
% S6RPRRRP9
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
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:26:04
% EndTime: 2019-03-09 06:26:04
% DurationCPUTime: 0.35s
% Computational Cost: add. (435->91), mult. (824->133), div. (0->0), fcn. (502->6), ass. (0->33)
t90 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t104 = t90 * mrSges(4,3);
t103 = qJD(1) / 0.2e1;
t100 = cos(qJ(3));
t97 = sin(qJ(3));
t83 = (pkin(3) * t97 - pkin(8) * t100 + qJ(2)) * qJD(1);
t84 = qJD(3) * pkin(8) + t97 * t90;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t75 = t99 * t83 - t84 * t96;
t102 = t100 * qJD(1);
t87 = qJD(3) * t96 + t99 * t102;
t91 = t97 * qJD(1) + qJD(4);
t72 = pkin(4) * t91 - pkin(9) * t87 + t75;
t76 = t96 * t83 + t99 * t84;
t86 = qJD(3) * t99 - t96 * t102;
t74 = pkin(9) * t86 + t76;
t95 = sin(qJ(5));
t98 = cos(qJ(5));
t68 = t95 * t72 + t98 * t74;
t67 = t98 * t72 - t74 * t95;
t85 = -qJD(3) * pkin(3) - t100 * t90;
t79 = -pkin(4) * t86 + t85;
t101 = qJD(1) ^ 2;
t94 = t101 * qJ(2) ^ 2;
t92 = -qJD(1) * pkin(1) + qJD(2);
t89 = qJD(5) + t91;
t78 = t86 * t95 + t87 * t98;
t77 = t86 * t98 - t87 * t95;
t69 = -pkin(5) * t77 + qJD(6) + t79;
t66 = qJ(6) * t77 + t68;
t65 = pkin(5) * t89 - qJ(6) * t78 + t67;
t1 = m(4) * (t94 + (t100 ^ 2 + t97 ^ 2) * t90 ^ 2) / 0.2e1 + m(3) * (t92 ^ 2 + t94) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t75 ^ 2 + t76 ^ 2 + t85 ^ 2) / 0.2e1 + m(7) * (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) / 0.2e1 + (t75 * mrSges(5,1) - t76 * mrSges(5,2) + Ifges(5,3) * t91 / 0.2e1) * t91 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t101 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t100 * mrSges(4,1) - t97 * mrSges(4,2)) * t90) * qJD(3) + (t92 * mrSges(3,2) + (qJ(2) * mrSges(4,1) * qJD(1) - Ifges(4,6) * qJD(3) + (Ifges(4,2) * t103 - t104) * t97) * t97 + (Ifges(4,5) * qJD(3) + (qJ(2) * mrSges(4,2) - Ifges(4,4) * t97) * qJD(1) + (Ifges(4,1) * t103 - t104) * t100) * t100) * qJD(1) + (t85 * mrSges(5,2) - t75 * mrSges(5,3) + Ifges(5,5) * t91 + Ifges(5,1) * t87 / 0.2e1) * t87 + (-t85 * mrSges(5,1) + t76 * mrSges(5,3) + Ifges(5,4) * t87 + Ifges(5,6) * t91 + Ifges(5,2) * t86 / 0.2e1) * t86 + (t67 * mrSges(6,1) + t65 * mrSges(7,1) - t68 * mrSges(6,2) - t66 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t89) * t89 + (t79 * mrSges(6,2) + t69 * mrSges(7,2) - t67 * mrSges(6,3) - t65 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t78 + (Ifges(6,5) + Ifges(7,5)) * t89) * t78 + (-t79 * mrSges(6,1) - t69 * mrSges(7,1) + t68 * mrSges(6,3) + t66 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t77 + (Ifges(6,6) + Ifges(7,6)) * t89 + (Ifges(6,4) + Ifges(7,4)) * t78) * t77;
T  = t1;
