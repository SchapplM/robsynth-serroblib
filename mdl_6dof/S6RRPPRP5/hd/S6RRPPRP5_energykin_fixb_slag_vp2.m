% Calculate kinetic energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:30
% EndTime: 2019-03-09 08:41:30
% DurationCPUTime: 0.50s
% Computational Cost: add. (456->104), mult. (941->139), div. (0->0), fcn. (554->6), ass. (0->34)
t107 = pkin(7) * mrSges(3,3);
t106 = cos(qJ(5));
t105 = -pkin(2) - qJ(4);
t97 = sin(qJ(2));
t104 = t97 * qJD(1);
t101 = -qJ(3) * t97 - pkin(1);
t98 = cos(qJ(2));
t80 = (t105 * t98 + t101) * qJD(1);
t102 = pkin(7) * t104 + qJD(3);
t81 = pkin(3) * t104 + t105 * qJD(2) + t102;
t94 = sin(pkin(9));
t95 = cos(pkin(9));
t72 = -t80 * t94 + t95 * t81;
t103 = t98 * qJD(1);
t86 = qJD(2) * t95 - t94 * t103;
t69 = pkin(4) * t104 - pkin(8) * t86 + t72;
t73 = t95 * t80 + t94 * t81;
t85 = -qJD(2) * t94 - t95 * t103;
t71 = pkin(8) * t85 + t73;
t96 = sin(qJ(5));
t66 = t106 * t71 + t96 * t69;
t88 = -pkin(7) * t103 - qJD(2) * qJ(3);
t83 = pkin(3) * t103 + qJD(4) - t88;
t65 = t106 * t69 - t96 * t71;
t76 = -pkin(4) * t85 + t83;
t89 = qJD(5) + t104;
t87 = -qJD(2) * pkin(2) + t102;
t84 = (-pkin(2) * t98 + t101) * qJD(1);
t75 = t106 * t86 + t96 * t85;
t74 = -t106 * t85 + t86 * t96;
t67 = pkin(5) * t74 - qJ(6) * t75 + t76;
t64 = qJ(6) * t89 + t66;
t63 = -t89 * pkin(5) + qJD(6) - t65;
t1 = m(4) * (t84 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t76 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + (t83 * mrSges(5,2) - t72 * mrSges(5,3) + Ifges(5,1) * t86 / 0.2e1) * t86 + (-t83 * mrSges(5,1) + t73 * mrSges(5,3) + Ifges(5,4) * t86 + Ifges(5,2) * t85 / 0.2e1) * t85 + (t87 * mrSges(4,2) - t88 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(2)) * qJD(2) + (t65 * mrSges(6,1) - t63 * mrSges(7,1) - t66 * mrSges(6,2) + t64 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t89) * t89 + (t76 * mrSges(6,2) + t63 * mrSges(7,2) - t65 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t75 + (Ifges(7,4) + Ifges(6,5)) * t89) * t75 + (t76 * mrSges(6,1) + t67 * mrSges(7,1) - t64 * mrSges(7,2) - t66 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t74 + (-Ifges(6,6) + Ifges(7,6)) * t89 + (-Ifges(6,4) + Ifges(7,5)) * t75) * t74 + ((-t88 * mrSges(4,1) + t84 * mrSges(4,2) + (-pkin(7) * mrSges(3,2) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t98 + (t87 * mrSges(4,1) + t72 * mrSges(5,1) - t73 * mrSges(5,2) - t84 * mrSges(4,3) + Ifges(5,5) * t86 + Ifges(5,6) * t85 + (-pkin(7) * mrSges(3,1) - Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t97 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t97 ^ 2 + t98 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + t107) * t98) * t98 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + t107) * t97 + (Ifges(3,4) + Ifges(4,6)) * t98) * t97) * qJD(1)) * qJD(1);
T  = t1;
