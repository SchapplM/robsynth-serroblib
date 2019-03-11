% Calculate kinetic energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:09
% EndTime: 2019-03-09 09:54:10
% DurationCPUTime: 0.45s
% Computational Cost: add. (522->104), mult. (1141->138), div. (0->0), fcn. (760->6), ass. (0->34)
t105 = pkin(7) * mrSges(3,3);
t104 = cos(qJ(4));
t103 = pkin(4) + qJ(6);
t96 = cos(qJ(2));
t102 = qJD(1) * t96;
t95 = sin(qJ(2));
t84 = (-pkin(2) * t96 - qJ(3) * t95 - pkin(1)) * qJD(1);
t89 = pkin(7) * t102 + qJD(2) * qJ(3);
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t78 = t93 * t84 - t89 * t92;
t101 = t95 * qJD(1);
t86 = qJD(2) * t92 + t93 * t101;
t72 = -pkin(3) * t102 - pkin(8) * t86 + t78;
t79 = t92 * t84 + t93 * t89;
t85 = qJD(2) * t93 - t92 * t101;
t75 = pkin(8) * t85 + t79;
t94 = sin(qJ(4));
t69 = t104 * t75 + t94 * t72;
t90 = qJD(4) - t102;
t67 = -qJ(5) * t90 - t69;
t68 = t104 * t72 - t94 * t75;
t88 = -qJD(2) * pkin(2) + pkin(7) * t101 + qJD(3);
t100 = qJD(5) - t68;
t80 = -pkin(3) * t85 + t88;
t77 = t104 * t86 + t94 * t85;
t99 = -qJ(5) * t77 + t80;
t76 = -t104 * t85 + t86 * t94;
t70 = pkin(4) * t76 + t99;
t66 = -t90 * pkin(4) + t100;
t65 = t103 * t76 + t99;
t64 = -pkin(5) * t76 + qJD(6) - t67;
t63 = t77 * pkin(5) - t103 * t90 + t100;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t78 ^ 2 + t79 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t80 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t66 ^ 2 + t67 ^ 2 + t70 ^ 2) / 0.2e1 + (t88 * mrSges(4,2) - t78 * mrSges(4,3) + Ifges(4,1) * t86 / 0.2e1) * t86 + (-Ifges(4,6) * t102 - t88 * mrSges(4,1) + t79 * mrSges(4,3) + Ifges(4,4) * t86 + Ifges(4,2) * t85 / 0.2e1) * t85 + (t68 * mrSges(5,1) - t69 * mrSges(5,2) + t66 * mrSges(6,2) + t64 * mrSges(7,2) - t67 * mrSges(6,3) - t63 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t90) * t90 + (t66 * mrSges(6,1) + t63 * mrSges(7,1) + t80 * mrSges(5,2) - t65 * mrSges(7,2) - t68 * mrSges(5,3) - t70 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t77 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t90) * t77 + (t80 * mrSges(5,1) + t67 * mrSges(6,1) - t64 * mrSges(7,1) - t70 * mrSges(6,2) - t69 * mrSges(5,3) + t65 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t76 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t90 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t77) * t76 + ((-pkin(7) * mrSges(3,1) + Ifges(3,5)) * qJD(2) * t95 + (-t78 * mrSges(4,1) + t79 * mrSges(4,2) - Ifges(4,5) * t86 + (-pkin(7) * mrSges(3,2) + Ifges(3,6)) * qJD(2)) * t96 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t95 ^ 2 + t96 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + t105) * t95) * t95 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t95 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t105) * t96) * t96) * qJD(1)) * qJD(1);
T  = t1;
