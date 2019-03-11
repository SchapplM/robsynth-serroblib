% Calculate kinetic energy for
% S6RRPPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-03-09 08:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:33:29
% EndTime: 2019-03-09 08:33:29
% DurationCPUTime: 0.37s
% Computational Cost: add. (282->101), mult. (553->125), div. (0->0), fcn. (246->4), ass. (0->29)
t101 = -pkin(2) - pkin(3);
t100 = pkin(7) * mrSges(3,3);
t92 = cos(qJ(2));
t98 = qJD(1) * t92;
t90 = sin(qJ(2));
t99 = qJD(1) * t90;
t75 = -qJD(1) * pkin(1) - pkin(2) * t98 - qJ(3) * t99;
t71 = pkin(3) * t98 + qJD(4) - t75;
t68 = (pkin(4) * t90 + pkin(8) * t92) * qJD(1) + t71;
t96 = qJ(4) * qJD(1);
t97 = pkin(7) * t99 + qJD(3);
t95 = -t90 * t96 + t97;
t70 = (-pkin(8) + t101) * qJD(2) + t95;
t89 = sin(qJ(5));
t91 = cos(qJ(5));
t64 = t89 * t68 + t91 * t70;
t79 = pkin(7) * t98 + qJD(2) * qJ(3);
t63 = t91 * t68 - t70 * t89;
t74 = t92 * t96 - t79;
t73 = qJD(2) * pkin(4) - t74;
t80 = qJD(5) + t99;
t78 = -qJD(2) * pkin(2) + t97;
t77 = -qJD(2) * t89 - t91 * t98;
t76 = -qJD(2) * t91 + t89 * t98;
t72 = qJD(2) * t101 + t95;
t67 = -pkin(5) * t76 + qJD(6) + t73;
t62 = qJ(6) * t76 + t64;
t61 = pkin(5) * t80 - qJ(6) * t77 + t63;
t1 = m(4) * (t75 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t72 ^ 2 + t74 ^ 2) / 0.2e1 + m(7) * (t61 ^ 2 + t62 ^ 2 + t67 ^ 2) / 0.2e1 + (t63 * mrSges(6,1) + t61 * mrSges(7,1) - t64 * mrSges(6,2) - t62 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t80) * t80 + (t73 * mrSges(6,2) + t67 * mrSges(7,2) - t63 * mrSges(6,3) - t61 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t77 + (Ifges(6,5) + Ifges(7,5)) * t80) * t77 + (-t73 * mrSges(6,1) - t67 * mrSges(7,1) + t64 * mrSges(6,3) + t62 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t76 + (Ifges(6,6) + Ifges(7,6)) * t80 + (Ifges(6,4) + Ifges(7,4)) * t77) * t76 + (-t78 * mrSges(4,1) - t74 * mrSges(5,1) + t72 * mrSges(5,2) + t79 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t75 * mrSges(4,1) + t79 * mrSges(4,2) - t71 * mrSges(5,2) + t74 * mrSges(5,3) + (-pkin(7) * mrSges(3,2) + Ifges(5,5) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t92 + (t71 * mrSges(5,1) + t78 * mrSges(4,2) - t75 * mrSges(4,3) - t72 * mrSges(5,3) + (-pkin(7) * mrSges(3,1) + Ifges(4,4) + Ifges(3,5) + Ifges(5,6)) * qJD(2)) * t90 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t90 ^ 2 + t92 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + t100) * t92) * t92 + (-pkin(1) * mrSges(3,2) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1 + t100) * t90 + (Ifges(3,4) + Ifges(5,4) - Ifges(4,5)) * t92) * t90) * qJD(1)) * qJD(1);
T  = t1;
