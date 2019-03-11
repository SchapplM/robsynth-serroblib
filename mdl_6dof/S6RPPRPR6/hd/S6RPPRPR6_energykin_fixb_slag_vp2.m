% Calculate kinetic energy for
% S6RPPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% Datum: 2019-03-09 01:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:50:33
% EndTime: 2019-03-09 01:50:33
% DurationCPUTime: 0.20s
% Computational Cost: add. (180->78), mult. (315->108), div. (0->0), fcn. (110->4), ass. (0->30)
t89 = -pkin(1) - qJ(3);
t71 = -t89 * qJD(1) - qJD(2);
t92 = t71 ^ 2;
t91 = m(3) / 0.2e1;
t75 = qJD(1) * qJ(2) + qJD(3);
t70 = -qJD(1) * pkin(7) + t75;
t90 = t70 * mrSges(5,3);
t79 = sin(qJ(4));
t88 = qJD(1) * t79;
t81 = cos(qJ(4));
t87 = qJD(1) * t81;
t86 = pkin(4) * t88 - qJD(2);
t85 = qJD(4) * qJ(5);
t84 = pkin(5) * qJD(1) - t70;
t83 = -qJ(5) * t81 - t89;
t80 = cos(qJ(6));
t78 = sin(qJ(6));
t76 = -qJD(1) * pkin(1) + qJD(2);
t73 = qJD(6) + t87;
t68 = qJD(4) * t80 + t78 * t88;
t67 = -qJD(4) * t78 + t80 * t88;
t66 = -t70 * t79 - t85;
t65 = -qJD(4) * pkin(4) - t70 * t81 + qJD(5);
t64 = t83 * qJD(1) + t86;
t63 = -t84 * t79 + t85;
t62 = qJD(5) + t84 * t81 + (-pkin(4) - pkin(8)) * qJD(4);
t61 = (pkin(8) * t79 + t83) * qJD(1) + t86;
t60 = t61 * t80 + t62 * t78;
t59 = -t61 * t78 + t62 * t80;
t1 = m(4) * (t75 ^ 2 + t92) / 0.2e1 + m(5) * (t92 + (t79 ^ 2 + t81 ^ 2) * t70 ^ 2) / 0.2e1 + t76 ^ 2 * t91 + m(7) * (t59 ^ 2 + t60 ^ 2 + t63 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + (Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t91 + mrSges(3,3)) * qJ(2)) * qJD(1) ^ 2 + (t59 * mrSges(7,1) - t60 * mrSges(7,2) + Ifges(7,3) * t73 / 0.2e1) * t73 + (t63 * mrSges(7,2) - t59 * mrSges(7,3) + Ifges(7,5) * t73 + Ifges(7,1) * t68 / 0.2e1) * t68 + (-t63 * mrSges(7,1) + t60 * mrSges(7,3) + Ifges(7,4) * t68 + Ifges(7,6) * t73 + Ifges(7,2) * t67 / 0.2e1) * t67 + (t65 * mrSges(6,2) - t66 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(4) + (t81 * mrSges(5,1) - t79 * mrSges(5,2)) * t70) * qJD(4) + (t76 * mrSges(3,2) + t75 * mrSges(4,2) + t71 * mrSges(4,3) + (t65 * mrSges(6,1) + t71 * mrSges(5,2) - t64 * mrSges(6,3) + (-t90 + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * qJD(1)) * t81 + (-Ifges(6,4) + Ifges(5,5)) * qJD(4)) * t81 + (t66 * mrSges(6,1) - t64 * mrSges(6,2) + t71 * mrSges(5,1) + (-t90 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(1)) * t79 + (-Ifges(5,4) - Ifges(6,6)) * t87 + (Ifges(6,5) - Ifges(5,6)) * qJD(4)) * t79) * qJD(1);
T  = t1;
