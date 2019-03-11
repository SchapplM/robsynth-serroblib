% Calculate kinetic energy for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:49
% EndTime: 2019-03-09 04:53:49
% DurationCPUTime: 0.27s
% Computational Cost: add. (285->88), mult. (522->116), div. (0->0), fcn. (260->4), ass. (0->29)
t90 = cos(qJ(4));
t75 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t89 = t75 * mrSges(4,3);
t88 = qJD(1) / 0.2e1;
t87 = pkin(4) + qJ(6);
t81 = sin(qJ(3));
t82 = cos(qJ(3));
t69 = (pkin(3) * t81 - pkin(8) * t82 + qJ(2)) * qJD(1);
t70 = qJD(3) * pkin(8) + t81 * t75;
t80 = sin(qJ(4));
t65 = t80 * t69 + t90 * t70;
t86 = t82 * qJD(1);
t76 = t81 * qJD(1) + qJD(4);
t62 = -qJ(5) * t76 - t65;
t64 = t90 * t69 - t80 * t70;
t71 = -qJD(3) * pkin(3) - t82 * t75;
t85 = qJD(5) - t64;
t73 = t80 * qJD(3) + t90 * t86;
t84 = -qJ(5) * t73 + t71;
t83 = qJD(1) ^ 2;
t79 = t83 * qJ(2) ^ 2;
t77 = -qJD(1) * pkin(1) + qJD(2);
t72 = -t90 * qJD(3) + t80 * t86;
t63 = pkin(4) * t72 + t84;
t61 = -t76 * pkin(4) + t85;
t60 = t87 * t72 + t84;
t59 = -pkin(5) * t72 + qJD(6) - t62;
t58 = t73 * pkin(5) - t87 * t76 + t85;
t1 = m(4) * (t79 + (t81 ^ 2 + t82 ^ 2) * t75 ^ 2) / 0.2e1 + m(3) * (t77 ^ 2 + t79) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t71 ^ 2) / 0.2e1 + m(7) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t83 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t82 * mrSges(4,1) - t81 * mrSges(4,2)) * t75) * qJD(3) + (t77 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t88 - t89) * t82) * t82 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t82) * qJD(1) + (Ifges(4,2) * t88 - t89) * t81) * t81) * qJD(1) + (t64 * mrSges(5,1) - t65 * mrSges(5,2) + t61 * mrSges(6,2) + t59 * mrSges(7,2) - t62 * mrSges(6,3) - t58 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * t76) * t76 + (t61 * mrSges(6,1) + t58 * mrSges(7,1) + t71 * mrSges(5,2) - t60 * mrSges(7,2) - t64 * mrSges(5,3) - t63 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t73 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t76) * t73 + (t71 * mrSges(5,1) + t62 * mrSges(6,1) - t59 * mrSges(7,1) - t63 * mrSges(6,2) - t65 * mrSges(5,3) + t60 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t72 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t76 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t73) * t72;
T  = t1;
