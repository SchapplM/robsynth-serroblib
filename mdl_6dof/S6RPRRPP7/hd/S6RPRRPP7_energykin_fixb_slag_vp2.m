% Calculate kinetic energy for
% S6RPRRPP7
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:12
% EndTime: 2019-03-09 04:50:12
% DurationCPUTime: 0.25s
% Computational Cost: add. (285->88), mult. (522->116), div. (0->0), fcn. (260->4), ass. (0->29)
t95 = -pkin(4) - pkin(5);
t94 = cos(qJ(4));
t80 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t93 = t80 * mrSges(4,3);
t92 = qJD(1) / 0.2e1;
t86 = sin(qJ(3));
t87 = cos(qJ(3));
t73 = (pkin(3) * t86 - pkin(8) * t87 + qJ(2)) * qJD(1);
t74 = qJD(3) * pkin(8) + t86 * t80;
t85 = sin(qJ(4));
t69 = t85 * t73 + t94 * t74;
t91 = qJD(1) * t87;
t81 = t86 * qJD(1) + qJD(4);
t66 = t81 * qJ(5) + t69;
t68 = t73 * t94 - t85 * t74;
t75 = -qJD(3) * pkin(3) - t87 * t80;
t90 = qJD(5) - t68;
t77 = t85 * qJD(3) + t91 * t94;
t89 = t77 * qJ(5) - t75;
t88 = qJD(1) ^ 2;
t84 = t88 * qJ(2) ^ 2;
t82 = -qJD(1) * pkin(1) + qJD(2);
t76 = -qJD(3) * t94 + t85 * t91;
t67 = t76 * pkin(4) - t89;
t65 = -t81 * pkin(4) + t90;
t64 = t76 * t95 + qJD(6) + t89;
t63 = t76 * qJ(6) + t66;
t62 = -t77 * qJ(6) + t81 * t95 + t90;
t1 = m(3) * (t82 ^ 2 + t84) / 0.2e1 + m(4) * (t84 + (t86 ^ 2 + t87 ^ 2) * t80 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t75 ^ 2) / 0.2e1 + m(7) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t88 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t87 * mrSges(4,1) - t86 * mrSges(4,2)) * t80) * qJD(3) + (t82 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t92 - t93) * t87) * t87 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t87) * qJD(1) + (Ifges(4,2) * t92 - t93) * t86) * t86) * qJD(1) + (t68 * mrSges(5,1) - t65 * mrSges(6,1) - t62 * mrSges(7,1) - t69 * mrSges(5,2) + t63 * mrSges(7,2) + t66 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t81) * t81 + (t75 * mrSges(5,2) + t65 * mrSges(6,2) + t64 * mrSges(7,2) - t68 * mrSges(5,3) - t67 * mrSges(6,3) - t62 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t77 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t81) * t77 + (t75 * mrSges(5,1) + t67 * mrSges(6,1) - t64 * mrSges(7,1) - t66 * mrSges(6,2) - t69 * mrSges(5,3) + t63 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t76 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t81 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t77) * t76;
T  = t1;
