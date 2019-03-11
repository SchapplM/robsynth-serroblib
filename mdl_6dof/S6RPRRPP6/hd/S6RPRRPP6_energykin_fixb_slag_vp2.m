% Calculate kinetic energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:17
% EndTime: 2019-03-09 04:46:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (421->91), mult. (810->131), div. (0->0), fcn. (488->6), ass. (0->32)
t87 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t100 = t87 * mrSges(4,3);
t99 = qJD(1) / 0.2e1;
t93 = sin(qJ(3));
t95 = cos(qJ(3));
t81 = (pkin(3) * t93 - pkin(8) * t95 + qJ(2)) * qJD(1);
t82 = qJD(3) * pkin(8) + t93 * t87;
t92 = sin(qJ(4));
t94 = cos(qJ(4));
t72 = t94 * t81 - t82 * t92;
t97 = t95 * qJD(1);
t85 = qJD(3) * t92 + t94 * t97;
t88 = t93 * qJD(1) + qJD(4);
t69 = pkin(4) * t88 - qJ(5) * t85 + t72;
t73 = t92 * t81 + t94 * t82;
t84 = qJD(3) * t94 - t92 * t97;
t71 = qJ(5) * t84 + t73;
t91 = sin(pkin(9));
t98 = cos(pkin(9));
t66 = t91 * t69 + t98 * t71;
t83 = -qJD(3) * pkin(3) - t95 * t87;
t65 = t69 * t98 - t91 * t71;
t76 = -pkin(4) * t84 + qJD(5) + t83;
t96 = qJD(1) ^ 2;
t90 = t96 * qJ(2) ^ 2;
t89 = -qJD(1) * pkin(1) + qJD(2);
t75 = t91 * t84 + t85 * t98;
t74 = -t84 * t98 + t85 * t91;
t67 = pkin(5) * t74 - qJ(6) * t75 + t76;
t64 = qJ(6) * t88 + t66;
t63 = -t88 * pkin(5) + qJD(6) - t65;
t1 = m(4) * (t90 + (t93 ^ 2 + t95 ^ 2) * t87 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t90) / 0.2e1 + m(5) * (t72 ^ 2 + t73 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t65 ^ 2 + t66 ^ 2 + t76 ^ 2) / 0.2e1 + m(7) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t96 + (t83 * mrSges(5,2) - t72 * mrSges(5,3) + Ifges(5,1) * t85 / 0.2e1) * t85 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t95 * mrSges(4,1) - t93 * mrSges(4,2)) * t87) * qJD(3) + (t89 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t99 - t100) * t95) * t95 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t95) * qJD(1) + (Ifges(4,2) * t99 - t100) * t93) * t93) * qJD(1) + (-t83 * mrSges(5,1) + t73 * mrSges(5,3) + Ifges(5,4) * t85 + Ifges(5,2) * t84 / 0.2e1) * t84 + (t76 * mrSges(6,2) + t63 * mrSges(7,2) - t65 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t75) * t75 + (t76 * mrSges(6,1) + t67 * mrSges(7,1) - t64 * mrSges(7,2) - t66 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t74 + (-Ifges(6,4) + Ifges(7,5)) * t75) * t74 + (t72 * mrSges(5,1) + t65 * mrSges(6,1) - t63 * mrSges(7,1) - t73 * mrSges(5,2) - t66 * mrSges(6,2) + t64 * mrSges(7,3) + Ifges(5,5) * t85 + Ifges(5,6) * t84 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t88 + (Ifges(7,4) + Ifges(6,5)) * t75 + (-Ifges(6,6) + Ifges(7,6)) * t74) * t88;
T  = t1;
