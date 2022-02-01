% Calculate kinetic energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:45
% EndTime: 2022-01-20 11:16:46
% DurationCPUTime: 0.28s
% Computational Cost: add. (321->64), mult. (465->115), div. (0->0), fcn. (254->8), ass. (0->34)
t87 = cos(pkin(9));
t84 = t87 ^ 2;
t85 = qJD(1) + qJD(2);
t90 = sin(qJ(2));
t98 = pkin(1) * qJD(1);
t80 = t85 * qJ(3) + t90 * t98;
t101 = t80 * t87;
t86 = sin(pkin(9));
t100 = t85 * t86;
t99 = t87 * t85;
t93 = cos(qJ(2));
t96 = -t93 * t98 + qJD(3);
t72 = (-pkin(3) * t87 - pkin(7) * t86 - pkin(2)) * t85 + t96;
t89 = sin(qJ(4));
t92 = cos(qJ(4));
t69 = t92 * t101 + t89 * t72;
t97 = pkin(8) * t100;
t82 = qJD(4) - t99;
t68 = -t101 * t89 + t92 * t72;
t91 = cos(qJ(5));
t88 = sin(qJ(5));
t83 = t86 ^ 2;
t81 = qJD(5) + t82;
t79 = t80 ^ 2;
t78 = -t85 * pkin(2) + t96;
t77 = t83 * t79;
t75 = (-t88 * t89 + t91 * t92) * t100;
t74 = (-t88 * t92 - t89 * t91) * t100;
t73 = (pkin(4) * t85 * t89 + t80) * t86;
t67 = -t89 * t97 + t69;
t66 = t82 * pkin(4) - t92 * t97 + t68;
t65 = t88 * t66 + t91 * t67;
t64 = t91 * t66 - t88 * t67;
t1 = m(4) * (t78 ^ 2 + t84 * t79 + t77) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t77) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t90 ^ 2 + t93 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t68 * mrSges(5,1) - t69 * mrSges(5,2) + Ifges(5,3) * t82 / 0.2e1) * t82 + (t64 * mrSges(6,1) - t65 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (t73 * mrSges(6,2) - t64 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t75 / 0.2e1) * t75 + (-t73 * mrSges(6,1) + t65 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,6) * t81 + Ifges(6,2) * t74 / 0.2e1) * t74 + (-t78 * mrSges(4,1) * t87 + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t84 / 0.2e1) * t85 + (mrSges(3,1) * t93 - mrSges(3,2) * t90) * t98 + (t83 + t84) * t80 * mrSges(4,3) + (t78 * mrSges(4,2) + Ifges(4,4) * t99 + (t80 * (mrSges(5,1) * t89 + mrSges(5,2) * t92) + (Ifges(5,1) * t92 ^ 2 / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t92 + Ifges(5,2) * t89 / 0.2e1) * t89) * t85) * t86 + (-t68 * t92 - t69 * t89) * mrSges(5,3) + t82 * (Ifges(5,5) * t92 - Ifges(5,6) * t89)) * t86) * t85;
T = t1;
