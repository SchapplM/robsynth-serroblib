% Calculate kinetic energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:07:36
% EndTime: 2022-01-20 12:07:37
% DurationCPUTime: 0.27s
% Computational Cost: add. (412->71), mult. (579->120), div. (0->0), fcn. (346->8), ass. (0->33)
t86 = qJD(1) + qJD(2);
t101 = t86 / 0.2e1;
t90 = sin(qJ(2));
t99 = pkin(1) * qJD(1);
t82 = t86 * pkin(7) + t90 * t99;
t100 = t82 * mrSges(4,3);
t89 = sin(qJ(3));
t97 = pkin(8) * t86 + t82;
t76 = qJD(3) * pkin(3) - t97 * t89;
t93 = cos(qJ(3));
t77 = t97 * t93;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t69 = t88 * t76 + t92 * t77;
t85 = qJD(3) + qJD(4);
t94 = cos(qJ(2));
t98 = t94 * t99;
t68 = t92 * t76 - t88 * t77;
t80 = -t98 + (-pkin(3) * t93 - pkin(2)) * t86;
t91 = cos(qJ(5));
t87 = sin(qJ(5));
t84 = qJD(5) + t85;
t83 = -t86 * pkin(2) - t98;
t79 = (t88 * t93 + t89 * t92) * t86;
t78 = (-t88 * t89 + t92 * t93) * t86;
t72 = -t78 * pkin(4) + t80;
t71 = t87 * t78 + t91 * t79;
t70 = t91 * t78 - t87 * t79;
t67 = t78 * pkin(9) + t69;
t66 = t85 * pkin(4) - t79 * pkin(9) + t68;
t65 = t87 * t66 + t91 * t67;
t64 = t91 * t66 - t87 * t67;
t1 = m(4) * (t83 ^ 2 + (t89 ^ 2 + t93 ^ 2) * t82 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t72 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t90 ^ 2 + t94 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t68 * mrSges(5,1) - t69 * mrSges(5,2) + Ifges(5,3) * t85 / 0.2e1) * t85 + (t64 * mrSges(6,1) - t65 * mrSges(6,2) + Ifges(6,3) * t84 / 0.2e1) * t84 + (Ifges(4,3) * qJD(3) / 0.2e1 + (-t89 * mrSges(4,1) - t93 * mrSges(4,2)) * t82) * qJD(3) + (t80 * mrSges(5,2) - t68 * mrSges(5,3) + Ifges(5,5) * t85 + Ifges(5,1) * t79 / 0.2e1) * t79 + (t72 * mrSges(6,2) - t64 * mrSges(6,3) + Ifges(6,5) * t84 + Ifges(6,1) * t71 / 0.2e1) * t71 + (-t80 * mrSges(5,1) + t69 * mrSges(5,3) + Ifges(5,4) * t79 + Ifges(5,6) * t85 + Ifges(5,2) * t78 / 0.2e1) * t78 + (-t72 * mrSges(6,1) + t65 * mrSges(6,3) + Ifges(6,4) * t71 + Ifges(6,6) * t84 + Ifges(6,2) * t70 / 0.2e1) * t70 + (Ifges(3,3) * t101 + (mrSges(3,1) * t94 - mrSges(3,2) * t90) * t99 + (-t83 * mrSges(4,1) + Ifges(4,6) * qJD(3) + (Ifges(4,2) * t101 + t100) * t93) * t93 + (Ifges(4,4) * t93 * t86 + t83 * mrSges(4,2) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t101 + t100) * t89) * t89) * t86;
T = t1;
