% Calculate kinetic energy for
% S5RRPRR5
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:08
% EndTime: 2022-01-20 11:02:09
% DurationCPUTime: 0.32s
% Computational Cost: add. (381->64), mult. (556->111), div. (0->0), fcn. (342->8), ass. (0->32)
t88 = cos(pkin(9));
t101 = t88 ^ 2;
t87 = sin(pkin(9));
t86 = qJD(1) + qJD(2);
t91 = sin(qJ(2));
t99 = pkin(1) * qJD(1);
t82 = qJ(3) * t86 + t91 * t99;
t98 = pkin(7) * t86 + t82;
t75 = t98 * t87;
t76 = t98 * t88;
t90 = sin(qJ(4));
t93 = cos(qJ(4));
t68 = -t90 * t75 + t93 * t76;
t100 = t87 ^ 2 + t101;
t67 = -t93 * t75 - t76 * t90;
t94 = cos(qJ(2));
t97 = -t94 * t99 + qJD(3);
t77 = (-pkin(3) * t88 - pkin(2)) * t86 + t97;
t92 = cos(qJ(5));
t89 = sin(qJ(5));
t85 = qJD(4) + qJD(5);
t80 = -pkin(2) * t86 + t97;
t79 = (t87 * t93 + t88 * t90) * t86;
t78 = (-t87 * t90 + t88 * t93) * t86;
t71 = -pkin(4) * t78 + t77;
t70 = t78 * t89 + t79 * t92;
t69 = t78 * t92 - t79 * t89;
t66 = pkin(8) * t78 + t68;
t65 = qJD(4) * pkin(4) - pkin(8) * t79 + t67;
t64 = t65 * t89 + t66 * t92;
t63 = t65 * t92 - t66 * t89;
t1 = m(4) * (t100 * t82 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t77 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t91 ^ 2 + t94 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t85 / 0.2e1) * t85 + (t77 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,1) * t79 / 0.2e1) * t79 + (-t77 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t79 + Ifges(5,2) * t78 / 0.2e1) * t78 + (t71 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t85 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t71 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t85 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,5) * t79 + Ifges(5,6) * t78 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t80 * (-mrSges(4,1) * t88 + mrSges(4,2) * t87) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t101 / 0.2e1 + (Ifges(4,4) * t88 + Ifges(4,1) * t87 / 0.2e1) * t87) * t86 + (mrSges(3,1) * t94 - mrSges(3,2) * t91) * t99 + t100 * t82 * mrSges(4,3)) * t86;
T = t1;
