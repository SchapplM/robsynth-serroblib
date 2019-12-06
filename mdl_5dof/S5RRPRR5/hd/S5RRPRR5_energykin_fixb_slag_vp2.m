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
% m_mdh [6x1]
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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:33:34
% EndTime: 2019-12-05 18:33:35
% DurationCPUTime: 0.27s
% Computational Cost: add. (381->64), mult. (556->111), div. (0->0), fcn. (342->8), ass. (0->32)
t86 = cos(pkin(9));
t99 = t86 ^ 2;
t85 = sin(pkin(9));
t84 = qJD(1) + qJD(2);
t89 = sin(qJ(2));
t97 = qJD(1) * pkin(1);
t80 = t84 * qJ(3) + t89 * t97;
t96 = pkin(7) * t84 + t80;
t73 = t96 * t85;
t74 = t96 * t86;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t66 = -t88 * t73 + t91 * t74;
t98 = t85 ^ 2 + t99;
t65 = -t91 * t73 - t88 * t74;
t92 = cos(qJ(2));
t95 = -t92 * t97 + qJD(3);
t75 = (-pkin(3) * t86 - pkin(2)) * t84 + t95;
t90 = cos(qJ(5));
t87 = sin(qJ(5));
t83 = qJD(4) + qJD(5);
t78 = -t84 * pkin(2) + t95;
t77 = (t85 * t91 + t86 * t88) * t84;
t76 = (-t85 * t88 + t86 * t91) * t84;
t69 = -t76 * pkin(4) + t75;
t68 = t87 * t76 + t90 * t77;
t67 = t90 * t76 - t87 * t77;
t64 = t76 * pkin(8) + t66;
t63 = qJD(4) * pkin(4) - t77 * pkin(8) + t65;
t62 = t87 * t63 + t90 * t64;
t61 = t90 * t63 - t87 * t64;
t1 = m(4) * (t98 * t80 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t69 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t75 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t89 ^ 2 + t92 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t83 / 0.2e1) * t83 + (t75 * mrSges(5,2) - t65 * mrSges(5,3) + Ifges(5,1) * t77 / 0.2e1) * t77 + (-t75 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,2) * t76 / 0.2e1) * t76 + (t69 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t83 + Ifges(6,1) * t68 / 0.2e1) * t68 + (-t69 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t68 + Ifges(6,6) * t83 + Ifges(6,2) * t67 / 0.2e1) * t67 + (t65 * mrSges(5,1) - t66 * mrSges(5,2) + Ifges(5,5) * t77 + Ifges(5,6) * t76 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t78 * (-mrSges(4,1) * t86 + mrSges(4,2) * t85) + (Ifges(4,2) * t99 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t86 + Ifges(4,1) * t85 / 0.2e1) * t85) * t84 + (mrSges(3,1) * t92 - mrSges(3,2) * t89) * t97 + t98 * t80 * mrSges(4,3)) * t84;
T = t1;
