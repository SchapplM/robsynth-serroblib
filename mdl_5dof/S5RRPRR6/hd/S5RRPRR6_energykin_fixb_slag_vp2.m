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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:35:38
% EndTime: 2019-12-05 18:35:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (321->64), mult. (465->115), div. (0->0), fcn. (254->8), ass. (0->34)
t85 = cos(pkin(9));
t82 = t85 ^ 2;
t83 = qJD(1) + qJD(2);
t88 = sin(qJ(2));
t96 = qJD(1) * pkin(1);
t78 = qJ(3) * t83 + t88 * t96;
t99 = t78 * t85;
t84 = sin(pkin(9));
t98 = t83 * t84;
t97 = t83 * t85;
t91 = cos(qJ(2));
t94 = -t91 * t96 + qJD(3);
t70 = (-pkin(3) * t85 - pkin(7) * t84 - pkin(2)) * t83 + t94;
t87 = sin(qJ(4));
t90 = cos(qJ(4));
t67 = t87 * t70 + t90 * t99;
t95 = pkin(8) * t98;
t80 = qJD(4) - t97;
t66 = t90 * t70 - t87 * t99;
t89 = cos(qJ(5));
t86 = sin(qJ(5));
t81 = t84 ^ 2;
t79 = qJD(5) + t80;
t77 = t78 ^ 2;
t76 = -pkin(2) * t83 + t94;
t75 = t81 * t77;
t73 = (-t86 * t87 + t89 * t90) * t98;
t72 = (-t86 * t90 - t87 * t89) * t98;
t71 = (pkin(4) * t83 * t87 + t78) * t84;
t65 = -t87 * t95 + t67;
t64 = pkin(4) * t80 - t90 * t95 + t66;
t63 = t64 * t86 + t65 * t89;
t62 = t64 * t89 - t65 * t86;
t1 = m(4) * (t76 ^ 2 + t77 * t82 + t75) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t75) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t88 ^ 2 + t91 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,3) * t80 / 0.2e1) * t80 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t79 / 0.2e1) * t79 + (t71 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t79 + Ifges(6,1) * t73 / 0.2e1) * t73 + (-t71 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t73 + Ifges(6,6) * t79 + Ifges(6,2) * t72 / 0.2e1) * t72 + (-t76 * mrSges(4,1) * t85 + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t82 / 0.2e1) * t83 + (mrSges(3,1) * t91 - mrSges(3,2) * t88) * t96 + (t81 + t82) * t78 * mrSges(4,3) + (t76 * mrSges(4,2) + Ifges(4,4) * t97 + (t78 * (mrSges(5,1) * t87 + mrSges(5,2) * t90) + (Ifges(5,1) * t90 ^ 2 / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t90 + Ifges(5,2) * t87 / 0.2e1) * t87) * t83) * t84 + (-t66 * t90 - t67 * t87) * mrSges(5,3) + t80 * (Ifges(5,5) * t90 - Ifges(5,6) * t87)) * t84) * t83;
T = t1;
