% Calculate kinetic energy for
% S5PRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:32
% EndTime: 2019-12-05 16:15:32
% DurationCPUTime: 0.18s
% Computational Cost: add. (158->48), mult. (240->84), div. (0->0), fcn. (120->6), ass. (0->25)
t76 = qJD(2) + qJD(3);
t88 = pkin(7) * t76;
t80 = sin(qJ(3));
t87 = pkin(2) * qJD(2);
t73 = t76 * qJ(4) + t80 * t87;
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t67 = t77 * qJD(1) + t78 * t73;
t82 = cos(qJ(3));
t86 = -t82 * t87 + qJD(4);
t84 = qJD(1) ^ 2;
t83 = qJD(2) ^ 2;
t81 = cos(qJ(5));
t79 = sin(qJ(5));
t75 = t78 * qJD(1);
t72 = -t76 * pkin(3) + t86;
t70 = (t77 * t81 + t78 * t79) * t76;
t69 = (-t77 * t79 + t78 * t81) * t76;
t68 = (-pkin(4) * t78 - pkin(3)) * t76 + t86;
t66 = -t77 * t73 + t75;
t65 = t78 * t88 + t67;
t64 = t75 + (-t73 - t88) * t77;
t63 = t79 * t64 + t81 * t65;
t62 = t81 * t64 - t79 * t65;
t1 = t83 * Ifges(3,3) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t72 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t84 + (t80 ^ 2 + t82 ^ 2) * pkin(2) ^ 2 * t83) / 0.2e1 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,5) * t70 + Ifges(6,6) * t69 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t72 * (-mrSges(5,1) * t78 + mrSges(5,2) * t77) + (Ifges(5,2) * t78 ^ 2 / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t78 + Ifges(5,1) * t77 / 0.2e1) * t77) * t76 + (mrSges(4,1) * t82 - mrSges(4,2) * t80) * t87 + (-t66 * t77 + t67 * t78) * mrSges(5,3)) * t76 + (m(3) + m(2)) * t84 / 0.2e1;
T = t1;
