% Calculate kinetic energy for
% S5PRRPR2
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:16
% EndTime: 2019-12-05 16:17:16
% DurationCPUTime: 0.20s
% Computational Cost: add. (142->45), mult. (214->77), div. (0->0), fcn. (96->6), ass. (0->22)
t71 = qJD(2) + qJD(3);
t75 = sin(qJ(3));
t82 = pkin(2) * qJD(2);
t68 = qJ(4) * t71 + t75 * t82;
t72 = sin(pkin(9));
t73 = cos(pkin(9));
t64 = -t73 * qJD(1) + t68 * t72;
t85 = t64 ^ 2;
t83 = t71 * t73;
t77 = cos(qJ(3));
t81 = -t77 * t82 + qJD(4);
t79 = qJD(1) ^ 2;
t78 = qJD(2) ^ 2;
t76 = cos(qJ(5));
t74 = sin(qJ(5));
t69 = qJD(5) - t83;
t67 = -pkin(3) * t71 + t81;
t66 = qJD(1) * t72 + t68 * t73;
t63 = (-pkin(4) * t73 - pkin(7) * t72 - pkin(3)) * t71 + t81;
t62 = t63 * t74 + t66 * t76;
t61 = t63 * t76 - t66 * t74;
t1 = t78 * Ifges(3,3) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t85) / 0.2e1 + m(4) * (t79 + (t75 ^ 2 + t77 ^ 2) * pkin(2) ^ 2 * t78) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t85) / 0.2e1 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t69 / 0.2e1) * t69 + (Ifges(4,3) * t71 / 0.2e1 + (-t67 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,2) * t83 / 0.2e1) * t73 + (mrSges(4,1) * t77 - mrSges(4,2) * t75) * t82 + (t67 * mrSges(5,2) + (Ifges(5,4) * t73 + (Ifges(5,1) / 0.2e1 + Ifges(6,1) * t76 ^ 2 / 0.2e1 + (-Ifges(6,4) * t76 + Ifges(6,2) * t74 / 0.2e1) * t74) * t72) * t71 + (mrSges(6,1) * t74 + mrSges(6,2) * t76 + mrSges(5,3)) * t64 + (-t61 * t76 - t62 * t74) * mrSges(6,3) + t69 * (Ifges(6,5) * t76 - Ifges(6,6) * t74)) * t72) * t71 + (m(3) + m(2)) * t79 / 0.2e1;
T = t1;
