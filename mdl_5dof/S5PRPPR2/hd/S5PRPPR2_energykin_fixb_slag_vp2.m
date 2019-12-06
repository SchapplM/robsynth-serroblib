% Calculate kinetic energy for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:53
% EndTime: 2019-12-05 15:23:53
% DurationCPUTime: 0.20s
% Computational Cost: add. (138->53), mult. (307->90), div. (0->0), fcn. (184->8), ass. (0->27)
t88 = cos(qJ(2));
t76 = qJD(2) * pkin(2) + t88 * qJD(1);
t82 = sin(pkin(8));
t84 = cos(pkin(8));
t86 = sin(qJ(2));
t92 = qJD(1) * t86;
t72 = t82 * t76 + t84 * t92;
t70 = qJD(2) * qJ(4) + t72;
t81 = sin(pkin(9));
t83 = cos(pkin(9));
t66 = t81 * qJD(3) + t83 * t70;
t93 = pkin(6) * qJD(2);
t71 = t84 * t76 - t82 * t92;
t91 = qJD(4) - t71;
t87 = cos(qJ(5));
t85 = sin(qJ(5));
t80 = t83 * qJD(3);
t74 = (t81 * t87 + t83 * t85) * qJD(2);
t73 = (-t81 * t85 + t83 * t87) * qJD(2);
t69 = -qJD(2) * pkin(3) + t91;
t67 = (-pkin(4) * t83 - pkin(3)) * qJD(2) + t91;
t65 = -t81 * t70 + t80;
t64 = t83 * t93 + t66;
t63 = t80 + (-t70 - t93) * t81;
t62 = t85 * t63 + t87 * t64;
t61 = t87 * t63 - t85 * t64;
t1 = m(5) * (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + (m(3) * (t86 ^ 2 + t88 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t67 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,1) * t74 / 0.2e1) * t74 + (-t67 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t74 + Ifges(6,2) * t73 / 0.2e1) * t73 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,5) * t74 + Ifges(6,6) * t73 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t72 * mrSges(4,2) + t71 * mrSges(4,1) + t69 * (-mrSges(5,1) * t83 + mrSges(5,2) * t81) + (t88 * mrSges(3,1) - t86 * mrSges(3,2)) * qJD(1) + (-t65 * t81 + t66 * t83) * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) * t83 ^ 2 / 0.2e1 + (Ifges(5,4) * t83 + Ifges(5,1) * t81 / 0.2e1) * t81) * qJD(2)) * qJD(2);
T = t1;
