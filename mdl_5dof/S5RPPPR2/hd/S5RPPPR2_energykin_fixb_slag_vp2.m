% Calculate kinetic energy for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:30:54
% EndTime: 2019-12-05 17:30:54
% DurationCPUTime: 0.42s
% Computational Cost: add. (285->79), mult. (764->131), div. (0->0), fcn. (500->8), ass. (0->36)
t92 = sin(pkin(7));
t94 = cos(pkin(8));
t107 = t92 * t94;
t90 = sin(pkin(9));
t93 = cos(pkin(9));
t95 = cos(pkin(7));
t79 = (t90 * t107 + t93 * t95) * qJD(1);
t108 = m(3) / 0.2e1;
t104 = qJD(1) * t95;
t103 = qJD(1) * qJ(2);
t101 = t95 * t103;
t81 = qJD(2) + (-pkin(2) * t95 - qJ(3) * t92 - pkin(1)) * qJD(1);
t91 = sin(pkin(8));
t74 = t94 * t101 + t91 * t81;
t70 = -qJ(4) * t104 + t74;
t105 = qJD(1) * t92;
t84 = t92 * t103 + qJD(3);
t76 = (pkin(3) * t91 - qJ(4) * t94) * t105 + t84;
t67 = t93 * t70 + t90 * t76;
t102 = t91 * t105;
t73 = -t91 * t101 + t81 * t94;
t66 = -t70 * t90 + t76 * t93;
t69 = pkin(3) * t104 + qJD(4) - t73;
t97 = cos(qJ(5));
t96 = sin(qJ(5));
t87 = -qJD(1) * pkin(1) + qJD(2);
t80 = (t93 * t107 - t90 * t95) * qJD(1);
t77 = qJD(5) + t79;
t72 = t96 * t102 + t80 * t97;
t71 = t97 * t102 - t80 * t96;
t65 = pkin(6) * t102 + t67;
t64 = -pkin(4) * t102 - t66;
t63 = pkin(4) * t79 - pkin(6) * t80 + t69;
t62 = t63 * t96 + t65 * t97;
t61 = t63 * t97 - t65 * t96;
t1 = t87 ^ 2 * t108 + m(4) * (t73 ^ 2 + t74 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t64 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t69 ^ 2) / 0.2e1 + (t69 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,1) * t80 / 0.2e1) * t80 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t77 / 0.2e1) * t77 - (-t69 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t80 - Ifges(5,2) * t79 / 0.2e1) * t79 + (t64 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t77 + Ifges(6,1) * t72 / 0.2e1) * t72 + (-t64 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t72 + Ifges(6,6) * t77 + Ifges(6,2) * t71 / 0.2e1) * t71 + ((-t87 * mrSges(3,1) - t73 * mrSges(4,1) + t74 * mrSges(4,2) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t104) * t95 + (t87 * mrSges(3,2) + (t84 * mrSges(4,2) - t73 * mrSges(4,3)) * t94 + (t84 * mrSges(4,1) + t66 * mrSges(5,1) - t67 * mrSges(5,2) - t74 * mrSges(4,3) + Ifges(5,5) * t80 - Ifges(5,6) * t79) * t91) * t92 + (Ifges(2,3) / 0.2e1 + (qJ(2) * t108 + mrSges(3,3)) * (t92 ^ 2 + t95 ^ 2) * qJ(2) + ((Ifges(3,1) / 0.2e1 + Ifges(4,1) * t94 ^ 2 / 0.2e1) * t92 + (-Ifges(4,5) * t94 + Ifges(3,4)) * t95 + (Ifges(4,6) * t95 + (-Ifges(4,4) * t94 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t91) * t92) * t91) * t92) * qJD(1)) * qJD(1);
T = t1;
