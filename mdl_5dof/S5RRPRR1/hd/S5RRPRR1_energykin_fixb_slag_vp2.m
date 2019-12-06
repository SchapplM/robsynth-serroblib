% Calculate kinetic energy for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:47
% EndTime: 2019-12-05 18:23:47
% DurationCPUTime: 0.28s
% Computational Cost: add. (215->70), mult. (461->109), div. (0->0), fcn. (266->6), ass. (0->30)
t77 = sin(qJ(4));
t78 = sin(qJ(2));
t80 = cos(qJ(4));
t81 = cos(qJ(2));
t68 = (t77 * t78 - t80 * t81) * qJD(1);
t85 = qJD(1) * t81;
t87 = pkin(3) + qJ(3);
t71 = t87 * t85;
t75 = qJD(2) * pkin(1);
t86 = qJD(1) * t78;
t83 = qJD(2) * pkin(2) - t87 * t86 + t75;
t58 = t71 * t77 - t80 * t83;
t88 = t58 ^ 2;
t60 = t80 * t71 + t77 * t83;
t70 = qJD(3) + (-pkin(1) - pkin(2)) * t85;
t82 = qJD(1) ^ 2;
t79 = cos(qJ(5));
t76 = sin(qJ(5));
t74 = qJD(2) + qJD(4);
t73 = -pkin(1) * t85 + qJD(3);
t72 = -qJ(3) * t86 + t75;
t69 = (t77 * t81 + t78 * t80) * qJD(1);
t66 = qJD(5) + t68;
t63 = t69 * t79 + t74 * t76;
t62 = -t69 * t76 + t74 * t79;
t61 = -pkin(4) * t69 + t70;
t57 = pkin(4) * t74 + t60;
t56 = t57 * t79 + t61 * t76;
t55 = -t57 * t76 + t61 * t79;
t1 = t82 * Ifges(2,3) / 0.2e1 + m(4) * (qJ(3) ^ 2 * t81 ^ 2 * t82 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t60 ^ 2 + t70 ^ 2 + t88) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t88) / 0.2e1 + (-t58 * mrSges(5,1) - t60 * mrSges(5,2) + Ifges(5,3) * t74 / 0.2e1) * t74 + (t55 * mrSges(6,1) - t56 * mrSges(6,2) + Ifges(6,3) * t66 / 0.2e1) * t66 + (t70 * mrSges(5,2) + t58 * mrSges(5,3) + Ifges(5,5) * t74 + Ifges(5,1) * t69 / 0.2e1) * t69 + (t58 * mrSges(6,2) - t55 * mrSges(6,3) + Ifges(6,5) * t66 + Ifges(6,1) * t63 / 0.2e1) * t63 - (-t70 * mrSges(5,1) + t60 * mrSges(5,3) + Ifges(5,4) * t69 + Ifges(5,6) * t74 - Ifges(5,2) * t68 / 0.2e1) * t68 + (-t58 * mrSges(6,1) + t56 * mrSges(6,3) + Ifges(6,4) * t63 + Ifges(6,6) * t66 + Ifges(6,2) * t62 / 0.2e1) * t62 + ((t73 * mrSges(4,2) - t72 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t86) * t78 + (-t73 * mrSges(4,1) + ((Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1 + qJ(3) * mrSges(4,3)) * t81 + (Ifges(3,4) + Ifges(4,4)) * t78) * qJD(1)) * t81) * qJD(1) + (t72 * mrSges(4,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2) + ((Ifges(3,5) + Ifges(4,5)) * t78 + (-mrSges(4,2) * qJ(3) + Ifges(3,6) + Ifges(4,6)) * t81) * qJD(1)) * qJD(2);
T = t1;
