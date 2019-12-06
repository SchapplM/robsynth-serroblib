% Calculate kinetic energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:34
% EndTime: 2019-12-05 15:12:34
% DurationCPUTime: 0.11s
% Computational Cost: add. (108->40), mult. (217->68), div. (0->0), fcn. (138->8), ass. (0->22)
t73 = qJD(3) + qJD(4);
t84 = t73 / 0.2e1;
t74 = sin(pkin(9));
t75 = cos(pkin(9));
t78 = sin(qJ(3));
t81 = cos(qJ(3));
t70 = (-t74 * t78 + t75 * t81) * qJD(1);
t69 = qJD(3) * pkin(3) + t70;
t71 = (t74 * t81 + t75 * t78) * qJD(1);
t77 = sin(qJ(4));
t80 = cos(qJ(4));
t66 = t77 * t69 + t80 * t71;
t65 = t69 * t80 - t71 * t77;
t83 = qJD(1) ^ 2;
t82 = qJD(2) ^ 2;
t79 = cos(qJ(5));
t76 = sin(qJ(5));
t64 = pkin(7) * t73 + t66;
t63 = -pkin(4) * t73 - t65;
t62 = qJD(2) * t76 + t64 * t79;
t61 = qJD(2) * t79 - t64 * t76;
t1 = m(4) * (t70 ^ 2 + t71 ^ 2 + t82) / 0.2e1 + m(2) * t83 / 0.2e1 + m(3) * (t82 + (t74 ^ 2 + t75 ^ 2) * t83) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t82) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t70 * mrSges(4,1) - t71 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t66 * mrSges(5,2) + t65 * mrSges(5,1) + Ifges(5,3) * t84 + (Ifges(6,2) * t79 * t84 - t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t79 + (t63 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * qJD(5) + (Ifges(6,4) * t79 + Ifges(6,1) * t76 / 0.2e1) * t73) * t76) * t73;
T = t1;
