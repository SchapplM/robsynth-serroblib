% Calculate kinetic energy for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPPRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:45
% EndTime: 2019-12-05 14:57:45
% DurationCPUTime: 0.10s
% Computational Cost: add. (78->36), mult. (176->63), div. (0->0), fcn. (114->8), ass. (0->22)
t83 = qJD(4) / 0.2e1;
t73 = sin(pkin(9));
t75 = cos(pkin(9));
t74 = sin(pkin(8));
t82 = qJD(1) * t74;
t69 = t75 * qJD(2) - t73 * t82;
t70 = t73 * qJD(2) + t75 * t82;
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t66 = t78 * t69 + t80 * t70;
t65 = t80 * t69 - t78 * t70;
t81 = qJD(1) ^ 2;
t79 = cos(qJ(5));
t77 = sin(qJ(5));
t76 = cos(pkin(8));
t72 = -t76 * qJD(1) + qJD(3);
t71 = t72 ^ 2;
t64 = qJD(4) * pkin(6) + t66;
t63 = -qJD(4) * pkin(4) - t65;
t62 = t79 * t64 + t77 * t72;
t61 = -t77 * t64 + t79 * t72;
t1 = m(4) * (t69 ^ 2 + t70 ^ 2 + t71) / 0.2e1 + m(2) * t81 / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t71) / 0.2e1 + m(3) * (qJD(2) ^ 2 + (t74 ^ 2 + t76 ^ 2) * t81) / 0.2e1 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t66 * mrSges(5,2) + t65 * mrSges(5,1) + Ifges(5,3) * t83 + (Ifges(6,2) * t79 * t83 - t63 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,6) * qJD(5)) * t79 + (t63 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * qJD(5) + (Ifges(6,4) * t79 + Ifges(6,1) * t77 / 0.2e1) * qJD(4)) * t77) * qJD(4);
T = t1;
