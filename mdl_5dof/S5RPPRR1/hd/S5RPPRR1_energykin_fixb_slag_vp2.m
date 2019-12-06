% Calculate kinetic energy for
% S5RPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:37:57
% EndTime: 2019-12-05 17:37:57
% DurationCPUTime: 0.18s
% Computational Cost: add. (122->54), mult. (232->86), div. (0->0), fcn. (90->4), ass. (0->23)
t73 = -pkin(1) - qJ(3);
t61 = -t73 * qJD(1) - qJD(2);
t77 = t61 ^ 2;
t76 = m(3) / 0.2e1;
t63 = qJD(1) * qJ(2) + qJD(3);
t60 = -qJD(1) * pkin(6) + t63;
t75 = t60 * mrSges(5,3);
t74 = qJD(1) / 0.2e1;
t72 = -pkin(7) * qJD(1) + t60;
t70 = cos(qJ(4));
t69 = cos(qJ(5));
t68 = sin(qJ(4));
t67 = sin(qJ(5));
t65 = qJD(4) + qJD(5);
t64 = -qJD(1) * pkin(1) + qJD(2);
t58 = -qJD(2) + (pkin(4) * t68 - t73) * qJD(1);
t57 = (-t67 * t68 + t69 * t70) * qJD(1);
t56 = (-t67 * t70 - t68 * t69) * qJD(1);
t55 = t72 * t68;
t54 = qJD(4) * pkin(4) + t72 * t70;
t53 = t54 * t67 + t55 * t69;
t52 = t54 * t69 - t55 * t67;
t1 = t64 ^ 2 * t76 + m(6) * (t52 ^ 2 + t53 ^ 2 + t58 ^ 2) / 0.2e1 + m(4) * (t63 ^ 2 + t77) / 0.2e1 + m(5) * (t77 + (t68 ^ 2 + t70 ^ 2) * t60 ^ 2) / 0.2e1 + (t52 * mrSges(6,1) - t53 * mrSges(6,2) + Ifges(6,3) * t65 / 0.2e1) * t65 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t70 * mrSges(5,1) - t68 * mrSges(5,2)) * t60) * qJD(4) + (t58 * mrSges(6,2) - t52 * mrSges(6,3) + Ifges(6,5) * t65 + Ifges(6,1) * t57 / 0.2e1) * t57 + (-t58 * mrSges(6,1) + t53 * mrSges(6,3) + Ifges(6,4) * t57 + Ifges(6,6) * t65 + Ifges(6,2) * t56 / 0.2e1) * t56 + (t64 * mrSges(3,2) + t63 * mrSges(4,2) + t61 * mrSges(4,3) + (t61 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t74 - t75) * t70) * t70 + (t61 * mrSges(5,1) - Ifges(5,6) * qJD(4) + (Ifges(5,2) * t74 - t75) * t68) * t68 + (Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t76 + mrSges(3,3)) * qJ(2) - Ifges(5,4) * t70 * t68) * qJD(1)) * qJD(1);
T = t1;
