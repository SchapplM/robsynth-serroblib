% Calculate kinetic energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:32
% EndTime: 2019-12-05 17:35:32
% DurationCPUTime: 0.27s
% Computational Cost: add. (153->63), mult. (355->93), div. (0->0), fcn. (178->6), ass. (0->25)
t70 = sin(pkin(7));
t66 = (pkin(1) * t70 + qJ(3)) * qJD(1);
t71 = cos(pkin(8));
t68 = t71 * qJD(2);
t69 = sin(pkin(8));
t62 = t69 * t66 - t68;
t82 = t62 ^ 2;
t81 = m(3) / 0.2e1;
t72 = cos(pkin(7));
t78 = -pkin(1) * t72 - pkin(2);
t60 = qJD(3) + (-pkin(3) * t71 - pkin(6) * t69 + t78) * qJD(1);
t64 = t69 * qJD(2) + t71 * t66;
t73 = sin(qJ(4));
t74 = cos(qJ(4));
t56 = t73 * t60 + t74 * t64;
t80 = t69 * qJD(1);
t79 = t71 * qJD(1);
t77 = qJ(5) * t80;
t55 = t74 * t60 - t73 * t64;
t67 = qJD(4) - t79;
t65 = t78 * qJD(1) + qJD(3);
t57 = qJD(5) - t68 + (pkin(4) * qJD(1) * t73 + t66) * t69;
t54 = -t73 * t77 + t56;
t53 = t67 * pkin(4) - t74 * t77 + t55;
t1 = qJD(2) ^ 2 * t81 + m(6) * (t53 ^ 2 + t54 ^ 2 + t57 ^ 2) / 0.2e1 + m(5) * (t55 ^ 2 + t56 ^ 2 + t82) / 0.2e1 + m(4) * (t64 ^ 2 + t65 ^ 2 + t82) / 0.2e1 + (t55 * mrSges(5,1) + t53 * mrSges(6,1) - t56 * mrSges(5,2) - t54 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t67) * t67 + ((-t65 * mrSges(4,1) + t64 * mrSges(4,3) + Ifges(4,2) * t79 / 0.2e1) * t71 + (t65 * mrSges(4,2) + t62 * mrSges(4,3) + (t62 * mrSges(5,2) + t57 * mrSges(6,2) - t55 * mrSges(5,3) - t53 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t74 * t80 + (Ifges(5,5) + Ifges(6,5)) * t67) * t74 + (t62 * mrSges(5,1) + t57 * mrSges(6,1) - t56 * mrSges(5,3) - t54 * mrSges(6,3) + (-Ifges(5,6) - Ifges(6,6)) * t67 + ((Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t73 + (-Ifges(5,4) - Ifges(6,4)) * t74) * t80) * t73) * t69 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t72 * mrSges(3,1) - t70 * mrSges(3,2) + (t70 ^ 2 + t72 ^ 2) * t81 * pkin(1)) * pkin(1) + (Ifges(4,4) * t71 + Ifges(4,1) * t69 / 0.2e1) * t69) * qJD(1)) * qJD(1);
T = t1;
