% Calculate kinetic energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRP10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:20
% EndTime: 2019-12-31 20:09:20
% DurationCPUTime: 0.30s
% Computational Cost: add. (208->80), mult. (432->106), div. (0->0), fcn. (204->4), ass. (0->25)
t83 = -pkin(2) - pkin(7);
t82 = pkin(6) * mrSges(3,3);
t75 = cos(qJ(2));
t73 = sin(qJ(2));
t78 = -qJ(3) * t73 - pkin(1);
t59 = (t75 * t83 + t78) * qJD(1);
t80 = t73 * qJD(1);
t79 = pkin(6) * t80 + qJD(3);
t60 = pkin(3) * t80 + qJD(2) * t83 + t79;
t72 = sin(qJ(4));
t74 = cos(qJ(4));
t54 = t74 * t59 + t72 * t60;
t81 = qJD(1) * t75;
t66 = -pkin(6) * t81 - qJD(2) * qJ(3);
t61 = pkin(3) * t81 - t66;
t53 = -t72 * t59 + t74 * t60;
t67 = qJD(4) + t80;
t65 = -qJD(2) * pkin(2) + t79;
t64 = t74 * qJD(2) - t72 * t81;
t63 = -t72 * qJD(2) - t74 * t81;
t62 = (-pkin(2) * t75 + t78) * qJD(1);
t55 = -t63 * pkin(4) + qJD(5) + t61;
t52 = t63 * qJ(5) + t54;
t51 = t67 * pkin(4) - t64 * qJ(5) + t53;
t1 = m(4) * (t62 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t51 ^ 2 + t52 ^ 2 + t55 ^ 2) / 0.2e1 + m(5) * (t53 ^ 2 + t54 ^ 2 + t61 ^ 2) / 0.2e1 + (t65 * mrSges(4,2) - t66 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(2)) * qJD(2) + (t53 * mrSges(5,1) + t51 * mrSges(6,1) - t54 * mrSges(5,2) - t52 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t67) * t67 + (t61 * mrSges(5,2) + t55 * mrSges(6,2) - t53 * mrSges(5,3) - t51 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t64 + (Ifges(5,5) + Ifges(6,5)) * t67) * t64 + (-t61 * mrSges(5,1) - t55 * mrSges(6,1) + t54 * mrSges(5,3) + t52 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t63 + (Ifges(5,6) + Ifges(6,6)) * t67 + (Ifges(5,4) + Ifges(6,4)) * t64) * t63 + ((-t66 * mrSges(4,1) + t62 * mrSges(4,2) + (-pkin(6) * mrSges(3,2) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t75 + (t65 * mrSges(4,1) - t62 * mrSges(4,3) + (-pkin(6) * mrSges(3,1) - Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t73 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t73 ^ 2 + t75 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t82 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t75) * t75 + (-pkin(1) * mrSges(3,2) + (t82 + Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1) * t73 + (Ifges(3,4) + Ifges(4,6)) * t75) * t73) * qJD(1)) * qJD(1);
T = t1;
