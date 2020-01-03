% Calculate kinetic energy for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR3_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR3_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:15
% EndTime: 2019-12-31 17:24:15
% DurationCPUTime: 0.26s
% Computational Cost: add. (204->62), mult. (491->102), div. (0->0), fcn. (308->6), ass. (0->25)
t80 = qJD(1) * (-pkin(6) - pkin(5));
t78 = pkin(5) * mrSges(3,3);
t72 = sin(qJ(2));
t65 = qJD(2) * pkin(2) + t72 * t80;
t75 = cos(qJ(2));
t66 = t75 * t80;
t71 = sin(qJ(3));
t74 = cos(qJ(3));
t58 = t71 * t65 - t74 * t66;
t69 = qJD(2) + qJD(3);
t57 = t74 * t65 + t71 * t66;
t67 = (-pkin(2) * t75 - pkin(1)) * qJD(1);
t73 = cos(qJ(4));
t70 = sin(qJ(4));
t68 = qJD(4) + t69;
t63 = (t71 * t75 + t72 * t74) * qJD(1);
t62 = (-t71 * t72 + t74 * t75) * qJD(1);
t59 = -t62 * pkin(3) + t67;
t56 = t70 * t62 + t73 * t63;
t55 = t73 * t62 - t70 * t63;
t54 = t62 * pkin(7) + t58;
t53 = t69 * pkin(3) - t63 * pkin(7) + t57;
t52 = t70 * t53 + t73 * t54;
t51 = t73 * t53 - t70 * t54;
t1 = m(4) * (t57 ^ 2 + t58 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t51 ^ 2 + t52 ^ 2 + t59 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + (t57 * mrSges(4,1) - t58 * mrSges(4,2) + Ifges(4,3) * t69 / 0.2e1) * t69 + (t51 * mrSges(5,1) - t52 * mrSges(5,2) + Ifges(5,3) * t68 / 0.2e1) * t68 + (t67 * mrSges(4,2) - t57 * mrSges(4,3) + Ifges(4,5) * t69 + Ifges(4,1) * t63 / 0.2e1) * t63 + (t59 * mrSges(5,2) - t51 * mrSges(5,3) + Ifges(5,5) * t68 + Ifges(5,1) * t56 / 0.2e1) * t56 + (-t67 * mrSges(4,1) + t58 * mrSges(4,3) + Ifges(4,4) * t63 + Ifges(4,6) * t69 + Ifges(4,2) * t62 / 0.2e1) * t62 + (-t59 * mrSges(5,1) + t52 * mrSges(5,3) + Ifges(5,4) * t56 + Ifges(5,6) * t68 + Ifges(5,2) * t55 / 0.2e1) * t55 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t72 ^ 2 + t75 ^ 2) * pkin(5) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t78 + Ifges(3,2) / 0.2e1) * t75) * t75 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t75 + (t78 + Ifges(3,1) / 0.2e1) * t72) * t72) * qJD(1) + ((-pkin(5) * mrSges(3,2) + Ifges(3,6)) * t75 + (-pkin(5) * mrSges(3,1) + Ifges(3,5)) * t72) * qJD(2)) * qJD(1);
T = t1;
