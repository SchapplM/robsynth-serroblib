% Calculate kinetic energy for
% S4RRRR2
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:06
% EndTime: 2019-12-31 17:23:06
% DurationCPUTime: 0.14s
% Computational Cost: add. (158->47), mult. (238->85), div. (0->0), fcn. (110->6), ass. (0->23)
t62 = qJD(1) + qJD(2);
t75 = t62 / 0.2e1;
t65 = sin(qJ(2));
t73 = pkin(1) * qJD(1);
t59 = t62 * pkin(6) + t65 * t73;
t74 = t59 * mrSges(4,3);
t68 = cos(qJ(2));
t72 = t68 * t73;
t71 = pkin(7) * t62 + t59;
t67 = cos(qJ(3));
t66 = cos(qJ(4));
t64 = sin(qJ(3));
t63 = sin(qJ(4));
t61 = qJD(3) + qJD(4);
t60 = -t62 * pkin(2) - t72;
t57 = -t72 + (-pkin(3) * t67 - pkin(2)) * t62;
t56 = (t63 * t67 + t64 * t66) * t62;
t55 = (-t63 * t64 + t66 * t67) * t62;
t54 = t71 * t67;
t53 = qJD(3) * pkin(3) - t71 * t64;
t52 = t63 * t53 + t66 * t54;
t51 = t66 * t53 - t63 * t54;
t1 = m(4) * (t60 ^ 2 + (t64 ^ 2 + t67 ^ 2) * t59 ^ 2) / 0.2e1 + m(5) * (t51 ^ 2 + t52 ^ 2 + t57 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t65 ^ 2 + t68 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t51 * mrSges(5,1) - t52 * mrSges(5,2) + Ifges(5,3) * t61 / 0.2e1) * t61 + (Ifges(4,3) * qJD(3) / 0.2e1 + (-t64 * mrSges(4,1) - t67 * mrSges(4,2)) * t59) * qJD(3) + (t57 * mrSges(5,2) - t51 * mrSges(5,3) + Ifges(5,5) * t61 + Ifges(5,1) * t56 / 0.2e1) * t56 + (-t57 * mrSges(5,1) + t52 * mrSges(5,3) + Ifges(5,4) * t56 + Ifges(5,6) * t61 + Ifges(5,2) * t55 / 0.2e1) * t55 + (Ifges(3,3) * t75 + (mrSges(3,1) * t68 - mrSges(3,2) * t65) * t73 + (-t60 * mrSges(4,1) + Ifges(4,6) * qJD(3) + (Ifges(4,2) * t75 + t74) * t67) * t67 + (Ifges(4,4) * t67 * t62 + t60 * mrSges(4,2) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t75 + t74) * t64) * t64) * t62;
T = t1;
