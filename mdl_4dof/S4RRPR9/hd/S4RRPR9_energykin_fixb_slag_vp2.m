% Calculate kinetic energy for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:18
% EndTime: 2019-12-31 17:09:18
% DurationCPUTime: 0.24s
% Computational Cost: add. (192->63), mult. (449->102), div. (0->0), fcn. (266->6), ass. (0->25)
t79 = pkin(5) * mrSges(3,3);
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t61 = (-pkin(2) * t74 - qJ(3) * t72 - pkin(1)) * qJD(1);
t77 = t74 * qJD(1);
t66 = pkin(5) * t77 + qJD(2) * qJ(3);
t69 = sin(pkin(7));
t70 = cos(pkin(7));
t57 = t69 * t61 + t70 * t66;
t78 = t72 * qJD(1);
t56 = t70 * t61 - t66 * t69;
t65 = -qJD(2) * pkin(2) + pkin(5) * t78 + qJD(3);
t73 = cos(qJ(4));
t71 = sin(qJ(4));
t67 = qJD(4) - t77;
t63 = qJD(2) * t69 + t70 * t78;
t62 = qJD(2) * t70 - t69 * t78;
t58 = -pkin(3) * t62 + t65;
t55 = t62 * t71 + t63 * t73;
t54 = t62 * t73 - t63 * t71;
t53 = pkin(6) * t62 + t57;
t52 = -pkin(3) * t77 - pkin(6) * t63 + t56;
t51 = t52 * t71 + t53 * t73;
t50 = t52 * t73 - t53 * t71;
t1 = m(5) * (t50 ^ 2 + t51 ^ 2 + t58 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t56 ^ 2 + t57 ^ 2 + t65 ^ 2) / 0.2e1 + (t50 * mrSges(5,1) - t51 * mrSges(5,2) + Ifges(5,3) * t67 / 0.2e1) * t67 + (t65 * mrSges(4,2) - t56 * mrSges(4,3) + Ifges(4,1) * t63 / 0.2e1) * t63 + (t58 * mrSges(5,2) - t50 * mrSges(5,3) + Ifges(5,5) * t67 + Ifges(5,1) * t55 / 0.2e1) * t55 + (-Ifges(4,6) * t77 - t65 * mrSges(4,1) + t57 * mrSges(4,3) + Ifges(4,4) * t63 + Ifges(4,2) * t62 / 0.2e1) * t62 + (-t58 * mrSges(5,1) + t51 * mrSges(5,3) + Ifges(5,4) * t55 + Ifges(5,6) * t67 + Ifges(5,2) * t54 / 0.2e1) * t54 + ((-pkin(5) * mrSges(3,1) + Ifges(3,5)) * qJD(2) * t72 + (-t56 * mrSges(4,1) + t57 * mrSges(4,2) - Ifges(4,5) * t63 + (-pkin(5) * mrSges(3,2) + Ifges(3,6)) * qJD(2)) * t74 + (m(3) * (pkin(1) ^ 2 + (t72 ^ 2 + t74 ^ 2) * pkin(5) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (t79 + Ifges(3,1) / 0.2e1) * t72) * t72 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t72 + (t79 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t74) * t74) * qJD(1)) * qJD(1);
T = t1;
