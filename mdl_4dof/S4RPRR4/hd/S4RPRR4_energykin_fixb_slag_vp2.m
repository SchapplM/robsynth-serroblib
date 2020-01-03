% Calculate kinetic energy for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:12
% EndTime: 2019-12-31 16:50:12
% DurationCPUTime: 0.19s
% Computational Cost: add. (104->49), mult. (237->86), div. (0->0), fcn. (110->6), ass. (0->22)
t72 = m(3) / 0.2e1;
t62 = sin(pkin(7));
t58 = (pkin(1) * t62 + pkin(5)) * qJD(1);
t65 = sin(qJ(3));
t67 = cos(qJ(3));
t54 = t65 * qJD(2) + t67 * t58;
t71 = qJD(1) * t65;
t63 = cos(pkin(7));
t70 = -pkin(1) * t63 - pkin(2);
t53 = qJD(2) * t67 - t58 * t65;
t66 = cos(qJ(4));
t64 = sin(qJ(4));
t60 = -qJD(1) * t67 + qJD(4);
t59 = t70 * qJD(1);
t57 = qJD(3) * t64 + t66 * t71;
t56 = qJD(3) * t66 - t64 * t71;
t52 = (-pkin(3) * t67 - pkin(6) * t65 + t70) * qJD(1);
t51 = qJD(3) * pkin(6) + t54;
t50 = -qJD(3) * pkin(3) - t53;
t49 = t51 * t66 + t52 * t64;
t48 = -t51 * t64 + t52 * t66;
t1 = m(5) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + qJD(2) ^ 2 * t72 + m(4) * (t53 ^ 2 + t54 ^ 2 + t59 ^ 2) / 0.2e1 + (t48 * mrSges(5,1) - t49 * mrSges(5,2) + Ifges(5,3) * t60 / 0.2e1) * t60 + (t50 * mrSges(5,2) - t48 * mrSges(5,3) + Ifges(5,5) * t60 + Ifges(5,1) * t57 / 0.2e1) * t57 + (-t50 * mrSges(5,1) + t49 * mrSges(5,3) + Ifges(5,4) * t57 + Ifges(5,6) * t60 + Ifges(5,2) * t56 / 0.2e1) * t56 + (t53 * mrSges(4,1) - t54 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t59 * (-mrSges(4,1) * t67 + mrSges(4,2) * t65) + (-t53 * t65 + t54 * t67) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t65 + Ifges(4,6) * t67) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (mrSges(3,1) * t63 - mrSges(3,2) * t62 + (t62 ^ 2 + t63 ^ 2) * t72 * pkin(1)) * pkin(1) + Ifges(4,2) * t67 ^ 2 / 0.2e1 + (Ifges(4,4) * t67 + Ifges(4,1) * t65 / 0.2e1) * t65) * qJD(1)) * qJD(1);
T = t1;
