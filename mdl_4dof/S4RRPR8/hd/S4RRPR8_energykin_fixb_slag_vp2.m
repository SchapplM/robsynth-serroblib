% Calculate kinetic energy for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:45
% EndTime: 2019-12-31 17:07:45
% DurationCPUTime: 0.24s
% Computational Cost: add. (116->60), mult. (265->90), div. (0->0), fcn. (114->4), ass. (0->22)
t73 = pkin(2) + pkin(3);
t72 = pkin(5) * mrSges(3,3);
t65 = cos(qJ(2));
t71 = qJD(1) * t65;
t57 = pkin(5) * t71 + qJD(2) * qJ(3);
t63 = sin(qJ(2));
t70 = t63 * qJD(1);
t69 = pkin(5) * t70 + qJD(3);
t68 = qJ(3) * t63 + pkin(1);
t64 = cos(qJ(4));
t62 = sin(qJ(4));
t60 = -qJD(2) + qJD(4);
t56 = -qJD(2) * pkin(2) + t69;
t55 = (-pkin(2) * t65 - t68) * qJD(1);
t54 = -pkin(6) * t71 + t57;
t53 = (-t62 * t65 + t63 * t64) * qJD(1);
t52 = (-t62 * t63 - t64 * t65) * qJD(1);
t51 = -pkin(6) * t70 - qJD(2) * t73 + t69;
t50 = (t65 * t73 + t68) * qJD(1);
t49 = t51 * t62 + t54 * t64;
t48 = t51 * t64 - t54 * t62;
t1 = m(4) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(5) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + (t48 * mrSges(5,1) - t49 * mrSges(5,2) + Ifges(5,3) * t60 / 0.2e1) * t60 + (t50 * mrSges(5,2) - t48 * mrSges(5,3) + Ifges(5,5) * t60 + Ifges(5,1) * t53 / 0.2e1) * t53 + (-t50 * mrSges(5,1) + t49 * mrSges(5,3) + Ifges(5,4) * t53 + Ifges(5,6) * t60 + Ifges(5,2) * t52 / 0.2e1) * t52 + (-t56 * mrSges(4,1) + t57 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2)) * qJD(2) + ((-t55 * mrSges(4,1) + t57 * mrSges(4,2) + (-mrSges(3,2) * pkin(5) + Ifges(3,6) - Ifges(4,6)) * qJD(2)) * t65 + (t56 * mrSges(4,2) - t55 * mrSges(4,3) + (-mrSges(3,1) * pkin(5) + Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t63 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t63 ^ 2 + t65 ^ 2) * pkin(5) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t72 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t65) * t65 + (-pkin(1) * mrSges(3,2) + (t72 + Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t63 + (Ifges(3,4) - Ifges(4,5)) * t65) * t63) * qJD(1)) * qJD(1);
T = t1;
