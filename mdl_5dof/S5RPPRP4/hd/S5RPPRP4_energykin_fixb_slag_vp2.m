% Calculate kinetic energy for
% S5RPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:52:00
% EndTime: 2019-12-31 17:52:00
% DurationCPUTime: 0.13s
% Computational Cost: add. (127->56), mult. (235->75), div. (0->0), fcn. (78->4), ass. (0->21)
t76 = m(3) / 0.2e1;
t61 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t67 = sin(pkin(7));
t68 = cos(pkin(7));
t74 = qJ(2) * qJD(1);
t59 = t67 * t61 + t68 * t74;
t57 = -qJD(1) * pkin(6) + t59;
t70 = sin(qJ(4));
t71 = cos(qJ(4));
t53 = t70 * qJD(3) + t71 * t57;
t75 = qJD(1) * t71;
t73 = qJ(5) * qJD(1);
t58 = t68 * t61 - t67 * t74;
t56 = qJD(1) * pkin(3) - t58;
t66 = t71 * qJD(3);
t64 = -qJD(1) * pkin(1) + qJD(2);
t54 = pkin(4) * t75 + qJD(5) + t56;
t52 = -t70 * t57 + t66;
t51 = -t71 * t73 + t53;
t50 = qJD(4) * pkin(4) + t66 + (-t57 + t73) * t70;
t1 = m(6) * (t50 ^ 2 + t51 ^ 2 + t54 ^ 2) / 0.2e1 + t64 ^ 2 * t76 + m(4) * (qJD(3) ^ 2 + t58 ^ 2 + t59 ^ 2) / 0.2e1 + m(5) * (t52 ^ 2 + t53 ^ 2 + t56 ^ 2) / 0.2e1 + (t52 * mrSges(5,1) + t50 * mrSges(6,1) - t53 * mrSges(5,2) - t51 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(4)) * qJD(4) + (-t64 * mrSges(3,1) - t58 * mrSges(4,1) + t59 * mrSges(4,2) + (t56 * mrSges(5,1) + t54 * mrSges(6,1) - t53 * mrSges(5,3) - t51 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t75 + (-Ifges(5,6) - Ifges(6,6)) * qJD(4)) * t71 + (-t56 * mrSges(5,2) - t54 * mrSges(6,2) + t52 * mrSges(5,3) + t50 * mrSges(6,3) + (-Ifges(5,5) - Ifges(6,5)) * qJD(4)) * t70 + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (qJ(2) * t76 + mrSges(3,3)) * qJ(2) + ((Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t70 + (Ifges(5,4) + Ifges(6,4)) * t71) * t70) * qJD(1)) * qJD(1);
T = t1;
