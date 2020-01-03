% Calculate kinetic energy for
% S4RPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-31 16:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR3_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR3_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR3_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:37:49
% EndTime: 2019-12-31 16:37:49
% DurationCPUTime: 0.16s
% Computational Cost: add. (93->45), mult. (234->78), div. (0->0), fcn. (120->6), ass. (0->22)
t73 = m(3) / 0.2e1;
t64 = sin(pkin(6));
t60 = (pkin(1) * t64 + qJ(3)) * qJD(1);
t63 = sin(pkin(7));
t65 = cos(pkin(7));
t54 = t63 * qJD(2) + t65 * t60;
t72 = pkin(5) * qJD(1);
t66 = cos(pkin(6));
t71 = -pkin(1) * t66 - pkin(2);
t68 = cos(qJ(4));
t67 = sin(qJ(4));
t62 = t65 * qJD(2);
t59 = qJD(1) * t71 + qJD(3);
t57 = (t63 * t68 + t65 * t67) * qJD(1);
t56 = (-t63 * t67 + t65 * t68) * qJD(1);
t55 = qJD(3) + (-pkin(3) * t65 + t71) * qJD(1);
t53 = -t63 * t60 + t62;
t52 = t65 * t72 + t54;
t51 = t62 + (-t60 - t72) * t63;
t50 = t67 * t51 + t68 * t52;
t49 = t68 * t51 - t67 * t52;
t1 = qJD(2) ^ 2 * t73 + m(4) * (t53 ^ 2 + t54 ^ 2 + t59 ^ 2) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t55 ^ 2) / 0.2e1 + (t55 * mrSges(5,2) - t49 * mrSges(5,3) + Ifges(5,1) * t57 / 0.2e1) * t57 + (-t55 * mrSges(5,1) + t50 * mrSges(5,3) + Ifges(5,4) * t57 + Ifges(5,2) * t56 / 0.2e1) * t56 + (t49 * mrSges(5,1) - t50 * mrSges(5,2) + Ifges(5,5) * t57 + Ifges(5,6) * t56 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t59 * (-mrSges(4,1) * t65 + mrSges(4,2) * t63) + (-t53 * t63 + t54 * t65) * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t66 * mrSges(3,1) - t64 * mrSges(3,2) + (t64 ^ 2 + t66 ^ 2) * t73 * pkin(1)) * pkin(1) + Ifges(4,2) * t65 ^ 2 / 0.2e1 + (Ifges(4,4) * t65 + Ifges(4,1) * t63 / 0.2e1) * t63) * qJD(1)) * qJD(1);
T = t1;
