% Calculate kinetic energy for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:34
% EndTime: 2019-12-31 16:41:34
% DurationCPUTime: 0.12s
% Computational Cost: add. (90->39), mult. (196->69), div. (0->0), fcn. (86->4), ass. (0->19)
t59 = cos(pkin(6));
t66 = t59 ^ 2;
t65 = m(3) / 0.2e1;
t58 = sin(pkin(6));
t64 = t58 ^ 2 + t66;
t53 = qJD(1) * qJ(2) + qJD(3);
t52 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t63 = -pkin(5) * qJD(1) + t52;
t61 = cos(qJ(4));
t60 = sin(qJ(4));
t54 = -qJD(1) * pkin(1) + qJD(2);
t50 = t58 * qJD(1) * pkin(3) + t53;
t49 = (-t58 * t60 + t59 * t61) * qJD(1);
t48 = (-t58 * t61 - t59 * t60) * qJD(1);
t47 = t63 * t59;
t46 = t63 * t58;
t45 = t61 * t46 + t60 * t47;
t44 = -t60 * t46 + t61 * t47;
t1 = m(5) * (t44 ^ 2 + t45 ^ 2 + t50 ^ 2) / 0.2e1 + t54 ^ 2 * t65 + m(4) * (t64 * t52 ^ 2 + t53 ^ 2) / 0.2e1 + (t50 * mrSges(5,2) - t44 * mrSges(5,3) + Ifges(5,1) * t49 / 0.2e1) * t49 + (-t50 * mrSges(5,1) + t45 * mrSges(5,3) + Ifges(5,4) * t49 + Ifges(5,2) * t48 / 0.2e1) * t48 + (t44 * mrSges(5,1) - t45 * mrSges(5,2) + Ifges(5,5) * t49 + Ifges(5,6) * t48 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t54 * mrSges(3,2) + t53 * (mrSges(4,1) * t58 + mrSges(4,2) * t59) - t64 * t52 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t65 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t66 / 0.2e1 + (-Ifges(4,4) * t59 + Ifges(4,2) * t58 / 0.2e1) * t58) * qJD(1)) * qJD(1);
T = t1;
