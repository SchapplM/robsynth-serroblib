% Calculate kinetic energy for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:06
% EndTime: 2019-12-31 17:18:07
% DurationCPUTime: 0.23s
% Computational Cost: add. (144->59), mult. (321->86), div. (0->0), fcn. (164->4), ass. (0->19)
t66 = pkin(5) * mrSges(3,3);
t59 = sin(qJ(2));
t61 = cos(qJ(2));
t50 = (-pkin(2) * t61 - pkin(6) * t59 - pkin(1)) * qJD(1);
t64 = t61 * qJD(1);
t55 = pkin(5) * t64 + qJD(2) * pkin(6);
t58 = sin(qJ(3));
t60 = cos(qJ(3));
t46 = t58 * t50 + t60 * t55;
t65 = t59 * qJD(1);
t45 = t60 * t50 - t58 * t55;
t54 = -qJD(2) * pkin(2) + pkin(5) * t65;
t56 = qJD(3) - t64;
t52 = t58 * qJD(2) + t60 * t65;
t51 = t60 * qJD(2) - t58 * t65;
t47 = -t51 * pkin(3) + qJD(4) + t54;
t44 = t51 * qJ(4) + t46;
t43 = t56 * pkin(3) - t52 * qJ(4) + t45;
t1 = m(4) * (t45 ^ 2 + t46 ^ 2 + t54 ^ 2) / 0.2e1 + m(5) * (t43 ^ 2 + t44 ^ 2 + t47 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + (t45 * mrSges(4,1) + t43 * mrSges(5,1) - t46 * mrSges(4,2) - t44 * mrSges(5,2) + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t56) * t56 + (t54 * mrSges(4,2) + t47 * mrSges(5,2) - t45 * mrSges(4,3) - t43 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t52 + (Ifges(4,5) + Ifges(5,5)) * t56) * t52 + (-t54 * mrSges(4,1) - t47 * mrSges(5,1) + t46 * mrSges(4,3) + t44 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,2) / 0.2e1) * t51 + (Ifges(4,6) + Ifges(5,6)) * t56 + (Ifges(4,4) + Ifges(5,4)) * t52) * t51 + ((m(3) * (pkin(1) ^ 2 + (t59 ^ 2 + t61 ^ 2) * pkin(5) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t66 + Ifges(3,2) / 0.2e1) * t61) * t61 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t61 + (t66 + Ifges(3,1) / 0.2e1) * t59) * t59) * qJD(1) + ((-pkin(5) * mrSges(3,2) + Ifges(3,6)) * t61 + (-pkin(5) * mrSges(3,1) + Ifges(3,5)) * t59) * qJD(2)) * qJD(1);
T = t1;
