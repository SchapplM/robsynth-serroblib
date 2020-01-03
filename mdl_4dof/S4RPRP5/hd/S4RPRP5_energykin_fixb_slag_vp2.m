% Calculate kinetic energy for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:43
% EndTime: 2019-12-31 16:44:44
% DurationCPUTime: 0.19s
% Computational Cost: add. (119->52), mult. (317->77), div. (0->0), fcn. (176->4), ass. (0->21)
t60 = cos(pkin(6));
t70 = t60 ^ 2;
t69 = m(3) / 0.2e1;
t68 = cos(qJ(3));
t67 = pkin(5) + qJ(2);
t59 = sin(pkin(6));
t64 = t59 * qJD(1);
t52 = t67 * t64;
t65 = qJD(1) * t60;
t53 = t67 * t65;
t61 = sin(qJ(3));
t47 = -t61 * t52 + t68 * t53;
t46 = -t68 * t52 - t61 * t53;
t54 = qJD(2) + (-pkin(2) * t60 - pkin(1)) * qJD(1);
t56 = -qJD(1) * pkin(1) + qJD(2);
t51 = (t68 * t59 + t60 * t61) * qJD(1);
t50 = t61 * t64 - t68 * t65;
t45 = qJD(3) * qJ(4) + t47;
t44 = -qJD(3) * pkin(3) + qJD(4) - t46;
t43 = t50 * pkin(3) - t51 * qJ(4) + t54;
t1 = t56 ^ 2 * t69 + m(4) * (t46 ^ 2 + t47 ^ 2 + t54 ^ 2) / 0.2e1 + m(5) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + (t54 * mrSges(4,2) + t44 * mrSges(5,2) - t46 * mrSges(4,3) - t43 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t51) * t51 + (t54 * mrSges(4,1) + t43 * mrSges(5,1) - t45 * mrSges(5,2) - t47 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t50 + (-Ifges(4,4) + Ifges(5,5)) * t51) * t50 + (t46 * mrSges(4,1) - t44 * mrSges(5,1) - t47 * mrSges(4,2) + t45 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (Ifges(5,4) + Ifges(4,5)) * t51 + (-Ifges(4,6) + Ifges(5,6)) * t50) * qJD(3) + (t56 * (-mrSges(3,1) * t60 + mrSges(3,2) * t59) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t69 + mrSges(3,3)) * (t59 ^ 2 + t70) * qJ(2) + Ifges(3,2) * t70 / 0.2e1 + (Ifges(3,4) * t60 + Ifges(3,1) * t59 / 0.2e1) * t59) * qJD(1)) * qJD(1);
T = t1;
