% Calculate kinetic energy for
% S4RPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:00
% EndTime: 2019-12-31 16:47:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (69->43), mult. (144->62), div. (0->0), fcn. (36->2), ass. (0->11)
t42 = qJD(2) + (-pkin(1) - pkin(5)) * qJD(1);
t48 = t42 * mrSges(4,3);
t47 = qJD(1) ^ 2;
t46 = cos(qJ(3));
t45 = sin(qJ(3));
t44 = t47 * qJ(2) ^ 2;
t43 = -qJD(1) * pkin(1) + qJD(2);
t40 = qJD(3) * qJ(4) + t42 * t45;
t39 = (pkin(3) * t45 - qJ(4) * t46 + qJ(2)) * qJD(1);
t38 = -qJD(3) * pkin(3) - t42 * t46 + qJD(4);
t1 = m(4) * (t44 + (t45 ^ 2 + t46 ^ 2) * t42 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + m(3) * (t43 ^ 2 + t44) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t47 + (-t38 * mrSges(5,1) + t40 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(3) + (t46 * mrSges(4,1) - t45 * mrSges(4,2)) * t42) * qJD(3) + (t43 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + t38 * mrSges(5,2) - t39 * mrSges(5,3) + (-t48 + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(1)) * t46 + (Ifges(5,4) + Ifges(4,5)) * qJD(3)) * t46 + (-t40 * mrSges(5,2) + t39 * mrSges(5,1) + (qJ(2) * mrSges(4,1) + (-Ifges(4,4) + Ifges(5,5)) * t46) * qJD(1) + (-t48 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(1)) * t45 + (-Ifges(4,6) + Ifges(5,6)) * qJD(3)) * t45) * qJD(1);
T = t1;
