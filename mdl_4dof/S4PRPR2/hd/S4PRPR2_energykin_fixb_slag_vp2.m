% Calculate kinetic energy for
% S4PRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR2_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:21:48
% EndTime: 2019-03-08 18:21:48
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->24), mult. (116->42), div. (0->0), fcn. (60->6), ass. (0->16)
t51 = sin(qJ(2));
t57 = qJD(1) * t51;
t53 = cos(qJ(2));
t46 = qJD(2) * pkin(2) + t53 * qJD(1);
t48 = sin(pkin(6));
t49 = cos(pkin(6));
t43 = t49 * t46 - t48 * t57;
t54 = qJD(3) ^ 2;
t52 = cos(qJ(4));
t50 = sin(qJ(4));
t47 = qJD(2) + qJD(4);
t44 = t48 * t46 + t49 * t57;
t42 = qJD(2) * pkin(3) + t43;
t41 = t50 * t42 + t52 * t44;
t40 = t52 * t42 - t50 * t44;
t1 = m(4) * (t43 ^ 2 + t44 ^ 2 + t54) / 0.2e1 + m(5) * (t40 ^ 2 + t41 ^ 2 + t54) / 0.2e1 + (m(3) * (t51 ^ 2 + t53 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t40 * mrSges(5,1) - t41 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + (t43 * mrSges(4,1) - t44 * mrSges(4,2) + (mrSges(3,1) * t53 - mrSges(3,2) * t51) * qJD(1) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2)) * qJD(2);
T  = t1;
