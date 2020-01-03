% Calculate kinetic energy for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:23
% EndTime: 2019-12-31 16:18:23
% DurationCPUTime: 0.11s
% Computational Cost: add. (47->30), mult. (119->54), div. (0->0), fcn. (62->6), ass. (0->16)
t56 = sin(pkin(7));
t57 = cos(pkin(7));
t59 = sin(qJ(3));
t61 = cos(qJ(3));
t52 = (t56 * t59 - t57 * t61) * qJD(1);
t65 = qJD(3) / 0.2e1;
t53 = (t56 * t61 + t57 * t59) * qJD(1);
t63 = qJD(1) ^ 2;
t62 = qJD(2) ^ 2;
t60 = cos(qJ(4));
t58 = sin(qJ(4));
t51 = qJD(3) * pkin(5) + t53;
t50 = -qJD(3) * pkin(3) + t52;
t49 = t58 * qJD(2) + t60 * t51;
t48 = t60 * qJD(2) - t58 * t51;
t1 = m(5) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + m(3) * (t62 + (t56 ^ 2 + t57 ^ 2) * t63) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t62) / 0.2e1 + m(2) * t63 / 0.2e1 + (t48 * mrSges(5,1) - t49 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (-t52 * mrSges(4,1) - t53 * mrSges(4,2) + Ifges(4,3) * t65 + (Ifges(5,2) * t60 * t65 - t50 * mrSges(5,1) + t49 * mrSges(5,3) + Ifges(5,6) * qJD(4)) * t60 + (t50 * mrSges(5,2) - t48 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t60 + Ifges(5,1) * t58 / 0.2e1) * qJD(3)) * t58) * qJD(3);
T = t1;
