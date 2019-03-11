% Calculate kinetic energy for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-03-08 18:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR1_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR1_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR1_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:53
% EndTime: 2019-03-08 18:34:54
% DurationCPUTime: 0.07s
% Computational Cost: add. (77->25), mult. (146->45), div. (0->0), fcn. (60->6), ass. (0->18)
t69 = pkin(1) * qJD(1);
t58 = qJD(1) + qJD(2);
t62 = sin(qJ(2));
t68 = t62 * t69;
t64 = cos(qJ(2));
t55 = t58 * pkin(2) + t64 * t69;
t59 = sin(pkin(7));
t60 = cos(pkin(7));
t52 = t60 * t55 - t59 * t68;
t65 = qJD(3) ^ 2;
t63 = cos(qJ(4));
t61 = sin(qJ(4));
t56 = qJD(4) + t58;
t53 = t59 * t55 + t60 * t68;
t51 = t58 * pkin(3) + t52;
t50 = t61 * t51 + t63 * t53;
t49 = t63 * t51 - t61 * t53;
t1 = m(5) * (t49 ^ 2 + t50 ^ 2 + t65) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t65) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t62 ^ 2 + t64 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t49 * mrSges(5,1) - t50 * mrSges(5,2) + Ifges(5,3) * t56 / 0.2e1) * t56 + (t52 * mrSges(4,1) - t53 * mrSges(4,2) + (mrSges(3,1) * t64 - mrSges(3,2) * t62) * t69 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t58) * t58;
T  = t1;
