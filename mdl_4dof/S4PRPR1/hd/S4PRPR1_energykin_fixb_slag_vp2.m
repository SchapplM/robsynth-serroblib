% Calculate kinetic energy for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR1_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:11
% EndTime: 2018-11-14 13:42:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (31->20), mult. (58->29), div. (0->0), fcn. (8->2), ass. (0->11)
t56 = m(4) / 0.2e1;
t54 = qJD(2) * qJ(3);
t53 = qJD(1) ^ 2;
t51 = cos(qJ(4));
t50 = sin(qJ(4));
t49 = -qJD(2) + qJD(4);
t48 = -qJD(2) * pkin(2) + qJD(3);
t47 = qJD(3) + (-pkin(2) - pkin(3)) * qJD(2);
t46 = t50 * t47 + t51 * t54;
t45 = t51 * t47 - t50 * t54;
t1 = (t48 ^ 2 + t53) * t56 + m(5) * (t45 ^ 2 + t46 ^ 2 + t53) / 0.2e1 + (t45 * mrSges(5,1) - t46 * mrSges(5,2) + Ifges(5,3) * t49 / 0.2e1) * t49 + (m(3) + m(2)) * t53 / 0.2e1 + (-t48 * mrSges(4,1) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + (qJ(3) * t56 + mrSges(4,3)) * qJ(3)) * qJD(2)) * qJD(2);
T  = t1;
