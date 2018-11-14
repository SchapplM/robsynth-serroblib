% Calculate kinetic energy for
% S4RPPR2
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:25
% EndTime: 2018-11-14 13:47:25
% DurationCPUTime: 0.05s
% Computational Cost: add. (66->27), mult. (115->40), div. (0->0), fcn. (32->4), ass. (0->17)
t60 = m(3) / 0.2e1;
t59 = qJ(2) * qJD(1);
t57 = qJD(3) ^ 2;
t56 = cos(qJ(4));
t55 = sin(qJ(4));
t54 = cos(pkin(6));
t53 = sin(pkin(6));
t52 = -qJD(1) + qJD(4);
t51 = -qJD(1) * pkin(1) + qJD(2);
t50 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t49 = t54 * t50;
t48 = t53 * t50 + t54 * t59;
t47 = -t53 * t59 + t49;
t46 = t49 + (-qJ(2) * t53 - pkin(3)) * qJD(1);
t45 = t55 * t46 + t56 * t48;
t44 = t56 * t46 - t55 * t48;
t1 = m(4) * (t47 ^ 2 + t48 ^ 2 + t57) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t57) / 0.2e1 + t51 ^ 2 * t60 + (t44 * mrSges(5,1) - t45 * mrSges(5,2) + Ifges(5,3) * t52 / 0.2e1) * t52 + (-t51 * mrSges(3,1) - t47 * mrSges(4,1) + t48 * mrSges(4,2) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (qJ(2) * t60 + mrSges(3,3)) * qJ(2)) * qJD(1)) * qJD(1);
T  = t1;
