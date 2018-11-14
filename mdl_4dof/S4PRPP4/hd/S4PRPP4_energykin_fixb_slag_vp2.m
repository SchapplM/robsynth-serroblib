% Calculate kinetic energy for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2018-11-14 14:09
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP4_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP4_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP4_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:13
% EndTime: 2018-11-14 14:09:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (39->23), mult. (88->37), div. (0->0), fcn. (36->4), ass. (0->12)
t44 = cos(qJ(2));
t39 = qJD(2) * pkin(2) + t44 * qJD(1);
t41 = sin(pkin(5));
t42 = cos(pkin(5));
t43 = sin(qJ(2));
t48 = qJD(1) * t43;
t37 = t41 * t39 + t42 * t48;
t36 = t42 * t39 - t41 * t48;
t45 = qJD(3) ^ 2;
t35 = qJD(2) * qJ(4) + t37;
t34 = -qJD(2) * pkin(3) + qJD(4) - t36;
t1 = m(4) * (t36 ^ 2 + t37 ^ 2 + t45) / 0.2e1 + m(5) * (t34 ^ 2 + t35 ^ 2 + t45) / 0.2e1 + (m(2) / 0.2e1 + m(3) * (t43 ^ 2 + t44 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t36 * mrSges(4,1) - t34 * mrSges(5,1) - t37 * mrSges(4,2) + t35 * mrSges(5,3) + (mrSges(3,1) * t44 - mrSges(3,2) * t43) * qJD(1) + (Ifges(3,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2)) * qJD(2);
T  = t1;
