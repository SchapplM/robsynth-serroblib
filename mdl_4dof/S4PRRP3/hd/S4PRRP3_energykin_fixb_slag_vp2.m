% Calculate kinetic energy for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:18
% EndTime: 2018-11-14 14:12:18
% DurationCPUTime: 0.07s
% Computational Cost: add. (42->21), mult. (85->36), div. (0->0), fcn. (36->4), ass. (0->12)
t44 = sin(qJ(2));
t48 = qJD(1) * t44;
t46 = cos(qJ(2));
t40 = qJD(2) * pkin(2) + qJD(1) * t46;
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t37 = t40 * t45 - t43 * t48;
t42 = qJD(2) + qJD(3);
t38 = t40 * t43 + t45 * t48;
t36 = t38 ^ 2;
t35 = pkin(3) * t42 + t37;
t1 = m(4) * (t37 ^ 2 + t36) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t35 ^ 2 + t36) / 0.2e1 + (m(2) / 0.2e1 + m(3) * (t44 ^ 2 + t46 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t46 * mrSges(3,1) - t44 * mrSges(3,2)) * qJD(1)) * qJD(2) + (t37 * mrSges(4,1) + t35 * mrSges(5,1) + (-mrSges(4,2) - mrSges(5,2)) * t38 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t42) * t42;
T  = t1;
