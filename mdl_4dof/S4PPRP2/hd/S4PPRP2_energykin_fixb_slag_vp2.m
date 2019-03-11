% Calculate kinetic energy for
% S4PPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
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
% Datum: 2019-03-08 18:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:13:14
% EndTime: 2019-03-08 18:13:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->20), mult. (72->33), div. (0->0), fcn. (34->4), ass. (0->11)
t39 = sin(pkin(5));
t40 = cos(pkin(5));
t41 = sin(qJ(3));
t42 = cos(qJ(3));
t35 = (t39 * t41 - t40 * t42) * qJD(1);
t36 = (t39 * t42 + t40 * t41) * qJD(1);
t45 = qJD(1) ^ 2;
t44 = qJD(2) ^ 2;
t34 = qJD(3) * qJ(4) + t36;
t33 = -qJD(3) * pkin(3) + qJD(4) + t35;
t1 = m(5) * (t33 ^ 2 + t34 ^ 2 + t44) / 0.2e1 + m(2) * t45 / 0.2e1 + m(3) * (t44 + (t39 ^ 2 + t40 ^ 2) * t45) / 0.2e1 + m(4) * (t35 ^ 2 + t36 ^ 2 + t44) / 0.2e1 + (-t35 * mrSges(4,1) - t33 * mrSges(5,1) - t36 * mrSges(4,2) + t34 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3)) * qJD(3);
T  = t1;
