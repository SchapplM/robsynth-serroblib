% Calculate kinetic energy for
% S4PRRP2
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:58
% EndTime: 2019-03-08 18:23:58
% DurationCPUTime: 0.06s
% Computational Cost: add. (42->21), mult. (85->36), div. (0->0), fcn. (36->4), ass. (0->12)
t43 = sin(qJ(2));
t47 = qJD(1) * t43;
t45 = cos(qJ(2));
t39 = qJD(2) * pkin(2) + qJD(1) * t45;
t42 = sin(qJ(3));
t44 = cos(qJ(3));
t36 = t39 * t44 - t42 * t47;
t41 = qJD(2) + qJD(3);
t37 = t39 * t42 + t44 * t47;
t35 = t37 ^ 2;
t34 = pkin(3) * t41 + t36;
t1 = m(4) * (t36 ^ 2 + t35) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t34 ^ 2 + t35) / 0.2e1 + (m(3) * (t43 ^ 2 + t45 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t45 - mrSges(3,2) * t43) * qJD(1)) * qJD(2) + (t36 * mrSges(4,1) + t34 * mrSges(5,1) + (-mrSges(4,2) - mrSges(5,2)) * t37 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t41) * t41;
T  = t1;
