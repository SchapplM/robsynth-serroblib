% Calculate kinetic energy for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(2,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:21
% EndTime: 2019-07-18 13:27:21
% DurationCPUTime: 0.04s
% Computational Cost: add. (35->19), mult. (73->36), div. (0->0), fcn. (20->4), ass. (0->14)
t45 = qJD(2) * pkin(1);
t36 = qJD(2) + qJD(3);
t38 = sin(qJ(3));
t44 = t38 * t45;
t42 = qJD(1) ^ 2;
t41 = qJD(2) ^ 2;
t40 = cos(qJ(3));
t39 = cos(qJ(4));
t37 = sin(qJ(4));
t35 = qJD(4) + t36;
t34 = t36 * pkin(2) + t40 * t45;
t33 = t37 * t34 + t39 * t44;
t32 = t39 * t34 - t37 * t44;
t1 = t41 * Ifges(3,3) / 0.2e1 + m(4) * (t42 + (t38 ^ 2 + t40 ^ 2) * pkin(1) ^ 2 * t41) / 0.2e1 + m(5) * (t32 ^ 2 + t33 ^ 2 + t42) / 0.2e1 + (Ifges(4,3) * t36 / 0.2e1 + (mrSges(4,1) * t40 - mrSges(4,2) * t38) * t45) * t36 + (t32 * mrSges(5,1) - t33 * mrSges(5,2) + Ifges(5,3) * t35 / 0.2e1) * t35 + (m(3) + m(2)) * t42 / 0.2e1;
T  = t1;
