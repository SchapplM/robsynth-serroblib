% Calculate kinetic energy for
% S4PRPR3
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:11:14
% EndTime: 2018-11-14 14:11:14
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->24), mult. (116->42), div. (0->0), fcn. (60->6), ass. (0->16)
t52 = sin(qJ(2));
t58 = qJD(1) * t52;
t54 = cos(qJ(2));
t47 = qJD(2) * pkin(2) + t54 * qJD(1);
t49 = sin(pkin(6));
t50 = cos(pkin(6));
t44 = t50 * t47 - t49 * t58;
t55 = qJD(3) ^ 2;
t53 = cos(qJ(4));
t51 = sin(qJ(4));
t48 = qJD(2) + qJD(4);
t45 = t49 * t47 + t50 * t58;
t43 = qJD(2) * pkin(3) + t44;
t42 = t51 * t43 + t53 * t45;
t41 = t53 * t43 - t51 * t45;
t1 = m(4) * (t44 ^ 2 + t45 ^ 2 + t55) / 0.2e1 + m(5) * (t41 ^ 2 + t42 ^ 2 + t55) / 0.2e1 + (m(2) / 0.2e1 + m(3) * (t52 ^ 2 + t54 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t41 * mrSges(5,1) - t42 * mrSges(5,2) + Ifges(5,3) * t48 / 0.2e1) * t48 + (t44 * mrSges(4,1) - t45 * mrSges(4,2) + (mrSges(3,1) * t54 - mrSges(3,2) * t52) * qJD(1) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2)) * qJD(2);
T  = t1;
