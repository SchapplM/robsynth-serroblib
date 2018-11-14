% Calculate kinetic energy for
% S4PPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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

function T = S4PPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR3_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR3_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR3_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:08:15
% EndTime: 2018-11-14 14:08:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (38->21), mult. (96->38), div. (0->0), fcn. (58->6), ass. (0->15)
t48 = sin(pkin(6));
t49 = cos(pkin(6));
t51 = sin(qJ(3));
t53 = cos(qJ(3));
t44 = (-t48 * t51 + t49 * t53) * qJD(1);
t55 = qJD(1) ^ 2;
t54 = qJD(2) ^ 2;
t52 = cos(qJ(4));
t50 = sin(qJ(4));
t47 = qJD(3) + qJD(4);
t45 = (t48 * t53 + t49 * t51) * qJD(1);
t43 = qJD(3) * pkin(3) + t44;
t42 = t50 * t43 + t52 * t45;
t41 = t52 * t43 - t50 * t45;
t1 = m(2) * t55 / 0.2e1 + m(5) * (t41 ^ 2 + t42 ^ 2 + t54) / 0.2e1 + m(3) * (t54 + (t48 ^ 2 + t49 ^ 2) * t55) / 0.2e1 + m(4) * (t44 ^ 2 + t45 ^ 2 + t54) / 0.2e1 + (t41 * mrSges(5,1) - t42 * mrSges(5,2) + Ifges(5,3) * t47 / 0.2e1) * t47 + (t44 * mrSges(4,1) - t45 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3);
T  = t1;
