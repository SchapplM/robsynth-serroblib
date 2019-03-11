% Calculate kinetic energy for
% S4PPRR2
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:16:49
% EndTime: 2019-03-08 18:16:49
% DurationCPUTime: 0.04s
% Computational Cost: add. (38->21), mult. (96->38), div. (0->0), fcn. (58->6), ass. (0->15)
t47 = sin(pkin(6));
t48 = cos(pkin(6));
t50 = sin(qJ(3));
t52 = cos(qJ(3));
t43 = (-t47 * t50 + t48 * t52) * qJD(1);
t54 = qJD(1) ^ 2;
t53 = qJD(2) ^ 2;
t51 = cos(qJ(4));
t49 = sin(qJ(4));
t46 = qJD(3) + qJD(4);
t44 = (t47 * t52 + t48 * t50) * qJD(1);
t42 = qJD(3) * pkin(3) + t43;
t41 = t49 * t42 + t51 * t44;
t40 = t51 * t42 - t49 * t44;
t1 = m(5) * (t40 ^ 2 + t41 ^ 2 + t53) / 0.2e1 + m(3) * (t53 + (t47 ^ 2 + t48 ^ 2) * t54) / 0.2e1 + m(4) * (t43 ^ 2 + t44 ^ 2 + t53) / 0.2e1 + m(2) * t54 / 0.2e1 + (t40 * mrSges(5,1) - t41 * mrSges(5,2) + Ifges(5,3) * t46 / 0.2e1) * t46 + (t43 * mrSges(4,1) - t44 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3);
T  = t1;
