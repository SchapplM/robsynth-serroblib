% Calculate kinetic energy for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR2_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR2_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:29
% EndTime: 2019-07-18 18:16:29
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->24), mult. (96->41), div. (0->0), fcn. (24->4), ass. (0->14)
t58 = pkin(1) * qJD(1);
t50 = qJD(1) + qJD(2);
t54 = cos(qJ(2));
t57 = -t54 * t58 + qJD(3);
t53 = cos(qJ(4));
t52 = sin(qJ(2));
t51 = sin(qJ(4));
t48 = qJD(4) - t50;
t47 = t50 * qJ(3) + t52 * t58;
t46 = -t50 * pkin(2) + t57;
t45 = (-pkin(2) - pkin(3)) * t50 + t57;
t44 = t51 * t45 + t53 * t47;
t43 = t53 * t45 - t51 * t47;
t1 = m(4) * (t46 ^ 2 + t47 ^ 2) / 0.2e1 + m(5) * (t43 ^ 2 + t44 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t52 ^ 2 + t54 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t43 * mrSges(5,1) - t44 * mrSges(5,2) + Ifges(5,3) * t48 / 0.2e1) * t48 + (-t46 * mrSges(4,1) + t47 * mrSges(4,3) + (mrSges(3,1) * t54 - mrSges(3,2) * t52) * t58 + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t50) * t50;
T  = t1;
