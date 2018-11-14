% Calculate kinetic energy for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RRRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP1_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP1_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:19
% EndTime: 2018-11-14 13:54:19
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->22), mult. (107->39), div. (0->0), fcn. (36->4), ass. (0->14)
t59 = pkin(1) * qJD(1);
t51 = qJD(1) + qJD(2);
t53 = sin(qJ(2));
t58 = t53 * t59;
t55 = cos(qJ(2));
t48 = t51 * pkin(2) + t55 * t59;
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t45 = t54 * t48 - t52 * t58;
t50 = qJD(3) + t51;
t46 = t52 * t48 + t54 * t58;
t44 = t46 ^ 2;
t43 = t50 * pkin(3) + t45;
t1 = m(4) * (t45 ^ 2 + t44) / 0.2e1 + m(5) * (qJD(4) ^ 2 + t43 ^ 2 + t44) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t53 ^ 2 + t55 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * t51 / 0.2e1 + (mrSges(3,1) * t55 - mrSges(3,2) * t53) * t59) * t51 + (t45 * mrSges(4,1) + t43 * mrSges(5,1) + (-mrSges(4,2) - mrSges(5,2)) * t46 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * t50) * t50;
T  = t1;
