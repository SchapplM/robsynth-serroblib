% Calculate kinetic energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR6_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR6_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR6_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:35
% EndTime: 2019-12-31 16:40:35
% DurationCPUTime: 0.15s
% Computational Cost: add. (86->49), mult. (228->78), div. (0->0), fcn. (107->4), ass. (0->23)
t65 = qJD(1) ^ 2;
t70 = t65 * qJ(2) ^ 2;
t62 = cos(pkin(6));
t69 = qJD(1) * t62;
t61 = sin(pkin(6));
t68 = t61 * qJD(1);
t55 = qJ(2) * t68 + qJD(3);
t67 = qJ(3) * t61 + pkin(1);
t64 = cos(qJ(4));
t63 = sin(qJ(4));
t60 = t62 ^ 2;
t59 = t61 ^ 2;
t58 = -qJD(1) * pkin(1) + qJD(2);
t57 = t60 * t70;
t54 = (-pkin(5) + qJ(2)) * t69;
t53 = -pkin(5) * t68 + t55;
t52 = (t61 * t64 - t62 * t63) * qJD(1);
t51 = (-t61 * t63 - t62 * t64) * qJD(1);
t50 = qJD(2) + (-pkin(2) * t62 - t67) * qJD(1);
t49 = -qJD(2) + ((pkin(2) + pkin(3)) * t62 + t67) * qJD(1);
t48 = t63 * t53 + t64 * t54;
t47 = t64 * t53 - t63 * t54;
t1 = m(3) * (t58 ^ 2 + t59 * t70 + t57) / 0.2e1 + m(4) * (t50 ^ 2 + t55 ^ 2 + t57) / 0.2e1 + m(5) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + (t49 * mrSges(5,2) - t47 * mrSges(5,3) + Ifges(5,1) * t52 / 0.2e1) * t52 + (-t49 * mrSges(5,1) + t48 * mrSges(5,3) + Ifges(5,4) * t52 + Ifges(5,2) * t51 / 0.2e1) * t51 + (t47 * mrSges(5,1) - t48 * mrSges(5,2) + Ifges(5,5) * t52 + Ifges(5,6) * t51 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + ((-t58 * mrSges(3,1) - t50 * mrSges(4,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t69) * t62 + (t58 * mrSges(3,2) + t55 * mrSges(4,2) - t50 * mrSges(4,3) + ((Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t61 + (Ifges(3,4) - Ifges(4,5)) * t62) * qJD(1)) * t61) * qJD(1) + (Ifges(2,3) / 0.2e1 + (mrSges(4,2) * t60 + (t59 + t60) * mrSges(3,3)) * qJ(2)) * t65;
T = t1;
