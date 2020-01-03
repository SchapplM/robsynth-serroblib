% Calculate kinetic energy for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:24
% EndTime: 2019-12-31 17:34:24
% DurationCPUTime: 0.10s
% Computational Cost: add. (74->49), mult. (161->70), div. (0->0), fcn. (62->4), ass. (0->17)
t60 = cos(qJ(4));
t66 = t60 * qJD(1);
t61 = cos(qJ(3));
t65 = t61 * qJD(2);
t64 = qJ(5) * qJD(3);
t59 = sin(qJ(3));
t56 = qJD(3) * pkin(6) + t59 * qJD(2);
t58 = sin(qJ(4));
t53 = -t58 * qJD(1) + t60 * t56;
t63 = qJD(1) ^ 2;
t62 = qJD(2) ^ 2;
t57 = -qJD(3) * pkin(3) - t65;
t54 = -t65 + qJD(5) + (-pkin(4) * t60 - pkin(3)) * qJD(3);
t52 = -t58 * t56 - t66;
t51 = t60 * t64 + t53;
t50 = qJD(4) * pkin(4) - t66 + (-t56 - t64) * t58;
t1 = m(5) * (t52 ^ 2 + t53 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (t63 + (t59 ^ 2 + t61 ^ 2) * t62) / 0.2e1 + m(2) * t63 / 0.2e1 + m(3) * (t62 + t63) / 0.2e1 + m(6) * (t50 ^ 2 + t51 ^ 2 + t54 ^ 2) / 0.2e1 + (t52 * mrSges(5,1) + t50 * mrSges(6,1) - t53 * mrSges(5,2) - t51 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(4)) * qJD(4) + (Ifges(4,3) * qJD(3) / 0.2e1 + (t61 * mrSges(4,1) - t59 * mrSges(4,2)) * qJD(2) + (-t57 * mrSges(5,1) - t54 * mrSges(6,1) + t53 * mrSges(5,3) + t51 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * qJD(3) * t60 + (Ifges(5,6) + Ifges(6,6)) * qJD(4)) * t60 + (t57 * mrSges(5,2) + t54 * mrSges(6,2) - t52 * mrSges(5,3) - t50 * mrSges(6,3) + (Ifges(5,5) + Ifges(6,5)) * qJD(4) + ((Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t58 + (Ifges(5,4) + Ifges(6,4)) * t60) * qJD(3)) * t58) * qJD(3);
T = t1;
