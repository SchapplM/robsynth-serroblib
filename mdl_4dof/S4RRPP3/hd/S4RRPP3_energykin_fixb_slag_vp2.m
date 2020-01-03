% Calculate kinetic energy for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:20
% EndTime: 2019-12-31 16:57:20
% DurationCPUTime: 0.17s
% Computational Cost: add. (132->58), mult. (337->84), div. (0->0), fcn. (180->4), ass. (0->19)
t64 = pkin(5) * mrSges(3,3);
t63 = pkin(5) + qJ(3);
t56 = sin(qJ(2));
t60 = t56 * qJD(1);
t51 = qJD(2) * pkin(2) - t63 * t60;
t57 = cos(qJ(2));
t61 = qJD(1) * t57;
t52 = t63 * t61;
t55 = sin(pkin(6));
t62 = cos(pkin(6));
t46 = t55 * t51 + t62 * t52;
t45 = t62 * t51 - t55 * t52;
t53 = qJD(3) + (-pkin(2) * t57 - pkin(1)) * qJD(1);
t49 = (t55 * t57 + t62 * t56) * qJD(1);
t48 = t55 * t60 - t62 * t61;
t44 = qJD(2) * qJ(4) + t46;
t43 = -qJD(2) * pkin(3) + qJD(4) - t45;
t42 = t48 * pkin(3) - t49 * qJ(4) + t53;
t1 = m(4) * (t45 ^ 2 + t46 ^ 2 + t53 ^ 2) / 0.2e1 + m(5) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + (t53 * mrSges(4,2) + t43 * mrSges(5,2) - t45 * mrSges(4,3) - t42 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t49) * t49 + (t53 * mrSges(4,1) + t42 * mrSges(5,1) - t44 * mrSges(5,2) - t46 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t48 + (-Ifges(4,4) + Ifges(5,5)) * t49) * t48 + (t45 * mrSges(4,1) - t43 * mrSges(5,1) - t46 * mrSges(4,2) + t44 * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(5,4) + Ifges(4,5)) * t49 + (-Ifges(4,6) + Ifges(5,6)) * t48 + (Ifges(3,5) * t56 + Ifges(3,6) * t57 + (-mrSges(3,1) * t56 - mrSges(3,2) * t57) * pkin(5)) * qJD(1)) * qJD(2) + (m(3) * (pkin(1) ^ 2 + (t56 ^ 2 + t57 ^ 2) * pkin(5) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t64 + Ifges(3,2) / 0.2e1) * t57) * t57 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t57 + (t64 + Ifges(3,1) / 0.2e1) * t56) * t56) * qJD(1) ^ 2;
T = t1;
