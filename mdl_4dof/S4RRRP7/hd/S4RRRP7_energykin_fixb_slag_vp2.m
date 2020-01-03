% Calculate kinetic energy for
% S4RRRP7
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRP7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:06
% EndTime: 2019-12-31 17:20:06
% DurationCPUTime: 0.23s
% Computational Cost: add. (144->59), mult. (317->86), div. (0->0), fcn. (160->4), ass. (0->19)
t66 = pkin(5) * mrSges(3,3);
t65 = cos(qJ(3));
t59 = sin(qJ(2));
t60 = cos(qJ(2));
t49 = (-pkin(2) * t60 - pkin(6) * t59 - pkin(1)) * qJD(1);
t63 = t60 * qJD(1);
t54 = pkin(5) * t63 + qJD(2) * pkin(6);
t58 = sin(qJ(3));
t47 = t58 * t49 + t65 * t54;
t64 = t59 * qJD(1);
t53 = -qJD(2) * pkin(2) + pkin(5) * t64;
t46 = t65 * t49 - t58 * t54;
t55 = qJD(3) - t63;
t51 = t58 * qJD(2) + t65 * t64;
t50 = -t65 * qJD(2) + t58 * t64;
t45 = t50 * pkin(3) - t51 * qJ(4) + t53;
t44 = t55 * qJ(4) + t47;
t43 = -t55 * pkin(3) + qJD(4) - t46;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(5) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(4) * (t46 ^ 2 + t47 ^ 2 + t53 ^ 2) / 0.2e1 + (t46 * mrSges(4,1) - t43 * mrSges(5,1) - t47 * mrSges(4,2) + t44 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t55) * t55 + (t53 * mrSges(4,2) + t43 * mrSges(5,2) - t46 * mrSges(4,3) - t45 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t51 + (Ifges(5,4) + Ifges(4,5)) * t55) * t51 + (t53 * mrSges(4,1) + t45 * mrSges(5,1) - t44 * mrSges(5,2) - t47 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t50 + (-Ifges(4,6) + Ifges(5,6)) * t55 + (-Ifges(4,4) + Ifges(5,5)) * t51) * t50 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t59 ^ 2 + t60 ^ 2) * pkin(5) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t66) * t60) * t60 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t60 + (Ifges(3,1) / 0.2e1 + t66) * t59) * t59) * qJD(1) + ((-pkin(5) * mrSges(3,2) + Ifges(3,6)) * t60 + (-pkin(5) * mrSges(3,1) + Ifges(3,5)) * t59) * qJD(2)) * qJD(1);
T = t1;
