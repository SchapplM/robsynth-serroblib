% Calculate kinetic energy for
% S3RRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2,d3]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:16
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRR1_energykin_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRR1_energykin_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S3RRR1_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRR1_energykin_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRR1_energykin_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRR1_energykin_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:53
% EndTime: 2018-11-14 10:15:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (31->15), mult. (65->32), div. (0->0), fcn. (20->4), ass. (0->12)
t49 = pkin(1) * qJD(1);
t41 = qJD(1) + qJD(2);
t43 = sin(qJ(2));
t48 = t43 * t49;
t45 = cos(qJ(2));
t44 = cos(qJ(3));
t42 = sin(qJ(3));
t40 = qJD(3) + t41;
t39 = t41 * pkin(2) + t45 * t49;
t38 = t42 * t39 + t44 * t48;
t37 = t44 * t39 - t42 * t48;
t1 = m(4) * (t37 ^ 2 + t38 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t43 ^ 2 + t45 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * t41 / 0.2e1 + (mrSges(3,1) * t45 - mrSges(3,2) * t43) * t49) * t41 + (t37 * mrSges(4,1) - t38 * mrSges(4,2) + Ifges(4,3) * t40 / 0.2e1) * t40;
T  = t1;
