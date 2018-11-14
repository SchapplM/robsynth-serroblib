% Calculate kinetic energy for
% S3PRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
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
% Datum: 2018-11-14 10:13
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3PRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRR2_energykin_fixb_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRR2_energykin_fixb_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3PRR2_energykin_fixb_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRR2_energykin_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRR2_energykin_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRR2_energykin_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:12:44
% EndTime: 2018-11-14 10:12:44
% DurationCPUTime: 0.03s
% Computational Cost: add. (21->14), mult. (51->29), div. (0->0), fcn. (20->4), ass. (0->10)
t34 = sin(qJ(2));
t38 = qJD(1) * t34;
t36 = cos(qJ(2));
t35 = cos(qJ(3));
t33 = sin(qJ(3));
t32 = qJD(2) + qJD(3);
t31 = qJD(2) * pkin(2) + t36 * qJD(1);
t30 = t33 * t31 + t35 * t38;
t29 = t35 * t31 - t33 * t38;
t1 = m(4) * (t29 ^ 2 + t30 ^ 2) / 0.2e1 + (m(2) / 0.2e1 + m(3) * (t34 ^ 2 + t36 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t29 * mrSges(4,1) - t30 * mrSges(4,2) + Ifges(4,3) * t32 / 0.2e1) * t32 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t36 * mrSges(3,1) - t34 * mrSges(3,2)) * qJD(1)) * qJD(2);
T  = t1;
