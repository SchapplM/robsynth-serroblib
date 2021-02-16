% Return the minimum parameter vector for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% m [6x1]
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
% MPV [18x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:56
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPPRP1_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_convert_par2_MPV_fixb: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t9 = (-Ifges(5,2) - Ifges(6,2));
t8 = m(5) * pkin(6) + mrSges(5,3);
t2 = sin(pkin(7));
t4 = cos(pkin(7));
t7 = mrSges(3,1) * t4 - mrSges(3,2) * t2;
t6 = -2 * mrSges(5,3) * pkin(6) - Ifges(4,1) + t9;
t5 = pkin(6) ^ 2;
t3 = cos(pkin(8));
t1 = sin(pkin(8));
t10 = [(m(5) * t5) + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t1 * (t8 * pkin(3) + Ifges(4,4)) + (Ifges(4,2) + (pkin(3) ^ 2 - t5) * m(5) + t6) * t3) * t3 + 0.2e1 * t7 * pkin(1) - t6; mrSges(2,1) + t7; mrSges(3,1) * t2 + mrSges(3,2) * t4 + mrSges(2,2); m(3); ((m(5) * pkin(3) + mrSges(4,1)) * t3 + t1 * (-mrSges(4,2) + t8)) / t3; mrSges(4,3); m(4) + m(5); Ifges(5,1) + Ifges(6,1) + t9; Ifges(5,4) + Ifges(6,4); Ifges(5,5) + Ifges(6,5); Ifges(5,6) + Ifges(6,6); Ifges(5,3) + Ifges(6,3); mrSges(5,1); mrSges(5,2); mrSges(6,1); mrSges(6,2); mrSges(6,3); m(6);];
MPV = t10;
