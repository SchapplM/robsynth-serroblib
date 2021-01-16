% Return the minimum parameter vector for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% m [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRPRRP3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t4 = m(5) + m(6);
t11 = (-Ifges(6,2) - Ifges(7,2));
t10 = 2 * pkin(9) * mrSges(6,3) - t11;
t9 = pkin(9) * m(6) + mrSges(6,3);
t7 = (pkin(4) ^ 2);
t6 = (pkin(8) ^ 2);
t5 = pkin(9) ^ 2;
t3 = cos(pkin(11));
t2 = sin(pkin(11));
t1 = [m(2) + m(3); ((t6 + t7) * m(6)) + (m(5) * t6) + (2 * pkin(8) * mrSges(5,3)) + Ifges(5,2) + Ifges(3,3) + Ifges(4,1) + (0.2e1 * Ifges(4,4) * t2 + (t4 * pkin(3) ^ 2 - Ifges(4,1) + Ifges(4,2)) * t3) * t3; mrSges(3,1); mrSges(3,2); ((t4 * pkin(3) + mrSges(4,1)) * t3 - mrSges(4,2) * t2) / t3; t4 * pkin(8) + mrSges(4,3) + mrSges(5,3); m(4) + t4; Ifges(5,1) - Ifges(5,2) + ((t5 - t7) * m(6)) + t10; t9 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t5 + t7) * m(6) + t10; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2) - t9; Ifges(6,1) + Ifges(7,1) + t11; Ifges(6,4) + Ifges(7,4); Ifges(6,5) + Ifges(7,5); Ifges(6,6) + Ifges(7,6); Ifges(6,3) + Ifges(7,3); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV = t1;
