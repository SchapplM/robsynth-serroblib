% Return the minimum parameter vector for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% MPV [22x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:07
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRPR7_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t14 = pkin(7) ^ 2;
t15 = pkin(4) ^ 2;
t22 = Ifges(5,1) + Ifges(6,2);
t16 = (Ifges(5,2) - t22 + (-t14 + t15) * m(6));
t24 = (pkin(7) * mrSges(6,3));
t1 = -2 * t24 + t16;
t4 = pkin(7) * m(6) + mrSges(6,3);
t18 = t4 * pkin(4) + Ifges(5,4);
t10 = sin(pkin(9));
t12 = cos(pkin(9));
t23 = t10 * t12;
t8 = 2 * t24;
t7 = t12 ^ 2;
t21 = 0.2e1 * t7;
t20 = t18 * t23;
t11 = sin(pkin(8));
t13 = cos(pkin(8));
t3 = m(4) * pkin(6) - mrSges(3,2) + mrSges(4,3);
t5 = m(4) * pkin(2) + mrSges(3,1);
t17 = t3 * t11 + t5 * t13;
t2 = [t1 * t7 + 0.2e1 * t20 + m(6) * t14 + t8 + ((pkin(2) ^ 2 + pkin(6) ^ 2) * m(4)) + (2 * mrSges(4,3) * pkin(6)) + Ifges(4,2) + Ifges(2,3) + Ifges(3,3) + 0.2e1 * t17 * pkin(1) + t22; mrSges(2,1) + t17; t5 * t11 - t3 * t13 + mrSges(2,2); m(3) + m(4); (-t16 + t8) * t21 - 0.4e1 * t20 + Ifges(4,1) - Ifges(4,2) + t1; -t1 * t23 + t18 * t21 + Ifges(4,4) - t18; t12 * Ifges(5,5) - t10 * Ifges(5,6) + Ifges(4,5); t10 * Ifges(5,5) + t12 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + Ifges(6,2) + t8 + (t14 + t15) * m(6); mrSges(4,1); mrSges(4,2); pkin(4) * m(6) + mrSges(5,1); mrSges(5,2) - t4; mrSges(5,3); m(5) + m(6); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t2;
