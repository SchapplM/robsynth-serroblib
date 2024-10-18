% Return the minimum parameter vector for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
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
% MPV [25x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5PRRRR11_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t17 = 2 * pkin(9);
t3 = (m(5) + m(6));
t8 = (pkin(3) ^ 2);
t16 = (t3 * t8);
t7 = (pkin(4) ^ 2);
t15 = (t7 * m(6));
t2 = pkin(8) + pkin(9);
t14 = pkin(8) ^ 2 + t8;
t13 = (mrSges(5,3) + mrSges(6,3));
t12 = m(4) + t3;
t11 = -m(5) * pkin(8) - t13;
t10 = -t2 * m(6) + t11;
t6 = (pkin(7) ^ 2);
t4 = m(4) + m(5);
t1 = sin(pkin(5));
t5 = [m(2) + m(3) + t12; ((pkin(8) * t17 + pkin(9) ^ 2 + t14 + t6 + t7) * m(6) + t4 * t6 + mrSges(6,3) * t17 + Ifges(4,2) + Ifges(5,2) + Ifges(6,2) + t14 * m(5) + 2 * t13 * pkin(8) + 2 * (mrSges(4,3) - t10) * pkin(7)) * t1 ^ 2 + Ifges(3,3) + t12 * pkin(2) ^ 2; t12 * pkin(2) + mrSges(3,1); ((-pkin(7) - t2) * m(6) - t4 * pkin(7) - mrSges(4,3) + t11) * t1 + mrSges(3,2); Ifges(4,1) - Ifges(4,2) - t16; Ifges(4,4); t10 * pkin(3) + Ifges(4,5); Ifges(4,6); Ifges(4,3) + t16; t3 * pkin(3) + mrSges(4,1); mrSges(4,2); Ifges(5,1) - Ifges(5,2) - t15; Ifges(5,4); (-m(6) * pkin(9) - mrSges(6,3)) * pkin(4) + Ifges(5,5); Ifges(5,6); Ifges(5,3) + t15; pkin(4) * m(6) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t5;
