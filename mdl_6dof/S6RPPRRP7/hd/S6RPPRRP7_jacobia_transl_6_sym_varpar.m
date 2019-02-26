% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:36
% EndTime: 2019-02-26 20:33:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (97->30), mult. (112->37), div. (0->0), fcn. (127->7), ass. (0->24)
t26 = pkin(5) + r_i_i_C(1);
t23 = r_i_i_C(3) + qJ(6) + pkin(8);
t8 = pkin(9) + qJ(4);
t7 = cos(t8);
t25 = t23 * t7;
t12 = sin(qJ(5));
t14 = cos(qJ(5));
t5 = pkin(5) * t14 + pkin(4);
t17 = r_i_i_C(1) * t14 - r_i_i_C(2) * t12 + t5;
t6 = sin(t8);
t24 = t17 * t7 + t23 * t6;
t15 = cos(qJ(1));
t22 = t12 * t15;
t13 = sin(qJ(1));
t21 = t13 * t12;
t20 = t13 * t14;
t19 = t14 * t15;
t18 = pkin(5) * t12 + pkin(1) + pkin(7) + qJ(3);
t3 = t6 * t22 + t20;
t1 = -t6 * t21 + t19;
t16 = pkin(3) * sin(pkin(9)) - t25 + t6 * t5 + qJ(2);
t4 = t6 * t19 - t21;
t2 = t6 * t20 + t22;
t9 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t18 * t13 + t16 * t15, t13, t15, t24 * t13, -t2 * r_i_i_C(2) + t26 * t1, -t13 * t7; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t13 + t18 * t15, -t15, t13, -t24 * t15, t4 * r_i_i_C(2) + t26 * t3, t15 * t7; 0, 0, 0, -t17 * t6 + t25 (-r_i_i_C(2) * t14 - t26 * t12) * t7, t6;];
Ja_transl  = t9;
