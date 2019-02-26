% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR6_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR6_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:17:36
% EndTime: 2019-02-26 21:17:36
% DurationCPUTime: 0.09s
% Computational Cost: add. (71->24), mult. (85->33), div. (0->0), fcn. (95->7), ass. (0->22)
t22 = pkin(8) + r_i_i_C(3);
t8 = pkin(11) + qJ(3);
t6 = sin(t8);
t17 = t22 * t6;
t7 = cos(t8);
t23 = t17 + t7 * pkin(3) + cos(pkin(11)) * pkin(2) + pkin(1);
t10 = sin(qJ(4));
t13 = cos(qJ(1));
t21 = t10 * t13;
t11 = sin(qJ(1));
t20 = t11 * t10;
t12 = cos(qJ(4));
t19 = t11 * t12;
t18 = t12 * t13;
t15 = r_i_i_C(1) * t12 - r_i_i_C(2) * t10 + pkin(3);
t14 = -t15 * t6 + t22 * t7;
t9 = -pkin(7) - qJ(2);
t4 = t7 * t18 + t20;
t3 = -t7 * t21 + t19;
t2 = -t7 * t19 + t21;
t1 = t7 * t20 + t18;
t5 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t23 * t11 - t13 * t9, t11, t14 * t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t11 * t9 + t23 * t13, -t13, t14 * t11, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0, 0; 0, 0, t15 * t7 + t17 (-r_i_i_C(1) * t10 - r_i_i_C(2) * t12) * t6, 0, 0;];
Ja_transl  = t5;
