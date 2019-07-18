% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Ja_transl [3x5]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S5RPRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobia_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobia_transl_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:25
% EndTime: 2019-07-18 13:26:25
% DurationCPUTime: 0.11s
% Computational Cost: add. (54->30), mult. (150->59), div. (0->0), fcn. (184->8), ass. (0->30)
t10 = sin(qJ(4));
t31 = t10 * r_i_i_C(3);
t11 = sin(qJ(3));
t9 = sin(qJ(5));
t30 = t11 * t9;
t13 = cos(qJ(5));
t29 = t11 * t13;
t14 = cos(qJ(4));
t28 = t11 * t14;
t16 = cos(qJ(1));
t27 = t11 * t16;
t12 = sin(qJ(1));
t26 = t12 * t10;
t15 = cos(qJ(3));
t25 = t14 * t15;
t24 = t14 * t16;
t23 = t16 * t10;
t22 = t13 * r_i_i_C(1) - t9 * r_i_i_C(2);
t4 = t12 * t25 - t23;
t21 = -t12 * t30 - t13 * t4;
t20 = t12 * t29 - t4 * t9;
t19 = t13 * t15 + t9 * t28;
t18 = -t13 * t28 + t15 * t9;
t17 = t18 * r_i_i_C(1) + t19 * r_i_i_C(2) - t11 * t31;
t6 = t15 * t24 + t26;
t5 = -t12 * t14 + t15 * t23;
t3 = -t15 * t26 - t24;
t2 = t6 * t13 + t9 * t27;
t1 = t13 * t27 - t6 * t9;
t7 = [t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t3 * r_i_i_C(3) + t16 * qJ(2), t12, t17 * t16, r_i_i_C(3) * t6 - t22 * t5, r_i_i_C(1) * t1 - t2 * r_i_i_C(2); r_i_i_C(1) * t2 + r_i_i_C(2) * t1 + r_i_i_C(3) * t5 + qJ(2) * t12, -t16, t17 * t12, r_i_i_C(3) * t4 + t22 * t3, t20 * r_i_i_C(1) + t21 * r_i_i_C(2); 0, 0, (t13 * t25 + t30) * r_i_i_C(1) + (-t9 * t25 + t29) * r_i_i_C(2) + t15 * t31, (r_i_i_C(3) * t14 - t22 * t10) * t11, -t19 * r_i_i_C(1) + t18 * r_i_i_C(2);];
Ja_transl  = t7;
