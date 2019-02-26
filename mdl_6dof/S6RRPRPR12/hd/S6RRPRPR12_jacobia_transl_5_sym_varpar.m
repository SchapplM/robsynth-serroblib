% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR12_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR12_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:12
% EndTime: 2019-02-26 21:44:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (127->40), mult. (232->62), div. (0->0), fcn. (286->10), ass. (0->30)
t20 = cos(qJ(4));
t33 = t20 * pkin(4) + pkin(3) + pkin(8);
t14 = sin(pkin(6));
t19 = sin(qJ(1));
t32 = t14 * t19;
t21 = cos(qJ(2));
t31 = t14 * t21;
t22 = cos(qJ(1));
t30 = t14 * t22;
t18 = sin(qJ(2));
t29 = t19 * t18;
t28 = t19 * t21;
t27 = t22 * t18;
t26 = t22 * t21;
t25 = -r_i_i_C(3) - qJ(5) - pkin(9) - pkin(2);
t17 = sin(qJ(4));
t24 = t17 * pkin(4) + qJ(3);
t13 = qJ(4) + pkin(11);
t11 = sin(t13);
t12 = cos(t13);
t23 = t11 * r_i_i_C(1) + t12 * r_i_i_C(2) + t24;
t15 = cos(pkin(6));
t7 = t12 * t30;
t6 = -t15 * t29 + t26;
t5 = t15 * t28 + t27;
t4 = t15 * t27 + t28;
t3 = -t15 * t26 + t29;
t2 = t5 * t11 + t12 * t32;
t1 = -t11 * t32 + t5 * t12;
t8 = [-t19 * pkin(1) + t7 * r_i_i_C(1) - t23 * t3 + (-t11 * r_i_i_C(2) + t33) * t30 + t25 * t4, t23 * t6 + t25 * t5, t5, t1 * r_i_i_C(1) - t2 * r_i_i_C(2) + (-t17 * t32 + t20 * t5) * pkin(4), t6, 0; t22 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t24 * t5 - t25 * t6 + t33 * t32, t23 * t4 + t25 * t3, t3 (t11 * t30 + t3 * t12) * r_i_i_C(1) + (-t3 * t11 + t7) * r_i_i_C(2) + (t17 * t30 + t3 * t20) * pkin(4), t4, 0; 0 (t23 * t18 - t25 * t21) * t14, -t31 (-t15 * t11 - t12 * t31) * r_i_i_C(1) + (t11 * t31 - t15 * t12) * r_i_i_C(2) + (-t15 * t17 - t20 * t31) * pkin(4), t14 * t18, 0;];
Ja_transl  = t8;
