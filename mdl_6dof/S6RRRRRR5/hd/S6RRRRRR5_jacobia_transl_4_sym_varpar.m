% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR5_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR5_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobia_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:49:25
% EndTime: 2019-02-26 22:49:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (132->38), mult. (218->61), div. (0->0), fcn. (264->10), ass. (0->33)
t20 = cos(pkin(6));
t22 = sin(qJ(2));
t26 = cos(qJ(1));
t31 = t26 * t22;
t23 = sin(qJ(1));
t25 = cos(qJ(2));
t32 = t23 * t25;
t10 = t20 * t31 + t32;
t18 = qJ(3) + qJ(4);
t16 = sin(t18);
t19 = sin(pkin(6));
t34 = t19 * t26;
t13 = t16 * t34;
t17 = cos(t18);
t40 = (-t10 * t16 - t17 * t34) * r_i_i_C(1) + (-t10 * t17 + t13) * r_i_i_C(2);
t30 = t26 * t25;
t33 = t23 * t22;
t12 = -t20 * t33 + t30;
t35 = t19 * t23;
t5 = -t12 * t16 + t17 * t35;
t6 = t12 * t17 + t16 * t35;
t39 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t36 = t19 * t22;
t38 = (-t16 * t36 + t20 * t17) * r_i_i_C(1) + (-t20 * t16 - t17 * t36) * r_i_i_C(2);
t37 = r_i_i_C(3) + pkin(10) + pkin(9);
t21 = sin(qJ(3));
t29 = pkin(3) * t21 + pkin(8);
t24 = cos(qJ(3));
t15 = t24 * pkin(3) + pkin(2);
t28 = t17 * r_i_i_C(1) - t16 * r_i_i_C(2) + t15;
t11 = t20 * t32 + t31;
t9 = -t20 * t30 + t33;
t1 = [-t23 * pkin(1) + t13 * r_i_i_C(1) - t37 * t9 - t28 * t10 + (r_i_i_C(2) * t17 + t29) * t34, -t11 * t28 + t12 * t37 (-t12 * t21 + t24 * t35) * pkin(3) + t39, t39, 0, 0; t26 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t11 * t37 + t12 * t15 + t29 * t35, t10 * t37 - t28 * t9 (-t10 * t21 - t24 * t34) * pkin(3) + t40, t40, 0, 0; 0 (t22 * t37 + t25 * t28) * t19 (t20 * t24 - t21 * t36) * pkin(3) + t38, t38, 0, 0;];
Ja_transl  = t1;
