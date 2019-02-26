% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:36
% EndTime: 2019-02-26 20:18:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (194->36), mult. (218->59), div. (0->0), fcn. (262->12), ass. (0->27)
t22 = qJ(3) + qJ(4);
t19 = qJ(5) + t22;
t15 = sin(t19);
t16 = cos(t19);
t24 = sin(pkin(6));
t25 = cos(pkin(12));
t33 = t24 * t25;
t23 = sin(pkin(12));
t28 = cos(qJ(2));
t26 = cos(pkin(6));
t27 = sin(qJ(2));
t31 = t26 * t27;
t8 = t23 * t28 + t25 * t31;
t38 = (-t8 * t15 - t16 * t33) * r_i_i_C(1) + (t15 * t33 - t8 * t16) * r_i_i_C(2);
t10 = -t23 * t31 + t25 * t28;
t34 = t23 * t24;
t37 = (-t10 * t15 + t16 * t34) * r_i_i_C(1) + (-t10 * t16 - t15 * t34) * r_i_i_C(2);
t32 = t24 * t27;
t36 = (-t15 * t32 + t26 * t16) * r_i_i_C(1) + (-t26 * t15 - t16 * t32) * r_i_i_C(2);
t35 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t30 = t26 * t28;
t18 = cos(t22);
t13 = pkin(4) * t18 + cos(qJ(3)) * pkin(3);
t29 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + pkin(2) + t13;
t17 = sin(t22);
t12 = -pkin(4) * t17 - sin(qJ(3)) * pkin(3);
t1 = [0, t35 * t10 + t29 * (-t23 * t30 - t25 * t27) t10 * t12 + t13 * t34 + t37 (-t10 * t17 + t18 * t34) * pkin(4) + t37, t37, 0; 0, t35 * t8 + t29 * (-t23 * t27 + t25 * t30) t8 * t12 - t13 * t33 + t38 (-t17 * t8 - t18 * t33) * pkin(4) + t38, t38, 0; 1 (t35 * t27 + t29 * t28) * t24, t12 * t32 + t26 * t13 + t36 (-t17 * t32 + t18 * t26) * pkin(4) + t36, t36, 0;];
Ja_transl  = t1;
