% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:05
% EndTime: 2019-02-26 20:04:06
% DurationCPUTime: 0.10s
% Computational Cost: add. (146->31), mult. (168->51), div. (0->0), fcn. (205->12), ass. (0->27)
t21 = qJ(3) + pkin(12);
t18 = qJ(5) + t21;
t16 = sin(t18);
t17 = cos(t18);
t23 = sin(pkin(6));
t24 = cos(pkin(11));
t32 = t23 * t24;
t22 = sin(pkin(11));
t27 = cos(qJ(2));
t25 = cos(pkin(6));
t26 = sin(qJ(2));
t30 = t25 * t26;
t8 = t22 * t27 + t24 * t30;
t37 = (-t8 * t16 - t17 * t32) * r_i_i_C(1) + (t16 * t32 - t17 * t8) * r_i_i_C(2);
t10 = -t22 * t30 + t24 * t27;
t33 = t22 * t23;
t36 = (-t10 * t16 + t17 * t33) * r_i_i_C(1) + (-t10 * t17 - t16 * t33) * r_i_i_C(2);
t31 = t23 * t26;
t35 = (-t16 * t31 + t17 * t25) * r_i_i_C(1) + (-t16 * t25 - t17 * t31) * r_i_i_C(2);
t34 = r_i_i_C(3) + pkin(9) + qJ(4) + pkin(8);
t29 = t25 * t27;
t13 = pkin(4) * cos(t21) + cos(qJ(3)) * pkin(3);
t28 = r_i_i_C(1) * t17 - r_i_i_C(2) * t16 + pkin(2) + t13;
t12 = -pkin(4) * sin(t21) - sin(qJ(3)) * pkin(3);
t9 = t22 * t29 + t24 * t26;
t7 = t22 * t26 - t24 * t29;
t1 = [0, t34 * t10 - t28 * t9, t10 * t12 + t13 * t33 + t36, t9, t36, 0; 0, -t28 * t7 + t34 * t8, t12 * t8 - t13 * t32 + t37, t7, t37, 0; 1 (t34 * t26 + t28 * t27) * t23, t12 * t31 + t13 * t25 + t35, -t23 * t27, t35, 0;];
Ja_transl  = t1;
